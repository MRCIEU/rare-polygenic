# #!/bin/bash

# datadir="/mnt/storage/private/mrcieu/data/ukbiobank/phenotypic/applications/81499/released/2022-06-07/data"

# head -n 1 $datadir/data.51913.csv | tr ',' '\n' > ukb.datafields

# less ukb.datafields 

# grep -n "^\"129-" ukb.datafields | cut -d ":" -f 1 
# grep -n "^\"130-" ukb.datafields | cut -d ":" -f 1 


# cut -d "," -f 1,435-440 $datadir/data.51913.csv > ukb.coords
# cut -d "," -f 10005-10044 $datadir/data.51913.csv > ukb.pcs






library(dplyr)
library(data.table)
library(ggplot2)
library(here)
library(nnet)
library(caret)
library(parallel)



generate_conf <- function(coord1, coord2, centroid, sharpness) {
    # Calculate the distance from the centroid
    distance_from_centroid <- sqrt((coord1 - centroid[1])^2 + (coord2 - centroid[2])^2)
    # Frequency of the variant reduces with distance according to sharpness
    frequency <- distance_from_centroid - min(distance_from_centroid, na.rm=T)
    frequency <- frequency^sharpness
    frequency <- 1 - frequency / max(frequency, na.rm=T)    
}

generate_variant_based_on_location <- function(conf)  {
    variant <- rbinom(length(conf), size = 2, prob = conf)
    return(variant)
}

generate_phen <- function(conf, variant, b_conf, b_variant) {
    phen <- scale(conf) * b_conf + scale(variant) * b_variant + rnorm(length(conf), 0, sqrt(1 - b_conf^2 - b_variant^2))
    return(drop(phen))
}

dgm <- function(centroid, sharpness, b_conf, b_variant, coord1=dat$northing, coord2=dat$easting) {
    conf <- generate_conf(coord1, coord2, centroid, sharpness)
    variant <- generate_variant_based_on_location(conf)
    phen <- generate_phen(conf, variant, b_conf, b_variant)
    return(tibble(variant=variant, conf=conf, phen=phen, northing=coord1, easting=coord2))
}

hexbins_adjustment <- function(phendat) {
    bins <- hexbin(phendat$easting, phendat$northing, 
                   xbins = 30, 
                   IDs = TRUE)
    phendat$hex_id <- bins@cID
    lm(phen ~ as.factor(hex_id), data=phendat)$residuals -> phendat$phen_resid
    summary(lm(phen_resid ~ variant, data=phendat))
}

hexbins_adjustment(phendat)

estimation <- function(dat, phendat) {
    mod1 <- summary(lm(phen ~ scale(variant) + conf, data=phendat))
    mod2 <- summary(lm(phen ~ scale(variant), data=phendat))
    mod3 <- summary(lm(phen ~ scale(variant) + as.matrix(dat[,3:42]), data=phendat))
    mod4 <- summary(lm(phen ~ scale(variant) + northing + easting, data=phendat))
    coord_pred <- fit_genetic_nn(
        y = phendat$phen,
        genetic_pcs = phendat %>% select(northing, easting),
        hidden_units = 10,
        decay = 0.001,
        scale_data = TRUE,
        maxit=25
    )$predictions
    mod5 <- summary(lm(phen ~ scale(variant) + coord_pred, data=phendat))
    tibble(
        model = c("variant + conf", "variant", "variant + PCs", "variant + coords", "variant + coord_nn"),
        freq = mean(phendat$variant)/2,
        var_beta = c(mod1$coefficients[2,1], mod2$coefficients[2,1], mod3$coefficients[2,1], mod4$coefficients[2,1], mod5$coefficients[2,1]),
        var_se = c(mod1$coefficients[2,2], mod2$coefficients[2,2], mod3$coefficients[2,2], mod4$coefficients[2,2], mod5$coefficients[2,2]),
        var_p = c(mod1$coefficients[2,4], mod2$coefficients[2,4], mod3$coefficients[2,4], mod4$coefficients[2,4], mod5$coefficients[2,4]),
        r2 = c(mod1$r.squared, mod2$r.squared, mod3$r.squared, mod4$r.squared, mod5$r.squared),
        adj_r2 = c(mod1$adj.r.squared, mod2$adj.r.squared, mod3$adj.r.squared, mod4$adj.r.squared, mod5$adj.r.squared),
        n = nrow(phendat)
    )
}

run_sim <- function(centroid, sharpness, b_conf, b_variant, rep=1) {
    args <- tibble(centroid=paste(centroid, collapse=","), sharpness=sharpness, b_conf=b_conf, b_variant=b_variant, rep=rep)
    phendat <- dgm(centroid, sharpness, b_conf, b_variant)
    est <- estimation(dat, phendat)
    bind_cols(args, est)
}

fit_genetic_nn <- function(y, genetic_pcs, 
                          hidden_units = 10, 
                          decay = 0.001,
                          train_prop = 0.8,
                          scale_data = TRUE,
                          set_seed = 123, maxit=100) {
  
  # Set seed for reproducibility
  set.seed(set_seed)
  
  # Combine data
  data <- data.frame(y = y, genetic_pcs)
  
  # Remove any rows with missing values
  data <- na.omit(data)
  n_samples <- nrow(data)
  
#   cat("Sample size after removing NAs:", n_samples, "\n")
#   cat("Number of genetic PCs:", ncol(genetic_pcs), "\n")
  
  # Split into training and testing sets
  train_indices <- sample(1:n_samples, size = floor(train_prop * n_samples))
  train_data <- data[train_indices, ]
  test_data <- data[-train_indices, ]
  
  # Scale the predictors if requested
  if (scale_data) {
    # Calculate scaling parameters from training data
    pc_means <- apply(train_data[, -1], 2, mean)
    pc_sds <- apply(train_data[, -1], 2, sd)
    
    # Scale training data
    train_data[, -1] <- scale(train_data[, -1])
    
    # Scale test data using training parameters
    test_data[, -1] <- scale(test_data[, -1], center = pc_means, scale = pc_sds)
  }
  
  # Fit neural network
#   cat("Fitting neural network with", hidden_units, "hidden units...\n")
  
  nn_model <- nnet(y ~ ., 
                   data = train_data,
                   size = hidden_units,
                   decay = decay,
                   linout = TRUE,  # Linear output for regression
                   trace = FALSE,  # Suppress training output
                   maxit = maxit)   # Maximum iterations
  
  # Make predictions on all data
  if (scale_data) {
    # Scale full dataset for predictions
    full_data_scaled <- data
    full_data_scaled[, -1] <- scale(data[, -1], center = pc_means, scale = pc_sds)
    predictions <- predict(nn_model, full_data_scaled[, -1])
  } else {
    predictions <- predict(nn_model, data[, -1])
  }
  
  # Calculate residuals
  residuals <- data$y - predictions
  
  # Calculate performance metrics
  train_predictions <- predict(nn_model, train_data[, -1])
  test_predictions <- predict(nn_model, test_data[, -1])
  
  train_rmse <- sqrt(mean((train_data$y - train_predictions)^2))
  test_rmse <- sqrt(mean((test_data$y - test_predictions)^2))
  train_r2 <- cor(train_data$y, train_predictions)^2
  test_r2 <- cor(test_data$y, test_predictions)^2
  
  # Print performance
#   cat("\nModel Performance:\n")
#   cat("Training RMSE:", round(train_rmse, 4), "\n")
#   cat("Test RMSE:", round(test_rmse, 4), "\n")
#   cat("Training R²:", round(train_r2, 4), "\n")
#   cat("Test R²:", round(test_r2, 4), "\n")
  
  # Return results
  results <- list(
    model = nn_model,
    predictions = predictions,
    residuals = residuals,
    train_indices = train_indices,
    performance = list(
      train_rmse = train_rmse,
      test_rmse = test_rmse,
      train_r2 = train_r2,
      test_r2 = test_r2
    ),
    scaling_params = if (scale_data) list(means = pc_means, sds = pc_sds) else NULL
  )
  
  return(results)
}


coords <- fread(here("ukb.coords")) %>% filter(!is.na(`130-0.0`))
names(coords)[2] <- "northing"
names(coords)[5] <- "easting"

# pcs <- fread("/mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/derived/principal_components/data.pca1-40.plink.txt")
pcs <- fread(here("data.pca1-40.plink.txt"))

# linker <- fread("/mnt/storage/private/mrcieu/data/ukbiobank/phenotypic/applications/81499/released/2022-06-07/data/linker.81499.csv")
linker <- fread(here("linker.81499.csv"))
pcs <- inner_join(pcs, linker, by=c("V1"="ieu"))
pcs
dat <- inner_join(pcs, coords, by=c("app"="eid"))

centroids <- dat %>%
    filter(!is.na(northing) & northing > 0) %>%
    select(northing, easting) %>%
    slice_sample(n=100)

param <- expand.grid(
    centroid = 1:10,
    sharpness = c(0.001, 0.01, 0.1, 0.5),
    b_conf = c(0, 0.3),
    b_variant = c(0, 0.1),
    rep= 1:10
)
dim(param)

mclapply(1:nrow(param), function(i) {
    message(paste("Running row", i, "of", nrow(param)))
    tryCatch(
        run_sim(
            centroid = c(centroids$northing[param$centroid[i]], centroids$easting[param$centroid[i]]),
            sharpness = param$sharpness[i],
            b_conf = param$b_conf[i],
            b_variant = param$b_variant[i],
            rep = param$rep[i]
        ),
        error = function(e) {
            message(paste("Error in row", i, ":", e$message))
            return(NULL)
        }
    )   
}, mc.cores=60) %>% bind_rows() -> res



i <- sample(1:nrow(dat), 1)
phendat <- dgm(c(dat$northing[i], dat$easting[i]), 0.1, 0.3, 0.1)
estimation(dat, phendat)

run_sim(c(dat$northing[i], dat$easting[i]), 0.01, 0.5, 0.1)
run_sim(c(dat$northing[i], dat$easting[i]), 0.01, 0, 0)



saveRDS(res, here("results", "stratification_sim_1.rds"))

ggplot(res, aes(x=freq, y=var_beta - b_variant, colour=model)) +
    geom_point() +
    labs(x="Frequency of variant", y="Estimated beta - true beta") +
    scale_colour_brewer(type = "qual") +
    facet_grid(b_conf ~ b_variant, labeller = label_both) +
    theme_minimal()

p1 <- res %>%
ggplot(., aes(x=as.factor(sharpness), y=var_beta - b_variant, colour=model)) +
    geom_boxplot() +
    labs(x="Frequency of variant", y="Estimated beta - true beta") +
    scale_colour_brewer(type = "qual") +
    facet_grid(b_conf ~ b_variant, labeller = label_both) +
    theme_minimal()
ggsave(p1, file="temp.pdf", width=14, height=7)


i <- sample(1:nrow(dat), 1)
phendat <- dgm(c(dat$northing[i], dat$easting[i]), 0.1, 0.3, 0.1)
estimation(dat, phendat)
temp <- stat_summary_hex(aes(z=phen, x=easting, y=northing), fun = mean, data=phendat)
temp2 <- layer_data(ggplot(phendat) + temp, 1)

str(temp)
# Get bins from temp
str(temp$data)

bins <- hexbin(temp$data$easting, temp$data$northing, 
               xbins = 30, 
               IDs = TRUE)

str(bins)
length(table(bins@cID))





# Neural Network to predict spatially structured variable from genetic PCs
# and return residuals

# Load required libraries

# Function to fit neural network and return residuals


a <- fit_genetic_nn(
    y = phendat$phen,
    genetic_pcs = dat[, 3:42],
    hidden_units = 10,
    decay = 0.001,
    scale_data = TRUE
)


plot(a$predictions, a1$predictions)
plot(a$predictions, phendat$conf)
plot(a$predictions, phendat$phen)
plot(a1$predictions, phendat$phen)

summary(lm(a$predictions ~ phendat$phen))
summary(lm(a$predictions ~ phendat$phen))


str(a)
var(a$residuals)
var(a$predictions)

summary(lm(phendat$phen ~ as.matrix(dat[, 3:42])))


plot_uk_density <- function(phendat) {
    p1 <- ggplot(phendat, aes(x=easting, y=northing)) +
        stat_summary_hex(aes(z=phen), fun = mean) +
        labs(title="Spatial distribution of phenotypes", x="Easting", y="Northing") +
        theme_minimal()
    p2 <- ggplot(phendat, aes(x=easting, y=northing)) +
        stat_summary_hex(aes(z=variant), fun = mean) +
        labs(title="Spatial distribution of variants", x="Easting", y="Northing")
    p3 <- ggplot(phendat, aes(x=easting, y=northing)) +
        stat_summary_hex(aes(z=conf), fun = mean) +
        labs(title="Spatial distribution of confounders", x="Easting", y="Northing")
    gridExtra::grid.arrange(p1, p2, p3, ncol=3)

}

estimation(dat, phendat)

run_sim(c(dat$northing[i], dat$easting[i]), 2, 0.5, 0.1)


i <- sample(1:nrow(dat), 1)
phendat <- dgm(c(dat$northing[i], dat$easting[i]), 0.001, 0.3, 0.1)
plot_uk_density(phendat)

phendat <- dgm(c(dat$northing[i], dat$easting[i]), 0.01, 0.3, 0.1)
plot_uk_density(phendat)

estimation(dat, phendat)

phendat <- dgm(c(dat$northing[i], dat$easting[i]), 0.5, 0.3, 0.1)
plot_uk_density(phendat)

estimation(dat, phendat)


a2 <- fit_genetic_nn(
    y = phendat$phen,
    genetic_pcs = phendat %>% select(northing, easting),
    hidden_units = 10,
    decay = 0.001,
    scale_data = TRUE,
    maxit=200
)

summary(lm(phen ~ northing + easting, data=phendat))

summary(lm(a1$residuals ~ phendat$variant))










# Alternative function using cross-validation for model selection
fit_genetic_nn_cv <- function(y, genetic_pcs, 
                             hidden_units = c(5, 10, 15, 20),
                             decay_values = c(0.001, 0.01, 0.1),
                             cv_folds = 5,
                             scale_data = TRUE,
                             set_seed = 123) {
  
  set.seed(set_seed)
  
  # Prepare data
  data <- data.frame(y = y, genetic_pcs)
  data <- na.omit(data)
  
  # Create cross-validation folds
  folds <- createFolds(data$y, k = cv_folds, list = TRUE)
  
  best_rmse <- Inf
  best_params <- NULL
  
  cat("Performing cross-validation...\n")
  
  # Grid search over hyperparameters
  for (size in hidden_units) {
    for (decay in decay_values) {
      cv_rmse <- numeric(cv_folds)
      
      for (i in 1:cv_folds) {
        train_data <- data[-folds[[i]], ]
        val_data <- data[folds[[i]], ]
        
        # Scale data
        if (scale_data) {
          pc_means <- apply(train_data[, -1], 2, mean)
          pc_sds <- apply(train_data[, -1], 2, sd)
          train_data[, -1] <- scale(train_data[, -1])
          val_data[, -1] <- scale(val_data[, -1], center = pc_means, scale = pc_sds)
        }
        
        # Fit model
        nn <- nnet(y ~ ., data = train_data, size = size, decay = decay, 
                   linout = TRUE, trace = FALSE, maxit = 1000)
        
        # Validate
        pred <- predict(nn, val_data[, -1])
        cv_rmse[i] <- sqrt(mean((val_data$y - pred)^2))
      }
      
      mean_cv_rmse <- mean(cv_rmse)
      
      if (mean_cv_rmse < best_rmse) {
        best_rmse <- mean_cv_rmse
        best_params <- list(size = size, decay = decay)
      }
      
      cat("Size:", size, "Decay:", decay, "CV RMSE:", round(mean_cv_rmse, 4), "\n")
    }
  }
  
  cat("\nBest parameters - Size:", best_params$size, "Decay:", best_params$decay, "\n")
  cat("Best CV RMSE:", round(best_rmse, 4), "\n")
  
  # Fit final model with best parameters
  final_results <- fit_genetic_nn(y, genetic_pcs, 
                                 hidden_units = best_params$size,
                                 decay = best_params$decay,
                                 scale_data = scale_data,
                                 set_seed = set_seed)
  
  final_results$best_params <- best_params
  final_results$cv_rmse <- best_rmse
  
  return(final_results)
}

# Example usage (uncomment and modify as needed):
# 
# # Assuming you have:
# # y - your spatially structured response variable (vector)
# # genetic_pcs - matrix or data frame with 40 genetic principal components
# 
# # Basic usage:
# results <- fit_genetic_nn(y, genetic_pcs, hidden_units = 10, decay = 0.001)
# 
# # Extract residuals:
# y_residuals <- results$residuals
# 
# # With cross-validation for parameter tuning:
# results_cv <- fit_genetic_nn_cv(y, genetic_pcs)
# y_residuals_cv <- results_cv$residuals
# 
# # Plot residuals vs predictions
# plot(results$predictions, results$residuals, 
#      xlab = "Predicted values", ylab = "Residuals",
#      main = "Residuals vs Fitted")
# abline(h = 0, col = "red", lty = 2)
# 
# # Spatial plot of residuals (if you have coordinates)
# # plot(coordinates$x, coordinates$y, col = rainbow(100)[cut(y_residuals, 100)],
# #      pch = 16, main = "Spatial distribution of residuals")

cat("Functions loaded successfully!\n")
cat("Use fit_genetic_nn() for basic neural network fitting\n")
cat("Use fit_genetic_nn_cv() for cross-validated parameter selection\n")





