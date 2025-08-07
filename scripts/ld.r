library(ieugwasr)
library(pwr)
library(dplyr)
library(purrr)
library(ggplot2)

# Determine the maximum LD rsq between two variants based on their allele frequencies
calc_max_rsq <- function(p1, p2) {
    p12 <- min(p1, p2)
    D <- p12 - p1 * p2
    rsq <- D^2 / (p1 * (1-p1) * p2 * (1-p2))
    return(rsq)
}

# Determine the statistical power for asociation of a rare variant based on its max LD with a common causal variant
get_power <- function(rsid, rare_freq, n, maf, rsq, sig.level=5e-8) {
    mrsq <- calc_max_rsq(rare_freq, maf)
    total_rsq <- mrsq * rsq
    pow <- pwr.r.test(n=n, r=sqrt(total_rsq), sig.level=sig.level)
    return(tibble(
        rsid, rare_freq, n, maf, rsq, sig.level, total_rsq, mrsq, power=pow$power
    ))
}

# Get the GWAS hits for height from UK Biobank
# Add variance explained by each trait
a <- tophits("ukb-b-10787")
a$rsq <- a$beta^2 * 2 * a$eaf * (1-a$eaf)
a$maf <- a$eaf
a$maf[a$maf > 0.5] <- 1 - a$maf[a$maf > 0.5]

# Calculate power for false positives for a rare variant in LD with any of the other traits
param <- expand.grid(
    rsid=a$rsid, 
    rare_freq = c(0.01, 0.005, 0.001, 0.0005, 0.0001), 
    n=450000
) %>% 
    inner_join(., a %>% select(rsid, maf, rsq))


res <- purrr::pmap(param, get_power) %>% bind_rows()

res %>% group_by(rare_freq) %>%
summarise(power = sum(power > 0.8)/n())



