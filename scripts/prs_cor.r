library(here)
library(purrr)
source(here("scripts/functions.r"))

find_nth_gen_relatives <- function(kin, gen) {
    diag(kin) <- 0
    coef <- 0.5^gen
    coef_low <- coef * 0.7
    coef_high <- coef * 1.3
    ind <- which(kin > coef_low & kin < coef_high, arr.ind = TRUE)
    dat <- tibble(
        i = rownames(kin)[ind[, 1]],
        j = colnames(kin)[ind[, 2]],
        coef = kin[ind]
    ) %>%
        filter(i != j)
    return(dat)
}

get_fam_correlation <- function(dat, kin, gens=1:10) {
    lapply(gens, \(i) {
        rels <- find_nth_gen_relatives(kin, i)
        r <- cor(dat$ped$prs[match(rels$i, dat$ped$id)], dat$ped$prs[match(rels$j, dat$ped$id)])
        tibble(i, n=nrow(rels), r=r)
    }) %>% bind_rows()
}


sim <- function(nfam, h2, rho, nsnp) {
    args <- as.list(environment()) %>% as_tibble()
    betas <- rnorm(nsnp)
    dat <- make_founders(betas, nfam, h2, rho)
    for(i in 1:10){
        dat <- create_generation(dat, betas, h2 = h2, rho = rho, reset_rare = FALSE)
    }
    ped <- make_ped(dat)
    kin_all <- kinship(ped)

    n <- nrow(kin_all)
    kin <- kin_all[(n-nfam * 6 + 1):22000, (n-nfam * 10 + 1):22000]
    kin <- kin * 2
    famr <- get_fam_correlation(dat, kin, 1:10)
    famr %>% bind_cols(., args)
}


param <- expand.grid(
    h2 = c(0, 0.4, 0.8),
    rho = c(0, 0.4, 0.8),
    nfam = 1000,
    nsnp=100
)

res <- param %>% pmap(sim, .progress=TRUE) %>% bind_rows()
save(res, file=here("results", "prs_cor.RData"))

# Results from empirical analysis
emp <- tibble(
    i = 1:9,
    r=c(0.53, 0.33, 0.25, 0.077, 0.12, 0.044, 0.019, -0.0062, 0.05)
)

res %>%
    mutate(h2lab = paste("h^2 ==", h2), rholab = paste("rho ==", rho)) %>%
    ggplot(aes(i, r, colour)) +
    geom_line() +
    geom_line(data=emp, aes(i, r), linetype="dashed") +
    facet_grid(h2lab ~ rholab, labeller=label_parsed) +
    theme_minimal() +
    labs(x="Degree of relatedness", y="PRS correlation") +
    scale_x_continuous(breaks=1:9)
ggsave(file=here("results", "prs_cor.pdf"), width=8, height=8)
