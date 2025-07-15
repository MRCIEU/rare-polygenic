library(dplyr)
library(furrr)
library(here)

source(here("scripts/functions.r"))

#' Run the simulation to estimate correlations between variants and trait / prs
#'
#' @param nsnp Number of SNPs
#' @param nfam Number of families
#' @param h2 Heritability of the trait
#' @param rho Assortative mating coefficient
#' @param nrare Number of rare variants
#' @param rep Replicate number
#' @param ngen Number of generations
#'
#' @return A data frame containing the correlation between the PRS and the variants split by each type of variant
run_sim <- function(nsnp, nfam, h2, rho, nrare, rep=NULL, ngen=8, nrare_per_family=2, family_rank=1) {
    betas <- rnorm(nsnp)
    dat <- make_founders(betas = betas, nfam = nfam, h2 = h2, rho = rho, nrare = nrare, nrare_per_family = nrare_per_family, family_rank = family_rank)
    l <- list()
    for(i in 1:ngen) {
        dat <- create_generation(dat, betas, h2 = h2, rho = rho, reset_rare = FALSE)
        cdat <- collapse_dat(dat)
        l[[i]] <- bind_rows(
            tibble(gen=i, r=get_cors(cdat, wh="rare", out = "prs") %>% drop, h2=h2, rho=rho, rep=rep) %>% mutate(what="rare", out = "prs"),
            tibble(gen=i, r=get_cors(cdat, wh="common", out = "prs") %>% drop, h2=h2, rho=rho, rep=rep) %>% mutate(what="common", out = "prs"),
            tibble(gen=i, r=get_cors(cdat, wh="prs", out = "prs") %>% drop, h2=h2, rho=rho, rep=rep) %>% mutate(what="prs", out = "prs"),
            tibble(gen=i, r=get_cors(cdat, wh="rare", out = "y") %>% drop, h2=h2, rho=rho, rep=rep) %>% mutate(what="rare", out = "y"),
            tibble(gen=i, r=get_cors(cdat, wh="common", out = "y") %>% drop, h2=h2, rho=rho, rep=rep) %>% mutate(what="common", out = "y"),
            tibble(gen=i, r=get_cors(cdat, wh="prs", out = "y") %>% drop, h2=h2, rho=rho, rep=rep) %>% mutate(what="prs", out = "y")
        )
    }
    l <- bind_rows(l)
    return(l)
}

#' Estimate the correlation between the trait/PRS and the rare variants
#' @param datc Output from collapse_dat
#' @param wh Type of variant to use
#' @return The correlation between the PRS and the rare variants
get_cors <- function(datc, wh="rare", out="prs") {
    g <- datc$g
    prs <- datc$x[[out]]
    rare <- g[,subset(datc$map, what == wh)$snp]
    rare <- rare[, colSums(rare) > 0]
    cor(rare, prs)
}

####

## Example run
l1 <- run_sim(500, 2000, 0.8, 0.4, 1000, 1, nrare_per_family = 2, family_rank = 1)

# Plot results
p1 <- l1 %>% 
    group_by(gen, what, out) %>%
    summarise(r=mean(r^2)) %>%
    ggplot(., aes(x=as.factor(gen), y=r, colour=as.factor(what))) + geom_line(aes(group=what)) +
    facet_grid(out ~ ., labeller = label_both) + scale_colour_brewer(type = "qual")
ggplot2::ggsave(here("results/run_sim_example2.png"), p1, width=10, height=6)

l2 <- run_sim(500, 2000, 0.8, 0.4, 1000, 1, nrare_per_family = 1)
p1 <- l2 %>% 
    group_by(gen, what, out) %>%
    summarise(r=mean(r^2)) %>%
    ggplot(., aes(x=as.factor(gen), y=r, colour=as.factor(what))) + geom_line(aes(group=what)) +
    facet_grid(out ~ ., labeller = label_both) + scale_colour_brewer(type = "qual")
ggplot2::ggsave(here("results/run_sim_example1.png"), p1, width=10, height=6)

l3 <- run_sim(500, 2000, 0.8, 0.4, 1000, 1, nrare_per_family = 1, family_rank = 0)
p1 <- l3 %>% 
    group_by(gen, what, out) %>%
    summarise(r=mean(r^2)) %>%
    ggplot(., aes(x=as.factor(gen), y=r, colour=as.factor(what))) + geom_line(aes(group=what)) +
    facet_grid(out ~ ., labeller = label_both) + scale_colour_brewer(type = "qual")
ggplot2::ggsave(here("results/run_sim_example3.png"), p1, width=10, height=6)

# Full simulation
param <- 
    expand.grid(
        h2=c(0, 0.8, 0.4), 
        nsnp=500, 
        rho=c(0, 0.4, 0.8), 
        nfam=c(2000), 
        nrare=c(10000), 
        rep=1:10,
        family_rank=seq(0, 1, length.out=10),
        nrare_per_family=c(1)
    )
dim(param)
plan(multicore, workers=10)
options <- furrr_options(seed=TRUE)
l <- future_pmap(param, run_sim, .progress=TRUE, .options=options)
l <- bind_rows(l)
save(l, file=here("results/run_sim2.RData"))
