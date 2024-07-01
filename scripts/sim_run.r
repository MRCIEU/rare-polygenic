library(here)
source(here("scripts/sim_functions.r"))

# npoly <- 100
# ngen <- 10
# nfam <- 250
# h2 <- 0.8
# rho <- 0.4

param <- expand.grid(
    ngen = 15,
    nfam = 250,
    npoly = 100,
    h2 = c(0, 1),
    rho = c(0, 0.2, 0.4),
    sims = 1:50
)

o <- mclapply(1:nrow(param), \(i) {
    message(i)
    fam_sim(
        param$ngen[i],
        param$nfam[i],
        param$npoly[i],
        param$h2[i],
        param$rho[i]
    ) %>% mutate(sim=param$sims[i])
}, mc.cores=200) %>% bind_rows()

saveRDS(o, file=here("data", "res.rds"))

