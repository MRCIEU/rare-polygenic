# Simulations for polygenic confounding of rare variants

Members in extended families are likely to be more similar at loci influencing traits under assortative mating, which may lead to polygenic confounding of rare variants that tag those families.

These simulations aims to evaluate if this mechanism induces confounding of rare variant associations and if it does, how many generations that confounding lasts for.

## To run

In R:

```r
remv::restore()
```

Or optionally use renv with mamba:

```
mamba create -n renv r-renv
mamba activate renv
# In R:
renv::restore()
```

On ieu-p1:

```sh
Rscript scripts/sim_run.r
Rscript scripts/figs.r
```
