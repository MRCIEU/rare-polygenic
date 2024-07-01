rule all:
    input:

rule setup:
    input:
        "renv.lock"
    output:
        "scripts/sim_analysis.html"
    shell:
        "Rscript -e 'renv::restore()'"

rule run:
    input:
        "renv.lock",
        "scripts/sim_run.r",
        "scripts/sim_functions.r"
    output:
        "data/res.rds"
    script:
        "scripts/sim_run.r"

rule analysis:
    input:
        "data/res.rds",
        "scripts/sim_analysis.rmd"
    output:
        "scripts/sim_analysis.html"
    script:
        "scripts/sim_analysis.rmd"
