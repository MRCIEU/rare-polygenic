library(dplyr)
library(MASS)
library(ggplot2)
library(here)
# library(kinship2)

# Take a vector of trait values for males and females and pair them such that the correlation between m and f trait values = rho
generate_assortment <- function(m, f, rho) {
    stopifnot(length(m) == length(f))
    require(MASS)
    mvdat <- mvrnorm(n = length(m), mu=c(0,0), Sigma=matrix(c(1,rho,rho,1), 2,2))
    rm <- rank(mvdat[ , 1], ties.method = "first")
    rf <- rank(mvdat[ , 2], ties.method = "first")
    m_order <- order(m)
    f_order <- order(f)
    return(tibble(m = m_order[rm], f=f_order[rf]))
}

#' Generate founder population
#' 
#' Generate a set of males and females who each have a single polygenic trait constructed from a set of causal common variants
#' Each individual also has a set of null rare variants and a set of null common variants
#' The individuals are then paired such that the correlation between their trait values corresponds to the required assortative mating.
#' 
#' @param betas Vector of effect sizes for each SNP contributing to the PRS
#' @param nfam Number of families (2x this number of individuals)
#' @param h2 Heritability of the trait
#' @param rho Assortative mating coefficient
#' @param nrare Number of null rare variants
#' @param ncommon Number of null common variants
#' @param nrare_per_family Whether both parents or just one parent should have a rare variant
#' @param family_rank Which family should have the rare variant, 0 to 1 scale where 0 is the lowest family phenotype and 1 is the highest. Defauly=1
#'
#' @return A list containing the pedigree, genotype matrices for the mothers and fathers, the trait values for the mothers and fathers, the PRS for the mothers and fathers, and a map of the SNPs to their type
make_founders <- function(betas, nfam, h2, rho, nrare=1000, ncommon=1000, nrare_per_family=2, family_rank=1) {
    npoly <- length(betas)
    # Generate poly variants
    g_mother1 <- sapply(1:npoly, \(i) rbinom(nfam, 1, 0.5))
    g_mother2 <- sapply(1:npoly, \(i) rbinom(nfam, 1, 0.5))
    g_father1 <- sapply(1:npoly, \(i) rbinom(nfam, 1, 0.5))
    g_father2 <- sapply(1:npoly, \(i) rbinom(nfam, 1, 0.5))

    prs_mother <- (g_mother1[,1:ncol(g_mother1)] + g_mother2[,1:ncol(g_mother2)]) %*% betas
    prs_father <- (g_father1[,1:ncol(g_father1)] + g_father2[,1:ncol(g_father2)]) %*% betas
    y_mother <- scale(prs_mother) * sqrt(h2) + rnorm(nfam, 0, sqrt(1 - h2)) 
    y_father <- scale(prs_father) * sqrt(h2) + rnorm(nfam, 0, sqrt(1 - h2))

    m <- generate_assortment(y_mother, y_father, rho)
    g_mother1 <- g_mother1[m$m,]
    g_mother2 <- g_mother2[m$m,]
    g_father1 <- g_father1[m$f,]
    g_father2 <- g_father2[m$f,]
    y_mother <- y_mother[m$m]
    y_father <- y_father[m$f]
    prs_mother <- prs_mother[m$m]
    prs_father <- prs_father[m$f]


    # Generate rare variants
    g_mother1 <- cbind(g_mother1, matrix(0, nfam, nrare))
    g_mother2 <- cbind(g_mother2, matrix(0, nfam, nrare))
    g_father1 <- cbind(g_father1, matrix(0, nfam, nrare))
    g_father2 <- cbind(g_father2, matrix(0, nfam, nrare))

    # Find rank family

    y_fam <- y_mother + y_father
    quant <- quantile(y_fam, family_rank)
    i <- which.min(abs(y_fam - quant))

    g_mother1[i, (npoly+1):ncol(g_mother1)] <- 1
    g_father1[i, (npoly+1):ncol(g_father1)] <- ifelse(nrare_per_family > 1, 1, 0)

    # Common variants
    g_mother1 <- cbind(g_mother1, sapply(1:ncommon, \(i) rbinom(nfam, 1, 0.5)))
    g_mother2 <- cbind(g_mother2, sapply(1:ncommon, \(i) rbinom(nfam, 1, 0.5)))
    g_father1 <- cbind(g_father1, sapply(1:ncommon, \(i) rbinom(nfam, 1, 0.5)))
    g_father2 <- cbind(g_father2, sapply(1:ncommon, \(i) rbinom(nfam, 1, 0.5)))

    pedm <- tibble(
            generation=1,
            id=paste(generation, 1:nfam, "mother"),
            sex = 2,
            motherid="0",
            fatherid="0"
        )[m$m,]
    pedf <- tibble(
            generation=1,
            id=paste(generation, 1:nfam, "father"),
            sex = 1,
            motherid="0",
            fatherid="0"
        )[m$f,]

    ped <- bind_rows(pedm, pedf)
    ped$prs <- c(prs_mother, prs_father)
    ped$y <- c(y_mother, y_father)
    map <- tibble(snp=1:ncol(g_mother1), what=c(rep("prs", npoly), rep("rare", nrare), rep("common", ncommon)))
    return(list(ped = ped, g_mother1 = g_mother1, g_mother2 = g_mother2, g_father1 = g_father1, g_father2 = g_father2, x = tibble(id_mother = pedm$id, id_father = pedf$id, y_mother, y_father, prs_mother, prs_father), map = map))
}

#' Create a new generation from a previous generation
#'
#' @param dat output from make_founders or create_generation
#' @param betas Vector of effect sizes for each SNP contributing to the PRS. Can be the same as the betas used to generate the previous generation
#' @param h2 Heritability of the trait
#' @param rho Assortative mating coefficient
#' @param reset_rare If TRUE, reset the rare variants to be new mutations
#' @return A list containing the pedigree, genotype matrices for the mothers and fathers, the trait values for the mothers and fathers, the PRS for the mothers and fathers, and a map of the SNPs to their type
create_generation <- function(dat, betas, h2, rho, reset_rare = FALSE, h2_rare=0) {
    nfam <- nrow(dat$x)
    nsnp <- ncol(dat$g_mother1)

    if(reset_rare) {
        vind <- subset(dat$map, what == "rare")$snp
        fid <- which.max(dat$x$y_mother + dat$x$y_father)
        dat$g_mother1[1:nfam, vind] <- 0
        dat$g_mother2[1:nfam, vind] <- 0
        dat$g_father1[1:nfam, vind] <- 0
        dat$g_father2[1:nfam, vind] <- 0
        dat$g_mother1[fid, 1:nfam] <- 1
        dat$g_father1[fid, 1:nfam] <- 1
    }

    sib1_m <- matrix(0, nfam, nsnp)
    sib1_f <- matrix(0, nfam, nsnp)
    sib2_m <- matrix(0, nfam, nsnp)
    sib2_f <- matrix(0, nfam, nsnp)
    for(i in 1:nsnp) {
        ind <- sample(c(TRUE, FALSE), nfam, replace=TRUE)
        sib1_m[ind, i] <- dat$g_mother1[ind, i]
        sib1_m[!ind, i] <- dat$g_mother2[!ind, i]
        ind <- sample(c(TRUE, FALSE), nfam, replace=TRUE)
        sib1_f[ind, i] <- dat$g_father1[ind, i]
        sib1_f[!ind, i] <- dat$g_father2[!ind, i]
        ind <- sample(c(TRUE, FALSE), nfam, replace=TRUE)
        sib2_m[ind, i] <- dat$g_mother1[ind, i]
        sib2_m[!ind, i] <- dat$g_mother2[!ind, i]
        ind <- sample(c(TRUE, FALSE), nfam, replace=TRUE)
        sib2_f[ind, i] <- dat$g_father1[ind, i]
        sib2_f[!ind, i] <- dat$g_father2[!ind, i]
    }

    vind <- subset(dat$map, what == "prs")$snp
    vind_rare <- subset(dat$map, what == "rare")$snp
    beta_rare <- rnorm(length(vind_rare))
    prs_sib1 <- (sib1_m[,vind] + sib1_f[,vind]) %*% betas
    prs_sib2 <- (sib2_m[,vind] + sib2_f[,vind]) %*% betas
    bv_rare1 <- scale((sib1_m[,vind_rare] + sib1_f[,vind_rare]) %*% beta_rare) * sqrt(h2_rare)
    bv_rare2 <- scale((sib2_m[,vind_rare] + sib2_f[,vind_rare]) %*% beta_rare) * sqrt(h2_rare)
    y_sib1 <- scale(prs_sib1) * sqrt(h2) + rnorm(nfam, 0, sqrt(1 - h2 - h2_rare)) + bv_rare1
    y_sib2 <- scale(prs_sib2) * sqrt(h2) + rnorm(nfam, 0, sqrt(1 - h2 - h2_rare)) + bv_rare2


    # This doesn't exclude inbreeding
    m <- generate_assortment(y_sib1, y_sib2, rho)
    sib1_m <- sib1_m[m$m,]
    sib1_f <- sib1_f[m$m,]
    sib2_m <- sib2_m[m$f,]
    sib2_f <- sib2_f[m$f,]
    y_sib1 <- y_sib1[m$m]
    y_sib2 <- y_sib2[m$f]
    prs_sib1 <- prs_sib1[m$m]
    prs_sib2 <- prs_sib2[m$f]

    gen <- max(dat$ped$generation) + 1
    pedl <- subset(dat$ped, generation == gen-1)
    ped_m <- tibble(
            generation = gen,
            id = paste(generation, 1:nfam, "mother"),
            sex = 2,
            motherid = pedl$id[1:nfam],
            fatherid = pedl$id[(nfam+1):(2*nfam)]
        )[m$m,]
    ped_f <- tibble(
            generation = gen,
            id = paste(generation, 1:nfam, "father"),
            sex = 1,
            motherid = pedl$id[1:nfam],
            fatherid = pedl$id[(nfam+1):(2*nfam)]
        )[m$f,]
    ped_m$prs <- prs_sib1
    ped_f$prs <- prs_sib2
    ped_m$y <- y_sib1
    ped_f$y <- y_sib2
    ped <- bind_rows(dat$ped, ped_m, ped_f)

    return(list(ped = ped, g_mother1 = sib1_m, g_mother2 = sib1_f, g_father1 = sib2_m, g_father2 = sib2_f, x = tibble(id_mother = ped_m$id, id_father = ped_f$id, y_mother=y_sib1, y_father=y_sib2, prs_mother=prs_sib1, prs_father=prs_sib2), map = dat$map))
}

#' Collapse a dataset to a single generation ready for analysis
#' @param dat Output from make_founders or create_generation
collapse_dat <- function(dat) {
    g_mother <- dat$g_mother1 + dat$g_mother2
    g_father <- dat$g_father1 + dat$g_father2
    g <- rbind(g_mother, g_father)
    x <- tibble(id = c(dat$x$id_mother, dat$x$id_father), y = c(dat$x$y_mother, dat$x$y_father), prs = c(dat$x$prs_mother, dat$x$prs_father))
    dat$map$af <- colSums(g) / (2*nrow(g))
    return(list(g = g, x = x, map = dat$map))
}

make_ped <- function(dat, gen=NULL) {
    if(!is.null(gen)) {
        dat$ped <- subset(dat$ped, generation %in% gen)
        mingen <- min(gen)
        dat$ped$motherid[dat$ped$generation == mingen] <- 0
        dat$ped$fatherid[dat$ped$generation == mingen] <- 0
    }
    pedigree(id=dat$ped$id, dadid=dat$ped$fatherid, momid=dat$ped$motherid, sex=dat$ped$sex, missid="0")
}


###########






