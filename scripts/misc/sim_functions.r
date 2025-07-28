library(dplyr)
library(MASS)
library(ggplot2)
library(parallel)
# library(igraph)

# Sample pairings of two phenotype vectors to induce a specified correlation between the two phenotype vectors
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

# Create a 
make_founders <- function(betas, nfam, h2, rho) {
    npoly <- length(betas)
    g_mother1 <- sapply(1:npoly, \(i) rbinom(nfam, 1, 0.5))
    g_mother2 <- sapply(1:npoly, \(i) rbinom(nfam, 1, 0.5))
    g_father1 <- sapply(1:npoly, \(i) rbinom(nfam, 1, 0.5))
    g_father2 <- sapply(1:npoly, \(i) rbinom(nfam, 1, 0.5))
    g_mother1 <- cbind(diag(nfam), g_mother1)
    g_mother2 <- cbind(matrix(0, nfam, nfam), g_mother2)
    g_father1 <- cbind(matrix(0, nfam, nfam), g_father1)
    g_father2 <- cbind(matrix(0, nfam, nfam), g_father2)

    prs_mother <- (g_mother1[,(nfam+1):ncol(g_mother1)] + g_mother2[,(nfam+1):ncol(g_mother2)]) %*% betas
    prs_father <- (g_father1[,(nfam+1):ncol(g_father1)] + g_father2[,(nfam+1):ncol(g_father2)]) %*% betas
    y_mother <- scale(prs_mother) + rnorm(nfam, 0, sqrt(1 - h2))
    y_father <- scale(prs_father) + rnorm(nfam, 0, sqrt(1 - h2))

    m <- generate_assortment(y_mother, y_father, rho)
    g_mother1 <- g_mother1[m$m,]
    g_mother2 <- g_mother2[m$m,]
    g_father1 <- g_father1[m$f,]
    g_father2 <- g_father2[m$f,]
    y_mother <- y_mother[m$m]
    y_father <- y_father[m$f]
    prs_mother <- prs_mother[m$m]
    prs_father <- prs_father[m$f]
    return(list(g_mother1 = g_mother1, g_mother2 = g_mother2, g_father1 = g_father1, g_father2 = g_father2, x = tibble(y_mother, y_father, prs_mother, prs_father)))
}

create_child <- function(dat, betas, h2, rho) {
    nfam <- nrow(dat$x)
    nsnp <- ncol(dat$g_mother1)
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

    prs_sib1 <- (sib1_m[,(nfam+1):nsnp] + sib1_f[,(nfam+1):nsnp]) %*% betas
    prs_sib2 <- (sib2_m[,(nfam+1):nsnp] + sib2_f[,(nfam+1):nsnp]) %*% betas
    y_sib1 <- scale(prs_sib1) + rnorm(nfam, 0, sqrt(1 - h2))
    y_sib2 <- scale(prs_sib2) + rnorm(nfam, 0, sqrt(1 - h2))


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
    return(list(g_mother1 = sib1_m, g_mother2 = sib1_f, g_father1 = sib2_m, g_father2 = sib2_f, x = tibble(y_mother=y_sib1, y_father=y_sib2, prs_mother=prs_sib1, prs_father=prs_sib2)))
}

collapse_dat <- function(dat) {
    g_mother <- dat$g_mother1 + dat$g_mother2
    g_father <- dat$g_father1 + dat$g_father2
    g <- rbind(g_mother, g_father)
    x <- tibble(y = c(dat$x$y_mother, dat$x$y_father), prs = c(dat$x$prs_mother, dat$x$prs_father))
    return(list(g = g, x = x))
}

collapse_dat2 <- function(l) {
    g <- do.call(rbind, lapply(l, \(x) x$g))
    x <- do.call(rbind, lapply(1:length(l), \(i) l[[i]]$x %>% mutate(gen=i)))
    return(list(g = g, x = x))
}

fastAssoc <- function(y, x) {
	index <- is.finite(y) & is.finite(x)
	n <- sum(index)
	y <- y[index]
	x <- x[index]
	vx <- var(x)
	bhat <- cov(y, x) / vx
	ahat <- mean(y) - bhat * mean(x)
	# fitted <- ahat + x * bhat
	# residuals <- y - fitted
	# SSR <- sum((residuals - mean(residuals))^2)
	# SSF <- sum((fitted - mean(fitted))^2)

	rsq <- (bhat * vx)^2 / (vx * var(y))
	fval <- rsq * (n-2) / (1-rsq)
	tval <- sqrt(fval)
    af <- sum(x) / n / 2
	se <- abs(bhat / tval)

	# Fval <- (SSF) / (SSR/(n-2))
	# pval <- pf(Fval, 1, n-2, lowe=F)
	p <- pf(fval, 1, n-2, lowe=F)
	return(list(
		ahat=ahat, bhat=bhat, se=se, fval=fval, pval=p, af = af
	))
}

gwas <- function(y, g) {
	out <- matrix(0, ncol(g), 6)
	for(i in 1:ncol(g))
	{
		o <- fastAssoc(y, g[,i])
		out[i, ] <- unlist(o)
	}
	out <- as.data.frame(out)
	names(out) <- names(o)
	return(out)
}

genotype_means <- function(g, x) {
    lapply(1:ncol(g), \(i) {
        tibble(g = g[,i], x = x) %>% group_by(g) %>% summarise(snp=i, mean = mean(x), sd = sd(x), n = n())
    }) %>% bind_rows()
}

qqplot <- function(pvector) {
    pvector <- pvector[!is.na(pvector) & !is.nan(pvector) & !is.null(pvector) & is.finite(pvector) & pvector<1 & pvector>0]
    o <- -log10(sort(pvector, decreasing=FALSE))
    e <- -log10(ppoints(length(pvector)))
    cs <- qchisq(1-pvector, 1)
    lambda <- median(cs, na.rm=TRUE) / qchisq(0.5, 1)
    plot(e, o, xlab = "Expected", ylab = "Observed")
    abline(0, 1)
    return(lambda)
}

whole_sim <- function(betas, nfam, h2, rho, ngen) {
    dat <- make_founders(betas = betas, nfam = nfam, h2 = h2, rho = rho)
    l <- list()
    l2 <- list()
    dat2 <- collapse_dat(dat)
    l[[1]] <- gwas(dat2$x$prs, dat2$g) %>% mutate(gen = 1)
    l2[[1]] <- l[[1]]
    for(i in 1:ngen) {
        message(i)
        dat <- create_child(dat, betas, h2 = h2, rho = rho)
        print(cor(dat$x$y_mother, dat$x$y_father))
        dat3 <- collapse_dat(dat)
        dat2 <- collapse_dat2(list(dat2, dat3))
        l[[i+1]] <- gwas(dat2$x$prs, dat2$g) %>% mutate(gen = i+1)
        l2[[i+1]] <- gwas(dat3$x$prs, dat3$g) %>% mutate(gen = i+1)
    }
    res <- bind_rows(l)
    res2 <- bind_rows(l2)
    m <- genotype_means(dat2$g, dat2$x$prs)
    return(list(res, res2, m))
}

make_grm <- function(g) {
    f <- colMeans(g) / 2
    g <- scale(g)
    m <- matrix(0, nrow(g), nrow(g))
    npol <- sum(f > 0)
    for(i in 1:nrow(g)) {
        for(j in 1:i) {
            x <- sum(g[i,] * g[j,], na.rm=TRUE) / npol
            # print(x)
            m[i,j] <- x
            m[j,i] <- m[i,j]
        }
    }
    # hist(m[lower.tri(m)], breaks=100)
    return(m)
}

greedy_remove_relateds <- function(m) {
    n <- nrow(m)
    keep <- rep(TRUE, n)

    b <- tibble(
        id=1:n,
        rel=sapply(1:n, \(i) sum(m[i,] > 0.2) - 1)
    ) %>% arrange(desc(rel)) %>%
    filter(rel > 0)

    diag(m) <- 0
    while(nrow(b) > 0) {
        message(nrow(b))
        keep[b$id[1]] <- FALSE
        m[b$id[1],] <- 0
        m[,b$id[1]] <- 0
        b <- tibble(
            id=1:n,
            rel=sapply(1:n, \(i) sum(m[i,] > 0.2) - 1)
        ) %>% arrange(desc(rel)) %>%
        filter(rel > 0)
        print(head(b))
    }
    return(keep)
}

# find_families <- function(cdat) {
#     gens <- unique(cdat$x$gen)
#     lapply(gens, function(i) {
#         m <- make_grm(cdat$g[cdat$x$gen == i,])
#         nfam <- nrow(m) / 2
#         diag(m) <- 0
#         m[m < 0.2] <- 0
#         m[m >= 0.2] <- 1
#         gr <- igraph::graph_from_adjacency_matrix(m)
#         rw <- igraph::cluster_walktrap(gr)
#         rw$membership
#         x <- cdat$x[cdat$x$gen == i,] %>% mutate(fid = rw$membership)
#         xs <- x %>% 
#             group_by(fid) %>% 
#             summarise(mean_prs = mean(prs), n = n(), se_prs = sd(prs) / sqrt(n())) %>% 
#             arrange(desc(mean_prs)) %>% 
#             filter(n > 1) %>% 
#             mutate(gen=i)

#         # Get family - y association
#         l <- list()
#         fams <- unique(x$fid)
#         for(j in fams) {
#             temp <- as.numeric(x$fid == j)
#             a <- summary(lm(temp ~ x$y))$coefficients
#             l[[j]] <- tibble(fid=j, y_beta = a[2,1], y_se=a[2,2], y_pval=a[2,4])
#         }
#         l <- bind_rows(l) %>% mutate(y_padj = p.adjust(y_pval, "bonferroni"))
#         xs <- left_join(xs, l, by="fid")

#         # Get rare variant - family association
#         l <- list()
#         for(j in fams) {
#             temp <- as.numeric(x$fid == j)
#             a <- gwas(temp, cdat$g[cdat$x$gen == i,1:nfam])
#             l[[j]] <- a %>% arrange(desc(fval)) %>% mutate(padj = p.adjust(pval, "bonferroni")) %>% slice_head(n=1) %>% mutate(fid = j, nsig=sum(padj < 0.05))
#         }
#         l <- bind_rows(l)
#         xs <- left_join(xs, l, by="fid")

#         # Get rare variant - y assocs
#         l <- list()
#         for(j in fams) {
#             temp <- as.numeric(x$fid == j)
#             a <- gwas(x$y, cdat$g[cdat$x$gen == i,1:nfam])
#             l[[j]] <- a %>% arrange(desc(fval)) %>% mutate(padj = p.adjust(pval, "bonferroni")) %>% slice_head(n=1) %>% mutate(fid = j, nsig2=sum(padj < 0.05))
#         }
#         l <- bind_rows(l)
#         xs <- left_join(xs, l, by="fid")

#         return(xs)
#     }) %>% bind_rows()
# }

fam_sim <- function(ngen, nfam, npoly, h2, rho) {
    betas <- rnorm(npoly)
    l <- list()
    l[[1]] <- make_founders(betas = betas, nfam = nfam, h2 = h2, rho = rho)
    for(i in 2:ngen) {
        message(i)
        l[[i]] <- create_child(l[[i-1]], betas, h2 = h2, rho = rho)
    }
    cdat <- collapse_dat2(lapply(l, collapse_dat))
    xs <- find_families(cdat) %>% mutate(npoly=npoly, nfam=nfam, h2=h2, rho=rho)
    return(xs)
}
