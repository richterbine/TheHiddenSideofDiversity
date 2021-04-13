# Model to simulate an abundance-based community, adapted from 
# adapted from SimComm function in AHMBook package (Kéry & Royle, 2016)

simCommNmix <- function (nsites = 300, nreps = 5, nspecies = 40,
                         ng = 2, mean.beta1 = 1, sd.beta1 = .2,
                         mean.beta2 = 1, sd.beta2 = .2, mu.beta3 = .02,
                         sd.beta3 = .3, mean.alpha1 = .5, sd.alpha1 = .3,
                         mean.alpha2 = .5, sd.alpha2 = .3, mu.alpha3 = 0,
                         sd.alpha3 = .8, mu.alpha4 = 1.2, sd.alpha4 = .3,
                         ng.re1 = 5, ng.re2 = 6, sd.re1 = .5, sd.re2 = .5) 
{
  nsites <- round(nsites[1])
  nreps <- round(nreps[1])
  nspecies <- round(nspecies[1])
  ng.re1 <- round(ng.re1[1])
  ng.re2 <- round(ng.re2[1])
  
  # create a error message
  stopifNegative(mean.beta1)
  stopifNegative(mean.beta2)
  stopifNegative(sd.beta1)
  stopifNegative(sd.beta2)
  stopifNegative(sd.beta3)
  stopifnotProbability(mean.alpha1)
  stopifnotProbability(mean.alpha2)
  stopifNegative(sd.alpha1)
  stopifNegative(sd.alpha2)
  stopifNegative(sd.alpha3)
  stopifNegative(sd.alpha4)
  stopifNegative(sd.re1)
  stopifNegative(sd.re2)
  
  # generate the data set
  y.all <- y.obs <- p <- array(NA, c(nsites, nreps, nspecies))
  dimnames(y.all) <- dimnames(y.obs) <- dimnames(p) <- list(paste("site", 1:nsites, sep = ""), 
                                                            paste("rep", 1:nreps, sep = ""), 
                                                            paste("sp", 1:nspecies, sep = ""))
  N <- lambda <- matrix(NA, nsites, nspecies)
  dimnames(N) <- dimnames(lambda) <- list(paste("site", 1:nsites, sep = ""), 
                                          paste("sp", 1:nspecies, sep = ""))
  detected.at.all <- rep(NA, nspecies)
  
  # creating the covariates and random effects data
  g <- c(rep(0, nsites/2), rep(1, nsites/2))
  cov.abn <- sort(rnorm(nsites), decreasing = T)
  cov1.det <- matrix(rnorm(nsites * nreps), ncol = nreps)
  cov2.det <- matrix(rnorm(nsites * nreps), ncol = nreps)
  g.re1 <- gl(n = ng.re1, k = ng * ng.re2 * 5)
  g.re2 <- gl(n = ng.re2, k = 1, length = nsites)
  
  # use the hyperparameters informed to generate the parameters 
  # that take into account the heterogeneity of the species
  # put the hyperparameters in the link-scale
  mu.beta1 <- log(mean.beta1)
  mu.beta2 <- log(mean.beta2)
  mu.alpha1 <- ifelse(mean.alpha1 == "1", 500, qlogis(mean.alpha1))
  mu.alpha2 <- ifelse(mean.alpha2 == "1", 500, qlogis(mean.alpha2))
  
  # for then generate the parameters by species
  beta1 <- rnorm(nspecies, mu.beta1, sd.beta1)
  beta2 <- rnorm(nspecies, mu.beta2, sd.beta2)
  beta3 <- rnorm(nspecies, mu.beta3, sd.beta3)
  
  alpha1 <- rnorm(nspecies, mu.alpha1, sd.alpha1)
  alpha2 <- rnorm(nspecies, mu.alpha2, sd.alpha2)
  alpha3 <- rnorm(nspecies, mu.alpha3, sd.alpha3)
  alpha4 <- rnorm(nspecies, mu.alpha4, sd.alpha4)
  
  # random effects
  re1 <- matrix(rnorm(ng.re1 * nspecies, 0, sd.re1), ncol = nspecies)
  re2 <- matrix(rnorm(ng.re2 * nspecies, 0, sd.re2), ncol = nspecies)
  
  for (k in 1:nspecies) {
    for (i in 1:nsites) {
      lambda[i, k] <- exp(beta1[k] * (1 - g[i]) + beta2[k] * g[i] + 
                            beta3[k] * cov.abn[i] + re1[g.re1[i], k] +
                            re2[g.re2[i], k])
      
      N[i, k] <- rpois(1,lambda[i, k])
      
      for (j in 1:nreps) {
        p[i, j, k] <- plogis(alpha1[k] * (1 - g[i]) + alpha2[k] * g[i] +
                               alpha3[k] * cov1.det[i,j] + 
                               alpha4[k] * cov2.det[i,j])
      }
    }
  }
  
  tmp <- apply(N, 2, sum)
  occurring.in.sample <- as.numeric(tmp > 0)
  
  for (k in 1:nspecies) {
    for (i in 1:nsites) {
      for (j in 1:nreps) {
        y.all[i, j, k] <- rbinom(1, N[i, k], p[i, j, k])
      }
    }
    detected.at.all[k] <- if (any(y.all[, , k] > 0)) 
      TRUE
    else FALSE
  }
  y.obs <- y.all[, , detected.at.all]
  detected.at.site <- apply(y.obs > 0, c(1, 3), any)
  ymax.obs <- apply(y.all, c(1, 3), max)
  Ntotal.fs <- sum(occurring.in.sample)
  Ntotal.obs <- sum(detected.at.all)
  
  return(list(Parametres = list(nsites = nsites, nreps = nreps, nspecies = nspecies, ng = ng,
                                mean.beta1 = mean.beta1, sd.beta1 = sd.beta1, mean.beta2 = mean.beta2, 
                                sd.beta2 = sd.beta2, mu.beta3 = mu.beta3, sd.beta3 = sd.beta3,
                                mean.alpha1 = mean.alpha1, sd.alpha1 = sd.alpha1,mean.alpha2 = mean.alpha2,
                                sd.alpha2 = sd.alpha2, mu.alpha3 = mu.alpha3, sd.alpha3 = sd.alpha3,
                                mu.alpha4 = mu.alpha4, sd.alpha4 = sd.alpha4, ng.re1 = ng.re1, ng.re2 = ng.re2,
                                sd.re1 = sd.re1, sd.re2 = sd.re2),
              Params.estimated = list(g = g, cov.abn = cov.abn, cov1.det = cov1.det, 
                                      cov2.det = cov2.det, g.re1 = g.re1, g.re2 = g.re2, mu.beta1 = mu.beta1,
                                      mu.beta2 = mu.beta2, mu.alpha1 = mu.alpha1, mu.alpha2 = mu.alpha2,
                                      beta1 = beta1, beta2 = beta2, beta3 = beta3,
                                      alpha1 = alpha1, alpha2 = alpha2, alpha3 = alpha3, 
                                      alpha4 = alpha4, lambda = lambda, p = p, N = N, 
                                      y.all = y.all, y.obs = y.obs, ymax.obs = ymax.obs, 
                                      Ntotal.fs = Ntotal.fs, Ntotal.obs = Ntotal.obs)))
  
}

stopifNegative <- function(arg, allowNA=FALSE, allowZero=TRUE) {
  name <- deparse(substitute(arg))
  if(allowNA && all(is.na(arg))) {  # An all-NA vector is logical, but ok.
    # do nothing
  } else {
    if(!allowNA && any(is.na(arg)))
      stop("Argument '", name, "' must not contain NA or NaN.", call.=FALSE)
    if(!is.numeric(arg))
      stop("Argument '", name, "' must be numeric.", call.=FALSE)
    if(allowZero) {
      if(any(arg < 0, na.rm=TRUE)) {
        if(allowNA) {
          stop("Argument '", name, "' must be non-negative, or NA.", call.=FALSE)
        } else {
          stop("Argument '", name, "' must be non-negative.", call.=FALSE)
        }
      }
    } else {
      if(any(arg <= 0, na.rm=TRUE)) {
        if(allowNA) {
          stop("Argument '", name, "' must be greater than 0, or NA.", call.=FALSE)
        } else {
          stop("Argument '", name, "' must be greater than 0.", call.=FALSE)
        }
      }
    }
  }
}

stopifnotProbability <- function(arg, allowNA=FALSE) {
  name <- deparse(substitute(arg))
  if(allowNA && all(is.na(arg))) {  # An all-NA vector is logical, but ok.
    # do nothing
  } else {
    if(!allowNA && any(is.na(arg)))
      stop("Argument '", name, "' must not contain NA or NaN.", call.=FALSE)
    if(!is.numeric(arg))
      stop("Argument '", name, "' must be numeric.", call.=FALSE)
    if(any(arg < 0 | arg > 1, na.rm=TRUE)) {
      if(allowNA) {
        stop("Argument '", name, "' must be a probability between 0 and 1, or NA.", call.=FALSE)
      } else {
        stop("Argument '", name, "' must be a probability between 0 and 1.", call.=FALSE)
      }
    }
  }
}
