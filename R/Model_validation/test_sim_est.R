# Model to simulate an abundance-based community, adapted from 
# adapted from SimComm function in AHMBook package (Kéry & Royle, 2016)
source(here::here("R", "functions", "simCommNmix_function.R"))

real.bfly.comm <- simCommNmix(nspecies = 35, mean.beta1 = exp(bfly.mod$mean$mu.beta.can),
                            mean.beta2 = exp(bfly.mod$mean$mu.beta.und), sd.beta1 = bfly.mod$mean$sd.alpha.can,
                            sd.beta2 = bfly.mod$mean$sd.beta.und, mu.beta3 = bfly.mod$mean$mu.beta1,
                            sd.beta3 = bfly.mod$mean$sd.beta1, mean.alpha1 = plogis(bfly.mod$mean$mu.alpha.can),
                            sd.alpha1 = bfly.mod$mean$sd.alpha.can, mean.alpha2 = plogis(bfly.mod$mean$mu.alpha.und),
                            sd.alpha2 = bfly.mod$mean$sd.alpha.und, mu.alpha3 = bfly.mod$mean$mu.alpha1,
                            sd.alpha3 = bfly.mod$mean$sd.alpha1, mu.alpha4 = bfly.mod$mean$mu.alpha2,
                            sd.alpha4 = bfly.mod$mean$sd.alpha2, sd.re1 = bfly.mod$mean$sd.month, 
                            sd.re2 = bfly.mod$mean$sd.area)

# Real parameters ---------------------------------------------------------

params <- c("beta1", "beta2", "beta3", "alpha1", "alpha2",
            "alpha3", "alpha4")

params.e <- c("beta.can", "beta.und", "beta1", "alpha.can", 
              "alpha.und", "alpha1", "alpha2")

df.bfly <- data.frame(values = c(unlist(real.bfly.comm$Params.estimated[params]),
                                 unlist(bfly.mod$mean[params.e])),
                      params = rep(params.e, each = 35, 2),
                      names = c(rep("Model", length(params)*35), rep("Sim", length(params)*35)))

library(ggplot2)
compars <- ggplot(data = df.bfly, aes(x = values, colour = names, fill = names)) +
  geom_density(alpha = 0.5) + facet_grid(~ params)

compars

bfly.fitjags <- readRDS(here::here("output", "fitjags.mod.bfly.rds"))
bfly.N <- as.matrix(bfly.fitjags)


nsamp <- dim(bfly.N)[1]  # 3000 MCMC samples

# extract the abundance for each species in each site
N <- array(NA, dim = c(300, 35, nsamp))
for(j in 1:nsamp){  
  # Fill z matrix by column (default)
  # j = 1
  cat(paste("\nMCMC sample", j, "\n"))
  N[,,j] <- bfly.N[j, 1:(300*35)]
}


N.r <- apply(real.bfly.comm$Params.estimated$N, 1, sum)
N.e <- apply(round(apply(N, c(1,2), mean), 2), 1, sum)

plot(N.r ~ N.e)

