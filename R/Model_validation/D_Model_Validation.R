# Script to validate the bayesian model, adopting different values for betas 
# and alphas

source(here::here("R", "functions", "simCommNmix_function.R"))

# Metacommunity simulation ------------------------------------------------

# creating a data frame with different values for target parameters
df.parameters <- data.frame(mean.beta1 = c(.5, 1.5, 2),
                            mean.beta2 = c(.5, 1.5, 2),
                            mean.alpha1 = c(.2, .6, .8),
                            mean.alpha2 = c(.2, .6, .8))

# simulate the communities with the target parameters (beta and alpha means)
sim.comm.beta1 <- sim.comm.beta2 <- sim.comm.alpha1 <- sim.comm.alpha2 <- list()

for (i in 1:nrow(df.parameters)) {
  sim.comm.beta1[[i]] <- simCommNmix(mean.beta1 = df.parameters[i, 1])
  sim.comm.beta2[[i]] <- simCommNmix(mean.beta2 = df.parameters[i, 2])
  sim.comm.alpha1[[i]] <- simCommNmix(mean.alpha1 = df.parameters[i, 3])
  sim.comm.alpha2[[i]] <- simCommNmix(mean.alpha2 = df.parameters[i, 4])
}
sim.comm.list <- list(sim.comm.beta1, sim.comm.beta2, sim.comm.alpha1, sim.comm.alpha2)
saveRDS(sim.comm.list, here::here("R", "Model_validation", "Comm_simulated.rds"))


# Model validation: testing the proposed N-mixture model ------------------

# parameters to monitor
params<- c("beta.can", "beta.und", "beta1","alpha.can", "alpha.und", "alpha1", 
           "alpha2", "mu.beta.can", "sd.beta.can", "mu.beta.und",
           "sd.beta.und", "mu.beta1", "sd.beta1", "mu.alpha.can", "sd.alpha.can", 
           "mu.alpha.und", "sd.alpha.und", "mu.alpha1", "sd.alpha1", "mu.alpha2", 
           "sd.alpha2", "sd.month", "sd.area", "mlambda.can", 
           "mlambda.und", "mp.can", "mp.und")


# create a list to store the results for the bayesian model to each parameter
sim.mod <-list(mean.beta1 = c(list(NA), list(NA), list(NA)),
               mean.beta2 = c(list(NA), list(NA), list(NA)),
               mean.alpha1 =c(list(NA), list(NA), list(NA)),
               mean.alpha2 = c(list(NA), list(NA), list(NA)))

# Bundle and summarize data set 
for (i in 1:length(sim.comm.list)) {
  for (j in 1:length(sim.comm.list[[i]])) {
    #i=1
    #j=1
    win.data <- list(yc = sim.comm.list[[i]][[j]]$Params.estimated$y.obs, 
                     nsite = dim(sim.comm.list[[i]][[j]]$Params.estimated$y.obs)[1], 
                     nrep = dim(sim.comm.list[[i]][[j]]$Params.estimated$y.obs)[2], 
                     nspec = dim(sim.comm.list[[i]][[j]]$Params.estimated$y.obs)[3], 
                     Strata = sim.comm.list[[i]][[j]]$Params.estimated$g,
                     Temp = sim.comm.list[[i]][[j]]$Params.estimated$cov.abn,
                     Date = sim.comm.list[[i]][[j]]$Params.estimated$cov1.det, 
                     Temp_det = sim.comm.list[[i]][[j]]$Params.estimated$cov2.det,
                     Area = as.numeric(sim.comm.list[[i]][[j]]$Params.estimated$g.re2), 
                     Month = as.numeric(sim.comm.list[[i]][[j]]$Params.estimated$g.re1))
    
    ## Initial values 
    Nst1 <- sim.comm.list[[i]][[j]]$Params.estimated$ymax.obs
    inits <- function()list(N = Nst1)
    
    # running the model
    sim.mod[[i]][[j]] <- jagsUI::jags(data = win.data, inits = inits, parameters.to.save = params, 
                                      n.chains = 3,
                                      n.iter = 10500, n.burnin = 500, n.thin = 10, 
                                      model.file = here::here("R", "Bayesian_models", "Static_Nmixture_P_RE.txt"))
  }
}

saveRDS(object = sim.mod, file = here::here("R", "Model_validation", "model_comm_simulated.rds"))
