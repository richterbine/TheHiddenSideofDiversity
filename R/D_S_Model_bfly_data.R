# run the model with the community real data
# FLONA bfly - 2016 to 2017
library(reshape)
library(jagsUI)

# read the bfly data
data.bfly <- read.table(here::here("output", "joint_bfly_data.txt"), header = T)

# Organizing the data by sampling day -------------------------------------
data.bfly$Date <- as.Date(as.vector(data.bfly$Date))
data.bfly$julian.date <- as.numeric(format(data.bfly$Date, "%j"))          # data in the julian date format

for (i in 1:ncol(data.bfly)) {
  if(class(data.bfly[,i]) == "character"){
    data.bfly[,i]<- as.factor(data.bfly[,i])
  }
}
str(data.bfly)

(species.list <- levels(data.bfly$Species))                               # alphabetic list
(spec.name.list <- tapply(data.bfly$Abundance, data.bfly$Species, sum))    # species ID 
(spec.id.list <- unique(data.bfly$Species))                                # ID list for unique species
(ordered.spec.name.list <- spec.name.list[order(spec.name.list)])          # ID-order list

# covariates for biological process --------------------------

cov.abn <- data.frame(Temperature = round(tapply(data.bfly$Temperature, data.bfly$Sites, mean, na.rm = T), 3),
                      Strata = c(rep(0, length(grep("C", unique(data.bfly$Sites)))), 
                                 rep(1, length(grep("U", unique(data.bfly$Sites)))))) 

cov.abn$Month <- as.numeric(as.factor(data.bfly[match(rownames(cov.abn), data.bfly$Sites), "Month"]))
cov.abn$SU <- as.numeric(as.factor(data.bfly[match(rownames(cov.abn), data.bfly$Sites), "SU"]))

# scaling the temperature covariate
cov.abn$Temperature <- scale(cov.abn$Temperature)


# Covariates for the observational process --------------------------------

cov.date <- melt(data.bfly, id.vars = c("Sites","Samp.day"), measure.vars = "julian.date")
cov.date <- cast(cov.date, Sites ~ Samp.day, mean)
cov.date[is.na(cov.date)] <- NA

# scaling the julian date
scaled <- matrix(NA, nrow = nrow(cov.date), ncol = ncol(cov.date[,-1]))
for (i in 1:ncol(scaled)) {
  scaled[,i] <- as.numeric(cov.date[,i+1])
}

colnames(scaled) <- c("rep1", "rep2", "rep3", "rep4", "rep5")
scaled.dates <- scale(x = scaled) # scaling the dates covariate
scaled.dates[is.na(scaled.dates)] <- 0
head(scaled.dates)

# temperature covariate
cov.temp <- melt(data.bfly, id.vars = c("Sites","Samp.day"), measure.vars = "Temperature")
cov.temp <- cast(cov.temp, Sites ~ Samp.day, mean, na.rm = T)
cov.temp[is.na(cov.temp)] <- NA

scaled <- matrix(NA, nrow = nrow(cov.temp), ncol = ncol(cov.temp[,-1]))
for (i in 1:ncol(scaled)) {
  scaled[,i] <- as.numeric(cov.temp[,i+1])
}

colnames(scaled) <- c("rep1", "rep2", "rep3", "rep4", "rep5")
scaled.temp <- scale(x = scaled) # scaling the temperature covariate
scaled.temp[is.na(scaled.temp)] <- 0
head(scaled.temp)


# put the observed data in a three dimensional array ----------------------

# the first dimension i is the sites; the second dimension j, is the daily
# repetitions; and the last dimension k, is the species. 

junk.melt <- melt(data.bfly, id.var = c("Species", "Sites", "Samp.day"), 
                  measure.var = "Abundance")
Ycount <- cast(junk.melt, Sites ~ Samp.day ~ Species)
yc <- Ycount[,,-36] # removing the NA that is not a species
dimnames(yc) <- list(rownames(cov.abn), NULL, names(spec.name.list))
dim(yc)

# adding NAs in the reps when sampling did not occur
posit<- which(is.na(cov.date) == T, arr.ind = T)
for (j in 1:dim(yc)[3]) {
  for (i in 1:nrow(posit)) {
    yc[,,j][posit[i,"row"], posit[i,"col"]-1] = NA
  }
}


# Bundle and summarize data set -------------------------------------------

str(win.data <- list(yc = yc, nsite = dim(yc)[1], nrep = dim(yc)[2], 
                     nspec = dim(yc)[3], Strata = as.vector(cov.abn$Strata),
                     Temp = as.vector(cov.abn$Temperature),
                     Date = scaled.dates, Temp_det = scaled.temp,
                     Area = as.vector(cov.abn$SU), Month = as.vector(cov.abn$Month)))

## Initial values 
Nst1 <- apply(yc, c(1,3), max, na.rm = T) + 3 #+ some.more
inits <- function()list(N = Nst1)

# Parameters monitored 
params <- c("beta.can", "mu.beta.can", "sd.beta.can",
           "beta.und", "mu.beta.und", "sd.beta.und",
           "beta1", "mu.beta1", "sd.beta1",
           "alpha.can", "mu.alpha.can", "sd.alpha.can",
           "alpha.und", "mu.alpha.und", "sd.alpha.und",
           "alpha1", "mu.alpha1", "sd.alpha1",
           "alpha2", "mu.alpha2", "sd.alpha2",
           "month", "sd.month", "area", "sd.area",
           "lambda", "mlambda.can", "mlambda.und", 
           "mp.can", "mp.und")

bfly.mod <- jags(data = win.data, inits = inits, parameters.to.save = params, 
                n.chains = 3, n.iter = 150000, n.burnin = 50000, n.thin = 100, 
                model.file = here::here("R", "Bayesian_models", "Static_Nmixture_P_RE.txt"))

saveRDS(bfly.mod, here::here("output", "out.mod.bfly.rds"))

# run basic jags model to take the chains
# define a new set of parameters to monitor
params.basic <- c("N", "beta.can", "beta.und", "beta1",
                  "alpha.can", "alpha.und", "alpha1", "alpha2")

bfly.mod.basic <- jags.basic(data = win.data, inits = inits, parameters.to.save = params.basic,
                             n.chains = 3,  n.iter = 150000, n.burnin = 50000, n.thin = 100,
                             parallel = T, model.file = here::here("R", "Bayesian_models",
                                                                   "Static_Nmixture_P_RE.txt"))

saveRDS(bfly.mod.basic, here::here("output", "fitjags.mod.bfly.rds"))


# ################################################################ --------
library(plyr)
library(ggplot2)

# read the outputs for N-mixture model
bfly.mod <- readRDS(here::here("output", "out.mod.bfly.rds"))
bfly.mod.basic <- readRDS(here::here("output",'fitjags.mod.bfly.rds'))

# traceplots for each parameters
#par(mfrow = c(1, 2))
traceplot(bfly.mod.basic[, 10501:10535]) # trace for alpha.can
traceplot(bfly.mod.basic[, 10536:10570]) # trace for alpha.und
traceplot(bfly.mod.basic[, 10571:10605]) # trace for alpha1
traceplot(bfly.mod.basic[, 10606:10640]) # trace for alpha2
traceplot(bfly.mod.basic[, 10641:10675]) # trace for beta.can
traceplot(bfly.mod.basic[, 10676:10710]) # trace for beta.und
traceplot(bfly.mod.basic[, 10711:10745]) # trace for beta1

params.est <- c("mu.beta.can", "sd.beta.can", "mu.beta.und", "sd.beta.und",
                "mu.beta1", "sd.beta1", "mu.alpha.can", "sd.alpha.can",
                "mu.alpha.und", "sd.alpha.und", "mu.alpha1", "sd.alpha1", 
                "mu.alpha2", "sd.alpha2", "sd.month", "sd.area")

df.est.params <- data.frame(row.names = params.est)

for (i in 1:length(params.est)) {
  df.est.params[i,1] <- bfly.mod$mean[[params.est[i]]]
  df.est.params[i,2] <- bfly.mod$q2.5[[params.est[i]]]
  df.est.params[i,3] <- bfly.mod$q97.5[[params.est[i]]]
}

colnames(df.est.params) <- c("Mean", "LowCRI", "UpperCRI")
df.est.params

### Community-mean abundance and detection
# for canopy
n.sample <- exp(rnorm(10^5, mean = bfly.mod$mean$mu.beta.can, sd = bfly.mod$mean$sd.beta.can))
p.sample <- plogis(rnorm(10^5, mean = bfly.mod$mean$mu.alpha.can, sd = bfly.mod$mean$sd.alpha.can))
summary(n.sample)   ;   summary(p.sample)

# for understory
nu.sample <- exp(rnorm(10^5, mean = bfly.mod$mean$mu.beta.und, sd = bfly.mod$mean$sd.beta.und))
pu.sample <- plogis(rnorm(10^5, mean = bfly.mod$mean$mu.alpha.und, sd = bfly.mod$mean$sd.alpha.und))
summary(nu.sample)   ;   summary(pu.sample)
sum(round(sort(nu.sample, decreasing = T), 2) < 100)-10^5

df <- data.frame(Strata = c(rep("Canopy", each = 10^5), rep("Understory", each = 10^5)),
                 Abundance = as.vector(c(n.sample, nu.sample)),
                 Detection = as.vector(c(p.sample, pu.sample)))

# Plotting the community-mean abundance
mean <- c(mean(n.sample), mean(nu.sample), mean(p.sample), mean(pu.sample))

n <- ggplot(df, aes(x = Abundance, color = Strata, fill = Strata)) + 
  geom_density(alpha = 0.5, position = "identity") +
  geom_vline(data = ddply(df, "Strata", summarise, grp.mean = mean(Abundance)), aes(xintercept = grp.mean, color = Strata),
             linetype = "dashed") + scale_x_sqrt(limits = c(0, 100)) +
  labs(tag = "a)") + ylab("Density") +
  scale_color_viridis_d(option = "A") + scale_fill_viridis_d(option = "A") +
  theme(legend.position = "none") + xlab(" Expected abundance") +
  annotate("text", x = 75, y = 4.75, parse = T, 
           label = sprintf("mu[can] == %0.3f", mean[1])) +
  annotate("text", x = 75, y = 4.4, parse = T, 
           label = sprintf("mu[und] == %0.3f", mean[2]))
n

# plotting the community-mean detection probability
p <- ggplot(df, aes(x = Detection, color = Strata, fill = Strata)) + 
  geom_density(alpha = 0.5, position = "identity") +
  labs(tag = "b)") + ylab(" ") + xlab("Detection probability") +
  geom_vline(data = ddply(df, "Strata", summarise, grp.mean = mean(Detection)), aes(xintercept = grp.mean, color = Strata),
             linetype = "dashed") + scale_x_continuous(trans = "sqrt") +
  scale_color_viridis_d(option = "A") + scale_fill_viridis_d(option = "A") +
  theme(legend.position = "none") +
  annotate("text", x = .75, y = 14.25, label = sprintf("mu[can] == %0.3f", mean[3]), 
           parse = T) +
  annotate("text", x = .75, y = 13.25, label = sprintf("mu[und] == %0.3f", mean[4]), 
           parse = T)
p

cowplot::plot_grid(n, p)
cowplot::save_plot(here::here("output", "figures", "Fig2_meanC.png"), 
                   cowplot::plot_grid(n, p), base_width = 8)

### Inter-specific variability - species heterogeneity
params <- params.est[c(2,4,6,8,10,12,14:16)]
tmp.sd <- matrix(NA, ncol = length(params), nrow = length(bfly.mod$sims.list$sd.beta.can))
colnames(tmp.sd) <- params

for (i in 1:length(params)) {
  tmp.sd[,i] <- as.vector(bfly.mod$sims.list[[params[i]]])   
}

df.sd <- data.frame(values = as.vector(tmp.sd),
                    params = rep(colnames(tmp.sd), each = nrow(tmp.sd)))
t <- gsub("[0-9]", "", unlist(lapply(strsplit(df.sd$params, "\\."), function(x) x[2])))

df.sd$proc.type <- ifelse(t == "beta", "Biological",
                          ifelse(t == "alpha", "Observational", 
                                 "Biological"))


sd.plot <- ggplot(data = df.sd, aes(x = values, colour = params, fill = params)) + 
  geom_density(alpha = 0.5) + 
  scale_color_viridis_d(option = "A", name = "Parameters", labels = c(expression(alpha["can"]),
                                                                      expression(alpha["und"]),
                                                                      expression(alpha[1]),
                                                                      expression(alpha[2]),
                                                                      expression(beta["can"]),
                                                                      expression(beta["und"]),
                                                                      expression(beta[1]),
                                                                      "SM", "SU")) +
  scale_fill_viridis_d(option = "A", name = "Parameters", labels = c(expression(alpha["can"]),
                                                                     expression(alpha["und"]),
                                                                     expression(alpha[1]),
                                                                     expression(alpha[2]),
                                                                     expression(beta["can"]),
                                                                     expression(beta["und"]),
                                                                     expression(beta[1]),
                                                                     "SM", "SU"))

mean.sd <- data.frame(mean = unlist(bfly.mod$mean[params]), 
                      params = names(unlist(bfly.mod$mean[params])),
                      proc.type = c(rep("Biological", 3), 
                                    rep("Observational", 4),
                                    rep("Biological", 2)))

sd.plot <- sd.plot + geom_vline(data = mean.sd, 
                                aes(xintercept = mean, color = params),
                                linetype = "dashed") +
  facet_wrap(~proc.type) + labs(x = "Estimated SD values", y = "Density")
sd.plot

cowplot::save_plot(here::here("output", "figures", "FigB1_Csd.png"),
                   sd.plot, base_width = 8)


# Mean species-specific parameters ----------------------------------------

params.spp <- c("beta.can", "beta.und", "beta1", 
                "alpha.can", "alpha.und", "alpha1", "alpha2")

tmp.spp <- matrix(NA, ncol = length(params.spp), nrow = length(bfly.mod$mean$beta.can))
colnames(tmp.spp) <- params.spp
list.spp <- list(Mean = tmp.spp,
                 q2.5 = tmp.spp,
                 q97.5 = tmp.spp)

for (i in 1:length(params.spp)) {
  list.spp[[1]][,i] <- bfly.mod$mean[[params.spp[i]]]
  list.spp[[2]][,i] <- bfly.mod$q2.5[[params.spp[i]]]
  list.spp[[3]][,i] <- bfly.mod$q97.5[[params.spp[i]]]
}

df.spp <- data.frame(Mean = as.vector(list.spp[[1]]),
                     q2.5 = as.vector(list.spp[[2]]),
                     q97.5 = as.vector(list.spp[[3]]),
                     params = rep(params.spp, each = nrow(tmp.spp)), 
                     Spp = rep(na.omit(sort(unique(data.bfly$Species)))))
df.spp$params <- factor(df.spp$params, labels = c(expression(alpha[can]),
                                                  expression(alpha[und]),
                                                  expression(alpha[1]),
                                                  expression(alpha[2]),
                                                  expression(beta[can]),
                                                  expression(beta[und]),
                                                  expression(beta[1])))

spp.plot <- ggplot(df.spp, aes(x = Mean, y = Spp, colour = Spp)) +
  geom_point() + geom_errorbar(aes(xmin = q2.5, xmax = q97.5), width = 0.1) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~params, labeller = label_parsed, ncol = 4, scales = "free_x") +
  theme(legend.position = "none") + labs(y = "Species") +
  scale_color_viridis_d(option = "A")
spp.plot

cowplot::save_plot(here::here("output", "figures", "FigB2_Smean.png"),
                   spp.plot, base_width = 10, base_height = 9)


# Effects of the covariates -----------------------------------------------
# NO EFFECTS!
