# performing the hidden diversity function to estimate the extent to which 
# imperfect detection affects biodiversity descriptors
library(jagsUI)
library(ggplot2)
library(reshape)

# read the raw data -------------------------------------

data.bfly <- read.table(here::here("output", "joint_bfly_data.txt"), header = TRUE) # data of butterflies for bayesian model

(spec.name.list <- tapply(data.bfly$Abundance, data.bfly$Species, sum))         # species ID 

# building the community matrix
junk.melt <- melt(data.bfly, id.var = c("Sites", "Species"), 
                  measure.var = "Abundance")
comm <- cast(junk.melt, Sites ~ Species)
comm.bfly.obs <- as.data.frame(comm[,c(-1, -37)]) # removing the NA that is not a species
dimnames(comm.bfly.obs) <- list(comm$Sites, names(spec.name.list))

# read the N matrices estimated by the bayesian model
# Plug MCMC samples for full N matrix into 3D array
bfly.fitjags <- readRDS(here::here("output", "fitjags.mod.bfly.rds"))
bfly.N <- as.matrix(bfly.fitjags)

nsite <- nrow(comm.bfly.obs)
nspec <- ncol(comm.bfly.obs)
nsamp <- dim(bfly.N)[1]  # 3000 MCMC samples

# extract the abundance for each species in each site
N <- array(NA, dim = c(nsite, nspec, nsamp))
for(j in 1:nsamp){  
  # Fill z matrix by column (default)
  # j = 1
  cat(paste("\nMCMC sample", j, "\n"))
  N[,,j] <- bfly.N[j, 1:(nsite*nspec)]
}

comm.bfly.est <- N[, , sample(1:dim(N)[3], 10)]
dimnames(comm.bfly.est) <- list(rownames(comm.bfly.obs), colnames(comm.bfly.obs), NULL)
str(comm.bfly.est)

# read the phylogenetic pruned tree for the species sampled at FLONA
tree.bfly <- ape::read.tree(here::here("output", "tree_bfly_FLONA.txt"))
str(tree.bfly)

# read the functional matrix by species
traits.bfly <- read.table(here::here("data", "processed", "Mean_traits_bfly.txt"), header = T)
str(traits.bfly)


# Run the hidden diversity parallel function  -------------------------------------

source(here::here("R", "functions", "hidden_diversity_parallel.R"))
# ape::is.rooted(tree.bfly)

hd.bfly <- hidden.diversity(comm = comm.bfly.obs, N = comm.bfly.est, phy = tree.bfly,
                            trait = traits.bfly, metrics = c("pd", "mpd"), binary = T, 
                            abundance.weighted = T, null.model = "taxa.labels", runs = 499,
                            parallel = 4)

head(hd.bfly$TD)
head(hd.bfly$Abund)
head(hd.bfly$sesPD)
head(hd.bfly$sesMPD)
head(hd.bfly$sesFD)
head(hd.bfly$sesMFD)

saveRDS(hd.bfly, here::here("output/HD_bfly_parallel.rds"))

##################################################################################
############  EXTRACTING THE VALUES AND PLOTTING THE RESULTS #####################
##################################################################################

hd.bfly <- readRDS(here::here("output/HD_bfly_parallel.rds"))
names(hd.bfly)
# creating a data frame to storage the results

HD.df <- data.frame(HD.values = c(hd.bfly$TD$HD.TD, hd.bfly$Abund$HD.N, hd.bfly$sesPD$HD.PD,
                                  hd.bfly$sesMPD$HD.MPD, hd.bfly$sesFD$HD.FD, hd.bfly$sesMFD$HD.MFD),
                    Metric = c(rep(colnames(hd.bfly$TD[4]), nsite), rep(colnames(hd.bfly$Abund[4]), nsite),
                               rep(colnames(hd.bfly$sesPD[4]), nsite), rep(colnames(hd.bfly$sesMPD[4]), nsite),
                               rep(colnames(hd.bfly$sesFD[4]), nsite), rep(colnames(hd.bfly$sesMFD[4]), nsite)),
                    Strata = c(rep(substr(rownames(hd.bfly$TD), 1, 1), length(hd.bfly))),
                    Month = c(rep(substr(rownames(hd.bfly$TD), 6, 6), length(hd.bfly))),
                    SU = c(rep(substr(rownames(hd.bfly$TD), 3, 5), length(hd.bfly))))


library(ggplot2)
library(viridis)

p.hd <- ggplot(HD.df, aes(x = Metric, y = HD.values, colour = Strata, fill = Strata)) +
  geom_boxplot(alpha = .5) + facet_wrap(~ Metric, nrow = 1, scales = "free_x") +
  scale_color_viridis_d(option = "A") + scale_fill_viridis_d(option = "A") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "firebrick1") +
  theme(legend.position = "none", 
        axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(y = "Hidden Diversity", tag = "a)", x = "Diversity measures")
p.hd


# testing if the HD differs between strata
library(lme4)
library(lmerTest)

mod.td <- lmer(HD.values ~ Strata + (1|Month)  + (1|SU), 
               data = subset(HD.df, Metric == "HD.TD"))
summary(mod.td)

mod.pd <- lmer(HD.values ~ Strata + (1|Month)  + (1|SU), 
               data = subset(HD.df, Metric == "HD.PD"))
summary(mod.pd)

mod.fd <- lmer(HD.values ~ Strata + (1|Month)  + (1|SU), 
               data = subset(HD.df, Metric == "HD.FD"))
summary(mod.fd)

mod.mpd <- lmer(HD.values ~ Strata + (1|Month)  + (1|SU), 
                data = subset(HD.df, Metric == "HD.MPD"))
summary(mod.mpd)

mod.mfd <- lmer(HD.values ~ Strata + (1|Month)  + (1|SU), 
                data = subset(HD.df, Metric == "HD.MFD"))
summary(mod.mfd)

out.lmer <- rbind(summary(mod.td)$coefficients, summary(mod.pd)$coefficients,
                  summary(mod.fd)$coefficients, summary(mod.mpd)$coefficients,
                  summary(mod.mfd)$coefficients)

# visualization of underestimation or overestimation of data

HD.df$Obs.values <- c(hd.bfly$TD[,1], hd.bfly$Abund[,1], hd.bfly$sesPD[,1], 
                      hd.bfly$sesMPD[,1], hd.bfly$sesFD[,1], hd.bfly$sesMFD[,1])
HD.df$Est.values <- c(hd.bfly$TD[,2], hd.bfly$Abund[,2], hd.bfly$sesPD[,2], 
                      hd.bfly$sesMPD[,2], hd.bfly$sesFD[,2], hd.bfly$sesMFD[,2])

str(HD.df)

p.bias_noise <- ggplot(HD.df, aes(x = Obs.values, y = Est.values, colour = Strata), alpha, .5, size = 4) +
  geom_point(alpha = .5) + scale_color_viridis_d(option = "A") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "lightblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "lightblue") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "firebrick1") +
  facet_wrap(~ Metric, nrow = 1, scale = "free") + theme(legend.position = "none") +
  labs(y = "Estimated Diversity", x = "Observed Diversity", tag = "b)")
p.bias_noise

plot1 <- cowplot::plot_grid(p.hd, p.bias_noise,  ncol = 1)
plot1

cowplot::save_plot(here::here("output", "figures", "Fig3_HD_bfly.png"),
                   plot1, base_width = 10, base_height = 6)



# visualizing the diversity pattern for observed and estimated data -------

p.div.pat.obs <- ggplot(HD.df, aes(x = Metric, y = Obs.values, colour = Strata, fill = Strata)) +
  geom_boxplot(alpha = 0.5) + ylim(-4, 12) + scale_color_viridis_d(option = "A") +
  scale_fill_viridis_d(option = "A") + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "firebrick1") +
  labs(y = "Observed diversity", x = " ", tag = "a)") +
  theme(legend.position = "none")
p.div.pat.obs

p.div.pat.est <- ggplot(HD.df, aes(x = Metric, y = Est.values, colour = Strata, fill = Strata)) +
  geom_boxplot(alpha = 0.5) + ylim(-4, 12) + scale_color_viridis_d(option = "A") +
  scale_fill_viridis_d(option = "A") + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "firebrick1") +
  labs(y = "Estimated diversity", x = " ", tag = "b)") +
  theme(legend.position = "none")
p.div.pat.est

