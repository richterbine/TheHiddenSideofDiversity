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

comm.bfly.est <- N[, , sample(1:dim(N)[3], 100)]
dimnames(comm.bfly.est) <- list(rownames(comm.bfly.obs), colnames(comm.bfly.obs), NULL)
str(comm.bfly.est)

# read the phylogenetic pruned tree for the species sampled at FLONA
tree.bfly <- ape::read.tree(here::here("output", "tree_bfly_FLONA.txt"))
str(tree.bfly)

# read the functional matrix by species
traits.bfly <- read.table(here::here("data", "processed", "Mean_traits_bfly.txt"), header = T)
str(traits.bfly)


# Run in a different way to be faster -------------------------------------

source(here::here("R", "functions", "hidden_diversity_function.R"))
ape::is.rooted(tree.bfly)

library(parallel)
ncor <- detectCores()
CL <- makeCluster(ncor)
clusterExport(CL, c("comm.bfly.obs", "comm.bfly.est", "tree.bfly",
                    "traits.bfly"))
clusterEvalQ(CL,library(vegan))
clusterEvalQ(CL,library(ape))
clusterEvalQ(CL,library(picante))
clusterEvalQ(CL,library(phytools))
clusterEvalQ(CL,library(ade4))

res_sesPD_samp <- parApply(cl = CL, MARGIN = 3, X = comm.bfly.est, 
                            function(x) {
                              picante::ses.pd(samp = x, tree = tree.bfly, 
                                              null.model = "taxa.labels",
                                              include.root = F,
                                              runs = 299)
                            }
                            
)
ncor <- detectCores()
CL <- makeCluster(ncor)
clusterExport(CL, c("comm.bfly.obs", "comm.bfly.est", "tree.bfly",
                    "traits.bfly"))
clusterEvalQ(CL,library(vegan))
clusterEvalQ(CL,library(ape))
clusterEvalQ(CL,library(picante))
clusterEvalQ(CL,library(phytools))
clusterEvalQ(CL,library(ade4))

res_sesMPD_samp <- parApply(cl = CL, MARGIN = 3, X = comm.bfly.est, 
                                        function(x) {
                                          picante::ses.mpd(samp = x, dis = cophenetic(tree.bfly), 
                                                           abundance.weighted = TRUE, 
                                                           null.model = "taxa.labels", runs = 499)
                                        }
  
                               )

bin <- vector()
  for(i in 1:ncol(traits.bfly)){
    bin[i] <- is.integer(traits.bfly[, i]) | is.factor(traits.bfly[, i])
  }
con.t <- which(bin == F)
bin.t <- which(bin == T)
t.dist <- ade4::dist.ktab(ade4::ktab.list.df(list(log(traits.bfly[, con.t]),
                                                    ade4::prep.binary(traits.bfly[, bin.t],
                                                                      col.blocks = ncol(traits.bfly[, bin.t])))),
                            type = c("Q", "B")) # create a dist matrix, considering mixed-variables
  

tree.func <- hclust(d = t.dist, method = "average") # clustering using UPGMA
tree.func <- ape::as.phylo(tree.func)
ape::write.tree(tree.func, here::here("output", "functional_dendrogram.txt"))

ncor <- detectCores()
CL <- makeCluster(ncor)
clusterExport(CL, c("comm.bfly.obs", "comm.bfly.est", "t.dist", "tree.func",
                    "traits.bfly"))
clusterEvalQ(CL,library(vegan))
clusterEvalQ(CL,library(ape))
clusterEvalQ(CL,library(picante))
clusterEvalQ(CL,library(phytools))
clusterEvalQ(CL,library(ade4))
res_sesFD_samp <- parApply(cl = CL, MARGIN = 3, X = comm.bfly.est, 
                           function(x) {
                             picante::ses.pd(samp = x, tree = tree.func, 
                                             include.root = F, 
                                              null.model = "taxa.labels", 
                                             runs = 299)
                           }
                           
)

ncor <- detectCores()
CL <- makeCluster(ncor)
clusterExport(CL, c("comm.bfly.obs", "comm.bfly.est", "t.dist", "tree.func",
                    "traits.bfly"))
clusterEvalQ(CL,library(vegan))
clusterEvalQ(CL,library(ape))
clusterEvalQ(CL,library(picante))
clusterEvalQ(CL,library(phytools))
clusterEvalQ(CL,library(ade4))
res_sesMFD_samp <- parApply(cl = CL, MARGIN = 3, X = comm.bfly.est, 
                           function(x) {
                             picante::ses.mpd(samp = x, dis = cophenetic(tree.func),
                                              abundance.weighted = T, 
                                             null.model = "taxa.labels", runs = 499)
                           }
                           
)

stopCluster(CL)

# observed data
PD_obs <- picante::ses.pd(samp = comm.bfly.obs, tree = tree.bfly, 
                          include.root = F, 
                          null.model = "taxa.labels", runs = 299)

MPD_obs <- picante::ses.mpd(samp = comm.bfly.obs, dis = cophenetic(tree.bfly), 
                            abundance.weighted = TRUE, 
                            null.model = "taxa.labels", runs = 499)

FD_obs <- picante::ses.pd(samp = comm.bfly.obs, tree = tree.func, 
                          include.root = F, 
                          null.model = "taxa.labels", runs = 299)

MFD_obs <- picante::ses.mpd(samp = comm.bfly.obs, dis = cophenetic(tree.func), 
                            abundance.weighted = TRUE, 
                            null.model = "taxa.labels", runs = 499)

HD.bfly <- list(SES.PDest = res_sesPD_samp, SES.MPDest = res_sesMPD_samp,SES.FDest = res_sesFD_samp, 
                SES.MFDest = res_sesMFD_samp,SES.PDobs = PD_obs, SES.MPDobs = MPD_obs, 
                SES.FDobs = FD_obs, SES.MFDobs = MFD_obs)
saveRDS(HD.bfly, here::here("output", "HD_bfly.rds"))


# end of run --------------------------------------------------------------


### calculating the hidden diversity
HD.bfly <- readRDS(here::here("output", "HD_bfly.rds"))

# calculate the mean and standar deviation for estimated communities
list_mean_sd <- list(sesPD.est = vector(mode = "list", length =2),
                     sesMPD.est = vector(mode = "list", length = 2),
                     sesFD.est = vector(mode = "list", length = 2),
                     sesMFD.est = vector(mode = "list", length = 2))
names(list_mean_sd$sesPD.est) <- c("mean", "sd")
names(list_mean_sd$sesMPD.est) <- c("mean", "sd")
names(list_mean_sd$sesFD.est) <- c("mean", "sd")
names(list_mean_sd$sesMFD.est) <- c("mean", "sd")

tmp.SES <- vector(mode = "list", length = 8)

for (j in 1:length(list_mean_sd)) {
  for(i in 1:length(tmp.SES)){
    tmp.SES[[i]] <- matrix(unlist(lapply(lapply(HD.bfly[[j]], function(x) apply(x, 1, function(y) y)), 
                       function(z){
                         z[i,]
                       } )), nrow = 300, ncol = 100, byrow = FALSE, 
         dimnames = list(rownames(HD.bfly[[j]][[1]]), 
                         paste("samp", 1:100, sep = ".")))
}

matrix_ses <- matrix(unlist(lapply(tmp.SES, function(x){
  apply(x, 1, function(y) mean(y, na.rm = TRUE))
})), nrow = 300, ncol = 8, dimnames = list(rownames(HD.bfly[[j]][[1]]),
                                           colnames(HD.bfly[[j]][[1]])
                                           ), byrow = FALSE
)
matrix_ses.sd <- matrix(unlist(lapply(tmp.SES, function(x){
  apply(x, 1, function(y) sd(y, na.rm = TRUE))
})), nrow = 300, ncol = 8, dimnames = list(rownames(HD.bfly[[j]][[1]]),
                                           colnames(HD.bfly[[j]][[1]])
), byrow = FALSE
)

list_mean_sd[[j]][[1]] <- matrix_ses
list_mean_sd[[j]][[2]] <- matrix_ses.sd

}

length(list_mean_sd$sesMFD.est$mean[,1])

pos <- which(is.na(HD.bfly$SES.PDobs$pd.obs.z))
HD.bfly$SES.PDobs[pos, "pd.obs.z"] <- 0

pos <- which(is.na(HD.bfly$SES.FDobs$pd.obs.z))
HD.bfly$SES.FDobs[pos, "pd.obs.z"] <- 0

pos <- which(is.na(HD.bfly$SES.MPDobs$mpd.obs.z))
HD.bfly$SES.MPDobs[pos, "mpd.obs.z"] <- 0

pos <- which(is.na(HD.bfly$SES.MFDobs$mpd.obs.z))
HD.bfly$SES.MFDobs[pos, "mpd.obs.z"] <- 0

res.HD.bfly <- data.frame(H.TD = (HD.bfly[[5]]$ntaxa - list_mean_sd$sesPD.est$mean[,"ntaxa"])/list_mean_sd$sesPD.est$sd[,"ntaxa"],
                          H.PD = (HD.bfly[[5]]$pd.obs.z - list_mean_sd$sesPD.est$mean[,"pd.obs.z"])/list_mean_sd$sesPD.est$sd[,"pd.obs.z"],
                          H.FD = (HD.bfly[[7]]$pd.obs.z - list_mean_sd$sesFD.est$mean[,"pd.obs.z"])/list_mean_sd$sesFD.est$sd[,"pd.obs.z"],
                          H.MPD = (HD.bfly[[6]]$mpd.obs.z - list_mean_sd$sesMPD.est$mean[,"mpd.obs.z"])/list_mean_sd$sesMPD.est$sd[,"mpd.obs.z"],
                          H.MFD = (HD.bfly[[8]]$mpd.obs.z - list_mean_sd$sesMFD.est$mean[,"mpd.obs.z"])/list_mean_sd$sesMFD.est$sd[,"mpd.obs.z"])
res.HD.bfly$Strata <- substr(rownames(res.HD.bfly), 1, 1)
res.HD.bfly$Stratum <- factor(res.HD.bfly$Strata, labels = c("Canopy", "Understory"))

df.hd <- data.frame(HD.values = c(res.HD.bfly$H.TD, res.HD.bfly$H.PD, res.HD.bfly$H.FD,
                                  res.HD.bfly$H.MPD, res.HD.bfly$H.MFD), 
                    Div.names = rep(colnames(res.HD.bfly[, 1:5]), each = nrow(res.HD.bfly),1),
                    Strata = rep(res.HD.bfly$Stratum), each = nrow(res.HD.bfly), 1)

library(ggplot2)
library(viridis)

p.hd <- ggplot(data = df.hd, aes(y = HD.values, x = Div.names,
                                 colour = Strata, fill = Strata)) +
  geom_boxplot() + scale_color_viridis_d(option = "A") +
  scale_fill_viridis_d(option = "A", alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = 2, color = "firebrick1")
  
p.hd
cowplot::save_plot(here::here("output", "figures", "Fig1_hdbfly.png"), p.hd,
                   base_width = 8)

# testing if the HD differs between strata
mod.td <- lm(res.HD.bfly$H.TD ~ res.HD.bfly$Strata)
summary(mod.td)

mod.pd <- lm(res.HD.bfly$H.PD ~ res.HD.bfly$Strata)
summary(mod.pd)

mod.fd <- lm(res.HD.bfly$H.FD ~ res.HD.bfly$Strata)
summary(mod.fd)

mod.mpd <- lm(res.HD.bfly$H.MPD ~ res.HD.bfly$Strata)
summary(mod.mpd)

mod.mfd <- lm(res.HD.bfly$H.MFD ~ res.HD.bfly$Strata)
summary(mod.mfd)

out.lm <- rbind(summary(mod.td)[[4]], summary(mod.pd)[[4]],
                summary(mod.fd)[[4]], summary(mod.mpd)[[4]],
                summary(mod.mfd)[[4]])


# correlation among td and pd/fd
mod.td.pd <- lm(res.HD.bfly$H.TD ~ res.HD.bfly$H.PD)
summary(mod.td.pd)
cor.test(res.HD.bfly$H.TD, res.HD.bfly$H.PD)
plot(H.TD ~ H.PD, data = res.HD.bfly)
abline(mod.td.pd)

mod.td.fd <- lm(res.HD.bfly$H.TD ~ res.HD.bfly$H.FD)
summary(mod.td.fd)
cor.test(res.HD.bfly$H.TD, res.HD.bfly$H.FD)
plot(res.HD.bfly$H.TD ~res.HD.bfly$H.FD)
abline(mod.td.fd)

out.cor <- rbind(summary(mod.td.pd)[[4]], summary(mod.td.fd)[[4]])



HD.bfly <- readRDS(here::here("output", "HD_bfly.rds"))

bfly.obs.mpd <- HD.bfly$SES.MPDobs

bfly.obs.mpd$strata <- substr(rownames(bfly.obs.mpd), 1, 1)
bfly.obs.mpd <- na.omit(bfly.obs.mpd)

summary(lm(mpd.obs.z ~ strata, data = bfly.obs.mpd))
boxplot(bfly.obs.mpd$mpd.obs.z ~ bfly.obs.mpd$strata)

colnames(HD.bfly$SES.MPDobs) <- colnames(HD.bfly$SES.MFDobs) <- colnames(HD.bfly$SES.PDobs)

df.obs.div <- rbind(HD.bfly$SES.PDobs, HD.bfly$SES.FDobs, 
                    HD.bfly$SES.MPDobs, HD.bfly$SES.MFDobs)

df.obs.div$strata <- substr(rownames(df.obs.div), 1, 1)
df.obs.div$div <- rep(names(HD.bfly)[5:8], each = nrow(HD.bfly$SES.PDobs))

p.div.obd <- ggplot(data = na.omit(df.obs.div), aes(y = pd.obs.z, x = div,
                                      colour = strata, fill = strata)) +
  geom_boxplot() + scale_color_viridis_d(option = "A") +
  scale_fill_viridis_d(option = "A", alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = 2, color = "firebrick1")
p.div.obd

cowplot::plot_grid(p.hd, p.div.obd)


length(na.omit(HD.bfly[[5]]$pd.obs.z - list_mean_sd$sesPD.est$mean[,"pd.obs.z"]))

a = rnorm(100)
b = rnorm(100)
ts <- as.data.frame(cbind(a, b, a-b, a < b, a < 0, b < 0))
ts$stat <- ifelse(ts$V5 == 0 & ts$V6 == 0, "A",
                  ifelse(ts$V5 == 0 & ts$V6 == 1, "B",
                         ifelse(ts$V5 == 1 & ts$V6 == 0, "C", 
                                ifelse(ts$V5 == 1 & ts$V6 == 1, "D", NA))))


ggplot(ts, aes(x = a, y = b, col = stat, shape = as.factor(V4)), size = 3) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "lightblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "lightblue") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "firebrick1")
ts

df.hd <- data.frame(Dobs = c(HD.bfly$SES.PDobs$pd.obs.z, HD.bfly$SES.MPDobs$mpd.obs.z,
                            HD.bfly$SES.FDobs$pd.obs.z, HD.bfly$SES.MFDobs$mpd.obs.z), 
                    Dest = c(list_mean_sd$sesPD.est$mean[,"pd.obs.z"], 
                            list_mean_sd$sesMPD.est$mean[,"mpd.obs.z"], 
                            list_mean_sd$sesFD.est$mean[,"pd.obs.z"],
                            list_mean_sd$sesMFD.est$mean[,"mpd.obs.z"]),
                    HD = c(res.HD.bfly$H.PD, res.HD.bfly$H.MPD, res.HD.bfly$H.FD,
                           res.HD.bfly$H.MFD),
                  Type = rep(c("SES.PD", "SES.MPD", "SES.FD", "SES.MFD"), 
                             each = nrow(res.HD.bfly)),
                  Strata = rep(c("Canopy", "Undestory"), each = 150))


ggplot(df.hd, aes(x = Dobs, y = Dest, colour = Type), alpha, .5, size = 4) +
  geom_point() + scale_color_viridis_d(option = "A") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "lightblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "lightblue") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "firebrick1") +
  facet_wrap(~Strata)


