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

# source(here::here("R", "functions", "hidden_diversity_function.R"))
# ape::is.rooted(tree.bfly)

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

# SES.PD
res_sesPD_samp <- parApply(cl = CL, MARGIN = 3, X = comm.bfly.est, 
                           function(x) {
                             picante::ses.pd(samp = x, tree = tree.bfly, 
                                             null.model = "taxa.labels",
                                             include.root = F,
                                             runs = 299)
                           }
)

# abundance-based SES.MPD
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


# incidence-based SES.MPD 
ncor <- detectCores()
CL <- makeCluster(ncor)
clusterExport(CL, c("comm.bfly.obs", "comm.bfly.est", "tree.bfly",
                    "traits.bfly"))
clusterEvalQ(CL,library(vegan))
clusterEvalQ(CL,library(ape))
clusterEvalQ(CL,library(picante))
clusterEvalQ(CL,library(phytools))
clusterEvalQ(CL,library(ade4))

res_sesMPDi_samp <- parApply(cl = CL, MARGIN = 3, X = comm.bfly.est, 
                             function(x) {
                               picante::ses.mpd(samp = x, dis = cophenetic(tree.bfly), 
                                                abundance.weighted = FALSE, 
                                                null.model = "taxa.labels", runs = 499)
                             }
)


# Funcitional diversity ---------------------------------------------------

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

# SES.FD
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

# abundance-based SES.MFD
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
                                               abundance.weighted = TRUE, 
                                               null.model = "taxa.labels", runs = 499)
                            }
                            
)

# incidence-based SES.MFD
ncor <- detectCores()
CL <- makeCluster(ncor)
clusterExport(CL, c("comm.bfly.obs", "comm.bfly.est", "t.dist", "tree.func",
                    "traits.bfly"))
clusterEvalQ(CL,library(vegan))
clusterEvalQ(CL,library(ape))
clusterEvalQ(CL,library(picante))
clusterEvalQ(CL,library(phytools))
clusterEvalQ(CL,library(ade4))

res_sesMFDi_samp <- parApply(cl = CL, MARGIN = 3, X = comm.bfly.est, 
                            function(x) {
                              picante::ses.mpd(samp = x, dis = cophenetic(tree.func),
                                               abundance.weighted = FALSE, 
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

MPDi_obs <- picante::ses.mpd(samp = comm.bfly.obs, dis = cophenetic(tree.bfly), 
                            abundance.weighted = FALSE, 
                            null.model = "taxa.labels", runs = 499)

FD_obs <- picante::ses.pd(samp = comm.bfly.obs, tree = tree.func, 
                          include.root = F, 
                          null.model = "taxa.labels", runs = 299)

MFD_obs <- picante::ses.mpd(samp = comm.bfly.obs, dis = cophenetic(tree.func), 
                            abundance.weighted = TRUE, 
                            null.model = "taxa.labels", runs = 499)

MFDi_obs <- picante::ses.mpd(samp = comm.bfly.obs, dis = cophenetic(tree.func), 
                            abundance.weighted = FALSE, 
                            null.model = "taxa.labels", runs = 499)

HD.bfly <- list(SES.PDest = res_sesPD_samp, SES.MPDest = res_sesMPD_samp, SES.MPDiest = res_sesMPDi_samp,
                SES.FDest = res_sesFD_samp, SES.MFDest = res_sesMFD_samp, SES.MFDiest = res_sesMFDi_samp,
                SES.PDobs = PD_obs, SES.MPDobs = MPD_obs, SES.MPDiobs = MPDi_obs,
                SES.FDobs = FD_obs, SES.MFDobs = MFD_obs, SES.MFDiobs = MFDi_obs)

saveRDS(HD.bfly, here::here("output", "HD_bfly.rds"))


# end of run --------------------------------------------------------------
