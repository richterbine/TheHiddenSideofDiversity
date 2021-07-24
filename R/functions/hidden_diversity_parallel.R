# comm = a community matrix or data frame with sampling units in rows and species in columns
# N = an array with three dimensions: sampling units in rows, species in the column by each posterior 
# sampling. N can be a incidence or abundance based sampling.
# phy = a phylogenetic tree with branch lengths, preferably ultrametrized. If this object is provided, 
# the user can choose the null model to perform a SES for PD or MPD, the number or permutations to randomize
# the null models, and if in cases where N was a abundance-based array, if the abundance should be weighted.
# trait = a traits matrix, with species in the rows and traits in the columns. Besides the adjust of same 
# parameters described in phy, the user can set the argument binary = TRUE, if there are binary traits in the matrix

# function default
  # phy = NULL
  # trait = NULL
  # metrics = c("pd", "mpd")
  # binary = FALSE
  # abundance.weighted = FALSE
  # null.model = "taxa.labels"
  # runs = 499
  # parallel = 3

hidden.diversity<- function(comm, N, phy = NULL, trait = NULL, metrics = c("pd", "mpd"), 
                            binary = FALSE, 
                            abundance.weighted = FALSE, 
                            null.model = "taxa.labels", runs = 499, parallel = 3) {
  
  n.site <- dim(N)[1] # n.site: the number of sampling sites
  n.samp <- dim(N)[3] # n.samp: the number of posterior sampling
  
  # transforming N in occurrence data (y)
  y <- N
  for (i in 1: dim(y)[3]) {
    b = which(y[,,i] > 0) 
    y[,,i][b] = 1  
    y[,,i][-b] = 0  
  }
  
  # calculating the observed and estimated richness (TD)
  TD.df <- data.frame(TD.obs = apply(vegan::decostand(x = comm, method = "pa"), 1, sum),
                      TD.est = apply(apply(y, c(1,3), sum), 1, mean),
                      TD.sd = apply(apply(y, c(1,3), sum), 1, sd))
  TD.df$HD.TD <- (TD.df$TD.obs - TD.df$TD.est)/TD.df$TD.sd
  
  # calculating the observed and estimated abundance
  N.df <- data.frame(N.obs = apply(comm, 1, sum),
                     N.est = apply(apply(N, c(1,3), sum), 1, mean), 
                     N.sd = apply(apply(N, c(1,3), sum), 1, sd))
  N.df$HD.N <- (N.df$N.obs - N.df$N.est) / N.df$N.sd
  
  div_measures <- c("rich", "abund", "pd", "mpd")
  hd_metric <- pmatch(metrics, div_measures)
  
  # if a phylogenetic tree is provided 
  if (!is.null(phy)){
    if(any(hd_metric == 3)){ # calculation of PD
      # observed data
      pd.obs <- picante::ses.pd(samp = comm, tree = phy, null.model = null.model, 
                                runs = runs, include.root = F)
      
      # estimated data
      if (is.numeric(parallel)) {
        CL1 <- parallel::makeCluster(parallel, type = "PSOCK")
        newClusters <- TRUE
      }
      if (!inherits(CL1, "cluster")) {
        pd.ses <- array(NA, dim = c(n.site, 2, n.samp))
        for (i in 1:n.samp){
          temp_pd <- picante::ses.pd(samp = N[,,i], tree = phy, null.model = null.model, 
                                  runs = runs, include.root = F)
          pd.ses[ , 1, i] <- cbind(temp_pd[,6])
          }
        PD.df <- data.frame(SES.PD.obs = pd.obs[ ,"pd.obs.z"],
                            SES.PD.est = apply(pd.ses[,1,], 1, mean, na.rm = T), 
                            SES.PD.sd = apply(pd.ses[,1,], 1, sd, na.rm = T))
        
      }
      else {
        res_sesPD_samp <- parallel::parApply(cl = CL1, MARGIN = 3, X = y, FUN = picante::ses.pd, 
                                             tree = phy, 
                                             null.model = null.model,
                                             include.root = F,
                                             runs = runs)
       
        HD.comm <- list(SES.PDest = res_sesPD_samp)
        PD_est <- matrix(unlist(lapply(HD.comm$SES.PDest, function(x) x$pd.obs.z)), 
                         nrow = nrow(comm), ncol = dim(y)[3], 
               dimnames = list(rownames(comm), paste("samp", 1:dim(y)[3], sep = "_")
                               )
               )
        matrix_mean_SES_PD <- data.frame(matrix(c(apply(PD_est, MARGIN = 1, mean), 
                                                  apply(PD_est, MARGIN = 1, sd)), 
                                                nrow = nrow(comm), ncol = 2, 
                                                dimnames = list(rownames(comm), 
                                                                c("mean_ses.pd", "sd_ses.pd")),
                                                byrow = FALSE))
        
        PD.df <- data.frame(SES.PD.obs = pd.obs$pd.obs.z, 
                            SES.PD.est = matrix_mean_SES_PD$mean_ses.pd,
                            SES.PD.sd = matrix_mean_SES_PD$sd_ses.pd)
        
      }
     
      # calculating the Hidden Diversity for SES.PD
      PD.df$HD.PD <- (PD.df$SES.PD.obs - PD.df$SES.PD.est) / PD.df$SES.PD.sd
    }
    
    if(any(hd_metric == 4)){ # calculation of MPD
      # observed data
      mpd.obs <- picante::ses.mpd(samp = comm, dis = cophenetic(x = phy), null.model = null.model, 
                                  runs = runs)
      
      # estimated data
      mpd.ses <- array(NA, dim = c(n.site, 2, n.samp))
      
      if (is.numeric(parallel)) {
        CL1 <- parallel::makeCluster(parallel, type = "PSOCK")
        newClusters <- TRUE
      }
      
      if (!inherits(CL1, "cluster")) {
        mpd.ses <- array(NA, dim = c(n.site, 1, n.samp))
        for (i in 1:n.samp){
          temp_mpd <- picante::ses.mpd(samp = N[,,i], dis = cophenetic(x = phy), null.model = null.model, 
                                       runs = runs)
          
          mpd.ses[,1,i] <- cbind(temp_mpd[,6])
        }
        
        MPD.df <- data.frame(SES.MPD.obs = mpd.obs[ , "mpd.obs.z"],
                             SES.MPD.est = apply(mpd.ses[,1,], 1, mean, na.rm = T), 
                             SES.MPD.sd = apply(mpd.ses[,1,], 1, sd, na.rm = T))
        
      } 
      else {
        res_sesMPD_samp <- parallel::parApply(cl = CL1, MARGIN = 3, X = y, FUN = picante::ses.mpd, 
                                              dis = cophenetic(phy), 
                                              abundance.weighted = abundance.weighted,
                                              null.model = null.model,
                                              runs = runs)
        HD.comm <- list(SES.MPDest = res_sesMPD_samp)
        
        MPD_est <- matrix(unlist(lapply(HD.comm$SES.MPDest, function(x) x$mpd.obs.z)), 
                          nrow = nrow(comm), ncol = dim(y)[3], 
                          dimnames = list(rownames(comm), 
                                          paste("samp", 1:dim(y)[3], sep = "_")
                          )
        )
        matrix_mean_SES_MPD <- data.frame(matrix(c(apply(MPD_est, MARGIN = 1, mean),
                                                   apply(MPD_est, MARGIN = 1, sd)), 
                                      nrow = nrow(comm), ncol = 2, 
                                      dimnames = list(rownames(comm), 
                                                      c("mean_ses.pd", "sd_ses.mpd")), 
                                      byrow = FALSE))
        
        MPD.df <- data.frame(SES.MPD.obs = mpd.obs$mpd.obs.z, 
                             SES.MPD.est = matrix_mean_SES_MPD$mean_ses.pd, 
                             SES.MPD.sd = matrix_mean_SES_MPD$sd_ses.mpd)
        
      }
      
    # hidden for mpd phylo
    MPD.df$HD.MPD <- (MPD.df$SES.MPD.obs - MPD.df$SES.MPD.est) / MPD.df$SES.MPD.sd
    }
  }

  # calculation for traits
  if (!is.null(trait)){
    if(binary == TRUE){
      bin <- vector()
      for(i in 1:ncol(trait)){
        bin[i] <- is.integer(trait[, i]) | is.factor(trait[, i])
      }
      con.t <- which(bin == F)
      bin.t <- which(bin == T)
      t.dist <- ade4::dist.ktab(ade4::ktab.list.df(list(log(trait[, con.t]),
                                                        ade4::prep.binary(trait[, bin.t],
                                                                          col.blocks = ncol(trait[, bin.t])))),
                                type = c("Q", "B")) # create a dist matrix, considering mixed-variables
    } else {
      t.dist <- ade4::dist.ktab(ade4::ktab.list.df(list(log(trait))), type = "Q")
    }
    tree.func <- hclust(d = t.dist, method = "average") # clustering using UPGMA
    tree.func <- ape::as.phylo(tree.func)
    
    if(any(hd_metric == 3)){ # calculation of FD
      # observed data
      fd.obs <- picante::ses.pd(samp = comm, tree = tree.func, null.model = null.model, 
                                runs = runs, include.root = F)
      
      # estimated data
      if (is.numeric(parallel)) {
        CL1 <- parallel::makeCluster(parallel, type = "PSOCK")
        newClusters <- TRUE
      }
      
      if (!inherits(CL1, "cluster")) {
        fd.ses <- array(NA, dim = c(n.site, 1, n.samp))
        for (i in 1:n.samp){
          temp_fd <- picante::ses.pd(samp = N[,,i], tree = tree.func, null.model = null.model, 
                                     runs = runs, include.root = F)
          fd.ses[ , 1, i] <- cbind(temp_fd[,6])
        }
        FD.df <- data.frame(SES.FD.obs = fd.obs[ , "pd.obs.z"],
                            SES.FD.est = apply(fd.ses[,1,], 1, mean, na.rm = T), 
                            SES.FD.sd = apply(fd.ses[,1,], 1, sd, na.rm = T))
      }
      else {
        res_sesFD_samp <- parallel::parApply(cl = CL1, MARGIN = 3, X = y, FUN = picante::ses.pd, 
                                             tree = tree.func, 
                                             null.model = null.model,
                                             include.root = F,
                                             runs = runs)
        
        HD.comm <- list(SES.FDest = res_sesFD_samp)
        FD_est <- matrix(unlist(lapply(HD.comm$SES.FDest, function(x) x$pd.obs.z)), 
                         nrow = nrow(comm), ncol = dim(y)[3], 
                         dimnames = list(rownames(comm), 
                                         paste("samp", 1:dim(y)[3], sep = "_")
                         )
        )
        matrix_mean_SES_FD <- data.frame(matrix(c(apply(FD_est, MARGIN = 1, mean), 
                                                  apply(FD_est, MARGIN = 1, sd)), 
                                     nrow = nrow(comm), ncol = 2, 
                                     dimnames = list(rownames(comm), 
                                                     c("mean_ses.fd", "sd_ses.fd")), 
                                     byrow = FALSE))
        
        FD.df <- data.frame(SES.FD.obs = fd.obs$pd.obs.z, 
                            SES.FD.est = matrix_mean_SES_FD$mean_ses.fd,
                            SES.FD.sd = matrix_mean_SES_FD$sd_ses.fd)
        
      }
      
      # hidden FD
      FD.df$HD.FD <- (FD.df$SES.FD.obs - FD.df$SES.FD.est) / FD.df$SES.FD.sd
    }
    
    if(any(hd_metric == 4)){
      dist.func <- cophenetic(x = tree.func)
      mfd.obs <- picante::ses.mpd(samp = comm, dis = dist.func, null.model = null.model, 
                                  runs = runs)
      mfd.ses <- array(NA, dim = c(n.site, 1, n.samp))
      
      if (is.numeric(parallel)) {
        CL1 <- parallel::makeCluster(parallel, type = "PSOCK")
        newClusters <- TRUE
      }
      if (!inherits(CL1, "cluster")) {
        for (i in 1:n.samp){
        temp_mfd <- picante::ses.mpd(samp = N[,,i], dis = dist.func, null.model = null.model, 
                                       runs = runs)
          
          mfd.ses[ , 1, i] <- cbind(temp_mfd[, 6])
        }
        MFD.df <- data.frame(SES.MFD.obs = mfd.obs[ , "mpd.obs.z"],
                             SES.MFD.est = apply(mfd.ses[,1,], 1, mean, na.rm = T), 
                             SES.MFD.sd = apply(mfd.ses[,1,], 1, sd, na.rm = T))
        
      } 
        else {
        res_sesMFD_samp <- parallel::parApply(cl = CL1, MARGIN = 3, X = y, FUN = picante::ses.mpd, 
                                              dis = cophenetic(tree.func), 
                                              abundance.weighted = abundance.weighted,
                                              null.model = null.model,
                                              runs = runs)
        HD.comm <- list(SES.MFDest = res_sesMFD_samp)
        
        MFD_est <- data.frame(matrix(unlist(lapply(HD.comm$SES.MFDest, function(x) x$mpd.obs.z)), 
                                     nrow = nrow(comm), ncol = dim(y)[3], 
                          dimnames = list(rownames(comm), 
                                          paste("samp", 1:dim(y)[3], sep = "_")
                          )
        ))
        
        matrix_mean_SES_MFD <- data.frame(matrix(c(apply(MFD_est, MARGIN = 1, mean), 
                                                   apply(MFD_est, MARGIN = 1, sd)), 
                                      nrow = nrow(comm), ncol = 2, 
                                      dimnames = list(rownames(comm), 
                                                      c("mean_ses.mfd", "sd_ses.mfd")),
                                      byrow = FALSE))
        
        MFD.df <- data.frame(SES.MFD.obs = mfd.obs$mpd.obs.z, 
                             SES.MFD.est = matrix_mean_SES_MFD$mean_ses.mfd, 
                             SES.MFD.sd = matrix_mean_SES_MFD$sd_ses.mfd)
            
      }
      # hidden mfd
      MFD.df$HD.MFD <- (MFD.df$SES.MFD.obs - MFD.df$SES.MFD.est) / MFD.df$SES.MFD.sd
    }
  }
  
  if (newClusters) {
    parallel::stopCluster(CL1)
  }
  
  if(!is.null(trait) & !is.null(phy)){
    if(all(hd_metric == c(3, 4))){
      list_res <- vector(mode = "list", length = 6)
      names(list_res) <- c("TD", "Abund", "sesPD", "sesMPD", "sesFD", "sesMFD")
      list_res$TD <- TD.df
      list_res$Abund <- N.df
      list_res$sesPD <- PD.df
      list_res$sesMPD <- MPD.df
      list_res$sesFD <- FD.df
      list_res$sesMFD <- MFD.df
      
      for (i in 3:length(list_res)) {
        pos_obs_na <- which(is.na(list_res[[i]][,1]) == TRUE & is.na(list_res[[i]][,2]) == FALSE)
        list_res[[i]][pos_obs_na, paste("HD", gsub("ses","", names(list_res)[i]), sep = ".")] <- list_res[[i]][pos_obs_na, 2]/list_res[[i]][pos_obs_na, 3]
      }
      
      return(list_res)
    }
      else {
      if(hd_metric == 3){
        list_res <- vector(mode = "list", length = 4)
        names(list_res) <- c("TD", "Abund", "sesPD", "sesFD")
        list_res$TD <- TD.df
        list_res$Abund <- N.df
        list_res$sesPD <- PD.df
        list_res$sesFD <- FD.df
        
        for (i in 3:length(list_res)) {
          pos_obs_na <- which(is.na(list_res[[i]][,1]) == TRUE & is.na(list_res[[i]][,2]) == FALSE)
          list_res[[i]][pos_obs_na, paste("HD", gsub("ses","", names(list_res)[i]), sep = ".")] <- list_res[[i]][pos_obs_na, 2]/list_res[[i]][pos_obs_na, 3]
        }
        
        return(list_res)
      }
        
      if(hd_metric == 4){
        list_res <- vector(mode = "list", length = 4)
        names(list_res) <- c("TD", "Abund", "sesMPD", "sesMFD")
        list_res$TD <- TD.df
        list_res$Abund <- N.df
        list_res$sesMPD <- MPD.df
        list_res$sesMFD <- MFD.df
        
        for (i in 3:length(list_res)) {
          pos_obs_na <- which(is.na(list_res[[i]][,1]) == TRUE & is.na(list_res[[i]][,2]) == FALSE)
          list_res[[i]][pos_obs_na, paste("HD", gsub("ses","", names(list_res)[i]), sep = ".")] <- list_res[[i]][pos_obs_na, 2]/list_res[[i]][pos_obs_na, 3]
        }
        
        return(list_res)
      }
    }
  }
  
  if(is.null(phy) & !is.null(trait)){
    if(all(hd_metric == c(3, 4))){
      list_res <- vector(mode = "list", length = 4)
      names(list_res) <- c("TD", "Abund", "sesFD", "sesMFD")
      list_res$TD <- TD.df
      list_res$Abund <- N.df
      list_res$sesFD <- FD.df
      list_res$sesMFD <- MFD.df
      
      for (i in 3:length(list_res)) {
        pos_obs_na <- which(is.na(list_res[[i]][,1]) == TRUE & is.na(list_res[[i]][,2]) == FALSE)
        list_res[[i]][pos_obs_na, paste("HD", gsub("ses","", names(list_res)[i]), sep = ".")] <- list_res[[i]][pos_obs_na, 2]/list_res[[i]][pos_obs_na, 3]
      }
      
      return(list_res)
    } 
      else{
      if(any(hd_metric == 3)){
        list_res <- vector(mode = "list", length = 3)
        names(list_res) <- c("TD", "Abund", "sesFD")
        list_res$TD <- TD.df
        list_res$Abund <- N.df
        list_res$sesFD <- FD.df
        
        for (i in 3:length(list_res)) {
          pos_obs_na <- which(is.na(list_res[[i]][,1]) == TRUE & is.na(list_res[[i]][,2]) == FALSE)
          list_res[[i]][pos_obs_na, paste("HD", gsub("ses","", names(list_res)[i]), sep = ".")] <- list_res[[i]][pos_obs_na, 2]/list_res[[i]][pos_obs_na, 3]
        }
        
        return(list_res)
      }
      if(any(hd_metric == 4)){
        list_res <- vector(mode = "list", length = 3)
        names(list_res) <- c("TD", "Abund", "sesMFD")
        list_res$TD <- TD.df
        list_res$Abund <- N.df
        list_res$sesMFD <- MFD.df
        
        for (i in 3:length(list_res)) {
          pos_obs_na <- which(is.na(list_res[[i]][,1]) == TRUE & is.na(list_res[[i]][,2]) == FALSE)
          list_res[[i]][pos_obs_na, paste("HD", gsub("ses","", names(list_res)[i]), sep = ".")] <- list_res[[i]][pos_obs_na, 2]/list_res[[i]][pos_obs_na, 3]
        }
        
        return(list_res)
      }
    }
  }
  if(!is.null(phy) & is.null(trait)){
    if(all(hd_metric == c(3, 4))){
      list_res <- vector(mode = "list", length = 4)
      names(list_res) <- c("TD", "Abund", "sesPD", "sesMPD")
      list_res$TD <- TD.df
      list_res$Abund <- N.df
      list_res$sesPD <- PD.df
      list_res$sesMPD <- MPD.df
      
      for (i in 3:length(list_res)) {
        pos_obs_na <- which(is.na(list_res[[i]][,1]) == TRUE & is.na(list_res[[i]][,2]) == FALSE)
        list_res[[i]][pos_obs_na, paste("HD", gsub("ses","", names(list_res)[i]), sep = ".")] <- list_res[[i]][pos_obs_na, 2]/list_res[[i]][pos_obs_na, 3]
      }
      
      return(list_res)
    } 
      else {
      if(hd_metric == 3){
        list_res <- vector(mode = "list", length = 3)
        names(list_res) <- c("TD", "Abund", "sesPD")
        list_res$TD <- TD.df
        list_res$Abund <- N.df
        list_res$sesPD <- PD.df
        for (i in 3:length(list_res)) {
          pos_obs_na <- which(is.na(list_res[[i]][,1]) == TRUE & is.na(list_res[[i]][,2]) == FALSE)
          list_res[[i]][pos_obs_na, paste("HD", gsub("ses","", names(list_res)[i]), sep = ".")] <- list_res[[i]][pos_obs_na, 2]/list_res[[i]][pos_obs_na, 3]
        }
        
        return(list_res)
      }
        
      if(hd_metric == 4){
        list_res <- vector(mode = "list", length = 3)
        names(list_res) <- c("TD", "Abund", "sesMPD")
        list_res$TD <- TD.df
        list_res$Abund <- N.df
        list_res$sesMPD <- MPD.df
        
        for (i in 3:length(list_res)) {
          pos_obs_na <- which(is.na(list_res[[i]][,1]) == TRUE & is.na(list_res[[i]][,2]) == FALSE)
          list_res[[i]][pos_obs_na, paste("HD", gsub("ses","", names(list_res)[i]), sep = ".")] <- list_res[[i]][pos_obs_na, 2]/list_res[[i]][pos_obs_na, 3]
        }
        
        return(list_res)
      }
    }
  }
}


