comm = comm.bfly.obs
N = comm.bfly.est
phy = tree.bfly
trait = traits.bfly
abundance.weighted = TRUE
metrics = c("pd", "mpd")
null.model = "taxa.labels"
runs = 10
parallel = 4

hidden.diversity<- function(comm, N, phy = NULL, trait = NULL, metrics = c("pd", "mpd"), binary = FALSE, 
                            abundance.weighted = FALSE, 
                            null.model = "taxa.labels", runs = 10, parallel = 4) {
  
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
                      TD.mean = apply(apply(y, c(1,3), sum), 1, mean),
                      TD.sd = apply(apply(y, c(1,3), sum), 1, sd))
  
  # calculating the observed and estimated abundance
  N.df <- data.frame(N.obs = apply(comm, 1, sum),
                     N.mean = apply(apply(N, c(1,3), sum), 1, mean), 
                     N.sd = apply(apply(N, c(1,3), sum), 1, sd))
  
  div_measures <- c("rich", "abund", "pd", "mpd")
  hd_metric <- pmatch(metrics, div_measures)
  # if only the phylogenetic tree is provided 
  if (!is.null(phy)){
    if(any(hd_metric == 3)){ # calculation of PD
      pd.obs <- picante::ses.pd(samp = comm, tree = phy, null.model = null.model, 
                                runs = runs, include.root = F)
      
      # estimated data
      if (is.numeric(parallel)) {
        parallel <- parallel::makeCluster(parallel, type = "PSOCK")
        newClusters <- TRUE
      }
      if (!inherits(parallel, "cluster")) {
        pd.ses <- array(NA, dim = c(n.site, 2, n.samp))
        for (i in 1:n.samp){
          temp_pd <- picante::ses.pd(samp = N[,,i], tree = phy, null.model = null.model, 
                                  runs = runs, include.root = F)
          pd.ses[ , 1, i] <- cbind(temp_pd[,2])
          }
        PD.df <- data.frame(PD.obs = pd.obs[ , "pd.obs"],
                            SES.PD.obs = pd.obs[ , "pd.obs.z"],
                            PD.est = apply(pd.ses[,1,], 1, mean, na.rm = T),
                            SES.PD.est = apply(pd.ses[,2,], 1, mean, na.rm = T), 
                            PD.sd = apply(pd.ses[,1,], 1, sd, na.rm = T), 
                            SES.PD.sd = apply(pd.ses[,2,], 1, sd, na.rm = T),)
        
      }
      else {
        res_sesPD_samp <- parallel::parApply(cl = parallel, MARGIN = 3, X = y, FUN = picante::ses.pd, 
                                             tree = phy, 
                                             null.model = null.model,
                                             include.root = F,
                                             runs = runs)
       
        HD.comm <- list(SES.PDest = res_sesPD_samp)
        PD_est <- matrix(unlist(lapply(HD.comm$SES.PDest, function(x) x$pd.obs.z)), nrow = nrow(comm), ncol = dim(y)[3], 
               dimnames = list(rownames(comm), paste("samp", 1:dim(y)[3], sep = "_")
                               )
               )
        matrix_mean_SES_PD <- data.frame(matrix(c(apply(PD_est, MARGIN = 1, mean), apply(PD_est, MARGIN = 1, sd)), 
               nrow = nrow(comm), ncol = 2, 
               dimnames = list(rownames(comm), c("mean_ses.pd", "sd_ses.pd")), byrow = FALSE))
        PD.df <- data.frame(SES.PD.obs = pd.obs$pd.obs.z, 
                            SES.PD.est = matrix_mean_SES_PD$mean_ses.pd,
                            SES.PD.sd = matrix_mean_SES_PD$sd_ses.pd)
        
      }
     
      # hidden pd
      PD.df$hidden.PD <- (PD.df$SES.PD.obs - PD.df$SES.PD.est) / PD.df$SES.PD.sd
    }
    if(any(hd_metric == 4)){
      dist.phylo <- cophenetic(x = phy)
      mpd.obs <- picante::ses.mpd(samp = comm, dis = dist.phylo, null.model = null.model, 
                                  runs = runs)
      mpd.ses <- array(NA, dim = c(n.site, 2, n.samp))
      if (is.numeric(parallel)) {
        parallel <- parallel::makeCluster(parallel, type = "PSOCK")
        newClusters <- TRUE
      }
      if (!inherits(parallel, "cluster")) {
        mpd.ses <- array(NA, dim = c(n.site, 2, n.samp))
        for (i in 1:n.samp){
          
          temp_mpd <- picante::ses.mpd(samp = N[,,i], dis = dist.phylo, null.model = null.model, 
                                       runs = runs, include.root = F)
          
          mpd.ses[ , 2, i] <- cbind(temp_mpd[,6])
        }
        MPD.df <- data.frame(SES.MPD.obs = mpd.obs[ , "mpd.obs.z"],
                             SES.MPD.est = apply(mpd.ses[,2,], 1, mean, na.rm = T), 
                             SES.MPD.sd = apply(mpd.ses[,2,], 1, sd, na.rm = T),)
        
      } 
      else {
        res_sesMPD_samp <- parallel::parApply(cl = parallel, MARGIN = 3, X = y, FUN = picante::ses.mpd, 
                                              dis = cophenetic(phy), 
                                              abundance.weighted = abundance.weighted,
                                              null.model = null.model,
                                              runs = runs)
        HD.comm <- list(SES.MPDest = res_sesMPD_samp)
        
        MPD_est <- matrix(unlist(lapply(HD.comm$SES.MPDest, function(x) x$mpd.obs.z)), nrow = nrow(comm), ncol = dim(y)[3], 
                          dimnames = list(rownames(comm), paste("samp", 1:dim(y)[3], sep = "_")
                          )
        )
        matrix_mean_SES_MPD <- data.frame(matrix(c(apply(MPD_est, MARGIN = 1, mean), apply(MPD_est, MARGIN = 1, sd)), 
                                      nrow = nrow(comm), ncol = 2, 
                                      dimnames = list(rownames(comm), c("mean_ses.pd", "sd_ses.mpd")), byrow = FALSE))
        
        MPD.df <- data.frame(SES.MPD.obs = mpd.obs$mpd.obs.z, 
                             SES.MPD.est = matrix_mean_SES_MPD$mean_ses.pd, 
                             SES.MPD.sd = matrix_mean_SES_MPD$sd_ses.mpd)
        
      }
    }
    # hidden for mpd phylo
    MPD.df$hidden.MPD <- (MPD.df$SES.MPD.obs - MPD.df$SES.MPD.est) / MPD.df$SES.MPD.sd
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
      fd.obs <- picante::ses.pd(samp = comm, tree = tree.func, null.model = null.model, 
                                runs = runs, include.root = F)
      
      # estimated data
      if (is.numeric(parallel)) {
        parallel <- parallel::makeCluster(parallel, type = "PSOCK")
        newClusters <- TRUE
      }
      if (!inherits(parallel, "cluster")) {
        fd.ses <- array(NA, dim = c(n.site, 2, n.samp))
        for (i in 1:n.samp){
          temp_fd <- picante::ses.pd(samp = N[,,i], tree = tree.func, null.model = null.model, 
                                     runs = runs, include.root = F)
          fd.ses[ , 1, i] <- cbind(temp_fd[,2])
        }
        FD.df <- data.frame(SES.FD.obs = fd.obs[ , "pd.obs.z"],
                            SES.FD.est = apply(fd.ses[,2,], 1, mean, na.rm = T), 
                            SES.FD.sd = apply(fd.ses[,2,], 1, sd, na.rm = T))
        
        
      }
      else {
        res_sesFD_samp <- parallel::parApply(cl = parallel, MARGIN = 3, X = y, FUN = picante::ses.pd, 
                                             tree = tree.func, 
                                             null.model = null.model,
                                             include.root = F,
                                             runs = runs)
        
        HD.comm <- list(SES.FDest = res_sesFD_samp)
        FD_est <- matrix(unlist(lapply(HD.comm$SES.FDest, function(x) x$pd.obs.z)), nrow = nrow(comm), ncol = dim(y)[3], 
                         dimnames = list(rownames(comm), paste("samp", 1:dim(y)[3], sep = "_")
                         )
        )
        matrix_mean_SES_FD <- data.frame(matrix(c(apply(FD_est, MARGIN = 1, mean), apply(FD_est, MARGIN = 1, sd)), 
                                     nrow = nrow(comm), ncol = 2, 
                                     dimnames = list(rownames(comm), c("mean_ses.fd", "sd_ses.fd")), byrow = FALSE))
        FD.df <- data.frame(SES.FD.obs = pd.obs$pd.obs.z, 
                            SES.FD.est = matrix_mean_SES_FD$mean_ses.fd,
                            SES.FD.sd = matrix_mean_SES_FD$sd_ses.fd
                            )
        
      }
      # hidden FD
      FD.df$hidden.FD <- (FD.df$SES.FD.obs - FD.df$SES.FD.est) / FD.df$SES.FD.sd
      
    }
    if(any(hd_metric == 4)){
      dist.func <- cophenetic(x = tree.func)
      mfd.obs <- picante::ses.mpd(samp = comm, dis = dist.func, null.model = null.model, 
                                  runs = runs)
      mfd.ses <- array(NA, dim = c(n.site, 2, n.samp))
      if (!inherits(parallel, "cluster")) {
        mfd.ses <- array(NA, dim = c(n.site, 2, n.samp))
        for (i in 1:n.samp){
          
          temp_mfd <- picante::ses.mpd(samp = N[,,i], dis = dist.func, null.model = null.model, 
                                       runs = runs, include.root = F)
          
          mfd.ses[ , 2, i] <- cbind(temp_mfd[,6])
        }
        MFD.df <- data.frame(SES.MFD.obs = mfd.obs[ , "mfd.obs.z"],
                             SES.MFD.est = apply(mfd.ses[,2,], 1, mean, na.rm = T), 
                             SES.MFD.sd = apply(mfd.ses[,2,], 1, sd, na.rm = T),)
        
      } else {
        res_sesMFD_samp <- parallel::parApply(cl = parallel, MARGIN = 3, X = y, FUN = picante::ses.mpd, 
                                              dis = cophenetic(tree.func), 
                                              abundance.weighted = abundance.weighted,
                                              null.model = null.model,
                                              runs = runs)
        HD.comm <- list(SES.MFDest = res_sesMFD_samp)
        
        MFD_est <- data.frame(matrix(unlist(lapply(HD.comm$SES.MFDest, function(x) x$mpd.obs.z)), nrow = nrow(comm), ncol = dim(y)[3], 
                          dimnames = list(rownames(comm), paste("samp", 1:dim(y)[3], sep = "_")
                          )
        ))
        matrix_mean_SES_MFD <- data.frame(matrix(c(apply(MFD_est, MARGIN = 1, mean), apply(MFD_est, MARGIN = 1, sd)), 
                                      nrow = nrow(comm), ncol = 2, 
                                      dimnames = list(rownames(comm), c("mean_ses.mfd", "sd_ses.mfd")), byrow = FALSE))
        
        MFD.df <- data.frame(SES.MFD.obs = mfd.obs$mpd.obs.z, 
                             SES.MFD.est = matrix_mean_SES_MFD$mean_ses.mfd, 
                             SES.MFD.sd = matrix_mean_SES_MFD$sd_ses.mfd)
            
      }
      # hidden mfd
      MFD.df$hidden.MFD <- (MFD.df$SES.MFD.obs - MFD.df$SES.MFD.est) / MFD.df$SES.MFD.sd
    }
  }
  if (newClusters) {
    parallel::stopCluster(parallel)
  }
  if(!is.null(trait) & !is.null(phy)){
    if(all(hd_metric == c(3, 4))){
      list_res <- vector(mode = "list", length = 6)
      names(list_res) <- c("TD", "Abund", "sesPD", "sesMPD", "sesFD", "ses.MFD")
      list_res$TD <- TD.df
      list_res$Abund <- N.df
      list_res$sesPD <- PD.df
      list_res$sesMPD <- MPD.df
      list_res$sesFD <- FD.df
      list_res$sesMFD <- MFD.df
      return(list_res)
    } else{
      if(hd_metric == 3){
        list_res <- vector(mode = "list", length = 4)
        names(list_res) <- c("TD", "Abund", "sesPD", "sesFD")
        list_res$TD <- TD.df
        list_res$Abund <- N.df
        list_res$sesPD <- PD.df
        list_res$sesFD <- FD.df
        return(list_res)
      }
      if(hd_metric == 4){
        list_res <- vector(mode = "list", length = 4)
        names(list_res) <- c("TD", "Abund", "sesMPD", "sesMFD")
        list_res$TD <- TD.df
        list_res$Abund <- N.df
        list_res$sesMPD <- MPD.df
        list_res$sesMFD <- MFD.df
        return(list_res)
      }
    }
  }
  if(is.null(phylo) & !is.null(trait)){
    if(all(hd_metric == c(3, 4))){
      list_res <- vector(mode = "list", length = 4)
      names(list_res) <- c("TD", "Abund", "sesFD", "sesMFD")
      list_res$TD <- TD.df
      list_res$Abund <- N.df
      list_res$sesFD <- FD.df
      list_res$sesMFD <- MFD.df
      return(list_res)
    } else{
      if(any(hd_metric == 3)){
        list_res <- vector(mode = "list", length = 3)
        names(list_res) <- c("TD", "Abund", "sesFD")
        list_res$TD <- TD.df
        list_res$Abund <- N.df
        list_res$sesFD <- FD.df
        return(list_res)
      }
      if(any(hd_metric == 4)){
        list_res <- vector(mode = "list", length = 3)
        names(list_res) <- c("TD", "Abund", "sesMFD")
        list_res$TD <- TD.df
        list_res$Abund <- N.df
        list_res$sesMFD <- MFD.df
        return(list_res)
      }
    }
  }
  if(!is.null(phylo) & is.null(trait)){
    if(all(hd_metric == c(3, 4))){
      list_res <- vector(mode = "list", length = 4)
      names(list_res) <- c("TD", "Abund", "sesPD", "sesMPD")
      list_res$TD <- TD.df
      list_res$Abund <- N.df
      list_res$sesPD <- PD.df
      list_res$sesMPD <- MPD.df
      pos_hd0_phy <- na.omit(match(which(is.na(list_res$sesMPD$SES.MPD.obs) == TRUE), 
                                   which(is.na(list_res$sesMPD$SES.MPD.est) == TRUE))
      )
      list_res$sesMPD[pos_hd0_phy, "hidden.MPD"] <- 0
      pos_hd_1_phy <- which(is.na(list_res$sesMPD$SES.MPD.obs) == TRUE & is.na(list_res$sesMPD$SES.MPD.est) == FALSE)
      list_res$sesMPD[pos_hd_1_phy, "hidden.MPD"] <- list_res$sesMPD[pos_hd_1, "SES.MPD.est"]/list_res$sesMPD[pos_hd_1, "SES.MPD.sd"]
      
      pos_hd0_phy <- na.omit(match(which(is.na(list_res$sesPD$SES.PD.obs) == TRUE), 
                                   which(is.na(list_res$sesPD$SES.PD.est) == TRUE))
      )
      list_res$sesPD[pos_hd0_phy, "hidden.PD"] <- 0
      pos_hd_1_phy <- which(is.na(list_res$sesPD$SES.PD.obs) == TRUE & is.na(list_res$sesPD$SES.PD.est) == FALSE)
      list_res$sesPD[pos_hd_1_phy, "hidden.PD"] <- list_res$sesPD[pos_hd_1, "SES.PD.est"]/list_res$sesPD[pos_hd_1, "SES.PD.sd"]
      
      return(list_res)
    } else{
      if(hd_metric == 3){
        list_res <- vector(mode = "list", length = 3)
        names(list_res) <- c("TD", "Abund", "sesPD")
        list_res$TD <- TD.df
        list_res$Abund <- N.df
        list_res$sesPD <- PD.df
        pos_hd0_phy <- na.omit(match(which(is.na(list_res$sesPD$SES.PD.obs) == TRUE), 
                                 which(is.na(list_res$sesPD$SES.PD.est) == TRUE))
        )
        list_res$sesPD[pos_hd0_phy, "hidden.PD"] <- 0
        pos_hd_1_phy <- which(is.na(list_res$sesPD$SES.PD.obs) == TRUE & is.na(list_res$sesPD$SES.PD.est) == FALSE)
        list_res$sesPD[pos_hd_1_phy, "hidden.PD"] <- list_res$sesPD[pos_hd_1, "SES.PD.est"]/list_res$sesPD[pos_hd_1, "SES.PD.sd"]
        list_res$sesPD[pos_hd_1_phy, "hidden.PD"] <- list_res$sesPD[pos_hd_1, "SES.PD.est"]/list_res$sesPD[pos_hd_1, "SES.PD.sd"]
        return(list_res)
      }
      if(hd_metric == 4){
        list_res <- vector(mode = "list", length = 3)
        names(list_res) <- c("TD", "Abund", "sesMPD")
        list_res$TD <- TD.df
        list_res$Abund <- N.df
        list_res$sesMPD <- MPD.df
        pos_hd0 <- na.omit(match(which(is.na(list_res$sesMFD$SES.MPD.obs) == TRUE), 
                                 which(is.na(list_res$sesMFD$SES.MPD.est) == TRUE))
                           )
        list_res$sesMPD[pos_hd0, "hidden.MPD"] <- 0
        pos_hd_1 <- which(is.na(list_res$sesMPD$SES.MPD.obs) == TRUE & is.na(list_res$sesMPD$SES.MPD.est) == FALSE)
        list_res$sesMPD[pos_hd_1, "hidden.MPD"] <- list_res$sesMPD[pos_hd_1, "SES.MPD.est"]/list_res$sesMPD[pos_hd_1, "SES.MPD.sd"]
        return(list_res)
      }
    }
  }
}

