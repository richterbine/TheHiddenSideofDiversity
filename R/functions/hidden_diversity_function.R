hidden.diversity<- function(comm, N, phy = NULL, trait = NULL, binary = FALSE, 
                            abundance.weighted = FALSE, 
                            null.model = "taxa.labels", runs = 999) {
  
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
  
  
  # if only the phylogenetic tree is provided 
  if (is.null(trait) == TRUE & is.null(phy) == FALSE){ 
    pd.ses <- array(NA, dim = c(n.site, 2, n.samp))
    mpd.ses <- array(NA, dim = c(n.site, 2, n.samp))
    
    dist.phylo <- cophenetic(x = phy)
    
    # calculating the phylogenetic diversity (SES.PD)
    # observed data
    pd.obs <- picante::ses.pd(samp = comm, tree = phy, null.model = null.model, 
                              runs = runs, include.root = F)
    
    # estimated data
    for (i in 1:n.samp){
      temp <- picante::ses.pd(samp = N[,,i], tree = phy, null.model = null.model, 
                              runs = runs, include.root = F)
      pd.ses[,1,i] <- cbind(temp[,2])
      pd.ses[,2,i] <- cbind(temp[,6])
    }
    
    PD.df <- data.frame(PD.obs = pd.obs[, "pd.obs"],
                        SES.PD.obs = pd.obs[, "pd.obs.z"],
                        PD.est = apply(pd.ses[,1,], 1, mean, na.rm = T),
                        SES.PD.est = apply(pd.ses[,2,], 1, mean, na.rm = T), 
                        PD.sd = apply(pd.ses[,1,], 1, sd, na.rm = T), 
                        SES.PD.sd = apply(pd.ses[,2,], 1, sd, na.rm = T))
    
    ## calculating the mean phylogenetic diversity (SES.MPD)
    # observed data
    mpd.obs <- picante::ses.mpd(samp = comm, dis = dist.phylo, abundance.weighted = abundance.weighted,
                                null.model =  null.model, runs = runs)
    
    # estimated data
    for (i in 1:n.samp){ 
      temp <- picante::ses.mpd(samp = N[,,i], dis = dist.phylo, abundance.weighted = abundance.weighted, 
                               null.model = null.model, runs=runs)
      mpd.ses[,1,i] <- cbind(temp[,2])
      mpd.ses[,2,i] <- cbind(temp[,6])
    }
    
    MPD.df <- data.frame(MPD.obs = pd.obs[, "pd.obs"],
                         SES.MPD.obs = pd.obs[, "pd.obs.z"],
                         MPD.est = apply(mpd.ses[,1,], 1, mean, na.rm = T),
                         SES.MPD.est = apply(mpd.ses[,2,], 1, mean, na.rm = T), 
                         MPD.sd = apply(mpd.ses[,1,], 1, sd, na.rm = T), 
                         SES.MPD.sd = apply(mpd.ses[,2,], 1, sd, na.rm = T))
    
  # calculating the hidden diversity
  TD.df$hidden.TD <- (TD.df$TD.obs - TD.df$TD.mean) / TD.df$TD.sd
  N.df$hidden.N <- (N.df$N.obs - N.df$N.mean) / N.df$N.sd
  PD.df$hidden.PD <- (PD.df$SES.PD.obs - PD.df$SES.PD.est) / PD.df$SES.PD.sd
  MPD.df$hidden.MPD <- (MPD.df$SES.MPD.obs - MPD.df$SES.MPD.est) / MPD.df$SES.MPD.sd
  
  return(list(TD.df, N.df, PD.df, MPD.df)) 
  }
  
  
  # if only a traits matrix is provided
  if (is.null(trait) == FALSE & is.null(phy) == TRUE){
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
    
    fd.ses <- array(NA, dim = c(n.site, 2, n.samp))
    mfd.ses <- array(NA, dim = c(n.site, 2, n.samp))
    
    dist.func <- cophenetic(tree.func)
    
    # calculating the functional diversity (SES.FD)
    # observed data
    fd.obs <- picante::ses.pd(samp = comm, tree = tree.func, null.model = null.model,
                              runs = runs, include.root = F)
    
    # estimated data
    for (i in 1:n.samp){ 
      temp <- picante::ses.pd(samp = N[,,i], tree = tree.func, null.model = null.model, 
                              runs = runs, include.root = F)
      fd.ses[,1,i] <- cbind(temp[,2])
      fd.ses[,2,i] <- cbind(temp[,6])
    }
    
    FD.df <- data.frame(FD.obs = fd.obs[, "pd.obs"],
                        SES.FD.obs = fd.obs[, "pd.obs.z"],
                        FD.est = apply(fd.ses[,1,], 1, mean, na.rm = T),
                        SES.FD.est = apply(fd.ses[,2,], 1, mean, na.rm = T), 
                        FD.sd = apply(fd.ses[,1,], 1, sd, na.rm = T), 
                        SES.FD.sd = apply(fd.ses[,2,], 1, sd, na.rm = T))
    
    
    ## calculating the mean functional diversity (SES.MFD)
    # observed data
    mfd.obs <- picante::ses.mpd(samp = comm, dis = dist.func, null.model = null.model,
                                abundance.weighted = abundance.weighted, 
                                runs = runs)
    
    # estimated data
    for (i in 1:n.samp){ 
      temp <- picante::ses.mpd(samp = N[,,i], dis = dist.func, null.model = null.model,
                               abundance.weighted = abundance.weighted, runs = runs)
      mfd.ses[,1,i] <- cbind(temp[,2])
      mfd.ses[,2,i] <- cbind(temp[,6])
    }
    
    MFD.df <- data.frame(MFD.obs = mfd.obs[, "mpd.obs"],
                         SES.MFD.obs = mfd.obs[, "mpd.obs.z"],
                         MFD.est = apply(mfd.ses[,1,], 1, mean, na.rm = T),
                         SES.MFD.est = apply(mfd.ses[,2,], 1, mean, na.rm = T), 
                         MFD.sd = apply(mfd.ses[,1,], 1, sd, na.rm = T), 
                         SES.MFD.sd = apply(mfd.ses[,2,], 1, sd, na.rm = T))
    
    # calculating the hidden diversity
    TD.df$hidden.TD <- (TD.df$TD.obs - TD.df$TD.mean) / TD.df$TD.sd
    N.df$hidden.N <- (N.df$N.obs - N.df$N.mean) / N.df$N.sd
    FD.df$hidden.FD <- (FD.df$SES.FD.obs - FD.df$SES.FD.est) / FD.df$SES.FD.sd
    MFD.df$hidden.MFD <- (MFD.df$SES.MFD.obs - MFD.df$SES.MFD.est) / MFD.df$SES.MFD.sd
    
    return(list(TD.df, N.df, FD.df, MFD.df))
  }
  
  
  # if both phylogenetic tree and traits matrix are provided
  if (is.null(trait) == FALSE & is.null(phy) == FALSE){
    pd.ses <- array(NA, dim = c(n.site, 2, n.samp))
    mpd.ses <- array(NA, dim = c(n.site, 2, n.samp))
    
    dist.phylo <- cophenetic(x = phy)
    
    # calculating the phylogenetic diversity (SES.PD)
    # observed data
    pd.obs <- picante::ses.pd(samp = comm, tree = phy, null.model = null.model, 
                              runs = runs, include.root = F)
    
    # estimated data
    for (i in 1:n.samp){
      temp <- picante::ses.pd(samp = N[,,i], tree = phy, null.model = null.model, 
                              runs = runs, include.root = F)
      pd.ses[,1,i] <- cbind(temp[,2])
      pd.ses[,2,i] <- cbind(temp[,6])
    }
    
    PD.df <- data.frame(PD.obs = pd.obs[, "pd.obs"],
                        SES.PD.obs = pd.obs[, "pd.obs.z"],
                        PD.est = apply(pd.ses[,1,], 1, mean, na.rm = T),
                        SES.PD.est = apply(pd.ses[,2,], 1, mean, na.rm = T), 
                        PD.sd = apply(pd.ses[,1,], 1, sd, na.rm = T), 
                        SES.PD.sd = apply(pd.ses[,2,], 1, sd, na.rm = T))
    
    ## calculating the mean phylogenetic diversity (SES.MPD)
    # observed data
    mpd.obs <- picante::ses.mpd(samp = comm, dis = dist.phylo, abundance.weighted = abundance.weighted,
                                null.model =  null.model, runs = runs)
    
    # estimated data
    for (i in 1:n.samp){ 
      temp <- picante::ses.mpd(samp = N[,,i], dis = dist.phylo, abundance.weighted = abundance.weighted, 
                               null.model = null.model, runs=runs)
      mpd.ses[,1,i] <- cbind(temp[,2])
      mpd.ses[,2,i] <- cbind(temp[,6])
    }
    
    MPD.df <- data.frame(MPD.obs = pd.obs[, "pd.obs"],
                         SES.MPD.obs = pd.obs[, "pd.obs.z"],
                         MPD.est = apply(mpd.ses[,1,], 1, mean, na.rm = T),
                         SES.MPD.est = apply(mpd.ses[,2,], 1, mean, na.rm = T), 
                         MPD.sd = apply(mpd.ses[,1,], 1, sd, na.rm = T), 
                         SES.MPD.sd = apply(mpd.ses[,2,], 1, sd, na.rm = T))
    
    # CALCULATING THE FUNCTIONAL DIVERSITY
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
    
    fd.ses <- array(NA, dim = c(n.site, 2, n.samp))
    mfd.ses <- array(NA, dim = c(n.site, 2, n.samp))
    
    dist.func <- cophenetic(tree.func)
    
    # calculating the functional diversity (SES.FD)
    # observed data
    fd.obs <- picante::ses.pd(samp = comm, tree = tree.func, null.model = null.model,
                              runs = runs, include.root = F)
    
    # estimated data
    for (i in 1:n.samp){ 
      temp <- picante::ses.pd(samp = N[,,i], tree = tree.func, null.model = null.model, 
                              runs = runs, include.root = F)
      fd.ses[,1,i] <- cbind(temp[,2])
      fd.ses[,2,i] <- cbind(temp[,6])
    }
    
    FD.df <- data.frame(FD.obs = fd.obs[, "pd.obs"],
                        SES.FD.obs = fd.obs[, "pd.obs.z"],
                        FD.est = apply(fd.ses[,1,], 1, mean, na.rm = T),
                        SES.FD.est = apply(fd.ses[,2,], 1, mean, na.rm = T), 
                        FD.sd = apply(fd.ses[,1,], 1, sd, na.rm = T), 
                        SES.FD.sd = apply(fd.ses[,2,], 1, sd, na.rm = T))
    
    
    ## calculating the mean functional diversity (SES.MFD)
    # observed data
    mfd.obs <- picante::ses.mpd(samp = comm, dis = dist.func, null.model = null.model,
                                abundance.weighted = abundance.weighted, 
                                runs = runs)
    
    # estimated data
    for (i in 1:n.samp){ 
      temp <- picante::ses.mpd(samp = N[,,i], dis = dist.func, null.model = null.model,
                               abundance.weighted = abundance.weighted, runs = runs)
      mfd.ses[,1,i] <- cbind(temp[,2])
      mfd.ses[,2,i] <- cbind(temp[,6])
    }
    
    MFD.df <- data.frame(MFD.obs = mfd.obs[, "mpd.obs"],
                         SES.MFD.obs = mfd.obs[, "mpd.obs.z"],
                         MFD.est = apply(mfd.ses[,1,], 1, mean, na.rm = T),
                         SES.MFD.est = apply(mfd.ses[,2,], 1, mean, na.rm = T), 
                         MFD.sd = apply(mfd.ses[,1,], 1, sd, na.rm = T), 
                         SES.MFD.sd = apply(mfd.ses[,2,], 1, sd, na.rm = T))
    
    ## calculating the hidden diversity
    # calculating the hidden diversity
    TD.df$hidden.TD <- (TD.df$TD.obs - TD.df$TD.mean) / TD.df$TD.sd
    N.df$hidden.N <- (N.df$N.obs - N.df$N.mean) / N.df$N.sd
    PD.df$hidden.PD <- (PD.df$SES.PD.obs - PD.df$SES.PD.est) / PD.df$SES.PD.sd
    MPD.df$hidden.MPD <- (MPD.df$SES.MPD.obs - MPD.df$SES.MPD.est) / MPD.df$SES.MPD.sd
    FD.df$hidden.FD <- (FD.df$SES.FD.obs - FD.df$SES.FD.est) / FD.df$SES.FD.sd
    MFD.df$hidden.MFD <- (MFD.df$SES.MFD.obs - MFD.df$SES.MFD.est) / MFD.df$SES.MFD.sd
    
    return(list(TD.df, N.df, PD.df, MPD.df, FD.df, MFD.df))
  }
  
  # if neither the phylogenetic tree nor the traits matrix is provided 
  else { 
    ## calculating the hidden diversity
    TD.df$hidden.TD <- (TD.df$TD.obs - TD.df$TD.mean) / TD.df$TD.sd
    N.df$hidden.N <- (N.df$N.obs - N.df$N.mean) / N.df$N.sd
    
    return(list(TD.df, N.df))
  }
}


  