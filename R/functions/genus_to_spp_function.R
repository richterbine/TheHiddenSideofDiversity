## function to include species at a genre-level tree
# tree = a genus-level tree
# genus.species.list = a matrix or data.frame with 2 columns, where the 1st is the genus names
 # and the 2nd is the species names

treedata_modif<- function (phy, data, sort = FALSE, warnings = TRUE) 
{
  dm = length(dim(data))
  if (is.vector(data)) {
    data <- as.matrix(data)
  }
  if (is.factor(data)) {
    data <- as.matrix(data)
  }
  if (is.array(data) & length(dim(data)) == 1) {
    data <- as.matrix(data)
  }
  if (is.null(rownames(data))) {
    stop("names for 'data' must be supplied")
  }
  else {
    data.names <- rownames(data)
  }
  nc <- geiger::name.check(phy, data)
  if (is.na(nc[[1]][1]) | nc[[1]][1] != "OK") {
    if (length(nc[[1]] != 0)) {
      phy = ape::drop.tip(phy, as.character(nc[[1]]))
      if (warnings) {
        warning(paste("The following tips were not found in 'data' and were dropped from 'phy':\n\t", 
                      paste(nc[[1]], collapse = "\n\t"), sep = ""))
      }
    }
    if (length(nc[[2]] != 0)) {
      m <- match(data.names, nc[[2]])
      data = as.matrix(data[is.na(m), ])
      data.names <- data.names[is.na(m)]
      if (warnings) {
        warning(paste("The following tips were not found in 'phy' and were dropped from 'data':\n\t", 
                      paste(nc[[2]], collapse = "\n\t"), sep = ""))
      }
    }
  }
  order <- match(data.names, phy$tip.label)
  rownames(data) <- phy$tip.label[order]
  if (sort) {
    index <- match(phy$tip.label, rownames(data))
    data <- as.matrix(data[index, ])
  }
  if (dm == 2) {
    data <- as.matrix(data)
  }
  phy$node.label = NULL
  return(list(phy = phy, data = data, nc= nc))
}

genus_spp_tree <- function(tree, genus.species.list){
  spec.id.list <- na.omit(unique(genus.species.list[, 2]))
  genus.id.list <- na.omit(unique(genus.species.list[, 1]))
  
  # attach names to lists
  spec.id.list <- setNames(spec.id.list, spec.id.list)
  genus.id.list <- setNames(genus.id.list, genus.id.list)
  
  
  ## creating a list with one species per genus
  spp.list <- data.frame(spp.insert = rep(NA, length(genus.id.list)))
  for (i in 1:length(genus.id.list)) {
    posit <- which((sub("_.*", "", spec.id.list)) == genus.id.list[i])
    spp.list[i,1] <- names(spec.id.list[posit[1]])
  }
  
  # cutting only the genus of interest
  genus.tree <- geiger::treedata(tree, genus.id.list)$phy
  
  # Changing the tip labels (genus to species level)
  for (i in 1:length(genus.id.list)) {
    genus.tree$tip.label[which(genus.tree$tip.label == genus.id.list[i])] <- spp.list[i,1]
  }
  
  # Insert the another species to the tree, as polytomies
  spp.insert <- treedata_modif(genus.tree, spec.id.list)$nc$data_not_tree
    for (i in 1:length(spp.insert)) {
    genus.tree <- phytools::add.species.to.genus(genus.tree, spp.insert[i])
  }
  
  #missing.spp <- geiger::name.check(genus.tree, spec.id.list)
  class(genus.tree) <- "phylo"
  
  ### return a species-level tree
  return(new.tree= genus.tree)
}
