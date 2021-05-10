################################################################################################
################## Model based in 5 moths of sampling in summmer of 2016/2017 ##################
############# Model reformulated of Kery and Royle 2016 (COMMUNITY N-MIXTURE (DRY) #############
################################################################################################

# improvements of the bayesian model --------------------------------------
# correction of priors for detection parameters, inclusion of sampling days 
# in the model (data matrix - StrataXtrapXday) and random effects for month and area

#### packages ####
library(jagsUI)
library(reshape)

##### data ####
data <- read.table(here::here("data", "processed", "Raw_data_bfly.txt"), header = TRUE) # data of butterflies for bayesian model
data$Date<- as.Date(as.vector(data$Date), "%m/%d/%y")
for (i in 1:ncol(data)) {
  if(class(data[,i]) == "character"){
    data[,i]<- as.factor(data[,i])
  }
}
str(data)

temps <- read.table(here::here("data", "processed", "raw_temp_bfly.txt"), header = T)
temps$Date <- as.Date(as.vector(temps$Date), "%m/%d/%y")
for (i in 1:ncol(temps)) {
  if(class(temps[,i]) == "character"){
    temps[,i]<- as.factor(temps[,i])
  }
}
temps$Month <- format(temps$Date, "%B")
pos <- which(temps$Month == "janeiro")
temps[pos, "Month"]<- "January"
pos <- which(temps$Month == "fevereiro")
temps[pos, "Month"]<- "February"
pos <- which(temps$Month == "março")
temps[pos, "Month"]<- "March"
pos <- which(temps$Month == "novembro")
temps[pos, "Month"]<- "November"
pos <- which(temps$Month == "dezembro")
temps[pos, "Month"]<- "December"

temps$Month <- as.factor(temps$Month)

# merging the data frames
tmp <-  merge(data, temps, by = "Code", all.y = T)
colnames(tmp)

data.bfly <- tmp[,c(1:5, 7, 13:19)]
colnames(data.bfly)

colnames(data.bfly) <- c("Code", "Subfamily", "Tribe", "Genus", "Species", "Abundance", 
                    "Strata", "Date", "SU", "Trap", "Temperature", "Samp.day", "Month")

pos <- which(is.na(data.bfly$Abundance) == TRUE)
data.bfly[pos, "Abundance"]<- 0

for (i in 1:nrow(data.bfly)) {
  tmp <- unlist(strsplit(as.character(data.bfly[i, "Code"]), "_"))
  data.bfly[i, "Sites"] <- paste(tmp[1], gsub("[a-z]","", as.character(data.bfly$Month[i])), 
                            sep = "")
}
write.table(data.bfly, here::here("output", "joint_bfly_data.txt"))


# Loading the Nymphalidae tree (Chazot et al. 2019) and the
# functional traits matrix (Mean traits by species)

source(here::here("R", "functions", "genus_to_spp_function.R"))

# Manipulating the proposed phylogenetic hypothesis
Nym.tree <- ape::read.tree(here::here("data", "Raw", "Nymphalidae_tree.txt"))

# Cutting the Nym.tree to my specific sampling pool
tree.bfly <- genus_spp_tree(tree = Nym.tree, 
                            genus.species.list = data.bfly[, c("Genus", "Species")])
tree.bfly
ape::write.tree(tree.bfly, here::here("output", "tree_bfly_FLONA.txt"))
