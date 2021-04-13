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

# Organizing the data by sampling day -------------------------------------

(species.list <- levels(data.bfly$Species))                                # alphabetic list
(spec.name.list <- tapply(data.bfly$Abundance, data.bfly$Species, sum))    # species ID 
(spec.id.list <- unique(data.bfly$Species))                                # ID list for unique species
(ordered.spec.name.list <- spec.name.list[order(spec.name.list)])          # ID-order list
data.bfly$julian.date <- as.numeric(format(data.bfly$Date, "%j"))          # data in the julian date format


# covariates for biological and sampling process --------------------------

## environmental data for biological process
cov.abn <- data.frame(Temperature = round(tapply(data.bfly$Temperature, data.bfly$Sites, mean, na.rm = T), 3),
                      Strata = c(rep(0, length(grep("C", unique(data.bfly$Sites)))), 
                                 rep(1, length(grep("U", unique(data.bfly$Sites)))))) 

cov.abn$Month <- as.numeric(as.factor(data.bfly[match(rownames(cov.abn), data.bfly$Sites), "Month"]))
cov.abn$SU <- as.numeric(data.bfly[match(rownames(cov.abn), data.bfly$Sites), "SU"])
str(cov.abn)

# scaling the temperature covariate
cov.abn$Temperature <- scale(cov.abn$Temperature)
head(cov.abn)


# Covariates for the observational process --------------------------------

## julian date and temperature for sampling process
cov.date <- melt(data.bfly, id.vars = c("Sites","Samp.day"), measure.vars = "julian.date")
cov.date <- cast(cov.date, Sites ~ Samp.day, mean)
cov.date[is.na(cov.date)] <- NA
summary(cov.date)

# scaling the temperature
scaled <- matrix(NA, nrow = nrow(cov.date), ncol = ncol(cov.date[,-1]))
for (i in 1:ncol(scaled)) {
  scaled[,i] <- as.numeric(cov.date[,i+1])
}

colnames(scaled) <- c("rep1", "rep2", "rep3", "rep4", "rep5")
scaled.dates <- scale(x = scaled) # scaling the dates covariate
scaled.dates[is.na(scaled.dates)] <- 0
str(scaled.dates)

# temperature covariate
cov.temp <- melt(data.bfly, id.vars = c("SU","Samp.day"), measure.vars = "Temperature")
cov.temp <- cast(cov.temp, SU ~ Samp.day, mean, na.rm = T)
cov.temp[is.na(cov.temp)] <- NA
cov.temp
scaled <- matrix(NA, nrow = nrow(cov.temp), ncol = ncol(cov.temp[,-1]))
for (i in 1:ncol(scaled)) {
  scaled[,i] <- as.numeric(cov.temp[,i+1])
}
colnames(scaled) <- c("rep1", "rep2", "rep3", "rep4", "rep5")
scaled.temp <- scale(x = scaled) # scaling the temperature covariate
scaled.temp[is.na(scaled.temp)] <- 0
head(scaled.temp)


# The data is reshaped into a three dimensional array yc where the first dimension
# i is the site; the second dimension, j, is the monthly repetitions;
# and the last dimension, k, is the species. 

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

# Loading the Nymphalidae tree (Chazot et al. 2019) and the
# functional traits matrix (Mean traits by species)

source(here::here("R", "functions", "genus_to_spp_function.R"))

# Manipulating the proposed phylogenetic hypothesis
Nym.tree <- ape::read.tree(here::here("data", "Raw", "Nymphalidae_tree.txt"))

# Cutting the Nym.tree to my specific sampling pool
tree.bfly <- genus_spp_tree(tree = Nym.tree, 
                            genus.species.list = data.bfly[, c("Genus", "Species")])
tree.bfly


traits.bfly <- read.table(here::here("data", "processed", "Mean_traits_bfly.txt"), header = T)
str(traits.bfly)

for (i in 1:ncol(traits.bfly)) {
  if(class(traits.bfly[,i]) == "integer"){
    traits.bfly[,i] <- as.factor(traits.bfly[,i])
  }
}


# Verifying if the detection probability is related with some functional trait
# for each stratum
# since the detection probability is truncated between 0 and 1, we need to employ
# a corrected model structure to fit the data
sim.det <- runif(n = dim(yc)[3])
det.can <- rbinom(n = dim(yc)[3], size = 1, prob = sim.det)
plot(det.can)
# it is necessary normalize and centralize the continuous predictor covariates
library(vegan)
std.traits<- vegan::decostand(traits.bfly[,1:7], MARGIN = 2, method = "standardize")
traits.bfly <- cbind(traits.bfly[, 8:12], std.traits)
str(traits.bfly)

# Canopy

mod.det.can <- glm(det.can ~., data = traits.bfly, family = binomial)

# removing multicolinear traits
traits.new <- traits.bfly
while (max(car::vif(mod.det.can)) > 3) {
  remove <- which(car::vif(mod.det.can) == max(car::vif(mod.det.can)))
  traits.new <- traits.new [, -remove]
  mod.det.can <- glm(det.can ~., data = traits.new, family = binomial)
}
summary(mod.det.can)
deviance(mod.det.can)/df.residual(mod.det.can)


mod.det.can <- glm(det.can ~., data = traits.bfly, family = quasibinomial)

# removing multicolinear traits
traits.new <- traits.bfly
while (max(car::vif(mod.det.can)) > 3) {
  remove <- which(car::vif(mod.det.can) == max(car::vif(mod.det.can)))
  traits.new <- traits.new [, -remove]
  mod.det.can <- glm(det.can ~., data = traits.new, family = quasibinomial)
}
obj <- summary(mod.det.can)
which(obj$coefficients[,4] < 0.5) 

host <- setNames(traits.bfly$Host, rownames(traits.bfly))
phytools::phylosig(tree.bfly, host, method = "K", test = T)

