# evaluation of the P_RE model output with simulated data

library(ggplot2)
library(jagsUI)
#detach("package:jagsUI", unload = T)

# read the simulated communities and the output for bayesian model
# Performed in the script D_Model_Validation.R
sim.comm.list <- readRDS(here::here("R", "Model_validation", "Comm_simulated.rds"))
sim.mod <- readRDS(here::here("R", "Model_validation", "model_comm_simulated.rds"))

# Real parameters ---------------------------------------------------------

params <- c("mean.beta1", "sd.beta1", "mean.beta2", "sd.beta2", "mu.beta3", "sd.beta3",
            "mean.alpha1", "sd.alpha1", "mean.alpha2", "sd.alpha2",
            "mu.alpha3", "sd.alpha3", "mu.alpha4", "sd.alpha4")

tmp <- list(b1 = vector(mode = "list", length = 3),
            b2 = vector(mode = "list", length = 3),
            a1 = vector(mode = "list", length = 3),
            a2 = vector(mode = "list", length = 3))
names(tmp$b1) <- c ("A", "B", "C")
names(tmp$b2) <- c ("D", "E", "F")
names(tmp$a1) <- c ("G", "H", "I")
names(tmp$a2) <- c ("J", "K", "L")

for (i in 1:length(sim.comm.list)) {
  for (j in 1:length(sim.comm.list[[i]])) {
    tmp[[i]][[j]] <- unlist(sim.comm.list[[i]][[j]]$Parametres[params])
  }
}

df.real.mean.sd <- data.frame(values = unlist(tmp))

df.real.mean.sd <-  cbind(df.real.mean.sd, matrix(unlist(strsplit(rownames(df.real.mean.sd), "\\.", " ")),
                                                  ncol = 4, nrow = nrow(df.real.mean.sd), byrow = T))
colnames(df.real.mean.sd) <- c("values", "code","set.code", "p.type", "params")

df.real.mean.sd$p.trans<- round(ifelse(df.real.mean.sd$p.type == "mean" & df.real.mean.sd$params == "alpha1", qlogis(df.real.mean.sd$values), 
                                 ifelse(df.real.mean.sd$p.type == "mean" & df.real.mean.sd$params == "alpha2", qlogis(df.real.mean.sd$values),
                                        ifelse(df.real.mean.sd$p.type == "mean" & df.real.mean.sd$params == "beta1", log(df.real.mean.sd$values),
                                               ifelse(df.real.mean.sd$p.type == "mean" & df.real.mean.sd$params == "beta2", log(df.real.mean.sd$values), df.real.mean.sd$values)))),
                                3)

df.real.mean.sd$params.2 <- factor(df.real.mean.sd$params, labels = c(expression(alpha[1]),
                                                                      expression(alpha[2]),
                                                                      expression(alpha[3]),
                                                                      expression(alpha[4]),
                                                                      expression(beta[1]),
                                                                      expression(beta[2]),
                                                                      expression(beta[3])))


# Parameters estimated by the N-mixture model -----------------------------

params.est <- c("mu.beta.can", "sd.beta.can", "mu.beta.und", "sd.beta.und",
                "mu.beta1", "sd.beta1", "mu.alpha.can", "sd.alpha.can",
                "mu.alpha.und", "sd.alpha.und", "mu.alpha1", "sd.alpha1", 
                "mu.alpha2", "sd.alpha2")


tmp <- list(b1 = vector(mode = "list", length = 3),
            b2 = vector(mode = "list", length = 3),
            a1 = vector(mode = "list", length = 3),
            a2 = vector(mode = "list", length = 3))

names(tmp$b1) <- c("A", "B", "C") 
names(tmp$b2) <- c("D", "E", "F") 
names(tmp$a1) <- c("G", "H", "I") 
names(tmp$a2) <- c("J", "K", "L") 


for (i in 1:length(sim.mod)) {
  for (j in 1:length(sim.mod[[i]])) {
    tmp[[i]][[j]] <- as.data.frame(sim.mod[[i]][[j]]$sims.list[params.est])
  }
}

df.est.mean.sd <- data.frame(values = c(unlist(tmp[1]), unlist(tmp[2]), 
                                        unlist(tmp[3]), unlist(tmp[4])))

t <- strsplit(rownames(df.est.mean.sd), "\\.", "")
df.est.mean.sd$code <- unlist(lapply(t, function(x) x[1]))
df.est.mean.sd$set.code <- unlist(lapply(t, function(x) x[2]))
df.est.mean.sd$p.type <- unlist(lapply(t, function(x) x[3]))

t1 <- unlist(lapply(t, function(x) x[4]))
t1 <- gsub("beta1", "beta1.", t1)
t1 <- gsub("alpha1", "alpha1.", t1)
t1 <- gsub("alpha2", "alpha2.", t1)
t2 <- strsplit(t1, "\\.")
t2 <- unlist(lapply(t2, function(x) x[1]))

df.est.mean.sd$params <- as.factor(paste(t2, ".", 
                               gsub("[0-9]", "", unlist(lapply(t, function(x) x[5]))), sep = ""))
df.est.mean.sd$params.2 <- factor(df.est.mean.sd$params, labels = c(expression(alpha[1]),
                                                                    expression(alpha[2]),
                                                                    expression(alpha[3]),
                                                                    expression(alpha[4]),
                                                                    expression(beta[1]),
                                                                    expression(beta[2]),
                                                                    expression(beta[3])))


# Preparing for plot ------------------------------------------------------
# Model setting: A to F

r.set <- df.real.mean.sd[which(df.real.mean.sd[,"code"] == "b1" | df.real.mean.sd[,"code"] == "b2"),]
r.set.m <- r.set[which(r.set[,"p.type"] == "mu" | r.set[,"p.type"] == "mean"),]
r.set.sd <- r.set[which(r.set[,"p.type"] == "sd"),]

e.set <- df.est.mean.sd[which(df.est.mean.sd[, "code"] == "b1" | df.est.mean.sd[, "code"] == "b2"),]
e.set.m <- e.set[which(e.set["p.type"] == "mu"),]
e.set.sd <- e.set[which(e.set["p.type"] == "sd"),]

# plotting the means
p.A_F <- ggplot(e.set.m, aes(x = values, colour = set.code, fill = set.code)) + 
  geom_density(alpha = 0.5, position = "identity") + 
  facet_wrap(~ params.2 , labeller = label_parsed, nrow = 1) +
  geom_vline(data = r.set.m, aes(xintercept = p.trans, colour = set.code), 
             linetype = "dashed") + theme(legend.position = "none") +
  scale_color_viridis_d(option = "A", name = "Setting code") +
  scale_fill_viridis_d(option = "A", name = "Setting code") +
  theme(legend.position = "none") + scale_x_continuous(breaks = c(-.5, .5, 1.5)) +
  labs(x = "Estimated mean values" , y = "Density", tag = "a)")
p.A_F

# plotting the standard deviation
p.A_F.sd <- ggplot(e.set.sd, aes(x = values, colour = set.code, fill = set.code)) + 
  geom_density(alpha = 0.5, position = "identity") + 
  facet_wrap(~ params.2 , labeller = label_parsed, nrow = 1) +
  geom_vline(data = r.set.sd, aes(xintercept = p.trans, colour = set.code), 
             linetype = "dashed") +  theme(legend.position = "none") +
  scale_color_viridis_d(option = "A", name = "Model code") +
  scale_fill_viridis_d(option = "A", name = "Model code") +
  labs(x = "Estimated SD values" , y = "Density", tag = "b)")
p.A_F.sd

legAF <- cowplot::get_legend(p.A_F + theme(legend.position = "bottom"))

SetAF <- cowplot::plot_grid(p.A_F, p.A_F.sd, legAF,
                            ncol = 1,
                            rel_heights = c(2, 2, 1))
cowplot::save_plot(here::here("output", "figures", "FigS3_hyperparams.png"),
                   SetAF, base_height = 6, base_width = 8)

# Model setting: G to L

r.set <- df.real.mean.sd[which(df.real.mean.sd[,"code"] == "a1" | df.real.mean.sd[,"code"] == "a2"),]
r.set.m <- r.set[which(r.set[,"p.type"] == "mu" | r.set[,"p.type"] == "mean"),]
r.set.sd <- r.set[which(r.set[,"p.type"] == "sd"),]

e.set <- df.est.mean.sd[which(df.est.mean.sd[, "code"] == "a1" | df.est.mean.sd[, "code"] == "a2"),]
e.set.m <- e.set[which(e.set["p.type"] == "mu"),]
e.set.sd <- e.set[which(e.set["p.type"] == "sd"),]

# plotting the means
p.G_L <- ggplot(e.set.m, aes(x = values, colour = set.code, fill = set.code)) + 
  geom_density(alpha = 0.5, position = "identity") + 
  facet_wrap(~ params.2 , labeller = label_parsed, nrow = 1) +
  geom_vline(data = r.set.m, aes(xintercept = p.trans, colour = set.code), 
             linetype = "dashed") + theme(legend.position = "none") +
  scale_color_viridis_d(option = "A", name = "Setting code") +
  scale_fill_viridis_d(option = "A", name = "Setting code") +
  labs(x = "Estimated mean values" , y = "Density", tag = "a)")
p.G_L

p.G_L.sd <- ggplot(e.set.sd, aes(x = values, colour = set.code, fill = set.code)) + 
  geom_density(alpha = 0.5, position = "identity") + 
  facet_wrap(~ params.2 , labeller = label_parsed, nrow = 1) +
  geom_vline(data = r.set.sd, aes(xintercept = p.trans, colour = set.code), 
             linetype = "dashed") + theme(legend.position = "none") +
  scale_color_viridis_d(option = "A", name = "Setting code") +
  scale_fill_viridis_d(option = "A", name = "Setting code") +
  labs(x = "Estimated SD values" , y = "Density", tag = "b)")
p.G_L.sd


legGL <- cowplot::get_legend(p.G_L + theme(legend.position = "bottom"))

SetGL <- cowplot::plot_grid(p.G_L, p.G_L.sd, legGL,
                             ncol = 1, rel_heights = c(2, 2, 1))

cowplot::save_plot(here::here("output", "figures", "FigS4_hyperparams.png"), 
                   SetGL, base_height = 6, base_width = 8)


# Community mean parameters -----------------------------------------------

# Comparison among real parameters simulated and parameters estimated by model ####
# for real values
betas <- c("mu.beta1", "mu.beta2", "sd.beta1", "sd.beta2")
alphas <- c("mu.alpha1", "mu.alpha2", "sd.alpha1", "sd.alpha2")

real.samples <- list(beta1 = vector(mode = "list", length = 3),
                     beta2 = vector(mode = "list", length = 3),
                     alpha1 = vector(mode = "list", length = 3),
                     alpha2 = vector(mode = "list", length = 3))
names(real.samples$beta1) <- c("A", "B", "C") 
names(real.samples$beta2) <- c("D", "E", "F") 
names(real.samples$alpha1) <- c("G", "H", "I") 
names(real.samples$alpha2) <- c("J", "K", "L") 

sample.tmp <- matrix(NA, ncol = length(sim.comm.list), nrow = 1000)
colnames(sample.tmp) <- c("abundance.b1", "abundance.b2", "detection.a1", "detection.a2")

for (i in 1:length(sim.comm.list)) {
  for (j in 1:length(sim.comm.list[[i]])) {
    for (n in 1:2) {
      sample.tmp[, n] <- exp(rnorm(1000, mean = sim.comm.list[[i]][[j]]$Params.estimated[[betas[n]]], 
                                   sd = sim.comm.list[[i]][[j]]$Parametres[[betas[n+2]]]))
      
      sample.tmp[, n+2] <- plogis(rnorm(1000, mean = sim.comm.list[[i]][[j]]$Params.estimated[[alphas[n]]],
                                        sd = sim.comm.list[[i]][[j]]$Parametres[[alphas[n+2]]]))
      
    }
    real.samples[[i]][[j]] <- sample.tmp
  }
}
tmp <- as.data.frame(real.samples)

df.real.sample <- data.frame(values = c(unlist(tmp[, c(1, 5, 9, 13, 17, 21, 25, 29, 33, 37,41, 45)]), 
                                        unlist(tmp[, c(2, 6, 10, 14, 18, 22, 26, 30, 34, 38, 42, 46)]),
                                        unlist(tmp[, c(3, 7, 11, 15, 19, 23, 27, 31, 35, 39, 43, 47)]),
                                        unlist(tmp[, c(4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48)])))

t <- strsplit(rownames(df.real.sample), "\\.")
t1 <- unlist(lapply(t, function(x) x[4]))
df.real.sample$parameters <- as.factor(substr(t1, 1, 2))
df.real.sample$code <- as.factor(gsub("[a-z|0-9.]", "", rownames(df.real.sample)))

df.real.sample$v.beta <- ifelse(df.real.sample[, "parameters"] == "b1" & df.real.sample[,"code"] == "A", 0.5,
                                ifelse(df.real.sample[, "parameters"] == "b1" & df.real.sample[,"code"] == "B", 1.5,
                                       ifelse(df.real.sample[, "parameters"] == "b1" & df.real.sample[,"code"] == "C", 2.0,
                                              ifelse(df.real.sample[, "parameters"] == "b2" & df.real.sample[,"code"] == "D", 0.5,
                                                     ifelse(df.real.sample[, "parameters"] == "b2" & df.real.sample[,"code"] == "E", 1.5,
                                                            ifelse(df.real.sample[, "parameters"] == "b2" & df.real.sample[,"code"] == "F", 2.0, 1.0))))))

df.real.sample$v.alpha <- ifelse(df.real.sample[, "parameters"] == "a1" & df.real.sample[,"code"] == "G", 0.2,
                                 ifelse(df.real.sample[, "parameters"] == "a1" & df.real.sample[,"code"] == "H", 0.6,
                                        ifelse(df.real.sample[, "parameters"] == "a1" & df.real.sample[,"code"] == "I", 0.8,
                                               ifelse(df.real.sample[, "parameters"] == "a2" & df.real.sample[,"code"] == "J", 0.2,
                                                      ifelse(df.real.sample[, "parameters"] == "a2" & df.real.sample[,"code"] == "K", 0.6,
                                                             ifelse(df.real.sample[, "parameters"] == "a2" & df.real.sample[,"code"] == "L", 0.8, 0.5))))))
df.real.sample$parameters.2 <- factor(df.real.sample$parameters, labels = c(expression(alpha[1]), 
                                                                            expression(alpha[2]),
                                                                            expression(beta[1]),
                                                                            expression(beta[2])))
betas.r <- df.real.sample[which(df.real.sample[, "parameters"] == "b1" | df.real.sample[, "parameters"] == "b2"),]

n.real <- ggplot(betas.r, aes(x = values, color = code, 
                            fill = code)) + 
  geom_density(alpha = 0.5, position = "identity") +
  scale_color_viridis_d(option = "A") + scale_fill_viridis_d(option = "A") +
  geom_vline(aes(xintercept = v.beta, color = code),
             linetype = "dashed") + labs(tag = "a)", y = "Density",
                                         x = "True mean abundance") +
  facet_wrap(~ parameters.2, labeller = label_parsed) + #+
  theme(legend.position = "none")
n.real

alphas.r <- df.real.sample[which(df.real.sample[, "parameters"] == "a1" | df.real.sample[, "parameters"] == "a2"),]

p.real <- ggplot(alphas.r, aes(x = values, color = code, 
                             fill = code)) + 
  geom_density(alpha = 0.5, position = "identity") +
  geom_vline(aes(xintercept = v.alpha, color = code),
             linetype = "dashed") + labs(tag = "b)", x = "True mean detection probability",
                                         y = "Density") +
  facet_wrap(~ parameters.2, labeller = label_parsed) +
  scale_color_viridis_d(option = "A") + scale_fill_viridis_d(option = "A")+
  theme(legend.position = "none")
p.real


# estimated parameters by N-mixture model ---------------------------------

betas <- c("mu.beta.can", "mu.beta.und", "sd.beta.can", "sd.beta.und")
alphas <- c("mu.alpha.can", "mu.alpha.und", "sd.alpha.can", "sd.alpha.und")

estimated.samples <- list(beta1 = vector(mode = "list", length = 3),
                          beta2 = vector(mode = "list", length = 3),
                          alpha1 = vector(mode = "list", length = 3),
                          alpha2 = vector(mode = "list", length = 3))

names(estimated.samples$beta1) <- c("A", "B", "C") 
names(estimated.samples$beta2) <- c("D", "E", "F") 
names(estimated.samples$alpha1) <- c("G", "H", "I") 
names(estimated.samples$alpha2) <- c("J", "K", "L") 

sample.tmp <- matrix(NA, ncol = length(sim.mod), nrow = 1000)
colnames(sample.tmp) <- c("abundance.b1", "abundance.b2", "detection.a1", "detection.a2")

for (i in 1:length(sim.mod)) {
  for (j in 1:length(sim.mod[[i]])) {
    for (n in 1:2) {
      sample.tmp[, n] <- exp(rnorm(1000, mean = sim.mod[[i]][[j]][[2]][[betas[n]]], 
                                   sd = sim.mod[[i]][[j]][[2]][[betas[n+2]]]))
      
      sample.tmp[, n+2] <- plogis(rnorm(1000, mean = sim.mod[[i]][[j]][[2]][[alphas[n]]],
                                        sd = sim.mod[[i]][[j]][[2]][[alphas[n+2]]]))
      
    }
    estimated.samples[[i]][[j]] <- sample.tmp
  }
}
tmp <- as.data.frame(estimated.samples)

df.estimated.sample <- data.frame(values = c(unlist(tmp[, c(1, 5, 9, 13, 17, 21, 25, 29, 33, 37,41, 45)]), 
                                             unlist(tmp[, c(2, 6, 10, 14, 18, 22, 26, 30, 34, 38, 42, 46)]),
                                             unlist(tmp[, c(3, 7, 11, 15, 19, 23, 27, 31, 35, 39, 43, 47)]),
                                             unlist(tmp[, c(4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48)])))

t <- strsplit(rownames(df.estimated.sample), "\\.")
t1 <- unlist(lapply(t, function(x) x[4]))
df.estimated.sample$parameters <- as.factor(substr(t1, 1, 2))
df.estimated.sample$code <- as.factor(gsub("[a-z|0-9.]", "", rownames(df.estimated.sample)))

df.estimated.sample$v.beta <- ifelse(df.estimated.sample[, "parameters"] == "b1" & df.estimated.sample[,"code"] == "A", 0.5,
                                     ifelse(df.estimated.sample[, "parameters"] == "b1" & df.estimated.sample[,"code"] == "B", 1.5,
                                            ifelse(df.estimated.sample[, "parameters"] == "b1" & df.estimated.sample[,"code"] == "C", 2.0,
                                                   ifelse(df.estimated.sample[, "parameters"] == "b2" & df.estimated.sample[,"code"] == "D", 0.5,
                                                          ifelse(df.estimated.sample[, "parameters"] == "b2" & df.estimated.sample[,"code"] == "E", 1.5,
                                                                 ifelse(df.estimated.sample[, "parameters"] == "b2" & df.estimated.sample[,"code"] == "F", 2.0, 1.0))))))

df.estimated.sample$v.alpha <- ifelse(df.estimated.sample[, "parameters"] == "a1" & df.estimated.sample[,"code"] == "G", 0.2,
                                      ifelse(df.estimated.sample[, "parameters"] == "a1" & df.estimated.sample[,"code"] == "H", 0.6,
                                             ifelse(df.estimated.sample[, "parameters"] == "a1" & df.estimated.sample[,"code"] == "I", 0.8,
                                                    ifelse(df.estimated.sample[, "parameters"] == "a2" & df.estimated.sample[,"code"] == "J", 0.2,
                                                           ifelse(df.estimated.sample[, "parameters"] == "a2" & df.estimated.sample[,"code"] == "K", 0.6,
                                                                  ifelse(df.estimated.sample[, "parameters"] == "a2" & df.estimated.sample[,"code"] == "L", 0.8, 0.5))))))
df.estimated.sample$parameters.2 <- factor(df.estimated.sample$parameters, labels = c(expression(alpha[can]), expression(alpha[und]),
                                                                                      expression(beta[can]), expression(beta[und])))
# preparing for plot

betas.e <- df.estimated.sample[which(df.estimated.sample[, "parameters"] == "b1" | df.estimated.sample[, "parameters"] == "b2"),]

n.est <- ggplot(betas.e, aes(x = values, color = code, 
                           fill = code)) + 
  geom_density(alpha = 0.5, position = "identity") +
  geom_vline(aes(xintercept = v.beta, color = code),
             linetype = "dashed") + labs(tag = "c)", y = "Density",
                                         x = "Estimated mean abundance") +
  facet_wrap(~ parameters.2, labeller = label_parsed) + 
  scale_color_viridis_d(option = "A") + scale_fill_viridis_d(option = "A") +
  theme(legend.position = "none")
n.est

alphas.e <- df.estimated.sample[which(df.estimated.sample[, "parameters"] == "a1" | df.estimated.sample[, "parameters"] == "a2"),]

p.est <- ggplot(alphas.e, aes(x = values, color = code, 
                            fill = code)) + 
  geom_density(alpha = 0.5, position = "identity") +
  scale_color_viridis_d(option = "A") + scale_fill_viridis_d(option = "A") +
  geom_vline(aes(xintercept = v.alpha, color = code),
             linetype = "dashed") + labs(tag = "d)", y = "Density",
                                         x = "Estimated mean detection probability") +
  facet_wrap(~ parameters.2, labeller = label_parsed) +
  theme(legend.position = "none")
p.est

legs <- cowplot::get_legend(p.est + theme(legend.position = "bottom"))

p.joint <- cowplot::plot_grid(n.real, p.real,
                              n.est, p.est,
                              ncol=2)

Mean_Community_params <- cowplot::plot_grid(p.joint, legs,
                                            rel_heights = c(3, 1),
                                            ncol = 1)
Mean_Community_params

cowplot::save_plot(here::here("output", "figures", "FigS1_meanCparams.png"), 
                   Mean_Community_params,
                   base_height = 8)
