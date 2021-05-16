# visual evaluation of hidden diversity (HD) to multiple dimension of biodiversity --------

### calculating the hidden diversity
HD.bfly <- readRDS(here::here("output", "HD_bfly.rds"))

# calculate the mean and standar deviation for estimated communities
list_mean_sd <- list(sesPD.est = vector(mode = "list", length =2),
                     sesMPD.est = vector(mode = "list", length = 2),
                     sesMPDi.est = vector(mode = "list", length = 2),
                     sesFD.est = vector(mode = "list", length = 2),
                     sesMFD.est = vector(mode = "list", length = 2),
                     sesMFDi.est = vector(mode = "list", length = 2))
names(list_mean_sd$sesPD.est) <- c("mean", "sd")
names(list_mean_sd$sesMPD.est) <- c("mean", "sd")
names(list_mean_sd$sesFD.est) <- c("mean", "sd")
names(list_mean_sd$sesMFD.est) <- c("mean", "sd")
names(list_mean_sd$sesMFDi.est) <- c("mean", "sd")
names(list_mean_sd$sesMPDi.est) <- c("mean", "sd")

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

# given low values to NAs for observed diversity
HD.bfly$SES.PDobs[which(is.na(HD.bfly$SES.PDobs$pd.obs.z)), "pd.obs.z"] <- 1e-05

HD.bfly$SES.FDobs[which(is.na(HD.bfly$SES.FDobs$pd.obs.z)), "pd.obs.z"] <- 1e-05

HD.bfly$SES.MPDobs[which(is.na(HD.bfly$SES.MPDobs$mpd.obs.z)), "mpd.obs.z"] <- 1e-05

HD.bfly$SES.MPDiobs[which(is.na(HD.bfly$SES.MPDiobs$mpd.obs.z)), "mpd.obs.z"] <- 1e-05

HD.bfly$SES.MFDobs[which(is.na(HD.bfly$SES.MFDobs$mpd.obs.z)), "mpd.obs.z"] <- 1e-05

HD.bfly$SES.MFDiobs[which(is.na(HD.bfly$SES.MFDiobs$mpd.obs.z)), "mpd.obs.z"] <- 1e-05

# calculating the hidden diversity (d.obs - d.est)/sd.est

res.HD.bfly <- data.frame(H.TD = (HD.bfly$SES.PDobs$ntaxa - list_mean_sd$sesPD.est$mean[,"ntaxa"])/list_mean_sd$sesPD.est$sd[,"ntaxa"],
                          H.PD = (HD.bfly$SES.PDobs$pd.obs.z - list_mean_sd$sesPD.est$mean[,"pd.obs.z"])/list_mean_sd$sesPD.est$sd[,"pd.obs.z"],
                          H.FD = (HD.bfly$SES.FDobs$pd.obs.z - list_mean_sd$sesFD.est$mean[,"pd.obs.z"])/list_mean_sd$sesFD.est$sd[,"pd.obs.z"],
                          H.MPDi = (HD.bfly$SES.MPDiobs$mpd.obs.z - list_mean_sd$sesMPDi.est$mean[,"mpd.obs.z"])/list_mean_sd$sesMPDi.est$sd[,"mpd.obs.z"],
                          H.MPD = (HD.bfly$SES.MPDobs$mpd.obs.z - list_mean_sd$sesMPD.est$mean[,"mpd.obs.z"])/list_mean_sd$sesMPD.est$sd[,"mpd.obs.z"],
                          H.MFDi = (HD.bfly$SES.MFDiobs$mpd.obs.z - list_mean_sd$sesMFDi.est$mean[,"mpd.obs.z"])/list_mean_sd$sesMFDi.est$sd[,"mpd.obs.z"],
                          H.MFD = (HD.bfly$SES.MFDobs$mpd.obs.z - list_mean_sd$sesMFD.est$mean[,"mpd.obs.z"])/list_mean_sd$sesMFD.est$sd[,"mpd.obs.z"])
res.HD.bfly$Strata <- ifelse(substr(rownames(res.HD.bfly), 1, 1) == "C", "Canopy", 
                              "Understory")
res.HD.bfly$Month <- ifelse(substr(rownames(res.HD.bfly), 6, 6) == "N", "November", 
                             ifelse(substr(rownames(res.HD.bfly), 6, 6) == "D", "December",
                                    ifelse(substr(rownames(res.HD.bfly), 6, 6) == "J", "January",
                                           ifelse(substr(rownames(res.HD.bfly), 6, 6) == "F", "February", "March"))))
res.HD.bfly$SU <- substr(rownames(res.HD.bfly), 5, 5) 

# converting in a data frame
df.hd <- data.frame(HD.values = c(res.HD.bfly$H.TD, res.HD.bfly$H.PD, res.HD.bfly$H.FD,
                                  res.HD.bfly$H.MPDi, res.HD.bfly$H.MPD, 
                                  res.HD.bfly$H.MFDi, res.HD.bfly$H.MFD), 
                    Div.names = rep(colnames(res.HD.bfly[, 1:7]), each = nrow(res.HD.bfly)),
                    Strata = rep(res.HD.bfly$Strata, 7))
df.hd$Div.names1 <- factor(df.hd$Div.names, levels = c("H.TD", "H.PD", "H.FD",
                                                       "H.MPDi", "H.MFDi", "H.MPD", "H.MFD"))

library(ggplot2)
library(viridis)

p.hd <- ggplot(data = df.hd, aes(y = HD.values, x = Div.names1,
                                 colour = Strata, fill = Strata)) +
  geom_boxplot() + scale_color_viridis_d(option = "A") +
  scale_fill_viridis_d(option = "A", alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = 2, color = "firebrick1") +
  labs(x = "Diversity Measures", y = "Hidden Diversity") +
  facet_wrap(~Div.names1)
p.hd

# testing if the HD differs between strata
library(lme4)
library(lmerTest)

mod.td <- lmer(H.TD ~ Strata + (1|Month)  + (1|SU), data = res.HD.bfly)
summary(mod.td)

mod.pd <- lmer(H.PD ~ Strata + (1|Month)  + (1|SU), data = res.HD.bfly)
summary(mod.pd)

mod.fd <- lmer(H.FD ~ Strata + (1|Month)  + (1|SU), data = res.HD.bfly)
summary(mod.fd)

mod.mpdi <- lmer(H.MPDi ~ Strata + (1|Month)  + (1|SU), data = res.HD.bfly)
summary(mod.mpdi)

mod.mfdi <- lmer(H.MFDi ~ Strata + (1|Month)  + (1|SU), data = res.HD.bfly)
summary(mod.mfdi)

mod.mpd <- lmer(H.MPD ~ Strata + (1|Month)  + (1|SU), data = res.HD.bfly)
summary(mod.mpd)

mod.mfd <- lmer(H.MFD ~ Strata + (1|Month)  + (1|SU), data = res.HD.bfly)
summary(mod.mfd)

out.lmer <- rbind(summary(mod.td)$coefficients, summary(mod.pd)$coefficients,
                  summary(mod.fd)$coefficients, summary(mod.mpdi)$coefficients,
                  summary(mod.mpd)$coefficients, summary(mod.mfdi)$coefficients,
                  summary(mod.mfd)$coefficients)


rownames(out.lmer) <- rep(colnames(res.HD.bfly[, 1:7]), each = 2)

round(out.lmer, 3)
write.table(out.lmer, here::here("output", "mod_lmer.txt"))


# visualization of the hidden diversity in another way --------------------
# hidden diversity framework #

df.HD <- data.frame(D.obs = c(HD.bfly$SES.PDobs$ntaxa, HD.bfly$SES.PDobs$pd.obs.z, HD.bfly$SES.FDobs$pd.obs.z,
                              HD.bfly$SES.MPDiobs$mpd.obs.z, HD.bfly$SES.MFDiobs$mpd.obs.z,
                              HD.bfly$SES.MPDobs$mpd.obs.z, HD.bfly$SES.MFDobs$mpd.obs.z),
                    D.est = c(list_mean_sd$sesPD.est$mean[,1], list_mean_sd$sesPD.est$mean[,6],
                              list_mean_sd$sesFD.est$mean[,6], list_mean_sd$sesMPDi.est$mean[,6],
                              list_mean_sd$sesMFDi.est$mean[,6], list_mean_sd$sesMPD.est$mean[,6],
                              list_mean_sd$sesMFD.est$mean[,6]), 
                    SD.est = c(list_mean_sd$sesPD.est$sd[,1], list_mean_sd$sesPD.est$sd[,6],
                               list_mean_sd$sesFD.est$sd[,6], list_mean_sd$sesMPDi.est$sd[,6],
                               list_mean_sd$sesMFDi.est$sd[,6], list_mean_sd$sesMPD.est$sd[,6],
                               list_mean_sd$sesMFD.est$sd[,6]),
                    Div.name = rep(c("TD","SES.PD", "SES.FD", "SES.MPDi", "SES.MFDi",
                                     "SES.MPD", "SES.MFD"), each = nrow(HD.bfly$SES.PDobs)), 
                    Strata = rep(substr(rownames(HD.bfly$SES.PDobs),1,1), 7))
df.HD$Strata <- ifelse(df.HD$Strata == "C", "Canopy", "Understory") 
df.HD$HD <- (df.HD$D.obs - df.HD$D.est)/df.HD$SD.est
df.HD$Div.name <- factor(df.HD$Div.name, levels = c("TD","SES.PD", "SES.MPDi", "SES.MPD", "SES.FD", 
                                                    "SES.MFDi", "SES.MFD"))

p.obs_est <- ggplot(df.HD, aes(x = D.obs, y = D.est, colour = Strata), alpha, .5, size = 4) +
  geom_point(alpha = .5) + scale_color_viridis_d(option = "A") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "lightblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "lightblue") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "firebrick1") +
  facet_wrap(~ Div.name, nrow = 1) + theme(legend.position = "none") +
  labs(y = "Estimated Diversity", x = "Observed Diversity", tag = "b)")
p.obs_est

library(ggplot2)
colnames(df.HD)
#HD.bfly <- readRDS(here::here("output", "HD_bfly.rds"))

p.div.pat <- ggplot(df.HD, aes(x = Div.name, y = D.obs, colour = Strata, fill = Strata)) +
  geom_boxplot(alpha = 0.5) + ylim(-4, 12) + scale_color_viridis_d(option = "A") +
  scale_fill_viridis_d(option = "A") + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "firebrick1") +
  labs(y = "Observed diversity", x = " ", tag = "a)") +
  theme(legend.position = "none")
p.div.pat

p.div.pat1 <- ggplot(df.HD, aes(x = Div.name, y = D.est, colour = Strata, fill = Strata)) +
  geom_boxplot(alpha = 0.5) + ylim(-4, 12) + scale_color_viridis_d(option = "A") +
  scale_fill_viridis_d(option = "A") + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "firebrick1") +
  labs(y = "Estimated diversity", x = " ", tag = "b)") +
  theme(legend.position = "none")
p.div.pat1

cowplot::save_plot(here::here("output","figures", "Fig.S3_obs_est.png"),
                   cowplot::plot_grid(p.div.pat, p.div.pat1),
                   base_width = 10)


p.hd <- ggplot(df.HD, aes(x = Div.name, y = HD, colour = Strata, fill = Strata)) +
  geom_boxplot(alpha = .5) + facet_wrap(~ Div.name, nrow = 1, scales = "free_x") +
  scale_color_viridis_d(option = "A") + scale_fill_viridis_d(option = "A") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "firebrick1") +
  theme(legend.position = "none") + labs(y = "Hidden Diversity", tag = "a)")
p.hd

plot1 <- cowplot::plot_grid(p.hd, p.obs_est,  ncol = 1)
plot1

cowplot::save_plot(here::here("output", "figures", "Fig3_HD.png"),
                   plot1, base_width = 10, base_height = 6)