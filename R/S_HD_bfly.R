
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

# converting in a data frame
df.hd <- data.frame(HD.values = c(res.HD.bfly$H.TD, res.HD.bfly$H.PD, res.HD.bfly$H.FD,
                                  res.HD.bfly$H.MPDi, res.HD.bfly$H.MPD, 
                                  res.HD.bfly$H.MFDi, res.HD.bfly$H.MFD), 
                    Div.names = rep(colnames(res.HD.bfly[, 1:7]), each = nrow(res.HD.bfly)),
                    Strata = rep(res.HD.bfly$Stratum, 7))
df.hd$Div.names1 <- factor(df.hd$Div.names, levels = c("H.TD", "H.PD", "H.FD",
                                                       "H.MPDi", "H.MFDi", "H.MPD", "H.MFD"))

library(ggplot2)
library(viridis)

p.hd <- ggplot(data = df.hd, aes(y = HD.values, x = Div.names1,
                                 colour = Strata, fill = Strata)) +
  geom_boxplot() + scale_color_viridis_d(option = "A") +
  scale_fill_viridis_d(option = "A", alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = 2, color = "firebrick1") +
  labs(x = "Diversity Measures", y = "Hidden Diversity") 
p.hd

cowplot::save_plot(here::here("output", "figures", "Fig1_hdbfly.png"), p.hd,
                   base_width = 8)

# testing if the HD differs between strata
mod.td <- lm(res.HD.bfly$H.TD ~ res.HD.bfly$Strata)
summary(mod.td)

mod.pd <- lm(res.HD.bfly$H.PD ~ res.HD.bfly$Strata)
summary(mod.pd)

mod.pd1 <- lm(HD.bfly$SES.PDobs$pd.obs.z ~ list_mean_sd$sesPD.est$mean[,6])
summary(mod.pd1)

plot(HD.bfly$SES.PDobs$pd.obs.z ~ list_mean_sd$sesPD.est$mean[,6])
abline(a = coef(mod.pd1)[1], b = coef(mod.pd1)[2])
abline(a = 0, b = 1, col = "red")

mod.fd <- lm(res.HD.bfly$H.FD ~ res.HD.bfly$Strata)
summary(mod.fd)

mod.mpdi <- lm(res.HD.bfly$H.MPDi ~ res.HD.bfly$Strata)
summary(mod.mpdi)

mod.mfdi <- lm(res.HD.bfly$H.MFDi ~ res.HD.bfly$Strata)
summary(mod.mfdi)

mod.mpd <- lm(res.HD.bfly$H.MPD ~ res.HD.bfly$Strata)
summary(mod.mpd)

mod.mfd <- lm(res.HD.bfly$H.MFD ~ res.HD.bfly$Strata)
summary(mod.mfd)

out.lm <- rbind(summary(mod.td)[[4]], summary(mod.pd)[[4]],
                summary(mod.fd)[[4]], summary(mod.mpdi)[[4]],
                summary(mod.mpd)[[4]], summary(mod.mfdi)[[4]],
                summary(mod.mfd)[[4]])
rownames(out.lm) <- rep(colnames(res.HD.bfly[, 1:7]), each = 2)

round(out.lm, 3)

# correlation between TD and SES.PD/FD
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



# visualization of the hidden diversity in another way --------------------

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
                    Div.name = rep(c("SES.TD","SES.PD", "SES.FD", "SES.MPDi", "SES.MFDi",
                                     "SES.MPD", "SES.MFD"), each = nrow(HD.bfly$SES.PDobs)), 
                    Strata = rep(substr(rownames(HD.bfly$SES.PDobs),1,1), 7))
df.HD$Strata <- ifelse(df.HD$Strata == "C", "Canopy", "Understory") 
df.HD$HD <- (df.HD$D.obs - df.HD$D.est)/df.HD$SD.est
df.HD$Div.name <- factor(df.HD$Div.name, levels = c("SES.TD","SES.PD", "SES.FD", "SES.MPDi", "SES.MFDi",
                                                    "SES.MPD", "SES.MFD"))

p.obs_est <- ggplot(df.HD, aes(x = D.est, y = D.obs, colour = Strata), alpha, .5, size = 4) +
  geom_point(alpha = .5) + scale_color_viridis_d(option = "A") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "lightblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "lightblue") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "firebrick1") +
  facet_wrap(~ Div.name, scales = "free") + theme(legend.position = "none")
p.obs_est

p.hd <- ggplot(df.HD, aes(x = Div.name, y = (HD*-1), colour = Strata, fill = Strata)) +
  geom_boxplot(alpha = .5) + facet_wrap(~ Div.name, nrow = 1, scales = "free_x") +
  scale_color_viridis_d(option = "A") + scale_fill_viridis_d(option = "A") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "firebrick1") +
  theme(legend.position = "none")

plot1 <- cowplot::plot_grid(p.obs_est, p.hd, ncol = 1)
cowplot::save_plot(here::here("output", "figures", "Figx_HD.png"),
                   plot1, base_width = 10, base_height = 6)



# Testing if the Estimated Diversity can predict the observed diversity --------

mod1 <- lm(HD.bfly$SES.PDobs$ntaxa ~ list_mean_sd$sesPD.est$mean[,1])
summary(mod1)
residuals(mod1)

mod2 <- lm(HD.bfly$SES.FDobs$pd.obs.z ~ list_mean_sd$sesFD.est$mean[,6])
summary(mod2)
residuals(mod2)
cor.test(residuals(mod2), res.HD.bfly$H.FD)

mod3 <- lm(HD.bfly$SES.MPDiobs$mpd.obs.z ~ list_mean_sd$sesMPDi.est$mean[,6])
summary(mod3)
cor.test(residuals(mod3), res.HD.bfly$H.MPDi)

mod4 <- lm(HD.bfly$SES.MPDobs$mpd.obs.z ~ list_mean_sd$sesMPD.est$mean[,6])
summary(mod4)
cor.test(residuals(mod4), res.HD.bfly$H.MPD)

mod5 <- lm(HD.bfly$SES.MFDiobs$mpd.obs.z ~ list_mean_sd$sesMFDi.est$mean[,6])
summary(mod5)
cor.test(residuals(mod5), res.HD.bfly$H.MFDi)

mod6 <- lm(HD.bfly$SES.MFDobs$mpd.obs.z ~ list_mean_sd$sesMFD.est$mean[,6])
summary(mod6)
cor.test(residuals(mod6), res.HD.bfly$H.MFD)
plot(residuals(mod6), res.HD.bfly$H.MFD)



obs.pd <- HD.bfly$SES.PDobs[,6]
est.pd <- list_mean_sd$sesPD.est$mean[,6]

sum(obs.pd > est.pd)
ts.o <- obs.pd[which(obs.pd < est.pd)]
ts.e <- est.pd[which(obs.pd < est.pd)]

plot(ts.e, ts.o)  
abline(a = 0, b = 1)
a = 0
b = 1
ts.e.pred <- (a + b*ts.e)
which(max(ts.o - ts.e.pred) == ts.o - ts.e.pred)

res.HD.bfly[which(rownames(res.HD.bfly) == "C2NF5D"),]
