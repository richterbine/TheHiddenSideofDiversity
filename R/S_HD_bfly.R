
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

list_mean_sd$sesPD.est$mean
list_mean_sd$sesMPD.est$mean
list_mean_sd$sesMPDi.est$mean
list_mean_sd$sesFD.est$mean
list_mean_sd$sesMFD.est$mean
list_mean_sd$sesMFDi.est$mean






HD.bfly$SES.PDobs[which(is.na(HD.bfly$SES.PDobs$pd.obs.z)), "pd.obs.z"] <- 0

HD.bfly$SES.FDobs[which(is.na(HD.bfly$SES.FDobs$pd.obs.z)), "pd.obs.z"] <- 0

HD.bfly$SES.MPDobs[which(is.na(HD.bfly$SES.MPDobs$mpd.obs.z)), "mpd.obs.z"] <- 0

HD.bfly$SES.MPDiest[which(is.na(HD.bfly$SES.MPDiest$mpd.obs.z)), "mpd.obs.z"] <- 0

HD.bfly$SES.MFDobs[which(is.na(HD.bfly$SES.MFDobs$mpd.obs.z)), "mpd.obs.z"] <- 0


res.HD.bfly <- data.frame(H.TD = (HD.bfly[[5]]$ntaxa - list_mean_sd$sesPD.est$mean[,"ntaxa"])/list_mean_sd$sesPD.est$sd[,"ntaxa"],
                          H.PD = (HD.bfly[[5]]$pd.obs.z - list_mean_sd$sesPD.est$mean[,"pd.obs.z"])/list_mean_sd$sesPD.est$sd[,"pd.obs.z"],
                          H.FD = (HD.bfly[[7]]$pd.obs.z - list_mean_sd$sesFD.est$mean[,"pd.obs.z"])/list_mean_sd$sesFD.est$sd[,"pd.obs.z"],
                          H.MPD = (HD.bfly[[6]]$mpd.obs.z - list_mean_sd$sesMPD.est$mean[,"mpd.obs.z"])/list_mean_sd$sesMPD.est$sd[,"mpd.obs.z"],
                          H.MFD = (HD.bfly[[8]]$mpd.obs.z - list_mean_sd$sesMFD.est$mean[,"mpd.obs.z"])/list_mean_sd$sesMFD.est$sd[,"mpd.obs.z"])
res.HD.bfly$Strata <- substr(rownames(res.HD.bfly), 1, 1)
res.HD.bfly$Stratum <- factor(res.HD.bfly$Strata, labels = c("Canopy", "Understory"))

df.hd <- data.frame(HD.values = c(res.HD.bfly$H.TD, res.HD.bfly$H.PD, res.HD.bfly$H.FD,
                                  res.HD.bfly$H.MPD, res.HD.bfly$H.MFD), 
                    Div.names = rep(colnames(res.HD.bfly[, 1:5]), each = nrow(res.HD.bfly),1),
                    Strata = rep(res.HD.bfly$Stratum), each = nrow(res.HD.bfly), 1)

library(ggplot2)
library(viridis)

p.hd <- ggplot(data = df.hd, aes(y = HD.values, x = Div.names,
                                 colour = Strata, fill = Strata)) +
  geom_boxplot() + scale_color_viridis_d(option = "A") +
  scale_fill_viridis_d(option = "A", alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = 2, color = "firebrick1")

p.hd
cowplot::save_plot(here::here("output", "figures", "Fig1_hdbfly.png"), p.hd,
                   base_width = 8)

# testing if the HD differs between strata
mod.td <- lm(res.HD.bfly$H.TD ~ res.HD.bfly$Strata)
summary(mod.td)

mod.pd <- lm(res.HD.bfly$H.PD ~ res.HD.bfly$Strata)
summary(mod.pd)

mod.fd <- lm(res.HD.bfly$H.FD ~ res.HD.bfly$Strata)
summary(mod.fd)

mod.mpd <- lm(res.HD.bfly$H.MPD ~ res.HD.bfly$Strata)
summary(mod.mpd)

mod.mfd <- lm(res.HD.bfly$H.MFD ~ res.HD.bfly$Strata)
summary(mod.mfd)

out.lm <- rbind(summary(mod.td)[[4]], summary(mod.pd)[[4]],
                summary(mod.fd)[[4]], summary(mod.mpd)[[4]],
                summary(mod.mfd)[[4]])


# correlation among td and pd/fd
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



HD.bfly <- readRDS(here::here("output", "HD_bfly.rds"))

bfly.obs.mpd <- HD.bfly$SES.MPDobs

bfly.obs.mpd$strata <- substr(rownames(bfly.obs.mpd), 1, 1)
bfly.obs.mpd <- na.omit(bfly.obs.mpd)

summary(lm(mpd.obs.z ~ strata, data = bfly.obs.mpd))
boxplot(bfly.obs.mpd$mpd.obs.z ~ bfly.obs.mpd$strata)

colnames(HD.bfly$SES.MPDobs) <- colnames(HD.bfly$SES.MFDobs) <- colnames(HD.bfly$SES.PDobs)

df.obs.div <- rbind(HD.bfly$SES.PDobs, HD.bfly$SES.FDobs, 
                    HD.bfly$SES.MPDobs, HD.bfly$SES.MFDobs)

df.obs.div$strata <- substr(rownames(df.obs.div), 1, 1)
df.obs.div$div <- rep(names(HD.bfly)[5:8], each = nrow(HD.bfly$SES.PDobs))

p.div.obd <- ggplot(data = na.omit(df.obs.div), aes(y = pd.obs.z, x = div,
                                                    colour = strata, fill = strata)) +
  geom_boxplot() + scale_color_viridis_d(option = "A") +
  scale_fill_viridis_d(option = "A", alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = 2, color = "firebrick1")
p.div.obd

cowplot::plot_grid(p.hd, p.div.obd)


length(na.omit(HD.bfly[[5]]$pd.obs.z - list_mean_sd$sesPD.est$mean[,"pd.obs.z"]))

a = rnorm(100)
b = rnorm(100)
ts <- as.data.frame(cbind(a, b, a-b, a < b, a < 0, b < 0))
ts$stat <- ifelse(ts$V5 == 0 & ts$V6 == 0, "pos-pos",
                  ifelse(ts$V5 == 0 & ts$V6 == 1, "pos-neg",
                         ifelse(ts$V5 == 1 & ts$V6 == 0, "neg-pos", 
                                ifelse(ts$V5 == 1 & ts$V6 == 1, "neg-neg", NA))))
ts$hd <- factor(ts$V4, labels = c("HD.pos", "HD.neg"))


ggplot(ts, aes(x = a, y = b, col = stat, shape = hd), size = 3) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "lightblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "lightblue") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "firebrick1")
ts

df.hd <- data.frame(Dobs = c(HD.bfly$SES.PDobs$pd.obs.z, HD.bfly$SES.MPDiobs$mpd.obs.z,
                             HD.bfly$SES.MPDobs$mpd.obs.z, HD.bfly$SES.FDobs$pd.obs.z, 
                             HD.bfly$SES.MFDiobs$mpd.obs.z, HD.bfly$SES.MFDobs$mpd.obs.z), 
                    Dest = c(list_mean_sd$sesPD.est$mean[,"pd.obs.z"], 
                             list_mean_sd$sesMPDi.est$mean[,"mpd.obs.z"],
                             list_mean_sd$sesMPD.est$mean[,"mpd.obs.z"], 
                             list_mean_sd$sesFD.est$mean[,"pd.obs.z"],
                             list_mean_sd$sesMFDi.est$mean[,"mpd.obs.z"],
                             list_mean_sd$sesMFD.est$mean[,"mpd.obs.z"]),
                    Type = rep(c("SES.PD", "SES.MPDi", "SES.MPD", "SES.FD", "SES.MFDi", "SES.MFD"), 
                               each = nrow(HD.bfly$SES.PDobs)),
                    Strata = rep(c("Canopy", "Undestory"), each = 150))


ggplot(df.hd, aes(x = Dobs, y = Dest, colour = Type), alpha, .5, size = 4) +
  geom_point() + scale_color_viridis_d(option = "A") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "lightblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "lightblue") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "firebrick1") +
  facet_wrap(~Strata)

ggplot(df.hd, aes(x = Dobs, y = Dest, colour = Strata), alpha, .5, size = 4) +
  geom_point() + scale_color_viridis_d(option = "A") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "lightblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "lightblue") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "firebrick1") +
  facet_wrap(~Type)


