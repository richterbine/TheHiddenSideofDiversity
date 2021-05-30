# testing if the HD differs between strata: Table C1 - Appendix C
# read files
HD.bfly <- readRDS(here::here("output", "HD_bfly.rds"))
list_mean_sd <- readRDS(here::here("output", "mean_estimated_diversity.rds"))

############
library(lme4)
library(lmerTest)

HD.bfly$SES.MPDobs

rownames(HD.bfly$SES.PDobs)
Strata <- substr(rownames(HD.bfly$SES.PDobs), 1,1)
Month <- substr(rownames(HD.bfly$SES.PDobs), 6, 6)
SU <- substr(rownames(HD.bfly$SES.PDobs), 3, 5)

mod.td <- lmer(SES.PDobs[,1] ~ Strata + (1|Month)  + (1|SU), data = HD.bfly)
summary(mod.td)

mod.pd <- lmer(SES.PDobs[,6] ~ Strata + (1|Month)  + (1|SU), data = HD.bfly)
summary(mod.pd)

mod.fd <- lmer(SES.FDobs[,6] ~ Strata + (1|Month)  + (1|SU), data = HD.bfly)
summary(mod.fd)

mod.mpdi <- lmer(SES.MPDiobs[,6] ~ Strata + (1|Month)  + (1|SU), data = HD.bfly)
summary(mod.mpdi)

mod.mfdi <- lmer(SES.MFDiobs[,6] ~ Strata + (1|Month)  + (1|SU), data = HD.bfly)
summary(mod.mfdi)

mod.mpd <- lmer(SES.MPDobs[,6] ~ Strata + (1|Month)  + (1|SU), data = HD.bfly)
summary(mod.mpd)

mod.mfd <- lmer(SES.MFDobs[,6] ~ Strata + (1|Month)  + (1|SU), data = HD.bfly)
summary(mod.mfd)

out.lmer <- rbind(summary(mod.td)$coefficients, summary(mod.pd)$coefficients,
                  summary(mod.fd)$coefficients, summary(mod.mpdi)$coefficients,
                  summary(mod.mpd)$coefficients, summary(mod.mfdi)$coefficients,
                  summary(mod.mfd)$coefficients)


rownames(out.lmer) <- rep(c("TD", "SES.PD", "SES.FD", "SES.MPDi",
                            "SES.MPD", "SES.MFDi", "SES.MFD"), each = 2)

names(HD.bfly)[7:12]

round(out.lmer, 3)
write.table(out.lmer, here::here("output", "mod_obs_lmer.txt"))


# model for estimated communities -----------------------------------------

list_mean_sd$sesMPDi.est$mean[,1]

mod.td <- lmer(sesPD.est$mean[,1] ~ Strata + (1|Month)  + (1|SU), data = list_mean_sd)
summary(mod.td)

mod.pd <- lmer(sesPD.est$mean[,6] ~ Strata + (1|Month)  + (1|SU), data = list_mean_sd)
summary(mod.pd)

mod.fd <- lmer(sesFD.est$mean[,6] ~ Strata + (1|Month)  + (1|SU), data = list_mean_sd)
summary(mod.fd)

mod.mpdi <- lmer(sesMPDi.est$mean[,6] ~ Strata + (1|Month)  + (1|SU), data = list_mean_sd)
summary(mod.mpdi)

mod.mfdi <- lmer(sesMFDi.est$mean[,6] ~ Strata + (1|Month)  + (1|SU), data = list_mean_sd)
summary(mod.mfdi)

mod.mpd <- lmer(sesMPD.est$mean[,6] ~ Strata + (1|Month)  + (1|SU), data = list_mean_sd)
summary(mod.mpd)

mod.mfd <- lmer(sesMFD.est$mean[,6] ~ Strata + (1|Month)  + (1|SU), data = list_mean_sd)
summary(mod.mfd)

out.lmer <- rbind(summary(mod.td)$coefficients, summary(mod.pd)$coefficients,
                  summary(mod.fd)$coefficients, summary(mod.mpdi)$coefficients,
                  summary(mod.mpd)$coefficients, summary(mod.mfdi)$coefficients,
                  summary(mod.mfd)$coefficients)


rownames(out.lmer) <- rep(c("TD", "SES.PD", "SES.FD", "SES.MPDi",
                            "SES.MPD", "SES.MFDi", "SES.MFD"), each = 2)

round(out.lmer, 3)
write.table(out.lmer, here::here("output", "mod_est_lmer.txt"))
