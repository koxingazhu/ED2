##########################################################################################################
library(rjags)
library(rlang)
library(R2jags)
library(car)
library(xtable)
library(broom)
library(ggplot2)
library(foreign)
library(reshape2)
library(devtools)
# install_github(repo = "johnbaums/jagstools")
library(jagstools)
library(gtools)
library(dplyr)
library(ggmcmc)
library(MASS)
library(lattice)
library(psych)
library(plyr)
library(readr)
library(lavaan)
library(RColorBrewer)
library(Metrics)
library(GGally)
library(wesanderson)
# install.packages('R2jags', dependencies = TRUE)
#
# library(remotes)
# remotes::install_github("communityecologist/densize")
library(densize)


setwd("E:\\")
#rm(list = ls())

################################################################################################################
# AGBS
################################################################################################################

AGB.data <- read.csv("E:\\data.csv")
names(AGB.data)
AGB.data <- AGB.data[AGB.data$Year>2050,]
plot(AGB.data$Density, AGB.data$Size)
hist(AGB.data$AGBS)
hist(AGB.data$AGBP)
plot(AGB.data$Year, AGB.data$AGBP)
AGB.data$lnAGBS <- log(AGB.data$AGBS)
AGB.data$lnAGBP <- log(AGB.data$AGBP)
hist(AGB.data$lnAGBS)
hist(AGB.data$lnAGBP)
#
DBH.SD.Plot <- tapply(AGB.data$DBH.SD, AGB.data$Plot, mean)
lnAGBP.Plot <- tapply(AGB.data$lnAGBP, AGB.data$Plot, mean)
cor.test(DBH.SD.Plot, lnAGBP.Plot)
#
# All
AGB.data$Density <- as.vector(scale(AGB.data$Density))
AGB.data$Size <- as.vector(scale(AGB.data$Size))
AGB.data$Simpson <- as.vector(scale(AGB.data$Simpson))
AGB.data$Shannon.Wiener <- as.vector(scale(AGB.data$Shannon.Wiener))
AGB.data$Pielou <- as.vector(scale(AGB.data$Pielou))
AGB.data$DBH.SD <- as.vector(scale(AGB.data$DBH.SD))
AGB.data$DBH.CV <- as.vector(scale(AGB.data$DBH.CV))
AGB.data$DBH.Gini <- as.vector(scale(AGB.data$DBH.Gini))
hist(AGB.data$lnAGBS)
summary(AGB.data)
dim(AGB.data)
AGB.data <- na.omit(AGB.data)
dim(AGB.data)
names(AGB.data)
plot(AGB.data$Simpson, AGB.data$AGBS)
plot(AGB.data$Simpson, AGB.data$lnAGBS)

#
plot(AGB.data$Size, AGB.data$Density)
#
AGB.data <- AGB.data[order(AGB.data$Year), ]
AGB.data$statecdid <- as.numeric(as.ordered(AGB.data$Year))
AGB.data.datjags <- as.list(AGB.data[, c("lnAGBS", "lnAGBP", "Simpson", "Shannon.Wiener", "Pielou",
                                             "AGBP", "AGBS", "Density", "Size", "DBH.SD", "DBH.CV", "DBH.Gini", "Year", "statecdid")])
AGB.data.datjags$YEAR <- as.vector(by(AGB.data.datjags$Year, AGB.data.datjags$statecdid, mean))
AGB.data.datjags$SIZE <- as.vector(by(AGB.data.datjags$Size, AGB.data.datjags$statecdid, mean))
AGB.data.datjags$DENSITY <- as.vector(by(AGB.data.datjags$Density, AGB.data.datjags$statecdid, mean))
AGB.data.datjags$N.resp <- 1750
AGB.data.datjags$N.unitcd <- length(unique(AGB.data.datjags$statecdid))
#
plot(AGB.data.datjags$YEAR, AGB.data.datjags$SIZE)
plot(AGB.data.datjags$YEAR, AGB.data.datjags$DENSITY)
plot(AGB.data.datjags$DENSITY, AGB.data.datjags$SIZE)
#
AGB.data.datjags.Y.D.S <- cbind(AGB.data.datjags$YEAR, AGB.data.datjags$DENSITY, AGB.data.datjags$SIZE)
as.data.frame(AGB.data.datjags.Y.D.S) -> df.AGB.data.datjags.Y.D.S
write.csv(df.AGB.data.datjags.Y.D.S, "E:\\BHM.Results\\...csv")
#
# Structural diversity
# 
##########################################################################################################
BUGS.DENSITY <- function(){
  for (i in 1:N.resp) {
    lnAGBP[i] ~ dnorm(mu[i], tau)
    mu[i] <- b.Simpson[statecdid[i]] * Simpson[i] + b.statecd[statecdid[i]]
  }
  
  for (j in 1:N.unitcd) {
    b.statecd[j] ~ dnorm(b.statecd.hat[j], tau.statecd)
    b.statecd.hat[j] <- b.statecd0 + b.DENSITY0 * DENSITY[j]
    b.Simpson[j] ~ dnorm(b.Simpson.hat[j], tau.Simpson)
    b.Simpson.hat[j] <- b.statecd1 + b.DENSITY1 * DENSITY[j]
  }
  #
  b.statecd0 ~ dnorm(0, 0.01)
  b.statecd1 ~ dnorm(0, 0.01)
  #
  b.DENSITY0 ~ dnorm(0, 0.01)
  b.DENSITY1 ~ dnorm(0, 0.01)
  #
  tau ~ dgamma(1, 1)
  sigma.statecd <- 1 / tau.statecd
  tau.statecd ~ dgamma(1, 1)
  tau.Simpson ~ dgamma(1,1)
}
#
BUGS.params <- c("b.DENSITY0", "b.DENSITY1",
                 "b.statecd0", "b.statecd1", "b.Simpson")
#
set.seed(123)
JAGS.Density.1 <- jags(data = AGB.data.datjags, inits = NULL, parameters.to.save = BUGS.params,
                       model.file = BUGS.DENSITY, n.chains = 3, n.iter = 20000, n.burnin = 5000)
#
options(max.print = 1000000)
JAGS.DENSITY.2 <- update(JAGS.Density.1, n.iter = 1000)
JAGS.DENSITY.2.DIC <- JAGS.DENSITY.2$BUGSoutput$DIC
JAGS.DENSITY.2.DIC
JAGS.DENSITY.2.mcmc <- as.mcmc(JAGS.DENSITY.2)
JAGS.DENSITY.AGBP <- JKmcmctab(JAGS.DENSITY.2.mcmc)
# Gelman, A and Rubin, DB (1992) Inference from iterative simulation using multiple sequences, Statistical Science, 7, 457-511.
JAGS.DENSITY.AGBP.diag <- gelman.diag(JAGS.DENSITY.2.mcmc)
JAGS.DENSITY.AGBP <- cbind(JAGS.DENSITY.AGBP, JAGS.DENSITY.AGBP.diag$psrf)
JAGS.DENSITY.AGBP
write.csv(JAGS.DENSITY.AGBP, "E:\\BHM.Results\\...csv")


#
BUGS.DENSITY <- function(){
  for (i in 1:N.resp) {
    lnAGBP[i] ~ dnorm(mu[i], tau)
    mu[i] <- b.Shannon.Wiener[statecdid[i]] * Shannon.Wiener[i] + b.statecd[statecdid[i]]
  }
  
  for (j in 1:N.unitcd) {
    b.statecd[j] ~ dnorm(b.statecd.hat[j], tau.statecd)
    b.statecd.hat[j] <- b.statecd0 + b.DENSITY0 * DENSITY[j]
    b.Shannon.Wiener[j] ~ dnorm(b.Shannon.Wiener.hat[j], tau.Shannon.Wiener)
    b.Shannon.Wiener.hat[j] <- b.statecd1 + b.DENSITY1 * DENSITY[j]
  }
  #
  b.statecd0 ~ dnorm(0, 0.01)
  b.statecd1 ~ dnorm(0, 0.01)
  #
  b.DENSITY0 ~ dnorm(0, 0.01)
  b.DENSITY1 ~ dnorm(0, 0.01)
  #
  tau ~ dgamma(1, 1)
  sigma.statecd <- 1 / tau.statecd
  tau.statecd ~ dgamma(1, 1)
  tau.Shannon.Wiener ~ dgamma(1,1)
}
#
BUGS.params <- c("b.DENSITY0", "b.DENSITY1",
                 "b.statecd0", "b.statecd1", "b.Shannon.Wiener")
#
set.seed(123)
JAGS.Density.1 <- jags(data = AGB.data.datjags, inits = NULL, parameters.to.save = BUGS.params,
                       model.file = BUGS.DENSITY, n.chains = 3, n.iter = 20000, n.burnin = 5000)
#
options(max.print = 1000000)
JAGS.DENSITY.2 <- update(JAGS.Density.1, n.iter = 1000)
JAGS.DENSITY.2.DIC <- JAGS.DENSITY.2$BUGSoutput$DIC
JAGS.DENSITY.2.DIC
JAGS.DENSITY.2.mcmc <- as.mcmc(JAGS.DENSITY.2)
JAGS.DENSITY.AGBP <- JKmcmctab(JAGS.DENSITY.2.mcmc)
# Gelman, A and Rubin, DB (1992) Inference from iterative simulation using multiple sequences, Statistical Science, 7, 457-511.
JAGS.DENSITY.AGBP.diag <- gelman.diag(JAGS.DENSITY.2.mcmc)
JAGS.DENSITY.AGBP <- cbind(JAGS.DENSITY.AGBP, JAGS.DENSITY.AGBP.diag$psrf)
JAGS.DENSITY.AGBP
write.csv(JAGS.DENSITY.AGBP, "E:\\BHM.Results\\...csv")


#
BUGS.DENSITY <- function(){
  for (i in 1:N.resp) {
    lnAGBP[i] ~ dnorm(mu[i], tau)
    mu[i] <- b.Pielou[statecdid[i]] * Pielou[i] + b.statecd[statecdid[i]]
  }
  
  for (j in 1:N.unitcd) {
    b.statecd[j] ~ dnorm(b.statecd.hat[j], tau.statecd)
    b.statecd.hat[j] <- b.statecd0 + b.DENSITY0 * DENSITY[j]
    b.Pielou[j] ~ dnorm(b.Pielou.hat[j], tau.Pielou)
    b.Pielou.hat[j] <- b.statecd1 + b.DENSITY1 * DENSITY[j]
  }
  #
  b.statecd0 ~ dnorm(0, 0.01)
  b.statecd1 ~ dnorm(0, 0.01)
  #
  b.DENSITY0 ~ dnorm(0, 0.01)
  b.DENSITY1 ~ dnorm(0, 0.01)
  #
  tau ~ dgamma(1, 1)
  sigma.statecd <- 1 / tau.statecd
  tau.statecd ~ dgamma(1, 1)
  tau.Pielou ~ dgamma(1,1)
}
#
BUGS.params <- c("b.DENSITY0", "b.DENSITY1",
                 "b.statecd0", "b.statecd1", "b.Pielou")
#
set.seed(123)
JAGS.Density.1 <- jags(data = AGB.data.datjags, inits = NULL, parameters.to.save = BUGS.params,
                       model.file = BUGS.DENSITY, n.chains = 3, n.iter = 20000, n.burnin = 5000)
#
options(max.print = 1000000)
JAGS.DENSITY.2 <- update(JAGS.Density.1, n.iter = 1000)
JAGS.DENSITY.2.DIC <- JAGS.DENSITY.2$BUGSoutput$DIC
JAGS.DENSITY.2.DIC
JAGS.DENSITY.2.mcmc <- as.mcmc(JAGS.DENSITY.2)
JAGS.DENSITY.AGBP <- JKmcmctab(JAGS.DENSITY.2.mcmc)
# Gelman, A and Rubin, DB (1992) Inference from iterative simulation using multiple sequences, Statistical Science, 7, 457-511.
JAGS.DENSITY.AGBP.diag <- gelman.diag(JAGS.DENSITY.2.mcmc)
JAGS.DENSITY.AGBP <- cbind(JAGS.DENSITY.AGBP, JAGS.DENSITY.AGBP.diag$psrf)
JAGS.DENSITY.AGBP
write.csv(JAGS.DENSITY.AGBP, "E:\\BHM.Results\\...csv")




BUGS.DENSITY <- function(){
  for (i in 1:N.resp) {
    lnAGBP[i] ~ dnorm(mu[i], tau)
    mu[i] <- b.DBH.SD[statecdid[i]] * DBH.SD[i] + b.statecd[statecdid[i]]
  }
  
  for (j in 1:N.unitcd) {
    b.statecd[j] ~ dnorm(b.statecd.hat[j], tau.statecd)
    b.statecd.hat[j] <- b.statecd0 + b.DENSITY0 * DENSITY[j]
    b.DBH.SD[j] ~ dnorm(b.DBH.SD.hat[j], tau.DBH.SD)
    b.DBH.SD.hat[j] <- b.statecd1 + b.DENSITY1 * DENSITY[j]
  }
  #
  b.statecd0 ~ dnorm(0, 0.01)
  b.statecd1 ~ dnorm(0, 0.01)
  #
  b.DENSITY0 ~ dnorm(0, 0.01)
  b.DENSITY1 ~ dnorm(0, 0.01)
  #
  tau ~ dgamma(1, 1)
  sigma.statecd <- 1 / tau.statecd
  tau.statecd ~ dgamma(1, 1)
  tau.DBH.SD ~ dgamma(1,1)
}
#
BUGS.params <- c("b.DENSITY0", "b.DENSITY1",
                 "b.statecd0", "b.statecd1", "b.DBH.SD")
#
set.seed(123)
JAGS.Density.1 <- jags(data = AGB.data.datjags, inits = NULL, parameters.to.save = BUGS.params,
                       model.file = BUGS.DENSITY, n.chains = 3, n.iter = 20000, n.burnin = 5000)
#
options(max.print = 1000000)
JAGS.DENSITY.2 <- update(JAGS.Density.1, n.iter = 1000)
JAGS.DENSITY.2.DIC <- JAGS.DENSITY.2$BUGSoutput$DIC
JAGS.DENSITY.2.DIC
JAGS.DENSITY.2.mcmc <- as.mcmc(JAGS.DENSITY.2)
JAGS.DENSITY.AGBP <- JKmcmctab(JAGS.DENSITY.2.mcmc)
# Gelman, A and Rubin, DB (1992) Inference from iterative simulation using multiple sequences, Statistical Science, 7, 457-511.
JAGS.DENSITY.AGBP.diag <- gelman.diag(JAGS.DENSITY.2.mcmc)
JAGS.DENSITY.AGBP <- cbind(JAGS.DENSITY.AGBP, JAGS.DENSITY.AGBP.diag$psrf)
JAGS.DENSITY.AGBP
write.csv(JAGS.DENSITY.AGBP, "E:\\BHM.Results\\...csv")


#
BUGS.DENSITY <- function(){
  for (i in 1:N.resp) {
    lnAGBP[i] ~ dnorm(mu[i], tau)
    mu[i] <- b.DBH.CV[statecdid[i]] * DBH.CV[i] + b.statecd[statecdid[i]]
  }
  
  for (j in 1:N.unitcd) {
    b.statecd[j] ~ dnorm(b.statecd.hat[j], tau.statecd)
    b.statecd.hat[j] <- b.statecd0 + b.DENSITY0 * DENSITY[j]
    b.DBH.CV[j] ~ dnorm(b.DBH.CV.hat[j], tau.DBH.CV)
    b.DBH.CV.hat[j] <- b.statecd1 + b.DENSITY1 * DENSITY[j]
  }
  #
  b.statecd0 ~ dnorm(0, 0.01)
  b.statecd1 ~ dnorm(0, 0.01)
  #
  b.DENSITY0 ~ dnorm(0, 0.01)
  b.DENSITY1 ~ dnorm(0, 0.01)
  #
  tau ~ dgamma(1, 1)
  sigma.statecd <- 1 / tau.statecd
  tau.statecd ~ dgamma(1, 1)
  tau.DBH.CV ~ dgamma(1,1)
}
#
BUGS.params <- c("b.DENSITY0", "b.DENSITY1",
                 "b.statecd0", "b.statecd1", "b.DBH.CV")
#
set.seed(123)
JAGS.Density.1 <- jags(data = AGB.data.datjags, inits = NULL, parameters.to.save = BUGS.params,
                       model.file = BUGS.DENSITY, n.chains = 3, n.iter = 20000, n.burnin = 5000)
#
options(max.print = 1000000)
JAGS.DENSITY.2 <- update(JAGS.Density.1, n.iter = 1000)
JAGS.DENSITY.2.DIC <- JAGS.DENSITY.2$BUGSoutput$DIC
JAGS.DENSITY.2.DIC
JAGS.DENSITY.2.mcmc <- as.mcmc(JAGS.DENSITY.2)
JAGS.DENSITY.AGBP <- JKmcmctab(JAGS.DENSITY.2.mcmc)
# Gelman, A and Rubin, DB (1992) Inference from iterative simulation using multiple sequences, Statistical Science, 7, 457-511.
JAGS.DENSITY.AGBP.diag <- gelman.diag(JAGS.DENSITY.2.mcmc)
JAGS.DENSITY.AGBP <- cbind(JAGS.DENSITY.AGBP, JAGS.DENSITY.AGBP.diag$psrf)
JAGS.DENSITY.AGBP
write.csv(JAGS.DENSITY.AGBP, "E:\\BHM.Results\\...csv")






BUGS.DENSITY <- function(){
  for (i in 1:N.resp) {
    lnAGBP[i] ~ dnorm(mu[i], tau)
    mu[i] <- b.DBH.Gini[statecdid[i]] * DBH.Gini[i] + b.statecd[statecdid[i]]
  }
  
  for (j in 1:N.unitcd) {
    b.statecd[j] ~ dnorm(b.statecd.hat[j], tau.statecd)
    b.statecd.hat[j] <- b.statecd0 + b.DENSITY0 * DENSITY[j]
    b.DBH.Gini[j] ~ dnorm(b.DBH.Gini.hat[j], tau.DBH.Gini)
    b.DBH.Gini.hat[j] <- b.statecd1 + b.DENSITY1 * DENSITY[j]
  }
  #
  b.statecd0 ~ dnorm(0, 0.01)
  b.statecd1 ~ dnorm(0, 0.01)
  #
  b.DENSITY0 ~ dnorm(0, 0.01)
  b.DENSITY1 ~ dnorm(0, 0.01)
  #
  tau ~ dgamma(1, 1)
  sigma.statecd <- 1 / tau.statecd
  tau.statecd ~ dgamma(1, 1)
  tau.DBH.Gini ~ dgamma(1,1)
}
#
BUGS.params <- c("b.DENSITY0", "b.DENSITY1",
                 "b.statecd0", "b.statecd1", "b.DBH.Gini")
#
set.seed(123)
JAGS.Density.1 <- jags(data = AGB.data.datjags, inits = NULL, parameters.to.save = BUGS.params,
                       model.file = BUGS.DENSITY, n.chains = 3, n.iter = 20000, n.burnin = 5000)
#
options(max.print = 1000000)
JAGS.DENSITY.2 <- update(JAGS.Density.1, n.iter = 1000)
JAGS.DENSITY.2.DIC <- JAGS.DENSITY.2$BUGSoutput$DIC
JAGS.DENSITY.2.DIC
JAGS.DENSITY.2.mcmc <- as.mcmc(JAGS.DENSITY.2)
JAGS.DENSITY.AGBP <- JKmcmctab(JAGS.DENSITY.2.mcmc)
# Gelman, A and Rubin, DB (1992) Inference from iterative simulation using multiple sequences, Statistical Science, 7, 457-511.
JAGS.DENSITY.AGBP.diag <- gelman.diag(JAGS.DENSITY.2.mcmc)
JAGS.DENSITY.AGBP <- cbind(JAGS.DENSITY.AGBP, JAGS.DENSITY.AGBP.diag$psrf)
JAGS.DENSITY.AGBP
write.csv(JAGS.DENSITY.AGBP, "E:\\BHM.Results\\...csv")



##############################################################################################################

BUGS.DENSITY <- function(){
  for (i in 1:N.resp) {
    AGBS[i] ~ dnorm(mu[i], tau)
    mu[i] <- b.Simpson[statecdid[i]] * Simpson[i] + b.statecd[statecdid[i]]
  }
  
  for (j in 1:N.unitcd) {
    b.statecd[j] ~ dnorm(b.statecd.hat[j], tau.statecd)
    b.statecd.hat[j] <- b.statecd0 + b.DENSITY0 * DENSITY[j]
    b.Simpson[j] ~ dnorm(b.Simpson.hat[j], tau.Simpson)
    b.Simpson.hat[j] <- b.statecd1 + b.DENSITY1 * DENSITY[j]
  }
  #
  b.statecd0 ~ dnorm(0, 0.01)
  b.statecd1 ~ dnorm(0, 0.01)
  #
  b.DENSITY0 ~ dnorm(0, 0.01)
  b.DENSITY1 ~ dnorm(0, 0.01)
  #
  tau ~ dgamma(1, 1)
  sigma.statecd <- 1 / tau.statecd
  tau.statecd ~ dgamma(1, 1)
  tau.Simpson ~ dgamma(1,1)
}
#
BUGS.params <- c("b.DENSITY0", "b.DENSITY1",
                 "b.statecd0", "b.statecd1", "b.Simpson")
#
set.seed(123)
JAGS.Density.1 <- jags(data = AGB.data.datjags, inits = NULL, parameters.to.save = BUGS.params,
                       model.file = BUGS.DENSITY, n.chains = 3, n.iter = 20000, n.burnin = 5000)
#
options(max.print = 1000000)
JAGS.DENSITY.2 <- update(JAGS.Density.1, n.iter = 1000)
JAGS.DENSITY.2.DIC <- JAGS.DENSITY.2$BUGSoutput$DIC
JAGS.DENSITY.2.DIC
JAGS.DENSITY.2.mcmc <- as.mcmc(JAGS.DENSITY.2)
JAGS.DENSITY.AGB <- JKmcmctab(JAGS.DENSITY.2.mcmc)
# Gelman, A and Rubin, DB (1992) Inference from iterative simulation using multiple sequences, Statistical Science, 7, 457-511.
JAGS.DENSITY.AGB.diag <- gelman.diag(JAGS.DENSITY.2.mcmc)
JAGS.DENSITY.AGB <- cbind(JAGS.DENSITY.AGB, JAGS.DENSITY.AGB.diag$psrf)
JAGS.DENSITY.AGB
write.csv(JAGS.DENSITY.AGB, "E:\\BHM.Results\\...csv")





#
BUGS.DENSITY <- function(){
  for (i in 1:N.resp) {
    AGBS[i] ~ dnorm(mu[i], tau)
    mu[i] <- b.Shannon.Wiener[statecdid[i]] * Shannon.Wiener[i] + b.statecd[statecdid[i]]
  }
  
  for (j in 1:N.unitcd) {
    b.statecd[j] ~ dnorm(b.statecd.hat[j], tau.statecd)
    b.statecd.hat[j] <- b.statecd0 + b.DENSITY0 * DENSITY[j]
    b.Shannon.Wiener[j] ~ dnorm(b.Shannon.Wiener.hat[j], tau.Shannon.Wiener)
    b.Shannon.Wiener.hat[j] <- b.statecd1 + b.DENSITY1 * DENSITY[j]
  }
  #
  b.statecd0 ~ dnorm(0, 0.01)
  b.statecd1 ~ dnorm(0, 0.01)
  #
  b.DENSITY0 ~ dnorm(0, 0.01)
  b.DENSITY1 ~ dnorm(0, 0.01)
  #
  tau ~ dgamma(1, 1)
  sigma.statecd <- 1 / tau.statecd
  tau.statecd ~ dgamma(1, 1)
  tau.Shannon.Wiener ~ dgamma(1,1)
}
#
BUGS.params <- c("b.DENSITY0", "b.DENSITY1",
                 "b.statecd0", "b.statecd1", "b.Shannon.Wiener")
#
set.seed(123)
JAGS.Density.1 <- jags(data = AGB.data.datjags, inits = NULL, parameters.to.save = BUGS.params,
                       model.file = BUGS.DENSITY, n.chains = 3, n.iter = 20000, n.burnin = 5000)
#
options(max.print = 1000000)
JAGS.DENSITY.2 <- update(JAGS.Density.1, n.iter = 1000)
JAGS.DENSITY.2.DIC <- JAGS.DENSITY.2$BUGSoutput$DIC
JAGS.DENSITY.2.DIC
JAGS.DENSITY.2.mcmc <- as.mcmc(JAGS.DENSITY.2)
JAGS.DENSITY.AGB <- JKmcmctab(JAGS.DENSITY.2.mcmc)
# Gelman, A and Rubin, DB (1992) Inference from iterative simulation using multiple sequences, Statistical Science, 7, 457-511.
JAGS.DENSITY.AGB.diag <- gelman.diag(JAGS.DENSITY.2.mcmc)
JAGS.DENSITY.AGB <- cbind(JAGS.DENSITY.AGB, JAGS.DENSITY.AGB.diag$psrf)
JAGS.DENSITY.AGB
write.csv(JAGS.DENSITY.AGB, "E:\\BHM.Results\\...csv")




#
BUGS.DENSITY <- function(){
  for (i in 1:N.resp) {
    AGBS[i] ~ dnorm(mu[i], tau)
    mu[i] <- b.Pielou[statecdid[i]] * Pielou[i] + b.statecd[statecdid[i]]
  }
  
  for (j in 1:N.unitcd) {
    b.statecd[j] ~ dnorm(b.statecd.hat[j], tau.statecd)
    b.statecd.hat[j] <- b.statecd0 + b.DENSITY0 * DENSITY[j]
    b.Pielou[j] ~ dnorm(b.Pielou.hat[j], tau.Pielou)
    b.Pielou.hat[j] <- b.statecd1 + b.DENSITY1 * DENSITY[j]
  }
  #
  b.statecd0 ~ dnorm(0, 0.01)
  b.statecd1 ~ dnorm(0, 0.01)
  #
  b.DENSITY0 ~ dnorm(0, 0.01)
  b.DENSITY1 ~ dnorm(0, 0.01)
  #
  tau ~ dgamma(1, 1)
  sigma.statecd <- 1 / tau.statecd
  tau.statecd ~ dgamma(1, 1)
  tau.Pielou ~ dgamma(1,1)
}
#
BUGS.params <- c("b.DENSITY0", "b.DENSITY1",
                 "b.statecd0", "b.statecd1", "b.Pielou")
#
set.seed(123)
JAGS.Density.1 <- jags(data = AGB.data.datjags, inits = NULL, parameters.to.save = BUGS.params,
                       model.file = BUGS.DENSITY, n.chains = 3, n.iter = 20000, n.burnin = 5000)
#
options(max.print = 1000000)
JAGS.DENSITY.2 <- update(JAGS.Density.1, n.iter = 1000)
JAGS.DENSITY.2.DIC <- JAGS.DENSITY.2$BUGSoutput$DIC
JAGS.DENSITY.2.DIC
JAGS.DENSITY.2.mcmc <- as.mcmc(JAGS.DENSITY.2)
JAGS.DENSITY.AGB <- JKmcmctab(JAGS.DENSITY.2.mcmc)
# Gelman, A and Rubin, DB (1992) Inference from iterative simulation using multiple sequences, Statistical Science, 7, 457-511.
JAGS.DENSITY.AGB.diag <- gelman.diag(JAGS.DENSITY.2.mcmc)
JAGS.DENSITY.AGB <- cbind(JAGS.DENSITY.AGB, JAGS.DENSITY.AGB.diag$psrf)
JAGS.DENSITY.AGB
write.csv(JAGS.DENSITY.AGB, "E:\\BHM.Results\\...csv")




BUGS.DENSITY <- function(){
  for (i in 1:N.resp) {
    AGBS[i] ~ dnorm(mu[i], tau)
    mu[i] <- b.DBH.SD[statecdid[i]] * DBH.SD[i] + b.statecd[statecdid[i]]
  }
  
  for (j in 1:N.unitcd) {
    b.statecd[j] ~ dnorm(b.statecd.hat[j], tau.statecd)
    b.statecd.hat[j] <- b.statecd0 + b.DENSITY0 * DENSITY[j]
    b.DBH.SD[j] ~ dnorm(b.DBH.SD.hat[j], tau.DBH.SD)
    b.DBH.SD.hat[j] <- b.statecd1 + b.DENSITY1 * DENSITY[j]
  }
  #
  b.statecd0 ~ dnorm(0, 0.01)
  b.statecd1 ~ dnorm(0, 0.01)
  #
  b.DENSITY0 ~ dnorm(0, 0.01)
  b.DENSITY1 ~ dnorm(0, 0.01)
  #
  tau ~ dgamma(1, 1)
  sigma.statecd <- 1 / tau.statecd
  tau.statecd ~ dgamma(1, 1)
  tau.DBH.SD ~ dgamma(1,1)
}
#
BUGS.params <- c("b.DENSITY0", "b.DENSITY1",
                 "b.statecd0", "b.statecd1", "b.DBH.SD")
#
set.seed(123)
JAGS.Density.1 <- jags(data = AGB.data.datjags, inits = NULL, parameters.to.save = BUGS.params,
                       model.file = BUGS.DENSITY, n.chains = 3, n.iter = 20000, n.burnin = 5000)
#
options(max.print = 1000000)
JAGS.DENSITY.2 <- update(JAGS.Density.1, n.iter = 1000)
JAGS.DENSITY.2.DIC <- JAGS.DENSITY.2$BUGSoutput$DIC
JAGS.DENSITY.2.DIC
JAGS.DENSITY.2.mcmc <- as.mcmc(JAGS.DENSITY.2)
JAGS.DENSITY.AGB <- JKmcmctab(JAGS.DENSITY.2.mcmc)
# Gelman, A and Rubin, DB (1992) Inference from iterative simulation using multiple sequences, Statistical Science, 7, 457-511.
JAGS.DENSITY.AGB.diag <- gelman.diag(JAGS.DENSITY.2.mcmc)
JAGS.DENSITY.AGB <- cbind(JAGS.DENSITY.AGB, JAGS.DENSITY.AGB.diag$psrf)
JAGS.DENSITY.AGB
write.csv(JAGS.DENSITY.AGB, "E:\\BHM.Results\\...csv")





BUGS.DENSITY <- function(){
  for (i in 1:N.resp) {
    AGBS[i] ~ dnorm(mu[i], tau)
    mu[i] <- b.DBH.CV[statecdid[i]] * DBH.CV[i] + b.statecd[statecdid[i]]
  }
  
  for (j in 1:N.unitcd) {
    b.statecd[j] ~ dnorm(b.statecd.hat[j], tau.statecd)
    b.statecd.hat[j] <- b.statecd0 + b.DENSITY0 * DENSITY[j]
    b.DBH.CV[j] ~ dnorm(b.DBH.CV.hat[j], tau.DBH.CV)
    b.DBH.CV.hat[j] <- b.statecd1 + b.DENSITY1 * DENSITY[j]
  }
  #
  b.statecd0 ~ dnorm(0, 0.01)
  b.statecd1 ~ dnorm(0, 0.01)
  #
  b.DENSITY0 ~ dnorm(0, 0.01)
  b.DENSITY1 ~ dnorm(0, 0.01)
  #
  tau ~ dgamma(1, 1)
  sigma.statecd <- 1 / tau.statecd
  tau.statecd ~ dgamma(1, 1)
  tau.DBH.CV ~ dgamma(1,1)
}
#
BUGS.params <- c("b.DENSITY0", "b.DENSITY1",
                 "b.statecd0", "b.statecd1", "b.DBH.CV")
#
set.seed(123)
JAGS.Density.1 <- jags(data = AGB.data.datjags, inits = NULL, parameters.to.save = BUGS.params,
                       model.file = BUGS.DENSITY, n.chains = 3, n.iter = 20000, n.burnin = 5000)
#
options(max.print = 1000000)
JAGS.DENSITY.2 <- update(JAGS.Density.1, n.iter = 1000)
JAGS.DENSITY.2.DIC <- JAGS.DENSITY.2$BUGSoutput$DIC
JAGS.DENSITY.2.DIC
JAGS.DENSITY.2.mcmc <- as.mcmc(JAGS.DENSITY.2)
JAGS.DENSITY.AGB <- JKmcmctab(JAGS.DENSITY.2.mcmc)
# Gelman, A and Rubin, DB (1992) Inference from iterative simulation using multiple sequences, Statistical Science, 7, 457-511.
JAGS.DENSITY.AGB.diag <- gelman.diag(JAGS.DENSITY.2.mcmc)
JAGS.DENSITY.AGB <- cbind(JAGS.DENSITY.AGB, JAGS.DENSITY.AGB.diag$psrf)
JAGS.DENSITY.AGB
write.csv(JAGS.DENSITY.AGB, "E:\\BHM.Results\\...csv")






BUGS.DENSITY <- function(){
  for (i in 1:N.resp) {
    AGBS[i] ~ dnorm(mu[i], tau)
    mu[i] <- b.DBH.Gini[statecdid[i]] * DBH.Gini[i] + b.statecd[statecdid[i]]
  }
  
  for (j in 1:N.unitcd) {
    b.statecd[j] ~ dnorm(b.statecd.hat[j], tau.statecd)
    b.statecd.hat[j] <- b.statecd0 + b.DENSITY0 * DENSITY[j]
    b.DBH.Gini[j] ~ dnorm(b.DBH.Gini.hat[j], tau.DBH.Gini)
    b.DBH.Gini.hat[j] <- b.statecd1 + b.DENSITY1 * DENSITY[j]
  }
  #
  b.statecd0 ~ dnorm(0, 0.01)
  b.statecd1 ~ dnorm(0, 0.01)
  #
  b.DENSITY0 ~ dnorm(0, 0.01)
  b.DENSITY1 ~ dnorm(0, 0.01)
  #
  tau ~ dgamma(1, 1)
  sigma.statecd <- 1 / tau.statecd
  tau.statecd ~ dgamma(1, 1)
  tau.DBH.Gini ~ dgamma(1,1)
}
#
BUGS.params <- c("b.DENSITY0", "b.DENSITY1",
                 "b.statecd0", "b.statecd1", "b.DBH.Gini")
#
set.seed(123)
JAGS.Density.1 <- jags(data = AGB.data.datjags, inits = NULL, parameters.to.save = BUGS.params,
                       model.file = BUGS.DENSITY, n.chains = 3, n.iter = 20000, n.burnin = 5000)
#
options(max.print = 1000000)
JAGS.DENSITY.2 <- update(JAGS.Density.1, n.iter = 1000)
JAGS.DENSITY.2.DIC <- JAGS.DENSITY.2$BUGSoutput$DIC
JAGS.DENSITY.2.DIC
JAGS.DENSITY.2.mcmc <- as.mcmc(JAGS.DENSITY.2)
JAGS.DENSITY.AGB <- JKmcmctab(JAGS.DENSITY.2.mcmc)
# Gelman, A and Rubin, DB (1992) Inference from iterative simulation using multiple sequences, Statistical Science, 7, 457-511.
JAGS.DENSITY.AGB.diag <- gelman.diag(JAGS.DENSITY.2.mcmc)
JAGS.DENSITY.AGB <- cbind(JAGS.DENSITY.AGB, JAGS.DENSITY.AGB.diag$psrf)
JAGS.DENSITY.AGB
write.csv(JAGS.DENSITY.AGB, "E:\\BHM.Results\\...csv")


#######################################################################################################


BUGS.SIZE <- function(){
  for (i in 1:N.resp) {
    lnAGBP[i] ~ dnorm(mu[i], tau)
    mu[i] <- b.Simpson[statecdid[i]] * Simpson[i] + b.statecd[statecdid[i]]
  }
  
  for (j in 1:N.unitcd) {
    b.statecd[j] ~ dnorm(b.statecd.hat[j], tau.statecd)
    b.statecd.hat[j] <- b.statecd0 + b.SIZE0 * SIZE[j]
    b.Simpson[j] ~ dnorm(b.Simpson.hat[j], tau.Simpson)
    b.Simpson.hat[j] <- b.statecd1 + b.SIZE1 * SIZE[j]
  }
  #
  b.statecd0 ~ dnorm(0, 0.01)
  b.statecd1 ~ dnorm(0, 0.01)
  #
  b.SIZE0 ~ dnorm(0, 0.01)
  b.SIZE1 ~ dnorm(0, 0.01)
  #
  tau ~ dgamma(1, 1)
  sigma.statecd <- 1 / tau.statecd
  tau.statecd ~ dgamma(1, 1)
  tau.Simpson ~ dgamma(1,1)
}
#
BUGS.params <- c("b.SIZE0", "b.SIZE1",
                 "b.statecd0", "b.statecd1", "b.Simpson")
#
set.seed(123)
JAGS.SIZE <- jags(data = AGB.data.datjags, inits = NULL, parameters.to.save = BUGS.params,
                  model.file = BUGS.SIZE, n.chains = 3, n.iter = 20000, n.burnin = 5000)
#
options(max.print = 1000000)
JAGS.SIZE.2 <- update(JAGS.SIZE, n.iter = 1000)
JAGS.SIZE.2.DIC <- JAGS.SIZE.2$BUGSoutput$DIC
JAGS.SIZE.2.DIC
JAGS.SIZE.2.mcmc <- as.mcmc(JAGS.SIZE.2)
JAGS.SIZE.AGBP <- JKmcmctab(JAGS.SIZE.2.mcmc)
JAGS.SIZE.AGBP.diag <- gelman.diag(JAGS.SIZE.2.mcmc)
JAGS.SIZE.AGBP <- cbind(JAGS.SIZE.AGBP, JAGS.SIZE.AGBP.diag$psrf)
JAGS.SIZE.AGBP
write.csv(JAGS.SIZE.AGBP, "E:\\BHM.Results\\...csv")




BUGS.SIZE <- function(){
  for (i in 1:N.resp) {
    lnAGBP[i] ~ dnorm(mu[i], tau)
    mu[i] <- b.Shannon.Wiener[statecdid[i]] * Shannon.Wiener[i] + b.statecd[statecdid[i]]
  }
  
  for (j in 1:N.unitcd) {
    b.statecd[j] ~ dnorm(b.statecd.hat[j], tau.statecd)
    b.statecd.hat[j] <- b.statecd0 + b.SIZE0 * SIZE[j]
    b.Shannon.Wiener[j] ~ dnorm(b.Shannon.Wiener.hat[j], tau.Shannon.Wiener)
    b.Shannon.Wiener.hat[j] <- b.statecd1 + b.SIZE1 * SIZE[j]
  }
  #
  b.statecd0 ~ dnorm(0, 0.01)
  b.statecd1 ~ dnorm(0, 0.01)
  #
  b.SIZE0 ~ dnorm(0, 0.01)
  b.SIZE1 ~ dnorm(0, 0.01)
  #
  tau ~ dgamma(1, 1)
  sigma.statecd <- 1 / tau.statecd
  tau.statecd ~ dgamma(1, 1)
  tau.Shannon.Wiener ~ dgamma(1,1)
}
#
BUGS.params <- c("b.SIZE0", "b.SIZE1",
                 "b.statecd0", "b.statecd1", "b.Shannon.Wiener")
#
set.seed(123)
JAGS.SIZE <- jags(data = AGB.data.datjags, inits = NULL, parameters.to.save = BUGS.params,
                  model.file = BUGS.SIZE, n.chains = 3, n.iter = 20000, n.burnin = 5000)
#
options(max.print = 1000000)
JAGS.SIZE.2 <- update(JAGS.SIZE, n.iter = 1000)
JAGS.SIZE.2.DIC <- JAGS.SIZE.2$BUGSoutput$DIC
JAGS.SIZE.2.DIC
JAGS.SIZE.2.mcmc <- as.mcmc(JAGS.SIZE.2)
JAGS.SIZE.AGBP <- JKmcmctab(JAGS.SIZE.2.mcmc)
JAGS.SIZE.AGBP.diag <- gelman.diag(JAGS.SIZE.2.mcmc)
JAGS.SIZE.AGBP <- cbind(JAGS.SIZE.AGBP, JAGS.SIZE.AGBP.diag$psrf)
JAGS.SIZE.AGBP
write.csv(JAGS.SIZE.AGBP, "E:\\BHM.Results\\...csv")






BUGS.SIZE <- function(){
  for (i in 1:N.resp) {
    lnAGBP[i] ~ dnorm(mu[i], tau)
    mu[i] <- b.Pielou[statecdid[i]] * Pielou[i] + b.statecd[statecdid[i]]
  }
  
  for (j in 1:N.unitcd) {
    b.statecd[j] ~ dnorm(b.statecd.hat[j], tau.statecd)
    b.statecd.hat[j] <- b.statecd0 + b.SIZE0 * SIZE[j]
    b.Pielou[j] ~ dnorm(b.Pielou.hat[j], tau.Pielou)
    b.Pielou.hat[j] <- b.statecd1 + b.SIZE1 * SIZE[j]
  }
  #
  b.statecd0 ~ dnorm(0, 0.01)
  b.statecd1 ~ dnorm(0, 0.01)
  #
  b.SIZE0 ~ dnorm(0, 0.01)
  b.SIZE1 ~ dnorm(0, 0.01)
  #
  tau ~ dgamma(1, 1)
  sigma.statecd <- 1 / tau.statecd
  tau.statecd ~ dgamma(1, 1)
  tau.Pielou ~ dgamma(1,1)
}
#
BUGS.params <- c("b.SIZE0", "b.SIZE1",
                 "b.statecd0", "b.statecd1", "b.Pielou")
#
set.seed(123)
JAGS.SIZE <- jags(data = AGB.data.datjags, inits = NULL, parameters.to.save = BUGS.params,
                  model.file = BUGS.SIZE, n.chains = 3, n.iter = 20000, n.burnin = 5000)
#
options(max.print = 1000000)
JAGS.SIZE.2 <- update(JAGS.SIZE, n.iter = 1000)
JAGS.SIZE.2.DIC <- JAGS.SIZE.2$BUGSoutput$DIC
JAGS.SIZE.2.DIC
JAGS.SIZE.2.mcmc <- as.mcmc(JAGS.SIZE.2)
JAGS.SIZE.AGBP <- JKmcmctab(JAGS.SIZE.2.mcmc)
JAGS.SIZE.AGBP.diag <- gelman.diag(JAGS.SIZE.2.mcmc)
JAGS.SIZE.AGBP <- cbind(JAGS.SIZE.AGBP, JAGS.SIZE.AGBP.diag$psrf)
JAGS.SIZE.AGBP
write.csv(JAGS.SIZE.AGBP, "E:\\BHM.Results\\...csv")





BUGS.SIZE <- function(){
  for (i in 1:N.resp) {
    lnAGBP[i] ~ dnorm(mu[i], tau)
    mu[i] <- b.DBH.SD[statecdid[i]] * DBH.SD[i] + b.statecd[statecdid[i]]
  }
  
  for (j in 1:N.unitcd) {
    b.statecd[j] ~ dnorm(b.statecd.hat[j], tau.statecd)
    b.statecd.hat[j] <- b.statecd0 + b.SIZE0 * SIZE[j]
    b.DBH.SD[j] ~ dnorm(b.DBH.SD.hat[j], tau.DBH.SD)
    b.DBH.SD.hat[j] <- b.statecd1 + b.SIZE1 * SIZE[j]
  }
  #
  b.statecd0 ~ dnorm(0, 0.01)
  b.statecd1 ~ dnorm(0, 0.01)
  #
  b.SIZE0 ~ dnorm(0, 0.01)
  b.SIZE1 ~ dnorm(0, 0.01)
  #
  tau ~ dgamma(1, 1)
  sigma.statecd <- 1 / tau.statecd
  tau.statecd ~ dgamma(1, 1)
  tau.DBH.SD ~ dgamma(1,1)
}
#
BUGS.params <- c("b.SIZE0", "b.SIZE1",
                 "b.statecd0", "b.statecd1", "b.DBH.SD")
#
set.seed(123)
JAGS.SIZE <- jags(data = AGB.data.datjags, inits = NULL, parameters.to.save = BUGS.params,
                  model.file = BUGS.SIZE, n.chains = 3, n.iter = 20000, n.burnin = 5000)
#
options(max.print = 1000000)
JAGS.SIZE.2 <- update(JAGS.SIZE, n.iter = 1000)
JAGS.SIZE.2.DIC <- JAGS.SIZE.2$BUGSoutput$DIC
JAGS.SIZE.2.DIC
JAGS.SIZE.2.mcmc <- as.mcmc(JAGS.SIZE.2)
JAGS.SIZE.AGBP <- JKmcmctab(JAGS.SIZE.2.mcmc)
JAGS.SIZE.AGBP.diag <- gelman.diag(JAGS.SIZE.2.mcmc)
JAGS.SIZE.AGBP <- cbind(JAGS.SIZE.AGBP, JAGS.SIZE.AGBP.diag$psrf)
JAGS.SIZE.AGBP
write.csv(JAGS.SIZE.AGBP, "E:\\BHM.Results\\...csv")



BUGS.SIZE <- function(){
  for (i in 1:N.resp) {
    lnAGBP[i] ~ dnorm(mu[i], tau)
    mu[i] <- b.DBH.CV[statecdid[i]] * DBH.CV[i] + b.statecd[statecdid[i]]
  }
  
  for (j in 1:N.unitcd) {
    b.statecd[j] ~ dnorm(b.statecd.hat[j], tau.statecd)
    b.statecd.hat[j] <- b.statecd0 + b.SIZE0 * SIZE[j]
    b.DBH.CV[j] ~ dnorm(b.DBH.CV.hat[j], tau.DBH.CV)
    b.DBH.CV.hat[j] <- b.statecd1 + b.SIZE1 * SIZE[j]
  }
  #
  b.statecd0 ~ dnorm(0, 0.01)
  b.statecd1 ~ dnorm(0, 0.01)
  #
  b.SIZE0 ~ dnorm(0, 0.01)
  b.SIZE1 ~ dnorm(0, 0.01)
  #
  tau ~ dgamma(1, 1)
  sigma.statecd <- 1 / tau.statecd
  tau.statecd ~ dgamma(1, 1)
  tau.DBH.CV ~ dgamma(1,1)
}
#
BUGS.params <- c("b.SIZE0", "b.SIZE1",
                 "b.statecd0", "b.statecd1", "b.DBH.CV")
#
set.seed(123)
JAGS.SIZE <- jags(data = AGB.data.datjags, inits = NULL, parameters.to.save = BUGS.params,
                  model.file = BUGS.SIZE, n.chains = 3, n.iter = 20000, n.burnin = 5000)
#
options(max.print = 1000000)
JAGS.SIZE.2 <- update(JAGS.SIZE, n.iter = 1000)
JAGS.SIZE.2.DIC <- JAGS.SIZE.2$BUGSoutput$DIC
JAGS.SIZE.2.DIC
JAGS.SIZE.2.mcmc <- as.mcmc(JAGS.SIZE.2)
JAGS.SIZE.AGBP <- JKmcmctab(JAGS.SIZE.2.mcmc)
JAGS.SIZE.AGBP.diag <- gelman.diag(JAGS.SIZE.2.mcmc)
JAGS.SIZE.AGBP <- cbind(JAGS.SIZE.AGBP, JAGS.SIZE.AGBP.diag$psrf)
JAGS.SIZE.AGBP
write.csv(JAGS.SIZE.AGBP, "E:\\BHM.Results\\...csv")



BUGS.SIZE <- function(){
  for (i in 1:N.resp) {
    lnAGBP[i] ~ dnorm(mu[i], tau)
    mu[i] <- b.DBH.Gini[statecdid[i]] * DBH.Gini[i] + b.statecd[statecdid[i]]
  }
  
  for (j in 1:N.unitcd) {
    b.statecd[j] ~ dnorm(b.statecd.hat[j], tau.statecd)
    b.statecd.hat[j] <- b.statecd0 + b.SIZE0 * SIZE[j]
    b.DBH.Gini[j] ~ dnorm(b.DBH.Gini.hat[j], tau.DBH.Gini)
    b.DBH.Gini.hat[j] <- b.statecd1 + b.SIZE1 * SIZE[j]
  }
  #
  b.statecd0 ~ dnorm(0, 0.01)
  b.statecd1 ~ dnorm(0, 0.01)
  #
  b.SIZE0 ~ dnorm(0, 0.01)
  b.SIZE1 ~ dnorm(0, 0.01)
  #
  tau ~ dgamma(1, 1)
  sigma.statecd <- 1 / tau.statecd
  tau.statecd ~ dgamma(1, 1)
  tau.DBH.Gini ~ dgamma(1,1)
}
#
BUGS.params <- c("b.SIZE0", "b.SIZE1",
                 "b.statecd0", "b.statecd1", "b.DBH.Gini")
#
set.seed(123)
JAGS.SIZE <- jags(data = AGB.data.datjags, inits = NULL, parameters.to.save = BUGS.params,
                  model.file = BUGS.SIZE, n.chains = 3, n.iter = 20000, n.burnin = 5000)
#
options(max.print = 1000000)
JAGS.SIZE.2 <- update(JAGS.SIZE, n.iter = 1000)
JAGS.SIZE.2.DIC <- JAGS.SIZE.2$BUGSoutput$DIC
JAGS.SIZE.2.DIC
JAGS.SIZE.2.mcmc <- as.mcmc(JAGS.SIZE.2)
JAGS.SIZE.AGBP <- JKmcmctab(JAGS.SIZE.2.mcmc)
JAGS.SIZE.AGBP.diag <- gelman.diag(JAGS.SIZE.2.mcmc)
JAGS.SIZE.AGBP <- cbind(JAGS.SIZE.AGBP, JAGS.SIZE.AGBP.diag$psrf)
JAGS.SIZE.AGBP
write.csv(JAGS.SIZE.AGBP, "E:\\BHM.Results\\...csv")




#######################################################################################################

BUGS.SIZE <- function(){
  for (i in 1:N.resp) {
    AGBS[i] ~ dnorm(mu[i], tau)
    mu[i] <- b.Simpson[statecdid[i]] * Simpson[i] + b.statecd[statecdid[i]]
  }
  
  for (j in 1:N.unitcd) {
    b.statecd[j] ~ dnorm(b.statecd.hat[j], tau.statecd)
    b.statecd.hat[j] <- b.statecd0 + b.SIZE0 * SIZE[j]
    b.Simpson[j] ~ dnorm(b.Simpson.hat[j], tau.Simpson)
    b.Simpson.hat[j] <- b.statecd1 + b.SIZE1 * SIZE[j]
  }
  #
  b.statecd0 ~ dnorm(0, 0.01)
  b.statecd1 ~ dnorm(0, 0.01)
  #
  b.SIZE0 ~ dnorm(0, 0.01)
  b.SIZE1 ~ dnorm(0, 0.01)
  #
  tau ~ dgamma(1, 1)
  sigma.statecd <- 1 / tau.statecd
  tau.statecd ~ dgamma(1, 1)
  tau.Simpson ~ dgamma(1,1)
}
#
BUGS.params <- c("b.SIZE0", "b.SIZE1",
                 "b.statecd0", "b.statecd1", "b.Simpson")
#
set.seed(123)
JAGS.SIZE <- jags(data = AGB.data.datjags, inits = NULL, parameters.to.save = BUGS.params,
                  model.file = BUGS.SIZE, n.chains = 3, n.iter = 20000, n.burnin = 5000)
#
options(max.print = 1000000)
JAGS.SIZE.2 <- update(JAGS.SIZE, n.iter = 1000)
JAGS.SIZE.2.DIC <- JAGS.SIZE.2$BUGSoutput$DIC
JAGS.SIZE.2.DIC
JAGS.SIZE.2.mcmc <- as.mcmc(JAGS.SIZE.2)
JAGS.SIZE.AGB <- JKmcmctab(JAGS.SIZE.2.mcmc)
JAGS.SIZE.AGB.diag <- gelman.diag(JAGS.SIZE.2.mcmc)
JAGS.SIZE.AGB <- cbind(JAGS.SIZE.AGB, JAGS.SIZE.AGB.diag$psrf)
JAGS.SIZE.AGB
write.csv(JAGS.SIZE.AGB, "E:\\BHM.Results\\...csv")




BUGS.SIZE <- function(){
  for (i in 1:N.resp) {
    AGBS[i] ~ dnorm(mu[i], tau)
    mu[i] <- b.Shannon.Wiener[statecdid[i]] * Shannon.Wiener[i] + b.statecd[statecdid[i]]
  }
  
  for (j in 1:N.unitcd) {
    b.statecd[j] ~ dnorm(b.statecd.hat[j], tau.statecd)
    b.statecd.hat[j] <- b.statecd0 + b.SIZE0 * SIZE[j]
    b.Shannon.Wiener[j] ~ dnorm(b.Shannon.Wiener.hat[j], tau.Shannon.Wiener)
    b.Shannon.Wiener.hat[j] <- b.statecd1 + b.SIZE1 * SIZE[j]
  }
  #
  b.statecd0 ~ dnorm(0, 0.01)
  b.statecd1 ~ dnorm(0, 0.01)
  #
  b.SIZE0 ~ dnorm(0, 0.01)
  b.SIZE1 ~ dnorm(0, 0.01)
  #
  tau ~ dgamma(1, 1)
  sigma.statecd <- 1 / tau.statecd
  tau.statecd ~ dgamma(1, 1)
  tau.Shannon.Wiener ~ dgamma(1,1)
}
#
BUGS.params <- c("b.SIZE0", "b.SIZE1",
                 "b.statecd0", "b.statecd1", "b.Shannon.Wiener")
#
set.seed(123)
JAGS.SIZE <- jags(data = AGB.data.datjags, inits = NULL, parameters.to.save = BUGS.params,
                  model.file = BUGS.SIZE, n.chains = 3, n.iter = 20000, n.burnin = 5000)
#
options(max.print = 1000000)
JAGS.SIZE.2 <- update(JAGS.SIZE, n.iter = 1000)
JAGS.SIZE.2.DIC <- JAGS.SIZE.2$BUGSoutput$DIC
JAGS.SIZE.2.DIC
JAGS.SIZE.2.mcmc <- as.mcmc(JAGS.SIZE.2)
JAGS.SIZE.AGB <- JKmcmctab(JAGS.SIZE.2.mcmc)
JAGS.SIZE.AGB.diag <- gelman.diag(JAGS.SIZE.2.mcmc)
JAGS.SIZE.AGB <- cbind(JAGS.SIZE.AGB, JAGS.SIZE.AGB.diag$psrf)
JAGS.SIZE.AGB
write.csv(JAGS.SIZE.AGB, "E:\\BHM.Results\\...csv")



BUGS.SIZE <- function(){
  for (i in 1:N.resp) {
    AGBS[i] ~ dnorm(mu[i], tau)
    mu[i] <- b.Pielou[statecdid[i]] * Pielou[i] + b.statecd[statecdid[i]]
  }
  
  for (j in 1:N.unitcd) {
    b.statecd[j] ~ dnorm(b.statecd.hat[j], tau.statecd)
    b.statecd.hat[j] <- b.statecd0 + b.SIZE0 * SIZE[j]
    b.Pielou[j] ~ dnorm(b.Pielou.hat[j], tau.Pielou)
    b.Pielou.hat[j] <- b.statecd1 + b.SIZE1 * SIZE[j]
  }
  #
  b.statecd0 ~ dnorm(0, 0.01)
  b.statecd1 ~ dnorm(0, 0.01)
  #
  b.SIZE0 ~ dnorm(0, 0.01)
  b.SIZE1 ~ dnorm(0, 0.01)
  #
  tau ~ dgamma(1, 1)
  sigma.statecd <- 1 / tau.statecd
  tau.statecd ~ dgamma(1, 1)
  tau.Pielou ~ dgamma(1,1)
}
#
BUGS.params <- c("b.SIZE0", "b.SIZE1",
                 "b.statecd0", "b.statecd1", "b.Pielou")
#
set.seed(123)
JAGS.SIZE <- jags(data = AGB.data.datjags, inits = NULL, parameters.to.save = BUGS.params,
                  model.file = BUGS.SIZE, n.chains = 3, n.iter = 20000, n.burnin = 5000)
#
options(max.print = 1000000)
JAGS.SIZE.2 <- update(JAGS.SIZE, n.iter = 1000)
JAGS.SIZE.2.DIC <- JAGS.SIZE.2$BUGSoutput$DIC
JAGS.SIZE.2.DIC
JAGS.SIZE.2.mcmc <- as.mcmc(JAGS.SIZE.2)
JAGS.SIZE.AGB <- JKmcmctab(JAGS.SIZE.2.mcmc)
JAGS.SIZE.AGB.diag <- gelman.diag(JAGS.SIZE.2.mcmc)
JAGS.SIZE.AGB <- cbind(JAGS.SIZE.AGB, JAGS.SIZE.AGB.diag$psrf)
JAGS.SIZE.AGB
write.csv(JAGS.SIZE.AGB, "E:\\BHM.Results\\...csv")






BUGS.SIZE <- function(){
  for (i in 1:N.resp) {
    AGBS[i] ~ dnorm(mu[i], tau)
    mu[i] <- b.DBH.SD[statecdid[i]] * DBH.SD[i] + b.statecd[statecdid[i]]
  }
  
  for (j in 1:N.unitcd) {
    b.statecd[j] ~ dnorm(b.statecd.hat[j], tau.statecd)
    b.statecd.hat[j] <- b.statecd0 + b.SIZE0 * SIZE[j]
    b.DBH.SD[j] ~ dnorm(b.DBH.SD.hat[j], tau.DBH.SD)
    b.DBH.SD.hat[j] <- b.statecd1 + b.SIZE1 * SIZE[j]
  }
  #
  b.statecd0 ~ dnorm(0, 0.01)
  b.statecd1 ~ dnorm(0, 0.01)
  #
  b.SIZE0 ~ dnorm(0, 0.01)
  b.SIZE1 ~ dnorm(0, 0.01)
  #
  tau ~ dgamma(1, 1)
  sigma.statecd <- 1 / tau.statecd
  tau.statecd ~ dgamma(1, 1)
  tau.DBH.SD ~ dgamma(1,1)
}
#
BUGS.params <- c("b.SIZE0", "b.SIZE1",
                 "b.statecd0", "b.statecd1", "b.DBH.SD")
#
set.seed(123)
JAGS.SIZE <- jags(data = AGB.data.datjags, inits = NULL, parameters.to.save = BUGS.params,
                  model.file = BUGS.SIZE, n.chains = 3, n.iter = 20000, n.burnin = 5000)
#
options(max.print = 1000000)
JAGS.SIZE.2 <- update(JAGS.SIZE, n.iter = 1000)
JAGS.SIZE.2.DIC <- JAGS.SIZE.2$BUGSoutput$DIC
JAGS.SIZE.2.DIC
JAGS.SIZE.2.mcmc <- as.mcmc(JAGS.SIZE.2)
JAGS.SIZE.AGB <- JKmcmctab(JAGS.SIZE.2.mcmc)
JAGS.SIZE.AGB.diag <- gelman.diag(JAGS.SIZE.2.mcmc)
JAGS.SIZE.AGB <- cbind(JAGS.SIZE.AGB, JAGS.SIZE.AGB.diag$psrf)
JAGS.SIZE.AGB
write.csv(JAGS.SIZE.AGB, "E:\\BHM.Results\\...csv")







BUGS.SIZE <- function(){
  for (i in 1:N.resp) {
    AGBS[i] ~ dnorm(mu[i], tau)
    mu[i] <- b.DBH.CV[statecdid[i]] * DBH.CV[i] + b.statecd[statecdid[i]]
  }
  
  for (j in 1:N.unitcd) {
    b.statecd[j] ~ dnorm(b.statecd.hat[j], tau.statecd)
    b.statecd.hat[j] <- b.statecd0 + b.SIZE0 * SIZE[j]
    b.DBH.CV[j] ~ dnorm(b.DBH.CV.hat[j], tau.DBH.CV)
    b.DBH.CV.hat[j] <- b.statecd1 + b.SIZE1 * SIZE[j]
  }
  #
  b.statecd0 ~ dnorm(0, 0.01)
  b.statecd1 ~ dnorm(0, 0.01)
  #
  b.SIZE0 ~ dnorm(0, 0.01)
  b.SIZE1 ~ dnorm(0, 0.01)
  #
  tau ~ dgamma(1, 1)
  sigma.statecd <- 1 / tau.statecd
  tau.statecd ~ dgamma(1, 1)
  tau.DBH.CV ~ dgamma(1,1)
}
#
BUGS.params <- c("b.SIZE0", "b.SIZE1",
                 "b.statecd0", "b.statecd1", "b.DBH.CV")
#
set.seed(123)
JAGS.SIZE <- jags(data = AGB.data.datjags, inits = NULL, parameters.to.save = BUGS.params,
                  model.file = BUGS.SIZE, n.chains = 3, n.iter = 20000, n.burnin = 5000)
#
options(max.print = 1000000)
JAGS.SIZE.2 <- update(JAGS.SIZE, n.iter = 1000)
JAGS.SIZE.2.DIC <- JAGS.SIZE.2$BUGSoutput$DIC
JAGS.SIZE.2.DIC
JAGS.SIZE.2.mcmc <- as.mcmc(JAGS.SIZE.2)
JAGS.SIZE.AGB <- JKmcmctab(JAGS.SIZE.2.mcmc)
JAGS.SIZE.AGB.diag <- gelman.diag(JAGS.SIZE.2.mcmc)
JAGS.SIZE.AGB <- cbind(JAGS.SIZE.AGB, JAGS.SIZE.AGB.diag$psrf)
JAGS.SIZE.AGB
write.csv(JAGS.SIZE.AGB, "E:\\BHM.Results\\...csv")






BUGS.SIZE <- function(){
  for (i in 1:N.resp) {
    AGBS[i] ~ dnorm(mu[i], tau)
    mu[i] <- b.DBH.Gini[statecdid[i]] * DBH.Gini[i] + b.statecd[statecdid[i]]
  }
  
  for (j in 1:N.unitcd) {
    b.statecd[j] ~ dnorm(b.statecd.hat[j], tau.statecd)
    b.statecd.hat[j] <- b.statecd0 + b.SIZE0 * SIZE[j]
    b.DBH.Gini[j] ~ dnorm(b.DBH.Gini.hat[j], tau.DBH.Gini)
    b.DBH.Gini.hat[j] <- b.statecd1 + b.SIZE1 * SIZE[j]
  }
  #
  b.statecd0 ~ dnorm(0, 0.01)
  b.statecd1 ~ dnorm(0, 0.01)
  #
  b.SIZE0 ~ dnorm(0, 0.01)
  b.SIZE1 ~ dnorm(0, 0.01)
  #
  tau ~ dgamma(1, 1)
  sigma.statecd <- 1 / tau.statecd
  tau.statecd ~ dgamma(1, 1)
  tau.DBH.Gini ~ dgamma(1,1)
}
#
BUGS.params <- c("b.SIZE0", "b.SIZE1",
                 "b.statecd0", "b.statecd1", "b.DBH.Gini")
#
set.seed(123)
JAGS.SIZE <- jags(data = AGB.data.datjags, inits = NULL, parameters.to.save = BUGS.params,
                  model.file = BUGS.SIZE, n.chains = 3, n.iter = 20000, n.burnin = 5000)
#
options(max.print = 1000000)
JAGS.SIZE.2 <- update(JAGS.SIZE, n.iter = 1000)
JAGS.SIZE.2.DIC <- JAGS.SIZE.2$BUGSoutput$DIC
JAGS.SIZE.2.DIC
JAGS.SIZE.2.mcmc <- as.mcmc(JAGS.SIZE.2)
JAGS.SIZE.AGB <- JKmcmctab(JAGS.SIZE.2.mcmc)
JAGS.SIZE.AGB.diag <- gelman.diag(JAGS.SIZE.2.mcmc)
JAGS.SIZE.AGB <- cbind(JAGS.SIZE.AGB, JAGS.SIZE.AGB.diag$psrf)
JAGS.SIZE.AGB
write.csv(JAGS.SIZE.AGB, "E:\\BHM.Results\\...csv")


