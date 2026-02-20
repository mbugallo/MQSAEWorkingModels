
##########################################################################
########    MAIN SCRIPT USED IN THE SIMULATIONS BASED ON REAL DATA     ###
##########################################################################

# Remove all objects from the current workspace to avoid conflicts
rm(list=ls())
set.seed(251998)

# Load required packages for mixed models and M-quantile regression
library(MASS)
library(rsae)
library(nlme)
library(dplyr)
library(numDeriv)

# Load required functions for M-quantile regression
source('AuxFunctions.R')
source('lossfunction.R')

###########################
###     Real dataset    ###
###########################

# Number of Monte Carlo replication
B = 500

# Use real auxiliary information and sampling structure from the Landsat data
data(landsat)

landsat <- landsat %>%  mutate(regioncode = as.numeric(factor(CountyName)))

landsat <- data.frame(landsat %>% group_by(regioncode) %>% 
                      mutate(nd = n(), mean.s = mean(HACorn), 
                      sum.s = sum(HACorn)) %>% ungroup() %>% arrange(regioncode))

landsat.mean <- unique(landsat[, c('regioncode', 'nd', 'SegmentsInCounty', 
                                   'MeanPixelsCorn', 'MeanPixelsSoybeans', 
                                   'mean.s', 'sum.s')])
colnames(landsat.mean)[3] <- 'Nd'

landsat <- landsat %>% dplyr::select(-c('SegementID', 'HASoybeans', 'nd', 'sum.s',
                                        'SegmentsInCounty', 'MeanPixelsCorn', 
                                        'MeanPixelsSoybeans', 'mean.s'))

y.s <- landsat$HACorn
x.s <- cbind(landsat$PixelsCorn, landsat$PixelsSoybeans)
regioncode.s <- landsat$regioncode

x.r <- cbind(1, landsat.mean$MeanPixelsCorn, landsat.mean$MeanPixelsSoybeans)
regioncode.r <- landsat.mean$regioncode

nd <- landsat.mean$nd
Nd <- landsat.mean$Nd
sum.s <- landsat.mean$sum.s
D <- length(unique(regioncode.s))
n <- sum(nd)

# Grid of candidate M-quantile orders used in the estimation step
tau <- sort(unique(c(seq(0.005,0.995,0.008), 
                     1- seq(0.005,0.995,0.008)), 0.5))

############################################
###   2. Nested error regression model   ###
############################################

modelo.NER <- lme(fixed = y.s~x.s, random=~1|regioncode.s)
beta <- modelo.NER$coefficients$fixed
beta.ud <- coef(modelo.NER)

# Estimate variance components from the nested error regression model
sigma.u <- as.numeric(VarCorr(modelo.NER)[1,2])
sigma.e <- as.numeric(VarCorr(modelo.NER)[2,2])

Yd <- pred.NER <- pred.MQ <- pred.BMQ <- pred.RUMQ <- 
  pred.RBUMQ <- mqo.save <- list()

# Monte Carlo simulation study for comparing small area predictors

for(b in 1:B){
  
  print(b)
  ##################################
  ###    0. Simulate the data    ###
  ##################################
  
  # Simulate area-level random effects and unit-level errors
  
  # Generate a new population and sample under the nested error model
  ud <- rnorm(D, mean=0, sd=sigma.u)
  edj.s <- sapply(nd, rnorm, mean=0, sd=sigma.e)
  edj.r <- rnorm(n, mean = 0, sd = sigma.e)
  
  y.s.b <- cbind(1, x.s)%*%beta + unlist(edj.s) + 
    unlist(sapply(1:D, function(d){ rep(ud[d], nd[d])  })) 
  sum.s.b <- aggregate(y.s.b, by=list(regioncode.s), sum)$V1
  
  sum.r.b <- sapply(1:D, function(d){Nd[d]*(x.r[regioncode.r==d, ]%*%beta
                                            + ud[d]) + sum(edj.r[[d]])})
  
  # Compute the true population mean for each area
  Yd[[b]] <- 1/(Nd+nd) * (sum.s.b + sum.r.b)

  # Compute NER and M-quantile-based predictors for each replication
  
  ############################################
  ###   2. Nested error regression model   ###
  ############################################
  
  modelo.NER <- lme(fixed = y.s.b~x.s, random=~1|regioncode.s)
  beta.ud <- coef(modelo.NER)
  
  pred.NER[[b]] <- sapply(1:D, function(d){
    1/(Nd[d] + nd[d]) * (sum.s.b[d] + Nd[d]*x.r[regioncode.r==d, ]%*% t(beta.ud[d, ]))} 
  )
  
  
  ##################################
  ###   3. Unit-level MQ model   ###
  ##################################
  
  mod <- QRLM(x=cbind(1,x.s), y=y.s.b, q=tau, maxit=35, k = 1.345)
  
  qo <- matrix(c(gridfitinter(y.s.b, mod$fitted.values,
                              mod$q.values)),nrow=n,ncol=1)
  qmat <- data.frame(qo, regioncode.s)
  mqo <- aggregate(qmat[,1], by=list(qmat[,2]), mean)[,2]
  mqo.save[[b]] <- mqo
  
  mod.SAE <- QRLM(x=cbind(1,x.s), y=y.s.b, q=mqo, maxit=35, k = 1.345)
  
  res.s <- lapply(1:D, function(d){ y.s.b-cbind(1,x.s)%*%mod.SAE$coef[,d]} )
  sd <- sapply(res.s, fun.MAD)
  
  # Plug-in type MQ predictor
  pred.MQ[[b]] <- sapply(1:D, function(d){
    1/(Nd[d] + nd[d]) * (sum.s.b[d] + Nd[d]*x.r[regioncode.r==d, ]%*% mod.SAE$coef[,d])} )
  
  # Plug-in type robust bias-corrected MQ predictor
  pred.BMQ[[b]] <- pred.MQ[[b]] + f.sum.BMQ(res.s, sd, regioncode.s, grid=2)[[1]]
  
  ############################################################
  ###   4. NEW PREDICTORS DERIVED FROM M-QUANTILE MODELS   ###
  ############################################################
  
  # Apply robust and adaptive bias-corrected M-quantile predictors
  
  # RUMQ predictor
  pred.RUMQ[[b]] <- calc_RUMQ(pred.MQ[[b]], res.s, sd, regioncode.s, mqo, Nd, nd, D)
  
  # RBUMQ predictor
  mqo.new     <- mqo.loss(x.s, y.s, mqo, Nd, nd, regioncode.s, binary = FALSE, ell = 1/2)
  mod.SAE.new <- QRLM(x=cbind(1,x.s), y=y.s.b, q=mqo.new, maxit=35, k = 1.345)
  
  res.s.new <- lapply(1:D, function(d){ y.s.b-cbind(1,x.s)%*%mod.SAE.new$coef[,d]} )
  sd.new <- sapply(res.s.new, fun.MAD)
  
  pred.RBUMQ[[b]] <- calc_RBUMQ(sum.s, x.r, mod.SAE.new, res.s.new, sd.new, 
                           mqo.new, regioncode.s, mqo, Nd, nd, D)
  }

#########################################################
###   PERFORMANCE OF THE POPULATION MEAN PREDICTORS   ###
#########################################################

Yd.mean <- Reduce('+', Yd)/B
mqo.mean <- Reduce('+', mqo.save)/B

# Compute the mean RBIAS and RRMSE
RBIAS.NER   <- 100 * apply(mapply('-', Yd, pred.NER),   1, mean)/Yd.mean
RBIAS.MQ    <- 100 * apply(mapply('-', Yd, pred.MQ),    1, mean)/Yd.mean
RBIAS.BMQ   <- 100 * apply(mapply('-', Yd, pred.BMQ),   1, mean)/Yd.mean
RBIAS.RUMQ  <- 100 * apply(mapply('-', Yd, pred.RUMQ),  1, mean)/Yd.mean
RBIAS.RBUMQ <- 100 * apply(mapply('-', Yd, pred.RBUMQ), 1, mean)/Yd.mean

RRMSE.NER   <- 100 * sqrt(apply(mapply('-', Yd, pred.NER)^2,   1, mean))/Yd.mean
RRMSE.MQ    <- 100 * sqrt(apply(mapply('-', Yd, pred.MQ)^2,    1, mean))/Yd.mean
RRMSE.BMQ   <- 100 * sqrt(apply(mapply('-', Yd, pred.BMQ)^2,   1, mean))/Yd.mean
RRMSE.RUMQ  <- 100 * sqrt(apply(mapply('-', Yd, pred.RUMQ)^2,  1, mean))/Yd.mean
RRMSE.RBUMQ <- 100 * sqrt(apply(mapply('-', Yd, pred.RBUMQ)^2, 1, mean))/Yd.mean

apply(abs(data.frame(RBIAS.NER, RBIAS.MQ, RBIAS.BMQ, 
                     RBIAS.RUMQ, RBIAS.RBUMQ)), 2, mean)
apply(data.frame(RRMSE.NER, RRMSE.MQ, RRMSE.BMQ, 
                 RRMSE.RUMQ, RRMSE.RBUMQ), 2, mean)


#############################
###     SOME BOXPLOTS     ###
#############################

# Graphical comparison of the finite-sample performance of the estimators

boxplot(RBIAS.NER, RBIAS.MQ, RBIAS.BMQ, RBIAS.RUMQ, RBIAS.RBUMQ, 
        las=2, ylab = "", xlab = "", pch=19, main='RBIAS (%)', cex.axis=1.4,
        cex.lab=1.4, cex.main=1.4, ylim=c(-1.1, 1.1), 
        names = c("ner", "mq", "mqbc", "rmqu", "rmqg"))
abline(h=0, col='red', lty=2, lwd=2)

boxplot(RRMSE.NER, RRMSE.MQ, RRMSE.BMQ, RRMSE.RUMQ, RRMSE.RBUMQ, 
        names = c("ner", "mq", "mqbc", "rmqu", "rmqg"), las=2, 
        ylab = "", xlab = "", pch=19, ylim=c(0, 17), main='RRMSE (%)',
        cex.axis=1.4, cex.lab=1.4, cex.main=1.4)



#############################
###      SOME  PLOTS      ###
#############################

cols <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854")
do.RRMSE <- 0

if( do.RRMSE==1 ){
  df_mean <- data.frame(nd, RRMSE.NER, RRMSE.MQ, RRMSE.BMQ, RRMSE.RUMQ, RRMSE.RBUMQ) %>%
    group_by(nd) %>% summarise(across(RRMSE.NER: RRMSE.RBUMQ, mean))
  
  plot(df_mean$nd, df_mean$RRMSE.NER, type='o', col=cols[1],
       xlab=expression(n[d]), ylab='Average RRMSE', lwd=2, pch=16, ylim=c(4,16))
  lines(df_mean$nd, df_mean$RRMSE.MQ, type='o', col=cols[2], lwd=2, pch=14)
  lines(df_mean$nd, df_mean$RRMSE.BMQ, type='o', col=cols[3], lwd=2, pch=17)
  lines(df_mean$nd, df_mean$RRMSE.RUMQ, type='o', col=cols[4], lwd=2, pch=18)
  lines(df_mean$nd, df_mean$RRMSE.RBUMQ, type='o', col=cols[5], lwd=2, pch=15)
  
} else {
  df_mean <- data.frame(nd, RBIAS.NER, RBIAS.MQ, RBIAS.BMQ, RBIAS.RUMQ, RBIAS.RBUMQ) %>%
    group_by(nd) %>% summarise(across(RBIAS.NER: RBIAS.RBUMQ, mean))

  plot(df_mean$nd, df_mean$RBIAS.NER, type='o', col=cols[1],
       xlab=expression(n[d]), ylab='Average ARBIAS', lwd=2, pch=16, ylim=c(-1, 1))
  abline(a=0, b=0, lwd=2)
  lines(df_mean$nd, df_mean$RBIAS.MQ, type='o', col=cols[2], lwd=2, pch=14)
  lines(df_mean$nd, df_mean$RBIAS.BMQ, type='o', col=cols[3], lwd=2, pch=17)
  lines(df_mean$nd, df_mean$RBIAS.RUMQ, type='o', col=cols[4], lwd=2, pch=18)
  lines(df_mean$nd, df_mean$RBIAS.RBUMQ, type='o', col=cols[5], lwd=2, pch=15)
}

legend("topright", legend=c("NER", "MQ", "MQBC", "RMQU", "RMQG"),
       col=cols, lty=1, lwd=2, pch=c(16,17,18,15), cex=.75)



