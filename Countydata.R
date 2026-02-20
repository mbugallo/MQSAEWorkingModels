
########################################################################
########    MAIN SCRIPT USED IN THE APPLICACION TO CORN DATA    ########
########################################################################

# Remove all objects from the current workspace to avoid conflicts
rm(list=ls())

set.seed(123)

# Load required packages for mixed models, M-quantile regression and visualization
library(MASS)
library(rsae)
library(nlme)
library(tidyr)
library(dplyr)

# Load required functions for M-quantile regression
source('AuxFunctions.R')
source('lossfunction.R')

###########################
###     Real dataset    ###
###########################

# Construct unit-level and area-level datasets with sample sizes and aggregates

data(landsat)

landsat <- landsat %>%  mutate(regioncode = as.numeric(factor(CountyName)))

landsat <- data.frame(landsat %>% group_by(regioncode) %>% 
        mutate(nd = n(), mean.s = mean(HACorn), 
               sum.s = sum(HACorn)) %>% ungroup() %>% arrange(regioncode))

landsat.mean <- unique(landsat[, c('regioncode', 'nd', 'SegmentsInCounty', 
                                   'MeanPixelsCorn', 'MeanPixelsSoybeans', 
                                   'mean.s', 'sum.s')])
colnames(landsat.mean)[3] <- 'Nd'

landsat <- landsat %>% dplyr::select(-c('SegementID', 'HASoybeans', 'nd', 
                                        'sum.s', 'SegmentsInCounty', 
                                        'MeanPixelsCorn', 'MeanPixelsSoybeans', 
                                        'mean.s'))

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

# Define a dense and symmetric grid of quantile orders
tau <- sort(unique(c(seq(0.005,0.995,0.008), 
                     1- seq(0.005,0.995,0.008)), 0.5))


############################################
###   2. Nested error regression model   ###
############################################

# Fit the nested error regression model as a benchmark predictor

modelo.NER <- lme(fixed = y.s~x.s, random=~1|regioncode.s)
beta.ud <- coef(modelo.NER)

# Compute area-level EBLUP-type predictions under the NER model
pred.NER <- sapply(1:D, function(d){
  1/(Nd[d] + nd[d]) * (sum.s[d] + Nd[d]*x.r[regioncode.r==d, ]%*% t(beta.ud[d, ]))} 
  )
  

##################################
###   3. Unit-level MQ model   ###
##################################

# Fit unit-level M-quantile regression models over the quantile grid
mod <- QRLM(x=cbind(1,x.s), y=y.s, q=tau, maxit=30, k = 1.345)

qo <- matrix(c(gridfitinter(y.s,mod$fitted.values,
                            mod$q.values)),nrow=n,ncol=1)
qmat <- data.frame(qo, regioncode.s)

# Estimate area-specific M-quantile orders by averaging unit-level solutions
mqo <- aggregate(qmat[,1],by=list(qmat[,2]),mean)[,2]

mod.SAE <- QRLM(x=cbind(1,x.s), y=y.s, q=mqo, maxit=30, k = 1.345)

res.s <- lapply(1:D, function(d){ y.s-cbind(1,x.s)%*%mod.SAE$coef[,d]} )
sd <- sapply(res.s, fun.MAD)

# Plug-in type MQ predictor
pred.MQ <- sapply(1:D, function(d){
  1/(Nd[d] + nd[d]) * (sum.s[d] + Nd[d]*x.r[regioncode.r==d, ]%*% mod.SAE$coef[,d])} )

# Plug-in type robust bias-corrected MQ predictor
pred.BMQ <- pred.MQ + f.sum.BMQ(res.s, sd, regioncode.s, grid=2)[[1]]

############################################################
###   4. NEW PREDICTORS DERIVED FROM M-QUANTILE MODELS   ###
############################################################

# Apply robust and bias-corrected M-quantile-based predictors

# RMQU predictor
pred.RMQU <- calc_RUMQ(pred.MQ, res.s, sd, regioncode.s, mqo, Nd, nd, D)
  
# RMQG predictor
mqo.new     <- mqo.loss(x.s, y.s, mqo, Nd, nd, regioncode.s, binary = FALSE, ell = 1/2)
mod.SAE.new <- QRLM(x=cbind(1,x.s), y=y.s, q=mqo.new, maxit=30, k = 1.345)

res.s.new <- lapply(1:D, function(d){ y.s-cbind(1,x.s)%*%mod.SAE.new$coef[,d]} )
sd.new <- sapply(res.s.new, fun.MAD)

pred.RMQG <- calc_RBUMQ(sum.s, x.r, mod.SAE.new, res.s.new, sd.new, 
                         mqo.new, regioncode.s, mqo, Nd, nd, D)


###############################################
###     SUMMARY AND PLOT OF THE RESULTS     ###
###############################################

# Summarize and compare predictors across areas using tables and diagnostic plots

results <- cbind(1:D, nd, round(cbind(pred.NER, pred.MQ, pred.BMQ,
                      pred.RMQU, pred.RMQG, mqo, mqo.new),3)); results


library(xtable)

results <- cbind(1:D, nd, round(cbind(pred.NER, pred.MQ, pred.BMQ,
                                      pred.RMQU, pred.RMQG, mqo, mqo.new), 3))

colnames(results) <- c("d","nd","pred.NER","pred.MQ","pred.BMQ", "pred.RMQU",
                       "pred.RMQG","mqo","mqo.new")

print(xtable::xtable(results, digits = 3), include.rownames = FALSE)
