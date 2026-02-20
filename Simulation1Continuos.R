
######################################
###    Simulation experiments      ###
######################################

# Monte Carlo simulation for unit-level SAE under nested error models with contamination

rm(list=ls())
set.seed(123)

library(MASS)
library(dplyr)
library(sae)
library(nlme)
library(numDeriv)

source('AuxFunctions.R')
source('lossfunction.R')

###########################
### Simulation scenario ###
###########################

B=200 # number of iterations
D=600  # number of areas

Nd=rep(100,D)
N=sum(Nd)
n=N*0.05 # p
regioncode=rep(1:D,Nd)

# Continuous Part
xdj1 <- rlnorm(N,  1, 0.5)

# Elevation factors
fdj <- runif(N, 1, 5)       
pdj <- n * fdj / sum(fdj)   
wdj <- 1/pdj

# Sample and non-sample subsets
s_list <- lapply(split(1:N, regioncode), function(idx) {
  sample(idx, 1, prob = pdj[idx]) # guarantee at least 1 from each region
})

# remaining units sampled as before
s <- sort(c(unlist(s_list), sample(setdiff(1:N, unlist(s_list)), 
            n - length(s_list), prob = pdj[setdiff(1:N, unlist(s_list))])))

x.s<-xdj1[s]
regioncode.s<-regioncode[s]
wdj.s<-wdj[s]

x.r<-xdj1[-s]
regioncode.r<-regioncode[-s]
wdj.r<-wdj[-s]

nd=as.numeric(table(regioncode.s))

# Quantiles on the grid
tau <- sort(unique(c(seq(0.005,0.995,0.008), 
                     1- seq(0.005,0.995,0.008)), 0.5))
grid <- seq(0, 5, by=0.01)

# Outliers 
mean.ui<-9
mean.e<-20
ki<-as.integer(.1 * D) 

res.s <- pred.Hajek <- pred.BMQ <- pred.NER <- pred.MQ <- C.MQ <- C.UMQ <-
  C.BUMQ <- pred.UMQ <- pred.BUMQ <- Yd <- asym.MSE.results <- 
  mqo <- mqo.new <- mqo.new.trim.50 <- list()


#############################
### Start the simulations ###
#############################

for (b in 1:B){
  set.seed(b)
  print(b)
  
  # Outliers 
  ud <- rnorm(D,0,sqrt(3))
  # % contamination in ui
  out.ui <- 1
  uui<-ud
  u1 <- rnorm(ki, mean.ui, sqrt(20))
  uui[(D-ki+1):D]<-u1
  out.ui <- rep(out.ui, D)
  ud <- ifelse(out.ui > 0, uui, ud)
  
  # % contamination in e
  out.e <- 0.03
  n1<-rbinom(N,1,out.e)
  edj <- (1-n1)*rnorm(N,0,sqrt(6))+n1*rnorm(N, mean.e, sqrt(150))
  
  ydj <- 100 + 5*xdj1 + rep(ud,Nd) + edj
  
  Yd[[b]] <- aggregate(ydj, by=list(regioncode), mean)[, 2]
  
  # Sample observations
  y.s <- ydj[s]
  data.s <- data.frame(regioncode.s, y.s, x.s, wdj.s)
  sum.s <- aggregate(y.s, by=list(regioncode.s), sum)$x
  
  #####################################
  ###   1. Hajek direct estimates   ###
  #####################################
  
  Y.d.Hajek <- as.data.frame(data.s %>% mutate(w = wdj.s) %>% 
               group_by(regioncode.s) %>% 
               summarise(Y.dir = sum(y.s * w) / sum(w),
                         Y.var = sum(w * (w - 1) * (y.s - Y.dir)^2) / (sum(w)^2),
                         Nd.hat =(sum(w)^2),  .groups = "drop")) %>%
               mutate(Y.var = ifelse(Y.var == 0, min(Y.var[Y.var != 0]), Y.var))
  pred.Hajek[[b]] <- Y.d.Hajek$Y.dir

  ############################################
  ###   2. Nested error regression model   ###
  ############################################
  
  # True population generated from a nested error regression model
  modelo.NER <- lme(fixed = y.s~x.s, random=~1|regioncode.s, data = data.s)
  
  # EBLUP-type predictor based on lme() (NER model)
  pred.NER[[b]] <- 1/Nd *(sum.s + sapply(1:D, function(d) 
        sum(as.matrix(cbind(1,x.r)[regioncode.r==d, ]) %*% 
              t(as.matrix(coef(modelo.NER)[d, ])))))
  
  ##################################
  ###   3. Unit-level MQ model   ###
  ##################################
  
  # Unit-level M-quantile regression fitted via QRLM()
  mod <- QRLM(x=cbind(1,x.s), y=y.s, q=tau, maxit=30, k = 1.345)
  qo <- matrix(c(gridfitinter(y.s,mod$fitted.values,
                              mod$q.values)),nrow=n,ncol=1)
  qmat <- data.frame(qo, regioncode.s)
  mqo[[b]] <- aggregate(qmat[,1],by=list(qmat[,2]),mean)[,2]
  
  mod.SAE <- QRLM(x=cbind(1,x.s), y=y.s, q=mqo[[b]], maxit=30, k = 1.345)
  res.s <- lapply(1:D, function(d){ y.s-cbind(1,x.s)%*%mod.SAE$coef[,d]} )
  sd <- sapply(res.s, fun.MAD)
  
  # Plug-in M-quantile predictor for area means (MQ)
  pred.MQ[[b]] <- 1/Nd *(sum.s + sapply(1:D, function(d) 
    sum(as.matrix(cbind(1,x.r)[regioncode.r==d, ]) %*% mod.SAE$coef[,d])))
  
  # Robust bias-corrected M-quantile predictor (BMQ)
  pred.BMQ[[b]] <- pred.MQ[[b]] + f.sum.BMQ(res.s, sd, regioncode.s, grid)[[1]]
  

  ############################################################
  ###   4. NEW PREDICTORS DERIVED FROM M-QUANTILE MODELS   ###
  ############################################################
  
  # NEW Plug-in type robust bias-corrected MQ predictor
  UMQ.aux <- f.sum.UMQ(res.s, sd, regioncode.s, mqo[[b]], grid, Nd, nd)
  pred.UMQ[[b]] <- pred.MQ[[b]] + UMQ.aux[[1]]
  
  mqo.new[[b]] <- mqo.loss(x.s, y.s, mqo[[b]], Nd, nd, regioncode.s, 
                           binary = FALSE, ell = 0.5)
  
  mod.SAE.new <- QRLM(x=cbind(1,x.s), y=y.s, q=mqo.new[[b]], maxit=30, k = 1.345)
  
  res.s.new <- lapply(1:D, function(d){ y.s-cbind(1,x.s)%*%mod.SAE.new$coef[,d]} )
  sd.new <- sapply(res.s.new, fun.MAD)
  
  # Robust unbiased M-quantile predictor (RUMQ)
  BUMQ.aux <- f.sum.UMQ(res.s.new, sd.new, regioncode.s, mqo.new[[b]], grid, Nd, nd)
  pred.BUMQ[[b]] <- 1/Nd *(sum.s + sapply(1:D, function(d) 
    sum(as.matrix(cbind(1,x.r)[regioncode.r==d, ]) %*% mod.SAE.new$coef[,d]))) + BUMQ.aux[[1]]

  ########################################
  ###   5. Asymptotic MSE estimation   ###
  ########################################
  
  asym.MSE.results[[b]] <- asym.MSE(x.s, y.s, mqo[[b]], x.r, Nd, nd, regioncode.s, 
                                    mod.SAE, res.s, sd, UMQ.aux[[2]])
  
  C.BUMQ[[b]] <- sum((Yd[[b]] > pred.BUMQ[[b]]-qnorm(1-0.05/2)*sqrt(asym.MSE.results[[b]]))&
                     (Yd[[b]] < pred.BUMQ[[b]]+qnorm(1-0.05/2)*sqrt(asym.MSE.results[[b]])))/D
  
  C.UMQ[[b]] <- sum((Yd[[b]] > pred.UMQ[[b]]-qnorm(1-0.05/2)*sqrt(asym.MSE.results[[b]]))&
         (Yd[[b]] < pred.UMQ[[b]]+qnorm(1-0.05/2)*sqrt(asym.MSE.results[[b]])))/D
  
  C.MQ[[b]] <- sum((Yd[[b]] > pred.MQ[[b]]-qnorm(1-0.05/2)*sqrt(asym.MSE.results[[b]]))&
                   (Yd[[b]] < pred.MQ[[b]]+qnorm(1-0.05/2)*sqrt(asym.MSE.results[[b]])))/D
  
}


#########################################################
###   PERFORMANCE OF THE POPULATION MEAN PREDICTORS   ###
#########################################################

# Coverage of the ICs
mean(unlist(C.MQ)); mean(unlist(C.UMQ)); mean(unlist(C.BUMQ))

Yd.mean <- Reduce('+', Yd)/B

rmse <- function(pred_list, Yd, B) {
  sq_errors <- sapply(1:B, function(b) (pred_list[[b]] - Yd[[b]])^2)
  sqrt(apply(sq_errors, 1, mean))
}

rrmse <- function(pred_list, Yd, Yd_mean, B) {
  sq_errors <- sapply(1:B, function(b) (pred_list[[b]] - Yd[[b]])^2)
  sqrt(apply(sq_errors, 1, median)) / Yd_mean
}

rbias <- function(pred_list, Yd, Yd_mean, B) {
  sq_bias <- sapply(1:B, function(b) (pred_list[[b]] - Yd[[b]]) )
  (apply(sq_bias, 1, median)) / Yd_mean
}

arbias <- function(pred_list, Yd, Yd_mean, B) {
  sq_bias <- sapply(1:B, function(b) abs(pred_list[[b]] - Yd[[b]]) )
  (apply(sq_bias, 1, median)) / Yd_mean
}


# Monte Carlo evaluation using relative bias and relative RMSE

rmse.MQ  <- rmse(pred.MQ,  Yd,  B)
rmse.UMQ <- rmse(pred.UMQ, Yd,  B)

rrmse.Hajek <- rrmse(pred.Hajek, Yd, Yd.mean, B)
rrmse.NER   <- rrmse(pred.NER, Yd, Yd.mean, B)
rrmse.MQ    <- rrmse(pred.MQ, Yd, Yd.mean, B)
rrmse.BMQ   <- rrmse(pred.BMQ, Yd, Yd.mean, B)
rrmse.UMQ   <- rrmse(pred.UMQ, Yd, Yd.mean, B)
rrmse.BUMQ  <- rrmse(pred.BUMQ, Yd, Yd.mean, B)
  
arbias.Hajek <- arbias(pred.Hajek, Yd, Yd.mean, B)
arbias.NER   <- arbias(pred.NER, Yd, Yd.mean, B)
arbias.MQ    <- arbias(pred.MQ, Yd, Yd.mean, B)
arbias.BMQ   <- arbias(pred.BMQ, Yd, Yd.mean, B)
arbias.UMQ   <- arbias(pred.UMQ, Yd, Yd.mean, B)
arbias.BUMQ  <- arbias(pred.BUMQ, Yd, Yd.mean, B)

rbias.Hajek <- rbias(pred.Hajek, Yd, Yd.mean, B)
rbias.NER   <- rbias(pred.NER, Yd, Yd.mean, B)
rbias.MQ    <- rbias(pred.MQ, Yd, Yd.mean, B)
rbias.BMQ   <- rbias(pred.BMQ, Yd, Yd.mean, B)
rbias.UMQ   <- rbias(pred.UMQ, Yd, Yd.mean, B)
rbias.BUMQ  <- rbias(pred.BUMQ, Yd, Yd.mean, B)

####################
###   Boxplots   ###
####################

box.rrmse  <- data.frame(rrmse.NER, rrmse.MQ, rrmse.BMQ, rrmse.UMQ, rrmse.BUMQ)
box.rbias <- data.frame(rbias.NER,  rbias.MQ, rbias.BMQ, rbias.UMQ, rbias.BUMQ)
box.arbias <- data.frame(arbias.NER, arbias.MQ, arbias.BMQ, arbias.UMQ, arbias.BUMQ)

round(colMeans(100*box.arbias),3);round(colMeans(100*box.rrmse),3)

boxplot(100*box.rrmse, names = c("ner", "mq", "bmq", "rmqu", "rmqg"),
        las=2, ylab = "", xlab = "", pch=4, ylim=c(0,1), 
        cex.axis = 1.35, cex.lab = 1.35, main='RRMSE (%) [0,0]')

boxplot(100*box.rbias, names = c("ner", "mq", "bmq", "rmqu", "rmqg"),
        las=2, ylab = "", xlab = "", pch=4,  main='RBIAS (%) [0,0]',
        cex.axis = 1.35, cex.lab = 1.35)
abline(h=0, col='red', lty=2, lwd=2)

# write.csv(box.rrmse,'results/boxAeu.csv')
# write.csv(box.rbias,'results/biasAeu.csv')

# write.csv(cbind(nd, box.rrmse), 'results/box15eu.csv')
# write.csv(cbind(nd, box.rbias),'results/bias15eu.csv')

########################
###   Other plots 1  ###
########################

p1  <- read.csv('results/box11eu.csv')
p2  <- read.csv('results/box12eu.csv')
p3  <- read.csv('results/box13eu.csv')
p4  <- read.csv('results/box14eu.csv')
p5  <- read.csv('results/box15eu.csv')

# Lista con los objetos
plist <- list(p1, p2, p3, p4, p5)
Nd_vals <- rep(100, 5)
D_vals <- c(50, 100, 200, 400, 600)
k_vals <-as.integer(.1 * D_vals)

# Configuración de la ventana gráfica
par(mfrow = c(3, 5),
    mar = c(3, 3, 2, 1),
    oma = c(0, 0, 0, 0),
    mgp = c(1.6, 0.5, 0))

# Bucle para dibujar los 5 plots
for (i in seq_along(plist)) {
  
  plot(plist[[i]]$nd[1:(D_vals[i]-k_vals[i])],
       (100*plist[[i]]$rrmse.MQ- 100*plist[[i]]$rrmse.BUMQ)[1:(D_vals[i]-k_vals[i])],
       cex.axis = 1.35, cex.lab = 1.35,   xlim=c(1, 12), ylim=c(-0.6, 0.7), #  ylim=c(-0.075, 0.8),
       main = bquote(D == .(D_vals[i])~ "," ~ N[d] == .(Nd_vals[i])),
       xlab = expression(n[d]), ylab = "MQ-RMQG")
     points(plist[[i]]$nd[(D_vals[i]-k_vals[i]+1):D_vals[i]], 
            0.6* (100*plist[[i]]$rrmse.MQ- 100*plist[[i]]$rrmse.BUMQ)[(D_vals[i]-k_vals[i]+1):D_vals[i]],
            col='blue', pch=4, cex=1.2)
  abline(h=0, col = "red", lty = 2, lwd=2)
}

for (i in seq_along(plist)) {
  
  plot(plist[[i]]$nd[1:(D_vals[i]-k_vals[i])], 
       (100*plist[[i]]$rrmse.BMQ-100*plist[[i]]$rrmse.BUMQ)[1:(D_vals[i]-k_vals[i])],
       cex.axis = 1.35, cex.lab = 1.35,  xlim=c(1, 12), ylim=c(-0.6, 0.7), # ylim=c(-0.075, 0.8),
       main = bquote(D == .(D_vals[i])~ "," ~ N[d] == .(Nd_vals[i])),
       xlab = expression(n[d]), ylab = "MQBC-RMQG")
     points(plist[[i]]$nd[(D_vals[i]-k_vals[i]+1):D_vals[i]], 
            0.6*(100*plist[[i]]$rrmse.BMQ- 100*plist[[i]]$rrmse.BUMQ)[(D_vals[i]-k_vals[i]+1):D_vals[i]],
          col='blue', pch=4, cex=1.2)
  abline(h=0, col = "red", lty = 2, lwd=2)
}

for (i in seq_along(plist)) {
  
  plot(plist[[i]]$nd[1:(D_vals[i]-k_vals[i])], 
       (100*plist[[i]]$rrmse.UMQ- 100*plist[[i]]$rrmse.BUMQ)[1:(D_vals[i]-k_vals[i])],
       cex.axis = 1.35, cex.lab = 1.35, xlim=c(1, 12), ylim=c(-0.6, 0.7), # ylim=c(-0.075, 0.8),
       main = bquote(D == .(D_vals[i])~ "," ~ N[d] == .(Nd_vals[i])),
       xlab =  expression(n[d]), ylab = "RMQU-RMQG")
   points(plist[[i]]$nd[(D_vals[i]-k_vals[i]+1):D_vals[i]], 
          0.6*(100*plist[[i]]$rrmse.UMQ- 100*plist[[i]]$rrmse.BUMQ)[(D_vals[i]-k_vals[i]+1):D_vals[i]],
          col='blue', pch=4, cex=1.2)
  abline(h=0, col = "red", lty = 2, lwd=2)
}

########################
###   Other plots 2  ###
########################

p1  <- read.csv('results/bias11eu.csv')
p2  <- read.csv('results/bias12eu.csv')
p3  <- read.csv('results/bias13eu.csv')
p4  <- read.csv('results/bias14eu.csv')
p5  <- read.csv('results/bias15eu.csv')

# Lista con los objetos
plist <- list(p1, p2, p3, p4, p5)

# Configuración de la ventana gráfica
par(mfrow = c(3, 5),
    mar = c(3, 3, 2, 1),
    oma = c(0, 0, 0, 0),
    mgp = c(1.6, 0.5, 0))

# Bucle para dibujar los 5 plots
for (i in seq_along(plist)) {
  
  plot(plist[[i]]$nd[1:(D_vals[i]-k_vals[i])], 
       (abs(100*plist[[i]]$rbias.MQ)- abs(100*plist[[i]]$rbias.BUMQ))[1:(D_vals[i]-k_vals[i])],
       cex.axis = 1.35, cex.lab = 1.35, xlim=c(1,12), ylim=c(-0.7, 0.5), # ylim=c(-0.12, 0.12),
       main = bquote(D == .(D_vals[i])~ "," ~ N[d] == .(Nd_vals[i])),
       xlab = expression(n[d]), ylab = "MQ-RMQG")
   points(plist[[i]]$nd[(D_vals[i]-k_vals[i]+1):D_vals[i]],
          (abs(100*plist[[i]]$rbias.MQ)- abs(100*plist[[i]]$rbias.BUMQ))[(D_vals[i]-k_vals[i]+1):D_vals[i]], 
           col='blue', pch=4, cex=1.2)
  abline(h=0, col = "red", lty = 2, lwd=2)
}

for (i in seq_along(plist)) {
  
  plot(plist[[i]]$nd[1:(D_vals[i]-k_vals[i])], 
       (abs(100*plist[[i]]$rbias.BMQ)- abs(100*plist[[i]]$rbias.BUMQ))[1:(D_vals[i]-k_vals[i])],
       cex.axis = 1.35, cex.lab = 1.35, xlim=c(1,12), ylim=c(-0.7, 0.5), # ylim=c(-0.12, 0.12),
       main = bquote(D == .(D_vals[i])~ "," ~ N[d] == .(Nd_vals[i])),
       xlab = expression(n[d]), ylab = "MQBC-RMQG")
    points(plist[[i]]$nd[(D_vals[i]-k_vals[i]+1):D_vals[i]],
           (abs(100*plist[[i]]$rbias.BMQ)- abs(100*plist[[i]]$rbias.BUMQ))[(D_vals[i]-k_vals[i]+1):D_vals[i]],
           col='blue', pch=4, cex=1.2)
  abline(h=0, col = "red", lty = 2, lwd=2)
}

for (i in seq_along(plist)) {
  
  plot(plist[[i]]$nd[1:(D_vals[i]-k_vals[i])], 
       (abs(100*plist[[i]]$rbias.UMQ)- abs(100*plist[[i]]$rbias.BUMQ))[1:(D_vals[i]-k_vals[i])],
       cex.axis = 1.35, cex.lab = 1.35, xlim=c(1,12), ylim=c(-0.7, 0.5), # ylim=c(-0.12, 0.12),
       main = bquote(D == .(D_vals[i])~ "," ~ N[d] == .(Nd_vals[i])),
       xlab = expression(n[d]), ylab = "RMQU-RMQG")
   points(plist[[i]]$nd[(D_vals[i]-k_vals[i]+1):D_vals[i]],
          (abs(100*plist[[i]]$rbias.UMQ)- abs(100*plist[[i]]$rbias.BUMQ))[(D_vals[i]-k_vals[i]+1):D_vals[i]], 
           col='blue', pch=4, cex=1.2)
  abline(h=0, col = "red", lty = 2, lwd=2)
}

########################
###   Other plots 3  ###
########################

# boxB00 <- read.csv('results/boxA00.csv')
# boxBe0 <- read.csv('results/boxAe0.csv')
# boxB0u <- read.csv('results/boxA0u.csv')
# boxBeu <- read.csv('results/boxAeu.csv')

boxB00 <- read.csv('results/biasA00.csv')
boxBe0 <- read.csv('results/biasAe0.csv')
boxB0u <- read.csv('results/biasA0u.csv')
boxBeu <- read.csv('results/biasAeu.csv')

round(apply(100*boxB00[, 2:6], 2, median), 3)

par(mfrow = c(1, 4),
    mar = c(4, 3, 2, 1),
    oma = c(0, 0, 0, 0),
    mgp = c(1.6, 0.5, 0))

boxplot(100*boxB00[, 2:6], names = c("ner", "mq", "bmq", "rumq", "rbumq"),
        las=2, ylab = "", xlab = "", pch=4, ylim=c(-1,0.5), cex.main=1.2,
        cex.axis = 1.4, cex.lab = 1.4, main='RBIAS (%) [0,0]')
abline(a=0, b=0, lwd=2)

boxplot(100*boxBe0[, 2:6], names = c("ner", "mq", "bmq", "rumq", "rbumq"),
        las=2, ylab = "", xlab = "", pch=4, ylim=c(-1,0.5),  cex.main=1.2,
        cex.axis = 1.4, cex.lab = 1.4, main='RBIAS (%) [e,0]')
abline(a=0, b=0, lwd=2)

boxplot(100*boxB0u[, 2:6], names = c("ner", "mq", "bmq", "rumq", "rbumq"),
        las=2, ylab = "", xlab = "", pch=4, ylim=c(-1,0.5), cex.main=1.2,
        cex.axis = 1.4, cex.lab = 1.4, main='RBIAS (%) [0,u]')
abline(a=0, b=0, lwd=2)

boxplot(100*boxBeu[, 2:6], names = c("ner", "mq", "bmq", "rumq", "rbumq"),
        las=2, ylab = "", xlab = "", pch=4, ylim=c(-1,0.5), cex.main=1.2,
        cex.axis = 1.4, cex.lab = 1.4, main='RBIAS (%) [e,u]') 
abline(a=0, b=0, lwd=2)

########################
###   Other plots 5  ###
########################

par(mfrow = c(1, 3),
    mar = c(4, 2, 2, 2),
    oma = c(0, 0, 0, 0),
    mgp = c(1.6, 0.5, 0))

grid <- seq(0.005, 0.995, 0.001)

plot(grid, xlab='q', main='Mean GALI(0,1,q)', ylim=c(-70,70),
     cex.axis = 1.35, cex.lab = 1.35, ylab='', lwd=2, lty=2,
     sapply(grid, function(q){(1-2*q)/(q*(1-q))}), type='l')
abline(h=0, lwd=2, col='red', lty=2)
lines(grid, col='blue', lwd=2, 
      sapply(grid, function(q){meanAH(q, m=0, s=1, k=1.345)}))

plot(grid, xlab='q', main='Variance GALI(0,1,q)',
     cex.axis = 1.35, cex.lab = 1.35, ylim=c(0,900),ylab='', lwd=2, lty=2,
     sapply(grid, function(q){((2*q^2-2*q+1)/(q^2*(1-q)^2))}), type='l')
lines(grid, col='blue', lwd=2,
      sapply(grid, function(q){varAH(q, m=0, s=1, k=1.345)$variance}))

plot(grid, xlab='q',main='CV GALI(0,1,q)',
     cex.axis = 1.35, cex.lab = 1.35,ylab='', lwd=2, lty=2,
     sapply(grid, function(q){sqrt((2*q^2-2*q+1)/(q^2*(1-q)^2))})/abs(sapply(grid, 
                function(q){(1-2*q)/(q*(1-q))})), type='l')
lines(grid, col='blue', lwd=2, sapply(grid, 
            function(q){sqrt(varAH(q, m=0, s=1, k=1.345)$variance)/abs(meanAH(q, m=0, s=1, k=1.345))}))

