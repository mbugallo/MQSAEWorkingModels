
###############################################################
###    Defining functions and loading from other sources    ###
###############################################################

# Monte Carlo simulation study for binary small area estimation using M-quantiles

# Clear workspace to ensure reproducibility
rm(list=ls())

# Load functions for binary M-quantile regression and robust score estimation
source('BinaryMQ/R/glm.mq.binom.R')
source('BinaryMQ/R/glm.rob.binom.R')
source('BinaryMQ/R/QLogit.Scores.R')
source('BinaryMQ/R/QSCORE.R')
source('lossfunction.R')
source('AuxFunctions.R')

# Suppress warnings to avoid excessive console output during simulations
options(warn = -1)

######################################
###    Simulation experiments      ###
######################################

set.seed(1234)

# Number of Monte Carlo replications, areas, and population sizes
B <- 500 # number of iterations
D <- 40  # number of areas

Nd <- rep(100,D)
N <- sum(Nd)
n <- N*0.05 # p
N.l.n <- N - n
regioncode <- rep(1:D,Nd)

# Generate population-level auxiliary covariates
Xpop  <- cbind(1, rnorm(N, mean = 0.5, sd = 0.5), rnorm(N, mean = 1, sd = 0.5))

# Unequal probability sampling with area-specific inclusion probabilities
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

X <- Xpop[s, ]
regioncode.s <- regioncode[s]
wdj.s <- wdj[s]

X.r <- Xpop[-s, ]
regioncode.r <- regioncode[-s]
wdj.r <- wdj[-s]

nd <- as.numeric(table(regioncode.s))
Nd <- as.numeric(table(regioncode))
Nd.l.nd <- Nd - nd

# Quantiles on the grid
tau <- sort(unique(c(seq(0.025,0.975,0.02), 
                     1- seq(0.025,0.975,0.02)), 0.5))

beta <- c(0.4, 0.8, -0.8)
area.out <- 0
m <- 4

p <- p.AG <- p.GLBE <- p.SAE.old <- p.SAE.new <- mqo1 <- mqo2 <- list()

for(b in 1:B){
  print(b)
  
  # Area-level random effects with optional outlying areas
  ud <- rnorm(D, 0, 1)
  if( area.out==1 ){ ud[((D-m+1):D)] <- rnorm(m, 4, 2) } # 3, 1
  
  # Sample data
  eta.s  <- X %*% beta + unlist(sapply(1:D, function(d) rep(ud[d], nd[d])))
  p.s    <- 1 / (1 + exp(-eta.s))
  
  # Generate binary responses from a logistic mixed model
  T.s    <- rbinom(n, 1, p.s)
  
  # Misclassification in T.s
  misprop <- 0.1
  mis.idx <- sample(1:n, round(n * misprop))
  T.s[mis.idx] <- 1 - T.s[mis.idx]
  
  # Population data
  eta.r  <- X.r %*% beta + unlist(sapply(1:D, function(d) rep(ud[d], Nd.l.nd[d]))) 
  
  # Compute true population probabilities
  p.r    <- 1 / (1 + exp(-eta.r))
  Y.r <- rbinom(N.l.n, 1, p.r)
  
  # True area-level proportions
  p.AG[[b]] <- aggregate(c(T.s, Y.r), by = list(c(regioncode.s, regioncode.r)),  mean)[,2]
  
  
  ##################################################
  ###    Logistic mixed model for binary data    ###
  ##################################################
  
  # Logistic mixed-effects model fitted as a benchmark SAE method
  data <- data.frame(T.s, x1 = X[, 2],  x2 = X[, 3], region = factor(regioncode.s))
  
  fit.GLBE <- lme4::glmer(T.s ~ x1 + x2 + (1 | region), data = data,
                          family = binomial(link = "logit"))
  
  # Area-level predictions from the logistic mixed model
  p.GLBE[[b]] <- 1/Nd * (aggregate(T.s, by=list(regioncode.s), sum)[,2] + 
                    aggregate(predict(fit.GLBE, data.frame(x1 = X.r[, 2],  
                      x2 = X.r[, 3], region = regioncode.r), type="response"), 
                                     by=list(regioncode.r), sum)[, 2])
  
  
  #######################################################
  ###    OLD SAE M-quantile models for binary data    ###
  #######################################################
  # Classical score-based M-quantile order estimation
  
  temp <- QLogit.Scores(y=T.s, x=X, qgrid = tau, k.value=1.345)
  scores <- (cbind(regioncode.s, temp$qscores))
  mqo1[[b]] <- aggregate(scores[,2], by=list(scores[,1]), mean)[,2]
  
  fit.SAE.old <- lapply(mqo1[[b]], function(q) try(glm.mq.binom(x=X, y=T.s, q = q), silent = T))
  eta.SAE.old <- unlist(lapply(1:D, function(d) 
    X[regioncode.s==d, ]%*%fit.SAE.old[[d]]$coefficients))
  
  eta.SAE.old.r <- unlist(lapply(1:D, function(d) 
    X.r[regioncode.r==d, ]%*%fit.SAE.old[[d]]$coefficients))
  p.SAE.old.r    <- 1 / (1 + exp(-eta.SAE.old.r))
  
  # Area-level predictions from the original M-quantile SAE approach
  p.SAE.old[[b]] <- 1/Nd * (aggregate(T.s, by=list(regioncode.s), sum)[,2] + 
                                 aggregate(p.SAE.old.r, by=list(regioncode.r), sum)[, 2])
  
  
  #######################################################
  ###    NEW SAE M-quantile models for binary data    ###
  #######################################################

  mqo2[[b]] <- mqo.loss(X[,-1], T.s, mqo1[[b]], Nd, nd, regioncode.s, binary=TRUE, ell=1/2)
  
  fit.SAE <- lapply(mqo2[[b]], function(q) try(glm.mq.binom(x=X, y=T.s, q = q), silent = T))
  eta.SAE <- unlist(lapply(1:D, function(d) 
    X[regioncode.s==d, ]%*%fit.SAE[[d]]$coefficients))
  
  eta.SAE.r <- unlist(lapply(1:D, function(d) 
    X.r[regioncode.r==d, ]%*%fit.SAE[[d]]$coefficients))
  p.SAE.r    <- 1 / (1 + exp(-eta.SAE.r))
  
  p.SAE.new[[b]] <- 1/Nd * (aggregate(T.s, by=list(regioncode.s), sum)[,2] + 
                           aggregate(p.SAE.r, by=list(regioncode.r), sum)[, 2])
  
}
 
#####################################
###    Average prediction of p    ###
#####################################  

rmse <- function(pred_list, Yd, B) {
  sq_errors <- sapply(1:B, function(b) (pred_list[[b]] - Yd[[b]])^2)
  sqrt(apply(sq_errors, 1, mean))
}


bias <- function(pred_list, Yd, B) {
  sq_bias <- sapply(1:B, function(b) (pred_list[[b]] - Yd[[b]]) )
  (apply(sq_bias, 1, mean)) 
}

# Monte Carlo estimation of area-level bias and RMSE

S1  <- bias(p.GLBE,  p.AG,  B)
S2  <- bias(p.SAE.old,  p.AG,  B)
S3  <- bias(p.SAE.new,  p.AG,  B)

V1  <- rmse(p.GLBE,  p.AG,  B)
V2  <- rmse(p.SAE.old,  p.AG,  B)
V3  <- rmse(p.SAE.new,  p.AG,  B)

box.bias  <- data.frame(S1, S2, S3)
box.rmse  <- data.frame(V1, V2, V3)

boxplot(box.bias); abline(a=0, b=0)
boxplot(box.rmse)

write.csv(box.rmse,'results/rmseGABEe0.csv')
write.csv(box.bias,'results/biasGABEe0.csv')

##############################
###   Some boxplots RBIAS  ### # 1200 * 300
##############################

# Graphical comparison of estimator performance across scenarios

ylim <- 100*c(-0.05, 0.05)

boxbiasB00 <- read.csv('results/biasGABE00.csv')[, -1]
boxbiasBe0 <- read.csv('results/biasGABEe0.csv')[, -1]
boxbiasB0u <- read.csv('results/biasGABE0u.csv')[, -1]
boxbiasBeu <- read.csv('results/biasGABEeu.csv')[, -1]


par(mfrow = c(1, 4),
    mar = c(4, 3, 2, 1),
    oma = c(0, 0, 0, 0),
    mgp = c(1.6, 0.5, 0))

# RBIAS
boxplot( 100*boxbiasB00, names = c("mlp ", "mq", "gmq"),
        las=2, ylab = "", xlab = "", pch=19,  ylim=ylim, 
        cex.axis = 1.45, cex.lab = 1.45, main='100 x BIAS  [0,0]')
abline(a=0, b=0, lwd=2)

boxplot( 100*boxbiasBe0, names = c("mlp ", "mq", "gmq"),
        las=2, ylab = "", xlab = "", pch=19,  ylim=ylim,
        cex.axis = 1.45, cex.lab = 1.45, main='100 x BIAS [e,0]')
abline(a=0, b=0, lwd=2)

boxplot( 100*boxbiasB0u, names = c("mlp ", "mq", "gmq"),
        las=2, ylab = "", xlab = "", pch=19,  ylim=ylim,
        cex.axis = 1.45, cex.lab = 1.45, main='100 x BIAS [0,u]')
abline(a=0, b=0, lwd=2)

boxplot( 100*boxbiasBeu, names = c("mlp ", "mq", "gmq"),
        las=2, ylab = "", xlab = "", pch=19,  ylim=ylim,
        cex.axis = 1.45, cex.lab = 1.45, main='100 x BIAS [e,u]')
abline(a=0, b=0, lwd=2)

round(100*apply(boxbiasB00, 2, median), 3)
round(100*apply(boxbiasBe0, 2, median), 3)
round(100*apply(boxbiasB0u[1:36,], 2, median), 3)
round(100*apply(boxbiasB0u[37:40,], 2, median), 3)
round(100*apply(boxbiasBeu[1:36,], 2, median), 3)
round(100*apply(boxbiasBeu[37:40,], 2, median), 3)


##############################
###   Some boxplots RRMSE  ### # 1200 * 300
##############################

boxrmseB00 <- read.csv('results/rmseGABE00.csv')[, -1]
boxrmseBe0 <- read.csv('results/rmseGABEe0.csv')[, -1]
boxrmseB0u <- read.csv('results/rmseGABE0u.csv')[, -1]
boxrmseBeu <- read.csv('results/rmseGABEeu.csv')[, -1]

par(mfrow = c(1, 4),
    mar = c(4, 3, 2, 1),
    oma = c(0, 0, 0, 0),
    mgp = c(1.6, 0.5, 0))

ylim <- 100 * c(0.11, 0.3)

# rmse
boxplot( 100*boxrmseB00, names = c("mlp ", "mq", "gmq"),
        las=2, ylab = "", xlab = "", pch=19, ylim=ylim,
        cex.axis = 1.45, cex.lab = 1.45, main='100 x RMSE [0,0]')
abline(a=0, b=0, lwd=2)

boxplot( 100*boxrmseBe0, names = c("mlp ", "mq", "gmq"),
        las=2, ylab = "", xlab = "", pch=19, ylim=ylim,
        cex.axis = 1.45, cex.lab = 1.45, main='100 x RMSE [e,0]')
abline(a=0, b=0, lwd=2)

boxplot( 100*boxrmseB0u, names = c("mlp ", "mq", "gmq"),
        las=2, ylab = "", xlab = "", pch=19, ylim=ylim,
        cex.axis = 1.45, cex.lab = 1.45, main='100 x RMSE [0,u]')
abline(a=0, b=0, lwd=2)

boxplot( 100*boxrmseBeu, names =c("mlp ", "mq", "gmq"),
        las=2, ylab = "", xlab = "", pch=19, ylim=ylim,
        cex.axis = 1.45, cex.lab = 1.45, main='100 x RMSE [e,u]')
abline(a=0, b=0, lwd=2)

round(100*apply(boxrmseB00, 2, median), 3)
round(100*apply(boxrmseBe0, 2, median), 3)
round(100*apply(boxrmseB0u[1:36,], 2, median), 3)
round(100*apply(boxrmseB0u[37:40,], 2, median), 3)
round(100*apply(boxrmseBeu[1:36,], 2, median), 3)
round(100*apply(boxrmseBeu[37:40,], 2, median), 3)

