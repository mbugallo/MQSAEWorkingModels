
###########################################################################
####     MEAN AND VARIANCE OF THE ASYMMETRIC HUBER MQ DISTRIBUTION     ####
###########################################################################

# NORMALIZING CONSTANT FOR BINARY DATA
Bq.binary <- function(q, k ,mu, s) {
  # Normalizing constant for the asymmetric Huber M-quantile distribution
  f <- exp(-rho_huber_asym_P4B( r = ( - mu ) / s, q = q, c_q = k )) + 
       exp(-rho_huber_asym_P4B( r = ( 1 - mu ) / s, q = q, c_q = k ))
  
  return(f)
}

# NORMALIZING CONSTANT FOR CONTINUOUS DATA
Bq.or <- function(q = q, k = k, s = 1) {
  # Normalizing constant for the asymmetric Huber M-quantile distribution
  s * (sqrt(pi / q)     * (pnorm(k * sqrt(2 * q)) - 0.5) +
         sqrt(pi / (1-q)) * (pnorm(k * sqrt(2 * (1 - q))) - 0.5) +
         (1 / (2 * k * q))      * exp(-q * k^2) +
         (1 / (2 * k * (1 - q))) * exp(-(1 - q) * k^2))
}

# MEAN OF A GALI RANDOM VARIABLE
meanAH <- function(q, m, s, k) {
  # Closed-form expression of the mean under the asymmetric Huber M-quantile model
  Bq <- function(q = q, k = k, s = 1) {
    s * (sqrt(pi / q)     * (pnorm(k * sqrt(2 * q)) - 0.5) +
           sqrt(pi / (1-q)) * (pnorm(k * sqrt(2 * (1 - q))) - 0.5) +
           (1 / (2 * k * q))      * exp(-q * k^2) +
           (1 / (2 * k * (1 - q))) * exp(-(1 - q) * k^2))
  }
  
  # Decomposition of the mean into left-tail, right-tail and central components
  part1 <- (1 / Bq(q = q, k = k, s = 1)) * (1 / (2 * q)) * (1 - exp(-q * k^2))
  part2 <- -(1 / Bq(q = q, k = k, s = 1)) * (1 / (2 * (1-q))) * (1 - exp(-(1 - q) * k^2))
  part3 <- -(1 / Bq(q = q, k = k, s = 1)) *
    exp(-k^2 * (1 - q)) * ((1 + 2 * k^2 * (1 - q)) / (4 * k^2 * (1 - q)^2))
  part4 <- (1 / Bq(q = q, k = k, s = 1)) *
    exp(-k^2 * q) * ((1 + 2 * k^2 * q) / (4 * k^2 * q^2))
  
  mean <- m + s * (part1 + part2 + part3 + part4)
  mean
}

# VARIANCE OF A GALI RANDOM VARIABLE
varAH <- function(q, m, s, k) {
  Bq <- function(q = q, k = k, s = 1) {
    s * (sqrt(pi / q)     * (pnorm(k * sqrt(2 * q)) - 0.5) +
           sqrt(pi / (1-q)) * (pnorm(k * sqrt(2 * (1 - q))) - 0.5) +
           (1 / (2 * k * q))      * exp(-q * k^2) +
           (1 / (2 * k * (1 - q))) * exp(-(1 - q) * k^2))
  }
  
  # First moment
  mean <- meanAH(q = q, m = m, s = s, k = k)
  
  # Components of second moment
  part1 <- (1 / Bq(q = q, k = k, s = 1)) *
    exp(-k^2 * (1 - q)) *
    ((1 + 2 * k^2 * (1 - q) + 2 * k^4 * (1 - q)^2) /
       (4 * k^3 * (1 - q)^3))
  
  part2 <- (1 / Bq(q = q, k = k, s = 1)) *
    exp(-k^2 * q) *
    ((1 + 2 * k^2 * q + 2 * k^4 * q^2) /
       (4 * k^3 * q^3))
  
  part3 <- (1 / (2 * Bq(q = q, k = k, s = 1))) *
    ((sqrt(pi) * (pnorm(k * sqrt(2 * q)) - 0.5)) / (q^(3/2)) -
       (k * exp(-k^2 * q)) / q)
  
  part4 <- (1 / (2 * Bq(q = q, k = k, s = 1))) *
    ((sqrt(pi) * (pnorm(k * sqrt(2 * (1 - q))) - 0.5)) /
       ((1 - q)^(3/2)) -
       (k * exp(-k^2 * (1 - q))) / (1 - q))
  
  # Second moment obtained by combining truncated and tail contributions
  moment2  <- (part1 + part2 + part3 + part4)
  
  # Analytical computation of the variance using first and second moments
  variance <- s^2 * (moment2 - mean^2)
  
  list(
    variance = variance,
    moment2  = moment2,
    moment1  = mean,
    part1    = part1,
    part2    = part2,
    part3    = part3,
    part4    = part4
  )
}


##########################################
###   MQ REGRESSION FITTING ALGORITHM  ###
##########################################

# Iteratively reweighted least squares algorithm for M-quantile regression
QRLM <-function (x, y, case.weights = rep(1, nrow(x)), k=1.345, 
                 var.weights = rep(1, nrow(x)), ..., w = rep(1, nrow(x)), 
                 init = "ls", psi = psi.huber, 
                 scale.est = c("MAD", "Huber", "proposal 2"), 
                 k2 = 1.345, method = c("M", "MM"), maxit = 20, 
                 acc = 1e-04, test.vec = "resid", q = 0.5)
{
  irls.delta <- function(old, new) sqrt(sum((old - new)^2)/max(1e-20, sum(old^2)))
  irls.rrxwr <- function(x, w, r) {
    w <- sqrt(w)
    max(abs((matrix(r*w,1,length(r)) %*% x)/
              sqrt(matrix(w,1,length(r)) %*% (x^2))))/sqrt(sum(w*r^2))
  }
  method <- match.arg(method)
  nmx <- deparse(substitute(x))
  if (is.null(dim(x))) {
    x <- as.matrix(x)
    colnames(x) <- nmx
  }
  else x <- as.matrix(x)
  if (is.null(colnames(x)))
    colnames(x) <- paste("X", seq(ncol(x)), sep = "")
  if (qr(x)$rank < ncol(x))
    stop("x is singular: singular fits are not implemented in rlm")
  if (!(any(test.vec == c("resid", "coef", "w", "NULL")) || is.null(test.vec)))
    stop("invalid testvec")
  if (length(var.weights) != nrow(x))
    stop("Length of var.weights must equal number of observations")
  if (any(var.weights < 0))
    stop("Negative var.weights value")
  if (length(case.weights) != nrow(x))
    stop("Length of case.weights must equal number of observations")
  w <- (w * case.weights)/var.weights
  if (method == "M") {
    scale.est <- match.arg(scale.est)
    if (!is.function(psi))
      psi <- get(psi, mode = "function")
    arguments <- list(...)
    if (length(arguments)) {
      pm <- pmatch(names(arguments), names(formals(psi)), nomatch = 0)
      if (any(pm == 0))
        warning(paste("some of ... do not match"))
      pm <- names(arguments)[pm > 0]
      formals(psi)[pm] <- unlist(arguments[pm])
    }
    if (is.character(init)) {
      if (init == "ls")
        temp <- lm.wfit(x, y, w, method = "qr")
      else if (init == "lts")
        temp <- lqs.default(x, y, intercept = FALSE, nsamp = 200)
      else stop("init method is unknown")
      coef <- temp$coef
      resid <- temp$resid
    }
    else {
      if (is.list(init))
        coef <- init$coef
      else coef <- init
      resid <- y - x %*% coef
    }
  }
  else if (method == "MM") {
    scale.est <- "MM"
    temp <- lqs.default(x, y, intercept = FALSE, method = "S", k0 = 1.548)
    coef <- temp$coef
    resid <- temp$resid
    psi <- psi.bisquare
    if (length(arguments <- list(...)))
      if (match("c", names(arguments), nomatch = FALSE)) {
        c0 <- arguments$c
        if (c0 > 1.548) {
          psi$c <- c0
        }
        else warning("c must be at least 1.548 and has been ignored")
      }
    scale <- temp$scale
  }
  else stop("method is unknown")
  done <- FALSE
  conv <- NULL
  n1 <- nrow(x) - ncol(x)
  if (scale.est != "MM")
    scale <- mad(resid/sqrt(var.weights), 0)
  
  qest <- matrix(0, nrow = ncol(x), ncol = length(q))
  qwt <- matrix(0, nrow = nrow(x), ncol = length(q))
  qfit <- matrix(0, nrow = nrow(x), ncol = length(q))
  qres <- matrix(0, nrow = nrow(x), ncol = length(q))
  qvar.matrix <- array(rep(0,ncol(x)*ncol(x)),dim=c(ncol(x),ncol(x),length(q)))
  qscale <- NULL
  
  # Separate IRLS fitting is performed for each quantile order
  for(i in 1:length(q)) {
    for (iiter in 1:maxit) {
      if (!is.null(test.vec))
        testpv <- get(test.vec)
      if (scale.est != "MM") {
        if (scale.est == "MAD")
          scale <- median(abs(resid/sqrt(var.weights)))/0.6745
        else {gamma<- 4*k2^2*(1-pnorm(k2))*((1-q[i])^2+q[i]^2) - 4*k2*dnorm(k2)*((1-q[i])^2+q[i]^2) + 
          4*(1-q[i])^2*(pnorm(0)-(1-pnorm(k2))) + 4*q[i]^2*(pnorm(k2)-pnorm(0))
        scale <- sqrt(sum(pmin(resid^2/var.weights,(k2*scale)^2))/(n1*gamma))
        }
        if (scale == 0) {
          done <- TRUE
          break
        }
      }
      w <- psi(resid/(scale * sqrt(var.weights)),k=k) * case.weights
      ww <- 2 * (1 - q[i]) * w
      ww[resid > 0] <- 2 * q[i] * w[resid > 0]
      w <- ww
      temp <- lm.wfit(x, y, w, method = "qr")
      coef <- temp$coef
      resid <- temp$residuals
      if (!is.null(test.vec))
        convi <- irls.delta(testpv, get(test.vec))
      else convi <- irls.rrxwr(x, wmod, resid)
      conv <- c(conv, convi)
      done <- (convi <= acc)
      if (done)
        break
    }
    if (!done)
      warning(paste("rlm failed to converge in", maxit, "steps at q = ", q[i]))
    qest[, i] <- coef
    qscale[i]<-scale
    qwt[, i] <- w
    qfit[, i] <- temp$fitted.values
    qres[,i] <- resid
    
    tmp.res.mq<-qres[,i]/qscale[i]
    Epsi2<-(sum((qwt[,i]*tmp.res.mq)^2)/(nrow(x)-ncol(x)))
    Epsi<-(1/qscale[i])*(sum(2*(q[i]*(0<=tmp.res.mq & tmp.res.mq<= k)+
                                  (1-q[i])*(-k <=tmp.res.mq & tmp.res.mq<0)))/nrow(x))
    qvar.matrix[,,i]<- (((Epsi2)/Epsi^2)*solve(t(x)%*%x))
    
  }
  list(fitted.values = qfit, residuals = qres, q.values = q, q.weights = qwt, coef= qest,
       qscale=qscale,var.beta=qvar.matrix)
}


###########################################
###   COMPUTING OF THE QUANTILE-ORDERS  ###
###########################################

# Linear interpolation to locate the zero-crossing of the estimating function
"zerovalinter"<-function(y, x)
{
  if(min(y) > 0) {
    xmin <- x[y == min(y)]
    if(length(xmin) > 0)
      xmin <- xmin[length(xmin)]
    xzero <- xmin
  }
  else {
    if(max(y) < 0) {
      xmin <- x[y == max(y)]
      if(length(xmin) > 0)
        xmin <- xmin[1]
      xzero <- xmin
    }
    else {
      y1 <- min(y[y > 0])
      if(length(y1) > 0)
        y1 <- y1[length(y1)]
      y2 <- max(y[y < 0])
      if(length(y2) > 0)
        y2 <- y2[1]
      x1 <- x[y == y1]
      if(length(x1) > 0)
        x1 <- x1[length(x1)]
      x2 <- x[y == y2]
      if(length(x2) > 0)
        x2 <- x2[1]
      xzero <- (x2 * y1 - x1 * y2)/(y1 - y2)
      xmin <- x1
      if(abs(y2) < y1)
        xmin <- x2
    }
  }
  resu <-  xzero
  resu
}


###########################################
###    LINEAR INTERPOLATION FUNCTION    ###
###########################################

"gridfitinter"<-function(y,expectile,Q)
{
  nq<-length(Q)
  diff <- y %*% t(as.matrix(rep(1, nq))) - expectile        
  vectordest <- apply(diff, 1, zerovalinter,Q)    
}


#######################################
###    MEDIAN ABSOLUTE DEVIATION    ###
#######################################

fun.MAD <- function(x){ median(abs( x-median(x)) )/0.6745 }


##################################
###       HUBER FUNCTIONS      ###
##################################

hub.psi <- function(x, k){ ifelse(abs(x) <= k, x, sign(x) * k) }
der.hub.psi <- function(x, k){ ifelse(abs(x) <= k, 1, 0) }


##########################################################
###     NEW WAY TO COMPUTE THE BEST PROBABILITIES      ###
##########################################################

# Selection of the M-quantile order in a novel way
mqo.loss <- function(x.s, y.s, mqo, Nd, nd, regioncode.s, binary, ell){
  
  tau.r <- seq(-0.495, 0.495, 0.02)
  eta <- 0.9 + 0.1 * (pmin(1, nd / Nd))^{ell}
  tau <- 0.5 + eta * tau.r
  
  D <- length(unique(regioncode.s))
  lambda.d <- (pmin(1, nd / Nd))^{ell}
  q.best.q <- vector("list", D)
  
  # DEVELOPMENTS FOR BINARY DATA  
  if(binary==TRUE){
    
    mod <- glm.mq.binom(y = y.s, x =  cbind(1, x.s), q = tau)
    
    res.s.q <- lapply(seq_along(tau), function(q){
      y.s - cbind(1, x.s) %*% mod$coefficients[, q]})
    
    sd.q <- sapply(res.s.q, fun.MAD)
    l.seq.tau <- length(tau)
    
    eta.q <- lapply(seq_along(tau), function(q){cbind(1, x.s) %*% mod$coefficients[, q]})
    
    for(d in 1:D){
      C1 <- median(abs(tau - mqo[d]), na.rm=T)
      C2 <- median(abs(sapply(1:l.seq.tau, function(q){
        sum(-rho_huber_asym_P4B(r = res.s.q[[q]][regioncode.s == d] / sd.q[[q]],
             q = tau[q], c_q = 1.345 ) - log(Bq.binary(q = tau[q], k = 1.345, 
              mu = eta.q[[q]][regioncode.s == d], s =  sd.q[[q]])))})), na.rm=T)
      
      arg.function <- sapply(1:l.seq.tau, function(q){
        mean((1-lambda.d[d]) * ( -rho_huber_asym_P4B( r = res.s.q[[q]][regioncode.s == d] / sd.q[[q]],
                      q = tau[q], c_q = 1.345 ) - log(Bq.binary(q = tau[q], 
                      k = 1.345, mu = eta.q[[q]][regioncode.s == d], s = sd.q[[q]]))
        )) / C2 - (lambda.d[d]) * abs(tau[q] - mqo[d]) / C1  })
      
      q.best.q[[d]] <- tau[which.min(abs(arg.function))]
    }
  
  # DEVELOPMENTS FOR CONTINUOUS DATA   
  } else {
    
    mod <- QRLM(x = cbind(1, x.s), y = y.s, q = tau, maxit = 30, k = 1.345)
    
    res.s.q <- lapply(seq_along(tau), function(q){
      y.s - cbind(1, x.s) %*% mod$coef[, q]})
    
    sd.q <- sapply(res.s.q, fun.MAD)
    l.seq.tau <- length(tau)
    
    for(d in 1:D){
      C1 <- median(abs(tau - mqo[d]), na.rm=T)
      C2 <- median(abs(sapply(1:l.seq.tau, function(q){
        sum(-rho_huber_asym_P4B(r = res.s.q[[q]][regioncode.s == d] / sd.q[[q]],
            q = tau[q], c_q = 1.345 ) - log(Bq.or(q = tau[q], 
                                        k = 1.345, s = sd.q[[q]])))})), na.rm=T)
      
      arg.function <- sapply(1:l.seq.tau, function(q){
        mean((1-lambda.d[d]) * ( -rho_huber_asym_P4B( r = res.s.q[[q]][regioncode.s == d] / sd.q[[q]],
                              q = tau[q], c_q = 1.345 ) - log(Bq.or(q = tau[q], 
                                                  k = 1.345, s = sd.q[[q]]))
        )) / C2 - (lambda.d[d]) * abs(tau[q] - mqo[d]) / C1  })
      
      q.best.q[[d]] <- tau[which.min(abs(arg.function))]
    }
  }
  
  return(unlist(q.best.q))
}


##########################################
###     ASYMPTOTIC MSE ESTIMATION      ###
##########################################

asym.MSE <- function(x.s, y.s, mqo, x.r, Nd, nd, regioncode.s, 
                     mod.SAE, res.s, sd, ck){
  
  n <- dim(cbind(1,x.s))[1]
  p <- dim(cbind(1,x.s))[2]
  D <- length(unique(regioncode.s))
  mse.UMQ <- var.beta.prod <- list()
  
  # other.term <- sapply(1:D, function(d){sd[d]^2/(Nd[d]-nd[d])*(2*mqo[d]^2-2*mqo[d]+1)/(mqo[d]^2*(1-mqo[d])^2)})
  # other.term  <- sapply(1:D, function(d) sum((res.s[[d]]^2)[regioncode.s==d]^2)/((Nd[d]-nd[d])*(nd[d]-1)))
  other.term <- sapply(1:D, function(d) sum(res.s[[d]]^2)/((Nd[d]-nd[d])*(n-D)))
  other.term2  <- sapply(1:D, function(d){ hub.psi( meanAH(q=mqo[d], m=0, s=1, k=1.345) , k=ck[d])^2 +
      (hub.psi( meanAH(q=mqo[d], m=0, s=1, k=1.345),k=ck[d])-sum(res.s[[d]][regioncode.s==d]/sd[d]) )^2})
  
  for(d in 1:D){
    xd.r.mean <-colMeans(cbind(1,x.r)[regioncode.r==d,])
    
    # Variance of beta
    var.beta <- n^2*sd[d]^2/(n-p)*sum(hub.psi(res.s[[d]]/sd[d], k=1.345)^2)/
        (sum(der.hub.psi(res.s[[d]]/sd[d], k=1.345))^2)*solve(t(cbind(1,x.s))%*%cbind(1,x.s))
      
    var.beta.prod[[d]] <- xd.r.mean%*%var.beta%*%xd.r.mean
    
    mse.UMQ[[d]] <- (1- nd[d]/Nd[d])^2*(var.beta.prod[[d]] + other.term[d]+ sd[d]^2*other.term2[d]) 
  }
  
  return(unlist(mse.UMQ))
}


###########################################
###     ASYMPTOTIC BIAS CORRECTION      ###
###########################################

f.sum.BMQ <- function(res.s, sd, regioncode.s, grid){
  
  res <- lapply(1:D, function(d) {sapply(grid, function(ck) {
    psi_val <- sum(hub.psi(res.s[[d]][regioncode.s==d]/sd[d], k = ck))
    term1 <- psi_val - sum(res.s[[d]][regioncode.s==d]/sd[d])
      
    term2 <- sum(hub.psi(res.s[[d]][regioncode.s==d]/sd[d], k = ck)^2)
    term1^2 + term2   })  })
  ck <- grid[sapply(res, which.min)]
  
  sum.BMQ <- sapply(1:D, function(d) 
    (Nd[d]-nd[d])/(Nd[d]*nd[d])*sd[d]*sum(hub.psi(res.s[[d]][regioncode.s==d]/sd[d], k = ck[d])) )
  
  return(list(sum.BMQ=sum.BMQ, ck=ck))
}


#################################################
###     BEST ASYMPTOTIC BIAS CORRECTION       ###
#################################################

f.sum.UMQ <- function(res.s, sd, regioncode.s, mqo, grid, Nd, nd){
  
  res <- lapply(1:D, function(d) {sapply(grid, function(ck) {
    psi_val <- hub.psi((1 - 2*mqo[d]) / (mqo[d] * (1-mqo[d])), k = ck)
    term1 <-  sd[d] * psi_val - sum(res.s[[d]][regioncode.s == d]) 
    term2 <-  sd[d]^2 * psi_val^2
    term1^2 + term2   })  })
  ck <- grid[sapply(res, which.min)]
  
  sum.UMQ <- sapply(1:D, function(d) 
    (Nd[d]-nd[d])/Nd[d]*sd[d]*hub.psi(meanAH(q=mqo[d], m=0, s=1, k=1.345), k=ck[d])) 
    
  return(list(sum.UMQ, ck))
  
}


############################################################
###     NEW PREDICTORS DERIVED FROM M-QUANTILE MODELS    ###
############################################################

# Function to compute the robust RUMQ predictor
calc_RUMQ <- function(pred.MQ, res.s, sd, regioncode.s, mqo, Nd, nd, D) {
  
  # Compute robust correction terms at two grid resolutions
  RUMQ.aux  <- f.sum.UMQ(res.s, sd, regioncode.s, mqo, grid = 0.5, Nd, nd)
  RUMQ.aux2 <- f.sum.UMQ(res.s, sd, regioncode.s, mqo, grid = 0.25, Nd, nd)
  
  # Initialize robust MQ predictor
  pred.RUMQ <- pred.MQ
  
  # Apply adaptive correction depending on the distance from the median (0.5)
  for(d in 1:D){
    if(abs(mqo[d] - 0.5) > 0.25){
      # Use coarser grid correction for extreme quantile orders
      pred.RUMQ[d] <- pred.RUMQ[d] + RUMQ.aux[[1]][d]
    } else {
      # Use finer grid correction for central quantile orders
      pred.RUMQ[d] <- pred.RUMQ[d] + RUMQ.aux2[[1]][d]
    }
  }
  
  return(pred.RUMQ)
}

# Function to compute the robust RBUMQ predictor
calc_RBUMQ <- function(sum.s, x.r, mod.SAE.new, res.s.new, sd.new, 
                       mqo.new, regioncode.s, mqo, Nd, nd, D) {
  
  # Compute robust correction terms at two grid resolutions
  RBUMQ.aux  <- f.sum.UMQ(res.s.new, sd.new, regioncode.s, mqo.new, grid=0.5, Nd, nd)
  RBUMQ.aux2 <- f.sum.UMQ(res.s.new, sd.new, regioncode.s, mqo.new, grid=0.25, Nd, nd)
  
  # Initialize robust MQ predictor
  pred.RBUMQ <- sapply(1:D, function(d){
    1/(Nd[d] + nd[d]) * (sum.s[d] + Nd[d]*x.r[regioncode.r==d, ]%*% mod.SAE.new$coef[,d])})
  
  # Apply adaptive correction depending on the distance from the median (0.5)
  for(d in 1:D){
    if(abs(mqo.new[d] - 0.5) > 0.25){
      # Use coarser grid correction for extreme quantile orders
      pred.RBUMQ[d] <- pred.RBUMQ[d] + RBUMQ.aux[[1]][d]
    } else {
      # Use finer grid correction for central quantile orders
      pred.RBUMQ[d] <- pred.RBUMQ[d] + RBUMQ.aux2[[1]][d]
    }
  }
  
  return(pred.RBUMQ)
}

