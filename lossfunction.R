###########################################################
#####   SOME REPRESENTATIONS FOR THE LOSS FUNCTIONS   #####
###########################################################

## -----------------------------
## Parameters
## -----------------------------
r   <- seq(-3, 3, length.out = 1000)

## -----------------------------
## 1. Quantile (check) loss
## -----------------------------
rho_quantile <- function(r, q) {
  r * (q - as.numeric(r < 0))
}


## -----------------------------
## 2. Asymmetric Huber loss 
## -----------------------------
rho_huber_asym_P4B <- function(r, q, c_q) {
  
  w <- abs(q - as.numeric(r <= 0))
  out <- numeric(length(r))
  
  idx1 <- abs(r) <= c_q
  idx2 <- abs(r) >  c_q
  
  out[idx1] <- 2 * (r[idx1]^2 / 2) * w[idx1]
  out[idx2] <- 2 * (c_q * abs(r[idx2]) - c_q^2 / 2) * w[idx2]
  
  out
}

## -----------------------------
## Visual comparison between check loss and asymmetric Huber loss
## -----------------------------
do.plots = FALSE
if(do.plots == TRUE )
{par(mfrow = c(2,3), mar = c(4, 4, 3, 1))

L1 <- rho_quantile(r, q = 0.25)
plot(r, L1, type = "l", lwd = 2,
     main = expression(paste("Check loss  ", rho[0.25](r))),
     xlab = "r", ylab = '', cex.axis=1.35, cex.lab=1.35, cex.main=1.35)
abline(v = 0, lty = 2, col = "gray", lwd = 2)
abline(h = 0, lty = 2, col = "gray", lwd = 2)

L1 <- rho_quantile(r, q = 0.50)
plot(r, L1, type = "l", lwd = 2,
     main = expression(paste("Check loss  ", rho[0.50](r))),
     xlab = "r", ylab = '', cex.axis=1.35, cex.lab=1.35, cex.main=1.35)
abline(v = 0, lty = 2, col = "gray", lwd = 2)
abline(h = 0, lty = 2, col = "gray", lwd = 2)

L1 <- rho_quantile(r, q = 0.90)
plot(r, L1, type = "l", lwd = 2,
     main = expression(paste("Check loss  ", rho[0.90](r))),
     xlab = "r", ylab = '', cex.axis=1.35, cex.lab=1.35, cex.main=1.35)
abline(v = 0, lty = 2, col = "gray", lwd = 2)
abline(h = 0, lty = 2, col = "gray", lwd = 2)

L2 <- rho_huber_asym_P4B(r, q=0.25, c_q=1.345)
plot(r, L2, type = "l", lwd = 2,
     main = expression(paste("Asymmetric Huber loss  ", rho[0.25](r))),
     xlab = "r", ylab = '', cex.axis=1.35, cex.lab=1.35, cex.main=1.35)
abline(v = 0, lty = 2, col = "gray", lwd = 2)
abline(h = 0, lty = 2, col = "gray", lwd = 2)

L2 <- rho_huber_asym_P4B(r, q=0.50, c_q=1.345)
plot(r, L2, type = "l", lwd = 2,
     main = expression(paste("Asymmetric Huber loss  ", rho[0.50](r))),
     xlab = "r", ylab = '', cex.axis=1.35, cex.lab=1.35, cex.main=1.35)
abline(v = 0, lty = 2, col = "gray", lwd = 2)
abline(h = 0, lty = 2, col = "gray", lwd = 2)

L2 <- rho_huber_asym_P4B(r, q=0.90, c_q=1.345)
plot(r, L2, type = "l", lwd = 2,
     main = expression(paste("Asymmetric Huber loss  ", rho[0.90](r))),
     xlab = "r", ylab = '', cex.axis=1.35, cex.lab=1.35, cex.main=1.35)
abline(v = 0, lty = 2, col = "gray", lwd = 2)
abline(h = 0, lty = 2, col = "gray", lwd = 2)
}


########################################
###     FINAL PLOTS (OTHER ONES)     ###
########################################

rho_q <- function(r, q) {
  ifelse(r >= 0, q * r, (q - 1) * r)  }

# Cumulative distribution function
F_u <- function(u, mu_q, sigma_q, q) {
  ifelse(u < mu_q,
         q * exp(-rho_q( r = ((u - mu_q) / sigma_q), q)),
         1 - (1 - q) * exp(-rho_q( r = ((u - mu_q) / sigma_q), q)))   }

# Quantile function
Q_p <- function(p, mu_q, sigma_q, q) {
  ifelse(p < q,
         mu_q + sigma_q/(1-q)*log(p/q),
         mu_q + sigma_q/q*log((1-q)/(1-p)))}

do.plots = FALSE
if(do.plots == TRUE )
{
# The plots 1
u <- seq(-50, 50, by = 0.1)
q <- c(0.1, 0.25, 0.5, 0.75, 0.9)
colors <- c("#1B3B6F", "#D55E00", "#009E73","#E69F00","#6A3D9A")

plot(u, F_u(u, mu_q = 0, sigma_q = 1, q = q[1]), 
     cex.axis = 1.55, cex.lab = 1.55, type = "l", main = "", lwd=3,
     col = colors[1], ylim = c(0, 1), ylab = "Distribution function", xlab = "u")

for(i in 2:length(q)) {
  lines(u, F_u(u, mu_q = 0, sigma_q = 1, q = q[i]), lwd=3, col = colors[i])}
legend("topleft", cex=1.55, legend = paste0("q = ", q), col = colors, lty = 1, lwd=3)


# The plots 2
p <- seq(0, 1, by = 0.001)
q <- c(0.1, 0.25, 0.5, 0.75, 0.9)
colors <- c("#1B3B6F", "#D55E00", "#009E73","#E69F00","#6A3D9A")

plot(p, Q_p(p, mu_q=0, sigma_q=1, q=q[1]), type = "l", 
     cex.axis = 1.55, cex.lab = 1.55, main = "", lwd=3,
     col = colors[1], ylab = "Quantile function", xlab = "p", ylim=c(-50, 50))

for(i in 2:length(q)) {
  lines(p, Q_p(p, mu_q = 0, sigma_q = 1, q = q[i]), lwd=3, col = colors[i])}
legend("topleft", cex=1.55, legend = paste0("q = ", q), col = colors, lty = 1, lwd=3)


###

par(mfrow = c(3, 4),
    mar = c(3, 3, 2, 1),
    oma = c(0, 0, 0, 0),
    mgp = c(1.6, 0.5, 0))
p <- seq(from = 0.05, to = 0.95, length.out = n)

r.squared <- list()

for(d in 1:D){ 
  
  q.emp <- quantile(res.s.new[[d]], p) 
  q.theo <- Q_p(p, mu_q = 0, sigma_q = 1, q = mqo.new[d])
  
  plot(q.emp, q.theo, cex.axis = 1.35, cex.lab = 1.35,
       main = bquote(hat(theta)[.(d)] == .(round(mqo[d], 3))))
  
  points(q.emp[regioncode.s==d], q.theo[regioncode.s==d], pch=19, col='blue', cex=1.1)
  fit <- lm(q.theo ~ q.emp)
  abline(fit, col = "red", lwd = 2)
  
  r.squared[[d]] <- summary(fit)$r.squared 
}
}

