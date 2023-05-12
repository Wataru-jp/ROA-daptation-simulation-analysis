## Numerical exercise of ROA model: climate change adaptation investment
## Author: Wataru Kodama
## Date: 15/3/2023
## Update: 5/5/2023

################################################################################
### 1. Define fundamental polynomial function and delta function
################################################################################
para <- c(0.03, 0.2, 0.25, 0.05, 0.2) # baseline (mu, sigma, eta, lambda, rho)
eps <- c(0, 0.5, 2, 5) # list of epsilon
# Delta function
delta_fn <- function(m, s, e, l, r, eps) { 
  r + m*(1-eps) + 1/2*eps*(1-eps)*s^2 + l- l*(1-e)^(1-eps) }
# Polynomial function
polynomial_fn <- function(x, m, s, e, l, r) {
  1/2*(x-1)*x*s^2 - m*x - (r+l) + l*(1-e)^x
  }
polynomial_fn(0, para[1], para[2], para[3], para[4], para[5]) # baseline when f(0)

# baseline value
delta_fn(para[1], para[2], para[3], 0.05, para[5], eps[1]) # baseline values eps=0
polynomial_fn(1-eps[1], para[1], para[2], para[3], para[4], para[5]) # f(1-eps)
delta_fn(para[1], para[2], para[3], 0.05, para[5], eps[2]) # baseline values eps=1/2
polynomial_fn(1-eps[2], para[1], para[2], para[3], para[4], para[5]) # f(1-eps)
delta_fn(para[1], para[2], para[3], 0.05, para[5], eps[3]) # baseline values eps=2
polynomial_fn(1-eps[3], para[1], para[2], para[3], para[4], para[5]) # f(1-eps)

# polynomial vs beta
beta <- function(x){polynomial_fn(x, para[1], para[2], para[3], para[4], para[5])}
plot(beta, -2, 0) # see plot


##############################################
### 2. Investment threshold under ROA and ENPV
##############################################
para <- c(0.03, 0.2, 0.25, 0.05, 0.2) # baseline (mu, sigma, eta, lambda, rho)
eps <- c(0, 0.5, 2) # list of epsilon
I <- 2

# Uniroot function
require(utils) # for str
root <- uniroot(polynomial_fn, c(-10, 0), tol = 0.0001, m=para[1], s=para[2], e=para[3], l=para[4], r=para[5])
root$root

# Utility function 
utility <- function(x, eps){(x^(1-eps))/(1-eps)}

# Threshold: ROA
delta <- delta_fn(para[1], para[2], para[3], para[4], para[5], eps)
beta <- uniroot(polynomial_fn, c(-10, 0), tol = 0.0001, m=para[1], 
                s=para[2], e=para[3], l=para[4], r=para[5])$root
thr <- ((delta/para[5])* beta/(beta-1+eps))^(1/(1-eps))* (1-I*0.1)
thr_EPV <- (delta/para[5])^(1/(1-eps))* (1-I*0.1)

rm(delta0, delta1, delta2, beta, hurdle0, hurdle1, hurdle2)
cat("The thresholds under ROA are: ", thr, "\n")
cat("The thresholds under ENPV are: ", thr_EPV, "\n")


########################################
### 3. Threshold vs different parameters
########################################
para <- c(0.03, 0.2, 0.25, 0.05, 0.2) # baseline (mu, sigma, eta, lambda, rho)
eps <- c(0, 0.5, 2) # list of epsilon
I <- 2
### Relationship with drift rate ### 
drift <- seq(0, 0.1, 0.01) # drift rate range

for (n in 0:2) {
  assign(paste0("delta_list", n), list())
  assign(paste0("hurdle_list", n), list())
  assign(paste0("thr_list", n, "a"), list())
} 

for (i in 1:length(drift)){
  delta_list0[[i]] <- delta_fn(drift[i], para[2], para[3], para[4], para[5], eps[1])
  delta_list1[[i]] <- delta_fn(drift[i], para[2], para[3], para[4], para[5], eps[2])
  delta_list2[[i]] <- delta_fn(drift[i], para[2], para[3], para[4], para[5], eps[3])
  beta_drift[[i]] <- uniroot(polynomial_fn, c(-10, 0), tol = 0.0001, m=drift[i], 
                           s=para[2], e=para[3], l=para[4], r=para[5])$root
  hurdle_list0[[i]] <- beta_drift[[i]]/(beta_drift[[i]]-1)
  hurdle_list1[[i]] <- (beta_drift[[i]]/(beta_drift[[i]]-0.5))^2
  hurdle_list2[[i]] <- (beta_drift[[i]]+1)/beta_drift[[i]]
  thr_list0a[[i]] <- hurdle_list0[[i]]* (delta_list0[[i]]/para[5])* (1 - I* 0.1)
  thr_list1a[[i]] <- hurdle_list1[[i]]* (delta_list1[[i]]/para[5])^2* (1 - I* 0.1)
  thr_list2a[[i]] <- hurdle_list2[[i]]* (para[5]/delta_list2[[i]])* (1 - I* 0.1)
  }

plot(drift, thr_list0a, type = "l")
plot(drift, thr_list1a, type = "l")
plot(drift, thr_list2a, type = "l")


### Relationship with volatility ### 
volatil <- seq(0, 0.5, 0.01) # volatility range

for (n in 0:2) {
  assign(paste0("delta_list", n), list())
  assign(paste0("hurdle_list", n), list())
  assign(paste0("thr_list", n, "b"), list())
} 
beta_volatil <- list()

for (i in 1:length(volatil)){
  delta_list0[[i]] <- delta_fn(para[1], volatil[i], para[3], para[4], para[5], eps[1])
  delta_list1[[i]] <- delta_fn(para[1], volatil[i], para[3], para[4], para[5], eps[2])
  delta_list2[[i]] <- delta_fn(para[1], volatil[i], para[3], para[4], para[5], eps[3])
  beta_volatil[[i]] <- uniroot(polynomial_fn, c(-10, 0), tol = 0.0001, m=para[1], 
                            s=volatil[i], e=para[3], l=para[4], r=para[5])$root
  hurdle_list0[[i]] <- beta_volatil[[i]]/(beta_volatil[[i]]-1)
  hurdle_list1[[i]] <- (beta_volatil[[i]]/(beta_volatil[[i]]-0.5))^2
  hurdle_list2[[i]] <- (beta_volatil[[i]]+1)/beta_volatil[[i]]
  thr_list0b[[i]] <- hurdle_list0[[i]]* (delta_list0[[i]]/para[5])* (1 - I* 0.1)
  thr_list1b[[i]] <- hurdle_list1[[i]]* (delta_list1[[i]]/para[5])^2* (1 - I* 0.1)
  thr_list2b[[i]] <- hurdle_list2[[i]]* (para[5]/delta_list2[[i]])* (1 - I* 0.1)
}

plot(volatil, thr_list0b, type='l')
plot(volatil, thr_list1b, type='l')
plot(volatil, thr_list2b, type='l')


### Relationship with magnitude of weather extremes ### 
etalist <- seq(0.0, 0.4, 0.01) # magnitude of shock range

for (n in 0:2) {
  assign(paste0("delta_list", n), list())
  assign(paste0("hurdle_list", n), list())
  assign(paste0("thr_list", n, "c"), list())
} 
beta_eta <- list()

for (i in 1:length(etalist)){
  delta_list0[[i]] <- delta_fn(para[1], para[2], etalist[i], para[4], para[5], eps[1])
  delta_list1[[i]] <- delta_fn(para[1], para[2], etalist[i], para[4], para[5], eps[2])
  delta_list2[[i]] <- delta_fn(para[1], para[2], etalist[i], para[4], para[5], eps[3])
  beta_eta[[i]] <- uniroot(polynomial_fn, c(-10, 0), tol = 0.0001, m=para[1], 
                            s=para[2], e=etalist[i], l=para[4], r=para[5])$root
  hurdle_list0[[i]] <- beta_eta[[i]]/(beta_eta[[i]]-1)
  hurdle_list1[[i]] <- (beta_eta[[i]]/(beta_eta[[i]]-0.5))^2
  hurdle_list2[[i]] <- (beta_eta[[i]]+1)/beta_eta[[i]]
  thr_list0c[[i]] <- hurdle_list0[[i]]* (delta_list0[[i]]/para[5])* (1 - I* 0.1)
  thr_list1c[[i]] <- hurdle_list1[[i]]* (delta_list1[[i]]/para[5])^2* (1 - I* 0.1)
  thr_list2c[[i]] <- hurdle_list2[[i]]* (para[5]/delta_list2[[i]])* (1 - I* 0.1)
}

plot(etalist, thr_list0c, type='l')
plot(etalist, thr_list1c, type='l')
plot(etalist, thr_list2c, type='l')


### Relationship with frequency of weather extreme ### 
lamlist <- seq(0, 0.15, 0.01) # magnitude of shock range

for (n in 0:2) {
  assign(paste0("delta_list", n), list())
  assign(paste0("hurdle_list", n), list())
  assign(paste0("thr_list", n, "d"), list())
} 
beta_lam <- list()

for (i in 1:length(lamlist)){
  delta_list0[[i]] <- delta_fn(para[1], para[2], para[3], lamlist[i], para[5], eps[1])
  delta_list1[[i]] <- delta_fn(para[1], para[2], para[3], lamlist[i], para[5], eps[2])
  delta_list2[[i]] <- delta_fn(para[1], para[2], para[3], lamlist[i], para[5], eps[3])
  beta_lam[[i]] <- uniroot(polynomial_fn, c(-10, 0), tol = 0.0001, m=para[1], 
                           s=para[2], e=para[3], l=lamlist[i], r=para[5])$root
  hurdle_list0[[i]] <- beta_lam[[i]]/(beta_lam[[i]]-1)
  hurdle_list1[[i]] <- (beta_lam[[i]]/(beta_lam[[i]]-0.5))^2
  hurdle_list2[[i]] <- (beta_lam[[i]]+1)/beta_lam[[i]]
  thr_list0d[[i]] <- hurdle_list0[[i]]* (delta_list0[[i]]/para[5])* (1 - I* 0.1)
  thr_list1d[[i]] <- hurdle_list1[[i]]* (delta_list1[[i]]/para[5])^2* (1 - I* 0.1)
  thr_list2d[[i]] <- hurdle_list2[[i]]* (para[5]/delta_list2[[i]])* (1 - I* 0.1)
}

plot(lamlist, thr_list0d, type='l')
plot(lamlist, thr_list1d, type='l')
plot(lamlist, thr_list2d, type='l')


### Graph ### 
setwd("C:/OwnCloud/A05/08_WP2_ROA/01_Irrigation_investment_numerical")
setwd("/Users/kodam1/ownCloud - wataru.kodama@uni-goettingen.de@owncloud.gwdg.de/A05/08_WP2_ROA/01_Irrigation_investment_numerical")
#par(mfrow = c(2, 2), cex.main = 0.8)
#par(mfrow=c(1,1), cex.main = 1)

# First plot
png("03_Programming/numerical/figure/threhold_mu.png", width = 5, height = 4, units = 'in', res = 500)
plot(drift, thr_list0a, type='l', lty=1, ylim=c(0.58,0.71), xlab='', ylab='')
lines(drift, thr_list1a, lty=2)
lines(drift, thr_list2a, lty=3)
legend("bottomright", cex = 0.8, lty=1:3,
       legend = c(expression(paste(epsilon," = 0")), 
                  expression(paste(epsilon," = 0.5 ")), 
                  expression(paste(epsilon," = 2"))))
dev.off()

# Second plot
png("03_Programming/numerical/figure/threhold_sigma.png", width = 5, height = 4, units = 'in', res = 500)
plot(volatil, thr_list0b, type='l', lty=1, ylim=c(0.43,0.8), xlab='', ylab='')
lines(volatil, thr_list1b, lty=2)
lines(volatil, thr_list2b, lty=3)
legend("topright", cex = 0.8, lty=1:3,
       legend = c(expression(paste(epsilon," = 0")), 
                  expression(paste(epsilon," = 0.5 ")), 
                  expression(paste(epsilon," = 2"))))
dev.off()

# Third plot
png("03_Programming/numerical/figure/threhold_lambda.png", width = 5, height = 4, units = 'in', res = 500)
plot(lamlist, thr_list0d, type='l', lty=1, ylim=c(0.625,0.665), xlab='', ylab='')
lines(lamlist, thr_list1d, lty=2)
lines(lamlist, thr_list2d, lty=3)
legend("bottomright", cex = 0.8, lty=1:3,
       legend = c(expression(paste(epsilon," = 0")), 
                  expression(paste(epsilon," = 0.5 ")), 
                  expression(paste(epsilon," = 2"))))
dev.off()

# Fourth plot
png("03_Programming/numerical/figure/threhold_eta.png", width = 5, height = 4, units = 'in', res = 500)
plot(etalist, thr_list0c, type='l', lty=1, ylim=c(0.52,0.67), xlab='', ylab='')
lines(etalist, thr_list1c, lty=2)
lines(etalist, thr_list2c, lty=3)
legend("topright", cex = 0.8, lty=1:3,
       legend = c(expression(paste(epsilon," = 0")), 
                expression(paste(epsilon," = 0.5 ")), 
                expression(paste(epsilon," = 2"))))
dev.off()



########################################
### 4. Relative importance of parameters
########################################
para <- c(0.03, 0.2, 0.25, 0.05, 0.2) # baseline (mu, sigma, eta, lambda, rho)
eps <- c(0, 0.5, 2) # list of epsilon
I <- 2
# Chnage in value
delta <- function(eps) { 
  delta_fn(para[1], para[2], para[3], para[4], para[5], eps) 
  }
value <- function(eps) { 
  utility(1-0.05, eps)/0.2 - utility(1, eps)/delta(eps) 
  }
dV <- function(eps) {
  value(eps)*0.01
}

# Change in parameters
delta_mu <- 0.01
delta_sigma <- function(eps) { 
  (delta_mu/eps + para[2]^2)^0.5 - para[2] }
delta_lambda <- function(eps) { 
  delta_mu*(1-eps)/ (1-(1-para[3])^(1-eps)) }
delta_eta <- function(eps) { 
  (delta_mu*(1-eps)/ para[4] - (1-para[3])^(1-eps))^ (1/(1-eps)) + (1-para[3]) }
cat("Change in parameters are: ", 
    delta_mu, delta_lambda(0), delta_eta(0), "\n")
#para2 <- c(para[1]+delta_mu, para[2], para[3]+delta_lambda(0), para[4]+delta_eta(0))
para2 <- c(para[1] + 0.01, para[2], para[3] + 0.20, para[4] + 0.04)

# Threshold: change in mu
delta_mu <- function(eps) { 
  delta_fn(para2[1], para[2], para[3], para[4], para[5], eps) 
}
beta0 <- uniroot(polynomial_fn, c(-10, 0), tol = 0.0001, m=para2[1], 
                s=para[2], e=para[3], l=para[4], r=para[5])$root
hurdle0 <- beta0/(beta0-1+eps[1])
hurdle1 <- (beta0/(beta0-1+eps[2]))^2
hurdle2 <- (beta0-1+eps[3])/beta0
thr0_mu <- hurdle0* (delta_mu(eps[1])/para[5])* (1 - I* 0.1)
thr1_mu <- hurdle1* (delta_mu(eps[2])/para[5])^2* (1 - I* 0.1)
thr2_mu <- hurdle2* (para[5]/delta_mu(eps[3]))* (1 - I* 0.1)

# Threshold: change in eta
delta_eta <- function(eps) { 
  delta_fn(para[1], para[2], para2[3], para[4], para[5], eps) 
}
beta0 <- uniroot(polynomial_fn, c(-10, 0), tol = 0.0001, m=para[1], 
                 s=para[2], e=para2[3], l=para[4], r=para[5])$root
hurdle0 <- beta0/(beta0-1+eps[1])
hurdle1 <- (beta0/(beta0-1+eps[2]))^2
hurdle2 <- (beta0-1+eps[3])/beta0
thr0_eta <- hurdle0* (delta_eta(eps[1])/para[5])* (1 - I* 0.1)
thr1_eta <- hurdle1* (delta_eta(eps[2])/para[5])^2* (1 - I* 0.1)
thr2_eta <- hurdle2* (para[5]/delta_eta(eps[3]))* (1 - I* 0.1)

# Threshold: change in lambda
delta_lambda <- function(eps) { 
  delta_fn(para[1], para[2], para[3], para2[4], para[5], eps) 
}
beta0 <- uniroot(polynomial_fn, c(-10, 0), tol = 0.0001, m=para[1], 
                 s=para[2], e=para[3], l=para2[4], r=para[5])$root
hurdle0 <- beta0/(beta0-1+eps[1])
hurdle1 <- (beta0/(beta0-1+eps[2]))^2
hurdle2 <- (beta0-1+eps[3])/beta0
thr0_lambda <- (delta_lambda(eps[1])/para[5])* beta0/(beta0-1+eps[1])* (1 - I* 0.1)
thr1_lambda <- ((delta_lambda(eps[2])/para[5])* beta0/(beta0-1+eps[2]))^2* (1 - I* 0.1)
thr2_lambda <- ((delta_lambda(eps[3])/para[5])* beta0/(beta0-1+eps[3]))^(-1)* (1 - I* 0.1)

# Results
thr0 <- thr[1]
thr1 <- thr[2]
thr2 <- thr[3]

cat("The thresholds under ROA are: ", thr0, thr1, thr2, "\n")

cat("Mu: thresholds are: ", thr0_mu, thr1_mu, thr2_mu, "\n")
cat("Mu: change in thresholds are: ", 
    (thr0_mu-thr0)/thr0*100, (thr1_mu-thr1)/thr1*100, (thr2_mu-thr2)/thr2*100, "\n")

cat("Lambda: thresholds are: ", thr0_lambda, thr1_lambda, thr2_lambda, "\n")
cat("Lambda: change in thresholds are: ", 
    (thr0_lambda-thr0)/thr0*100, (thr1_lambda-thr1)/thr1*100, (thr2_lambda-thr2)/thr2*100, "\n")

cat("Eta: thresholds are: ", thr0_eta, thr1_eta, thr2_eta, "\n")
cat("Eta: change in thresholds are: ", 
    (thr0_eta-thr0)/thr0*100, (thr1_eta-thr1)/thr1*100, (thr2_eta-thr2)/thr2*100, "\n")

