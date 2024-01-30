## Numerical exercise of ROA model: climate change adaptation investment
## Author: Wataru Kodama & Stefan Seifert
## Date: 15/3/2023
## Update: 04/01/2024

# Libraries
library(matrixStats) # for rowVars
require(utils) # for str
library(openxlsx)

############################################### ###
### 1. Investment threshold under ROA and EPV  ####
############################################### ###
# set baseline prarameter values
para <- c(0.0290, 0.1092, 0.300, 0.0600, 0.1000) #  (mu, sigma, eta, lambda, rho)
#para <- c(0.0309, 0.1233, 0.300, 0.0600, 0.1000) # old parameter set
eps <- c(0, 0.5, 2) # list of epsilon
I <- 2 # switching cost

# Delta function: convenience yield
delta_fn <- function(m, s, e, l, r, eps) { 
  r + m*(1-eps) + 1/2*eps*(1-eps)*s^2 + l - l*(1-e)^(1-eps) 
  }
# Delta function: convenience yield (approximated version)
delta_fn2 <- function(m, s, e, l, r, eps) { 
  r + (m+e*l)*(1-eps) + 1/2*eps*(1-eps)*s^2 
  }

# Polynomial function
polynomial_fn <- function(x, m, s, e, l, r, eps) {
  m2 = m*(1-eps) + (s^2)/2*eps*(1-eps)
  s2 = s*(1-eps)
  e2 = 1 - (1-e)^(1-eps)
  1/2*(x-1)*x*(s2)^2 - (m2)*x - l + l*(1-(e2))^x - r
  }

# ROA adaptation threshold function
ROA_fn <- function(m, s, e, l, r, eps) {
  d <- delta_fn(m, s, e, l, r, eps)
  # Find a positive root when eps > 1
  if (eps > 1) {
    b <- uniroot(polynomial_fn, c(0, 10), tol = 0.0001, m = m, 
                  s = s, e = e, l = l, r = r, eps = eps)$root
    ((d/r)* b/(b-1))^(1/(1-eps))* (1-I*r)
  }
  # Find a negative root when eps < 1
  else if (eps < 1) {
    b <- uniroot(polynomial_fn, c(-10, 0), tol = 0.0001, m = m, 
                  s = s, e = e, l = l, r = r, eps = eps)$root
    ((d/r)* b/(b-1))^(1/(1-eps))* (1-I*r)
  }
}

# ROA adaptation threshold function (closed form solution)
ROA_fn2 <- function(m, s, e, l, r, eps) {
  d <- delta_fn2(m, s, e, l, r, eps) # convenience yield
  a <- 1/2 + (m+e*l)/(s^2)
  b <- (a - sqrt(a^2 + 2*r/(s^2)))/(1-eps)
 ((d/r)* b/(b-1))^(1/(1-eps))* (1-I*r)
}

# EPV adaptation threshold function
EPV_fn <- function(m, s, e, l, r, eps) {
  d <- delta_fn(m, s, e, l, r, eps) # convenience yield
  ((d/r)* 1)^(1/(1-eps))* (1-I*r)
}


round(EPV_fn(para[1], para[2], para[3], para[4], para[5], eps), 3)

round(ROA_fn(para[1], para[2], para[3], para[4], para[5], eps[1]), 3)
round(ROA_fn(para[1], para[2], para[3], para[4], para[5], eps[2]), 3)
round(ROA_fn(para[1], para[2], para[3], para[4], para[5], eps[3]), 3)


############################ ### 
### 2. Adaptation threshold ####
############################ ###
# set baseline prarameter values
para <- c(0.0290, 0.1092, 0.300, 0.0600, 0.1000) #  (mu, sigma, eta, lambda, rho)
#para <- c(0.0309, 0.1233, 0.300, 0.0600, 0.1000) # old parameter set
eps <- c(0, 0.5, 2) # list of epsilon
I <- 2 # switching cost

### Relationship with drift rate ### 
drift <- seq(0, 0.05, 0.001) # drift rate range
threshold_list1 <- matrix(0, 3, length(drift)) # store thresholds in this matrix
for (i in 1:length(drift)){
  for (j in 1:3) {
    threshold_list1[j, i] <- ROA_fn(drift[i], para[2], para[3], para[4], para[5], eps[j])
  }}

### Relationship with volatility ### 
volatil <- seq(0, 0.25, 0.001) # volatility range
threshold_list2 <- matrix(0, 3, length(volatil)) # store thresholds in this matrix
for (i in 1:length(volatil)){
  for (j in 1:3) {
    threshold_list2[j, i] <- ROA_fn(para[1], volatil[i], para[3], para[4], para[5], eps[j])
  }}

### Relationship with magnitude of weather extremes ### 
etalist <- seq(0.0, 0.40, 0.001) # magnitude of shock range
threshold_list3 <- matrix(0, 3, length(etalist)) # store thresholds in this matrix
for (i in 1:length(etalist)){
  for (j in 1:3) {
    threshold_list3[j, i] <- ROA_fn(para[1], para[2], etalist[i], para[4], para[5], eps[j])
  }}
# closed-form solution
threshold_list3b <- matrix(0, 3, length(etalist)) # store thresholds in this matrix
for (i in 1:length(etalist)){
  for (j in 1:3) {
    threshold_list3b[j, i] <- ROA_fn2(para[1], para[2], etalist[i], para[4], para[5], eps[j])
  }}

### Relationship with frequency of weather extreme ### 
lamlist <- seq(0, 0.12, 0.001) # magnitude of shock range
threshold_list4 <- matrix(0, 3, length(lamlist))
for (i in 1:length(lamlist)){
  for (j in 1:3) {
    threshold_list4[j, i] <- ROA_fn(para[1], para[2], para[3], lamlist[i], para[5], eps[j])
  }}


### ** Figure #### 
## Path
setwd("/Users/kodam1/ownCloud - wataru.kodama@uni-goettingen.de@owncloud.gwdg.de/08_WP2_ROA/01_Irrigation_investment_numerical")
setwd("C:/OwnCloud/A05/08_WP2_ROA/01_Irrigation_investment_numerical")
setwd("D:/OC/DETECT/DETECT-SP/A05/08_WP2_ROA/01_Irrigation_investment_numerical")

## Combined figure: Adaptation threshold
png("03_Programming/simulation/figure/Figure1_v04.png", width = 5.5, height = 4, units = 'in', res = 500)
layout(matrix(c(1,2,3,4,5,5), nrow=3, byrow=TRUE), heights=c(5, 5, 0.75))
par(mai=c(0.45,0.5,0.3,0.1))

# drift
plot(drift, threshold_list1[1,], type='l', lty=1, xlab='', ylab='', ylim = c(0.63, 0.72),
     main = expression(bold(paste("Drift rate ", "\u03bc"))), cex.main = 1, yaxt="n")
axis(2, at=seq(0.64,0.72,0.02), las = 2)
lines(drift, threshold_list1[2,], lty=2)
lines(drift, threshold_list1[3,], lty=3)

# volatility
plot(volatil, threshold_list2[1,], type='l', lty=1, xlab='', ylab='', ylim = c(0.53, 0.76),
     main = expression(bold(paste("Volatility ", "\u03c3"))), cex.main = 1, yaxt="n")
axis(2, at=seq(0.55,0.75,0.05), labels=c("0.55","0.60","0.65","0.70","0.75"), las = 2)
lines(volatil, threshold_list2[2,], lty=2)
lines(volatil, threshold_list2[3,], lty=3)

# eta
plot(etalist, threshold_list3[1,], type='l', lty=1, xlab='', ylab='', ylim = c(0.66, 0.72),
     main = expression(bold(paste("Magnitude ", "\u03b7"))), cex.main = 1, yaxt="n", xaxt = "n")
axis(1, at=seq(0,0.40,0.10), labels=c("0","0.1","0.2","0.3","0.4"), las = 1)
axis(2, at=seq(0.66,0.72,0.02), las = 2)
lines(etalist, threshold_list3[2,], lty=2)
lines(etalist, threshold_list3[3,], lty=3)

# lambda
plot(lamlist, threshold_list4[1,], type='l', lty=1, xlab='', ylab='', ylim = c(0.66, 0.72),
     main = expression(bold(paste("Frequency ", "\u03bb"))), cex.main = 1, yaxt="n", xaxt = "n")
axis(1, at=seq(0,0.12,0.03), labels=c("0","0.03","0.06","0.09","0.12"), las = 1)
axis(2, at=seq(0.66,0.72,0.02), las = 2)
lines(lamlist, threshold_list4[2,], lty=2)
lines(lamlist, threshold_list4[3,], lty=3)

# legend
par(mai=c(0,0,0,0))
plot.new()
legend(x="center", ncol=3, lty=1:3,
       legend = c(expression(paste(epsilon," = 0")), 
                  expression(paste(epsilon," = 0.5")), 
                  expression(paste(epsilon," = 2"))),
       lwd = 1.5)

dev.off()


## Comparison: approximated closed form vs numerical solution
png("03_Programming/simulation/figure/FigureA1.png", width = 5, height = 4, units = 'in', res = 500)
plot(etalist, threshold_list3[1,], type='l', lty=1, xlab='', ylab='', ylim = c(0.655, 0.73),
     main = expression(bold(paste("Magnitude ", "\u03b7"))), cex.main = 1, yaxt="n")
axis(2, at=seq(0.66,0.73,0.02), labels=c("0.66","0.68","0.70","0.72"), las = 2)
lines(etalist, threshold_list3[2,], lty=2)
lines(etalist, threshold_list3[3,], lty=3)
lines(etalist, threshold_list3b[1,], lty=1, col="lightblue")
lines(etalist, threshold_list3b[2,], lty=2, col="lightblue")
lines(etalist, threshold_list3b[3,], lty=3, col="lightblue")
legend(x="top", ncol=3, lty=1:3,
       legend = c(expression(paste(epsilon," = 0")), 
                  expression(paste(epsilon," = 0.5")), 
                  expression(paste(epsilon," = 2"))),
       cex = 0.7)
dev.off()




############################## ###
### 3. Adaptation probability ####
############################## ###
para <- c(0.0290, 0.1092, 0.300, 0.0600, 0.1000) #  (mu, sigma, eta, lambda, rho)
#para <- c(0.0309, 0.1233, 0.300, 0.0600, 0.1000) # old parameter set
I <- 2

# Simulation parameters
n <- 100000                # Number of runs --> GO MUCH LARGER!
Tlength <- 5              # Time period
dt <- 1                    # Time step
t <- seq(0, Tlength, dt)  # steps
N <- length(t)            # Number of steps
x0 <- 1

# Simulation initialization
set.seed(1404)      # WK Birthday


### Relationship with drift rate ### 
drift <- seq(0, 0.05, 0.001) # drift rate range
## Probability list
prob_list_drift <- matrix(0, 3, length(drift)) # store probabilities in this matrix
adaptation <- list()
for (s in 1:length(drift)){
  # simulation
  x_c_t <- matrix(0, n, N)
  x_c_t[, 1] <- x0
  for (j in 1:n) {
    for (t in 1:(N-1)) {
      W <- rnorm(1, mean = 0, sd = sqrt(dt))
      J <- rpois(1, para[4] * dt)
      x_c_t[j, t+1] <- x_c_t[j, t] - drift[s] * x_c_t[j, t] + W * para[2] * x_c_t[j, t] - para[3] * J * x_c_t[j, t]
    }
    print(paste0("j=",j))}
  # Initialize a list of adaptation matrix
  adaptation[[s]] <- list() 
  for (e in 1:3) {
    # Create adaptation matrix
    adaptation[[s]][[e]] <- apply(x_c_t, 1, FUN = function(x) {x < threshold_list1[e, s]})
    for (j in 1:ncol(adaptation[[s]][[e]])) {
      if (sum(adaptation[[s]][[e]][, j]) > 0) {
        adaptation[[s]][[e]][, j][min(which(adaptation[[s]][[e]][, j] == TRUE)):length(adaptation[[s]][[e]][, j])] <- TRUE
      }}
    # Estimate prob
    prob_list_drift[e, s] <- rowSums(adaptation[[s]][[e]])[6] / ncol(adaptation[[s]][[e]]) * 100
    #print(paste0("j=",j))}
  }}


### Relationship with volatility ### 
volatil <- seq(0, 0.25, 0.001) # volatility range
## Probability list
prob_list_volatil <- matrix(0, 3, length(volatil)) # store probabilities in this matrix
adaptation <- list()
for (s in 1:length(volatil)){
  # simulation
  x_c_t <- matrix(0, n, N)
  x_c_t[, 1] <- x0
  for (j in 1:n) {
    for (t in 1:(N-1)) {
      W <- rnorm(1, mean = 0, sd = sqrt(dt))
      J <- rpois(1, para[4] * dt)
      x_c_t[j, t+1] <- x_c_t[j, t] - para[1] * x_c_t[j, t] + W * volatil[s] * x_c_t[j, t] - para[3] * J * x_c_t[j, t]
    }}
  # Initialize a list of adaptation matrix
  adaptation[[s]] <- list() 
  for (e in 1:3) {
    # Create adaptation matrix
    adaptation[[s]][[e]] <- apply(x_c_t, 1, FUN = function(x) {x < threshold_list2[e, s]})
    for (j in 1:ncol(adaptation[[s]][[e]])) {
      if (sum(adaptation[[s]][[e]][, j]) > 0) {
        adaptation[[s]][[e]][, j][min(which(adaptation[[s]][[e]][, j] == TRUE)):length(adaptation[[s]][[e]][, j])] <- TRUE
      }}
    # Estimate prob
    prob_list_volatil[e, s] <- rowSums(adaptation[[s]][[e]])[6] / ncol(adaptation[[s]][[e]]) * 100
  }}


### Relationship with magnitude of weather extremes ### 
etalist <- seq(0.0, 0.40, 0.001) # magnitude of shock range
## Probability list
prob_list_eta <- matrix(0, 3, length(etalist)) # store probabilities in this matrix
adaptation <- list()
for (s in 1:length(etalist)){
  # simulation
  x_c_t <- matrix(0, n, N)
  x_c_t[, 1] <- x0
  for (j in 1:n) {
    for (t in 1:(N-1)) {
      W <- rnorm(1, mean = 0, sd = sqrt(dt))
      J <- rpois(1, para[4] * dt)
      x_c_t[j, t+1] <- x_c_t[j, t] - para[1] * x_c_t[j, t] + W * para[2] * x_c_t[j, t] - etalist[s] * J * x_c_t[j, t]
    }}
  # Initialize a list of adaptation matrix
  adaptation[[s]] <- list() 
  for (e in 1:3) {
    # Create adaptation matrix
    adaptation[[s]][[e]] <- apply(x_c_t, 1, FUN = function(x) {x < threshold_list3[e, s]})
    for (j in 1:ncol(adaptation[[s]][[e]])) {
      if (sum(adaptation[[s]][[e]][, j]) > 0) {
        adaptation[[s]][[e]][, j][min(which(adaptation[[s]][[e]][, j] == TRUE)):length(adaptation[[s]][[e]][, j])] <- TRUE
      }}
    # Estimate prob
    prob_list_eta[e, s] <- rowSums(adaptation[[s]][[e]])[6] / ncol(adaptation[[s]][[e]]) * 100
  }}


### Relationship with frequency of weather extreme ### 
lamlist <- seq(0, 0.12, 0.001) # magnitude of shock range
## Probability list
prob_list_lambda <- matrix(0, 3, length(lamlist)) # store probabilities in this matrix
adaptation <- list()
for (s in 1:length(lamlist)){
  # simulation
  x_c_t <- matrix(0, n, N)
  x_c_t[, 1] <- x0
  for (j in 1:n) {
    for (t in 1:(N-1)) {
      W <- rnorm(1, mean = 0, sd = sqrt(dt))
      J <- rpois(1, lamlist[s] * dt)
      x_c_t[j, t+1] <- x_c_t[j, t] - para[1] * x_c_t[j, t] + W * para[2] * x_c_t[j, t] - para[3] * J * x_c_t[j, t]
    }}
  # Initialize a list of adaptation matrix
  adaptation[[s]] <- list() 
  for (e in 1:3) {
    # Create adaptation matrix
    adaptation[[s]][[e]] <- apply(x_c_t, 1, FUN = function(x) {x < threshold_list4[e, s]})
    for (j in 1:ncol(adaptation[[s]][[e]])) {
      if (sum(adaptation[[s]][[e]][, j]) > 0) {
        adaptation[[s]][[e]][, j][min(which(adaptation[[s]][[e]][, j] == TRUE)):length(adaptation[[s]][[e]][, j])] <- TRUE
      }}
    # Estimate prob
    prob_list_lambda[e, s] <- rowSums(adaptation[[s]][[e]])[6] / ncol(adaptation[[s]][[e]]) * 100
  }}


### ** Figure ####
## Path
setwd("/Users/kodam1/ownCloud - wataru.kodama@uni-goettingen.de@owncloud.gwdg.de/A05/08_WP2_ROA/01_Irrigation_investment_numerical")
setwd("C:/OwnCloud/A05/08_WP2_ROA/01_Irrigation_investment_numerical")
setwd("D:/OC/DETECT/DETECT-SP/A05/08_WP2_ROA/01_Irrigation_investment_numerical")

## Combined figure
png(paste0("03_Programming/simulation/figure/Figure2_v04_N=",n,".png"), width = 5, height = 4, units = 'in', res = 500)
#png(paste0("03_Programming/simulation/figure/Figure2_v02.png"), width = 5, height = 4, units = 'in', res = 500)
layout(matrix(c(1,2,3,4,5,5), nrow=3, byrow=TRUE), heights=c(5, 5, 1))
par(mai=c(0.5,0.5,0.3,0.1))

# drift
plot(drift, prob_list_drift[1, ], type='l', lty=1, xlab='', ylab='', ylim = c(15, 60),
     main = expression(bold(paste("Drift rate ", "\u03bc"))), cex.main = 1, yaxt="n")
axis(2, at=seq(20,60,10), labels=c(paste0(seq(20,60,10),"%")), las = 2)
lines(drift, prob_list_drift[2, ], lty=2)
lines(drift, prob_list_drift[3, ], lty=3)

# volatility
plot(volatil, prob_list_volatil[1, ], type='l', lty=1, xlab='', ylab='', ylim = c(22, 55),
     main = expression(bold(paste("Volatility ", "\u03c3"))), cex.main = 1, yaxt="n")
axis(2, at=seq(25,55,10), labels=c(paste0(seq(25,55,10),"%")), las = 2)
lines(volatil, prob_list_volatil[2, ], lty=2)
lines(volatil, prob_list_volatil[3, ], lty=3)

# eta
plot(etalist, prob_list_eta[1, ], type='l', lty=1, xlab='', ylab='', ylim = c(25, 45),
     main = expression(bold(paste("Magnitude ", "\u03b7"))), cex.main = 1, yaxt="n")
axis(2, at=seq(20,45,5), labels=c(paste0(seq(20,45,5),"%")), las = 2)
lines(etalist, prob_list_eta[2, ], lty=2)
lines(etalist, prob_list_eta[3, ], lty=3)

# lambda
plot(lamlist, prob_list_lambda[1, ], type='l', lty=1, xlab='', ylab='', ylim = c(25, 55),
     main = expression(bold(paste("Frequency ", "\u03bb"))), cex.main = 1, yaxt="n", xaxt="n")
axis(2, at=seq(25,55,10), labels=c(paste0(seq(25,55,10),"%")), las = 2)
axis(1, at=seq(0,0.12,0.03), labels=c("0","0.03","0.06","0.09","0.12"), las = 1)
lines(lamlist, prob_list_lambda[2, ], lty=2)
lines(lamlist, prob_list_lambda[3, ], lty=3)

# legend
par(mai=c(0,0,0,0))
plot.new()
legend(x="center", ncol=3, lty=1:3,
       legend = c(expression(paste(epsilon," = 0")), 
                  expression(paste(epsilon," = 0.5")), 
                  expression(paste(epsilon," = 2"))))
dev.off()




############################## ###
### 4. Austria case ####
############################## ###
## 4.1. Quantify the expected loss and GBM parameters ##
## Path
setwd("/Users/kodam1/ownCloud - wataru.kodama@uni-goettingen.de@owncloud.gwdg.de/08_WP2_ROA/01_Irrigation_investment_numerical")
setwd("C:/OwnCloud/A05/08_WP2_ROA/01_Irrigation_investment_numerical")
setwd("D:/OC/DETECT/DETECT-SP/A05/08_WP2_ROA/01_Irrigation_investment_numerical")

## ** Data ####
# Temperature mean for each model/ SSP
Austria_temp <- matrix(NA, 11, 5) # store mean temperature for each model
colnames(Austria_temp) <- c("SSP 1-1.9", "SSP 1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5")
for (i in 1:11) {
  file_path <- paste0("03_Programming/simulation/data/2040_model", i, ".csv")
  data <- read.csv(file_path)
  data <- data[c(3:9), -1] # Remove first row & Keep March - September temprature
  for (j in 1:5) {
    Austria_temp[i, j] <- colMeans(data)[[j]]
  }
}

# Rainfall mean for each model/ SSP
Austria_rain <- matrix(NA, 11, 5) # store mean rainfall for each model
colnames(Austria_rain) <- c("SSP 1-1.9", "SSP 1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5")
for (i in 1:11) {
  file_path <- paste0("03_Programming/simulation/data/2040_rain_model", i, ".csv")
  data <- read.csv(file_path)
  data <- data[c(3:9), -1] # Remove first row & Keep March - September temprature
  for (j in 1:5) {
    Austria_rain[i, j] <- colMeans(data)[[j]]
  }
}

## ** Expected loss ####
# Expected loss for each model/ SSP
x0 <- 815.62 # initial value
Austria_loss <- matrix(NA, 11, 5) # store expected loss for each model
Austria_x0 <- matrix(NA, 11, 5) # store expected loss % for each model
colnames(Austria_loss) <- c("SSP 1-1.9", "SSP 1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5")
colnames(Austria_x0) <- c("SSP 1-1.9", "SSP 1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5")
sign2 <- function(x) {ifelse(sign(x) == 1, 1, 0)}
for (i in 1:11) {
  for (j in 1:5) {
    x1 = Austria_temp[[i, j]]
    x2 = Austria_rain[[i, j]]
    loss <- 56.58*x1 + 0.13*(x1)^2 - 23.05*abs(x2)*sign2(x2) - 8.13*abs(x2)*(1 - sign2(x2)) + 9.54*(x2)^2*sign2(x2) + 6.25*(x2)^2*(1 - sign2(x2))
    Austria_loss[i, j] <- loss/x0*100
    Austria_x0[i, j] <- x0 - loss
  }
}

## drop SSP 1-2.6
Austria_loss <- Austria_loss[, -2]
Austria_x <- Austria_x0[, -2]
Austria_temp <- Austria_temp[, -2]
Austria_rain <- Austria_rain[, -2]

## Summary statistics: Table 2
T2 <- matrix(0, 11, 5, byrow = TRUE) # store mean, 10th and 90th percentile 
colnames(T2) <- c("SSP 1-1.9", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5", "All SSP")
rownames(T2) <- c("Mean temp", "T10th", "T90th", "Mean rain", "R10th", "R90th",
                       "Mean loss", "L10th", "L90th", "drift rate", "volatility")
# For each SSP
for (i in 1:4) {
  T2[1, i] <- round(mean(Austria_temp[, i]), 2) # mean
  T2[2, i] <- round(quantile(Austria_temp[, i], c(.10)), 2) # 10th percentile
  T2[3, i] <- round(quantile(Austria_temp[, i], c(.90)), 2) # 90th percentile
  T2[4, i] <- round(mean(Austria_rain[, i]), 2)
  T2[5, i] <- round(quantile(Austria_rain[, i], c(.10)), 2)
  T2[6, i] <- round(quantile(Austria_rain[, i], c(.90)), 2)
  T2[7, i] <- round(mean(Austria_loss[, i]), 2)
  T2[8, i] <- round(quantile(Austria_loss[, i], c(.10)), 2)
  T2[9, i] <- round(quantile(Austria_loss[, i], c(.90)), 2)
}
# All SSP
T2[1, 5] <- round(mean(Austria_temp), 2) # mean
T2[2, 5] <- round(quantile(Austria_temp, c(.10)), 2) # 10th percentile
T2[3, 5] <- round(quantile(Austria_temp, c(.90)), 2) # 90th percentile
T2[4, 5] <- round(mean(Austria_rain), 2)
T2[5, 5] <- round(quantile(Austria_rain, c(.10)), 2)
T2[6, 5] <- round(quantile(Austria_rain, c(.90)), 2)
T2[7, 5] <- round(mean(Austria_loss), 2)
T2[8, 5] <- round(quantile(Austria_loss, c(.10)), 2)
T2[9, 5] <- round(quantile(Austria_loss, c(.90)), 2)

## Drift rate and volatility: Table 2
years <- 20 # years from the reference period
# For each SSP
for (i in 1:4) {
  x1 <- quantile(Austria_x[, i], c(.10)) # 10th quantile of the expected payoff (x0 - loss)
  x2 <- quantile(Austria_x[, i], c(.90)) # 90th quantile of the expected payoff (x0 - loss)
  sigma <- log(x2/x1)/ (2*sqrt(years)*qnorm(0.90))
  mu <- (log(x2/x0) + log(x1/x0))/ (years*2) + sigma*sigma/2
  T2[10, i] <- round(mu, 4)
  T2[11, i] <- round(sigma, 4)
}
# All SSP
x1 <- quantile(Austria_x, c(.10)) # 10th quantile of the expected payoff (x0 - loss)
x2 <- quantile(Austria_x, c(.90)) # 90th quantile of the expected payoff (x0 - loss)
sigma <- log(x2/x1)/ (2*sqrt(years)*qnorm(0.90))
mu <- (log(x2/x0) + log(x1/x0))/ (years*2) + sigma*sigma/2
T2[10, 5] <- round(mu, 4)
T2[11, 5] <- round(sigma, 4)
# Table 2
T2


## ** Thresholds ####
## 4.2. Adaptation threshold and probability ##
para <- c(0.0290, 0.1092, 0.300, 0.0600, 0.1000) #  (mu, sigma, eta, lambda, rho)
#para <- c(0.0309, 0.1233, 0.300, 0.0600, 0.1000) # old parameter set
drift <- c(-T2[10, 1], -T2[10, 2], -T2[10, 3], -T2[10, 4], -T2[10, 5]) # for diff SSP
volat <- c(T2[11, 1], T2[11, 2], T2[11, 3], T2[11, 4], T2[11, 5]) # for diff SSP

## Thresholds
threshold_list <- matrix(0, 6, 5) # threshold list
colnames(threshold_list) <- c("SSP 1-1.9", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5", "All SSP")
rownames(threshold_list) <- c("ROA eps=0", "ROA eps=0.5", "ROA eps=2", 
                              "EPV eps=0", "EPV eps=0.5", "EPV eps=2")
# ROA thresholds
for (i in 1:5) {
  threshold_list[1, i] <- ROA_fn(drift[i], volat[i], para[3], para[4], para[5], eps[1]) # eps = 0
  threshold_list[2, i] <- ROA_fn(drift[i], volat[i], para[3], para[4], para[5], eps[2]) # eps = 0.5
  threshold_list[3, i] <- ROA_fn(drift[i], volat[i], para[3], para[4], para[5], eps[3]) # eps = 2
}
# EPV thresholds
for (i in 1:5) {
  threshold_list[4, i] <- EPV_fn(drift[i], volat[i], para[3], para[4], para[5], eps[1]) # eps = 0
  threshold_list[5, i] <- EPV_fn(drift[i], volat[i], para[3], para[4], para[5], eps[2]) # eps = 0.5
  threshold_list[6, i] <- EPV_fn(drift[i], volat[i], para[3], para[4], para[5], eps[3]) # eps = 2
}
threshold_list

### ** Simulation #### 
# Simulation parameters
n <- 100000               # Number of runs --> GO MUCH LARGER! @Stefan
Tlength <- 10              # Time period
dt <- 1                    # Time step
t <- seq(0, Tlength, dt)  # steps
N <- length(t)            # Number of steps
I <- 2
x0 <- 1

# Simulation initialization
set.seed(1404)      # WK Birthday
variables <- c("mu", "sigma", "lambda", "eta")

for (i in 0:5) {
  assign(paste0("xt", i), matrix(0, n, N))
}
xt <- list(xt1, xt2, xt3, xt4, xt5)
for (i in 1:5) {
  xt[[i]][, 1] <- x0
}

# Simulate GBM with Poisson Jump process
# STS: HERE I REPLACED X_C_T WITH xt[[s]][i, t]!!
# NOT SURE IF CORRECT
for (s in 1:5) {
  for (i in 1:n) {
    for (t in 1:(N-1)) {
      W <- rnorm(1, mean = 0, sd = sqrt(dt))
      J <- rpois(1, para[4] * dt)
      xt[[s]][i, t+1] <- xt[[s]][i, t] - drift[s] * xt[[s]][i, t] + W * volat[s] * xt[[s]][i, t] - para[3] * J * xt[[s]][i, t]
    }
  }
}


## Adaptation probability ## 
adaptation <- list()
for (s in 1:5) {
  adaptation[[s]] <- list()  # Initialize a list to store adaptation matrices for each 's'
  
  for (i in 1:3) {
    # Create adaptation matrix
    adaptation[[s]][[i]] <- apply(xt[[s]], 1, FUN = function(x) {x < threshold_list[i, s]})
    
    for (j in 1:ncol(adaptation[[s]][[i]])) {
      if (sum(adaptation[[s]][[i]][, j]) > 0) {
        adaptation[[s]][[i]][, j][min(which(adaptation[[s]][[i]][, j] == TRUE)):length(adaptation[[s]][[i]][, j])] <- TRUE
      }
    }
    
    # Create prob variable
    prob <- list()
    for (j in 1:N) {
      prob[[j]] <- rowSums(adaptation[[s]][[i]])[j] / ncol(adaptation[[s]][[i]]) * 100
    }
  
    # Assign prob variable to the corresponding epsilon
    assign(paste0("prob_", s, "_", i), prob)
  }
}

# Probability list
prob_list1 <- matrix(0, 3, 5) # 5-year
prob_list2 <- matrix(0, 3, 5) # 10-year
colnames(prob_list1) <- c("SSP 1-1.9", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5", "All SSP")
colnames(prob_list2) <- c("SSP 1-1.9", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5", "All SSP")
# overwrite the threshold values
for (s in 1:5) {
  for (i in 1:3) {
    prob_list1[i, s] <- get(paste0("prob_", s, "_", i))[[6]]
    prob_list2[i, s] <- get(paste0("prob_", s, "_", i))[[11]]
  }
}

## Table 3
T3 <- matrix(0, 12, 5) # threshold list
colnames(T3) <- c("SSP 1-1.9", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5", "All SSP")
rownames(T3) <- c("ROA eps=0", "EPV", "ROA eps=0.5", "EPV", "ROA eps=2", "EPV", 
                  "P5 eps=0", "P5 eps=0.5", "P5 eps=2",
                  "P10 eps=0", "P10 eps=0.5", "P10 eps=2")
T3[1, ] <- round(threshold_list[1, ], 3)
T3[2, ] <- round(threshold_list[4, ], 3)
T3[3, ] <- round(threshold_list[2, ], 3)
T3[4, ] <- round(threshold_list[5, ], 3)
T3[5, ] <- round(threshold_list[3, ], 3)
T3[6, ] <- round(threshold_list[6, ], 3)
T3[7, ] <- prob_list1[1, ]
T3[8, ] <- prob_list1[2, ]
T3[9, ] <- prob_list1[3, ]
T3[10, ] <- prob_list2[1, ]
T3[11, ] <- prob_list2[2, ]
T3[12, ] <- prob_list2[3, ]
T3
