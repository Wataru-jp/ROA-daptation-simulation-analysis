# Libraries
library(matrixStats) # for rowVars
require(utils) # for str
library(openxlsx)

# graph wd
setwd("C:/OwnCloud/A05/08_WP2_ROA/01_Irrigation_investment_numerical")
setwd("/Users/kodam1/ownCloud - wataru.kodama@uni-goettingen.de@owncloud.gwdg.de/A05/08_WP2_ROA/01_Irrigation_investment_numerical")

################################
### 0.1 Simulation Set Up
################################
para <- c(0.02, 0.2, 0.25, 0.05, 0.2) # baseline (mu, sigma, eta, lambda, rho)
eps <- c(0, 0.5, 2, 4) # list of epsilon
I <- 3
# Simulation parameters
n <- 10000                # Number of runs --> GO MUCH LARGER!
#Tlength <- 1              # Time period
#dt <- 0.01                # Time step
Tlength <- 30              # Time period
dt <- 1                    # Time step
t <- seq(0, Tlength, dt)  # steps
N <- length(t)            # Number of steps

# Model parameters
x0 <- 1                  # Initial profits (LeeZhao)/ payoffs / benefits (EisenackParsch)
x_a <- 1
r <- 0.1



################################
### 0.2 Threshold Set Up
################################
# Delta function
delta_fn <- function(m, s, e, l, r, eps) { 
  r + m*(1-eps) + 1/2*eps*(1-eps)*s^2 + l- l*(1-e)^(1-eps) }
delta <- function(eps) { 
  delta_fn(para[1], para[2], para[3], para[4], para[5], eps) 
}
# Polynomial function
polynomial_fn <- function(x, m, s, e, l, r) {
  1/2*(x-1)*x*s^2 - m*x - (r+l) + l*(1-e)^x
}

# Calculate threshold
beta <- uniroot(polynomial_fn, c(-10, 0), tol = 0.0001, m=para[1], 
                 s=para[2], e=para[3], l=para[4], r=para[5])$root
threshold <- ((delta(eps)/para[5])* beta/(beta-1+eps))^(1/(1-eps))* (1-I*r)
threshold_EPV <- ((delta(eps)/para[5]) )^(1/(1-eps))* (1-I*r)
threshold
threshold_EPV



################################
### 1. Simulation: Random draw
################################
para <- c(0.02, 0.2, 0.25, 0.05, 0.2) # baseline (mu, sigma, eta, lambda, rho)
# Initialize the pi matrix
x_c_t <- matrix(0, n, N)
x_c_t[, 1] <- x0

# Simulation initialization
set.seed(1404)      # WK Birthday

# Simulate GBM with Poisson Jump process
for (i in 1:n) {
  #W <- cumsum(rnorm(N - 1, mean = 0, sd = sqrt(dt)))    # Brownian motion
  #  J <- rpois(N - 1, para[4] * dt)                       # Number of jumps at each time step
  
  #jumps <- J * para[3] * x0                            # Jump sizes
  #cum_jumps <- cumsum(jumps)                            # Cumulative jumps
  
  # GBM with jumps
  #x_c_t[i, -1] <- x0 * exp((- para[1] - 0.5 * para[2]^2) * t[-1] + para[2] * W - cum_jumps)
  for (t in 1:(N-1)) {
    W <- rnorm(1, mean = 0, sd = sqrt(dt))
    J <- rpois(1, para[4] * dt)
    x_c_t[i, t+1] <- x_c_t[i, t] - para[1] * x_c_t[i, t] + W * para[2] * x_c_t[i, t] - para[3] * J * x_c_t[i, t]
  }
}


x_c_t <- matrix(0, n, N)
x_c_t[, 1] <- x0
for (i in 1:n) {
  for (t in 1:(N-1)) {
    W <- rnorm(1, mean = 0, sd = sqrt(dt))
    J <- rpois(1, para[4] * dt)
    x_c_t[i, t+1] <- x_c_t[i, t] - para[1]*x_c_t[i, t] + W*para[2]*x_c_t[i, t] - para[3] * J*x_c_t[i, t]
  }
}

# Results
rownames(x_c_t) <- paste0("Run_", 1:n)
colnames(x_c_t) <- t
t <- 1:N  


### Graph ###
setwd("C:/OwnCloud/A05/08_WP2_ROA/01_Irrigation_investment_numerical")
setwd("/Users/kodam1/ownCloud - wataru.kodama@uni-goettingen.de@owncloud.gwdg.de/A05/08_WP2_ROA/01_Irrigation_investment_numerical")

# Plot one draw
png("03_Programming/simulation/figure/random_draw.png", width = 5, height = 4, units = 'in', res = 500)
plot(t, x_c_t[1, ], type = "l", ylim = c(0, 1.5), 
     xlab = "Time", ylab = "Payoff")
abline(h=threshold_EPV[1], lty = 1)
abline(h=threshold[1], lty = 2)
abline(h=threshold[3], lty = 3)
legend("topright", cex = 0.7, lty=1:3,
       legend = c(expression(paste("EPV threshold (", epsilon," = 0)")), 
                  expression(paste("ROA threshold (", epsilon," = 0)")), 
                  expression(paste("ROA threshold (", epsilon," = 2)"))))
dev.off()

# Plot five draws
png("03_Programming/simulation/figure/random_draw2.png", width = 5, height = 4, units = 'in', res = 500)
plot(t, x_c_t[1, ], type = "l", ylim = c(0, 1.7), 
     xlab = "Time", ylab = "Payoff")
for (i in c(2, 3, 4)) {
  lines(t, x_c_t[i, ], lty = 1)
}
abline(h = threshold_EPV[1], lty = 1)
abline(h = threshold[1], lty = 2)
abline(h = threshold[3], lty = 3)
legend("topright", cex = 0.7, lty=1:3,
       legend = c(expression(paste("EPV threshold (", epsilon," = 0)")), 
                  expression(paste("ROA threshold (", epsilon," = 0)")), 
                  expression(paste("ROA threshold (", epsilon," = 2)"))))
dev.off()


############################################
### 2. Adaptation probability
############################################
para <- c(0.02, 0.2, 0.25, 0.05, 0.2) # baseline (mu, sigma, eta, lambda, rho)
eps <- c(0, 0.5, 2, 3) # list of epsilon

# EPV probability
for (e in 1:3) {
  assign(paste0("prob_EPV", e), list())
  adap_EPV <- apply(x_c_t, 1, FUN = function(x) {x < threshold_EPV[e]})
  for (i in 1:ncol(adap_EPV)) {
    if(sum(adap_EPV[,i])>0) {
      adap_EPV[,i][min(which(adap_EPV[,i]==TRUE)):length(adap_EPV[,i])] <- TRUE
    }
  }
  for (i in 1:N) {
    if (e == 1) {
      prob_EPV1[[i]] <- rowSums(adap_EPV)[i]/ncol(adap_EPV)*100
    } else if (e == 2) {
      prob_EPV2[[i]] <- rowSums(adap_EPV)[i]/ncol(adap_EPV)*100
    } else if (e == 3) {
      prob_EPV3[[i]] <- rowSums(adap_EPV)[i]/ncol(adap_EPV)*100
    }
  }
}

# ROA probability
for (e in 1:3) {
  assign(paste0("prob", e), list())
  adapt <- apply(x_c_t, 1, FUN = function(x) {x < threshold[e]})
  for (i in 1:ncol(adapt)) {
    if(sum(adapt[,i])>0) {
      adapt[,i][min(which(adapt[,i]==TRUE)):length(adapt[,i])] <- TRUE
    }
  }
  for (i in 1:N) {
    if (e == 1) {
      prob1[[i]] <- rowSums(adapt)[i]/ncol(adapt)*100
    } else if (e == 2) {
      prob2[[i]] <- rowSums(adapt)[i]/ncol(adapt)*100
    } else if (e == 3) {
      prob3[[i]] <- rowSums(adapt)[i]/ncol(adapt)*100
    } else if (e == 4) {
      prob4[[i]] <- rowSums(adapt)[i]/ncol(adapt)*100
    }
  }
}

plot(t, prob1, type='l')
plot(t, prob2, type='l')
plot(t, prob3, type='l')
plot(t, prob4, type='l')

### Graph ###
png("03_Programming/simulation/figure/probability.png", width = 5, height = 4, units = 'in', res = 500)
plot(t, prob_EPV1, type='l', lty=1, ylim=c(0, 100), 
     xlab='Time', ylab='Probability', col = "blue")
lines(t, prob1, lty=1)
lines(t, prob2, lty=2)
lines(t, prob3, lty=3)
legend("bottomright", cex = 0.8, lty = c(1, 1, 2, 3),
       col = c("blue", "black", "black", "black"),
       legend = c(expression(paste("EPV ", "(", epsilon," = 0)")),
                  expression(paste("ROA ", "(", epsilon," = 0)")),
                  expression(paste("ROA ", "(", epsilon," = 0.5)")),
                  expression(paste("ROA ", "(", epsilon," = 2)"))))
dev.off()


############################################
### 3. Sensitivity analysis
############################################
para <- c(0.02, 0.2, 0.25, 0.05, 0.2) # baseline (mu, sigma, eta, lambda, rho)
# expected loss
#para2 <- c(para[1] + para[3]*para[4], para[2] + para[3]/2, para[3] + para[3], para[4] + para[4])
# elasticity
#para2 <- c(para[1] + 0.025/para[4]*para[1], para[2] + 0.025/para[4]*para[2], para[3] + 0.025/para[4]*para[3], para[4] + 0.025) 
# expected loss & elasticity
para2 <- c(para[1] + 0.025*para[3], para[2] + 0.025/para[4]*para[2], para[3] + 0.025* para[3]/ para[4], para[4] + 0.025)
eps <- c(0, 0.5, 2) # list of epsilon

# Simulation parameters
n <- 10000                # Number of runs --> GO MUCH LARGER!
Tlength <- 30              # Time period
dt <- 1                    # Time step
t <- seq(0, Tlength, dt)  # steps
N <- length(t)            # Number of steps

# Simulation initialization
set.seed(1404)      # WK Birthday
variables <- c("mu", "sigma", "lambda", "eta")
for (var in variables) {
  assign(paste0("x_c_t_", var), matrix(x0, n, N))
}


# Simulation
for (i in 1:n) {
  #W <- cumsum(rnorm(N - 1, mean = 0, sd = sqrt(dt)))    # Brownian motion
  #J <- rpois(N - 1, para[4] * dt) * para[3] * x0       # Number of jumps and size at each time step
  #J_lambda <- rpois(N - 1, para2[4] * dt) * para[3] * x0 
  #J_eta <- rpois(N - 1, para[4] * dt) * para2[3] * x0       
  # GBM with jumps
  #x_c_t_mu[i, -1] <- x0 * exp((- para2[1] - 0.5 * (para[2])^2) * t[-1] + para[2] * W - cumsum(J))
  #x_c_t_sigma[i, -1] <- x0 * exp((- para[1] - 0.5 * (para2[2])^2) * t[-1] + para2[2] * W - cumsum(J))
  #x_c_t_lambda[i, -1] <- x0 * exp((- para[1] - 0.5 * (para[2])^2) * t[-1] + para[2] * W - cumsum(J_lambda))
  #x_c_t_eta[i, -1] <- x0 * exp((- para[1] - 0.5 * (para[2])^2) * t[-1] + para[2] * W - cumsum(J_eta))
  for (t in 1:(N-1)) {
    W <- rnorm(1, mean = 0, sd = sqrt(dt))
    J <- rpois(1, para[4] * dt)
    J2 <- rpois(1, para2[4] * dt)
    x_c_t_mu[i, t+1] <- x_c_t_mu[i, t] - para2[1] * x_c_t_mu[i, t] + W * para[2] * x_c_t_mu[i, t] - para[3] * J * x_c_t_mu[i, t]
    x_c_t_sigma[i, t+1] <- x_c_t_sigma[i, t] - para[1] * x_c_t_sigma[i, t] + W * para2[2] * x_c_t_sigma[i, t] - para[3] * J * x_c_t_sigma[i, t]
    x_c_t_lambda[i, t+1] <- x_c_t_lambda[i, t] - para[1] * x_c_t_lambda[i, t] + W * para[2] * x_c_t_lambda[i, t] - para[3] * J2 * x_c_t_lambda[i, t]
    x_c_t_eta[i, t+1] <- x_c_t_eta[i, t] - para[1] * x_c_t_eta[i, t] + W * para[2] * x_c_t_eta[i, t] - para2[3] * J * x_c_t_eta[i, t]
  }
}

#####change in mu #####
beta <- uniroot(polynomial_fn, c(-10, 0), tol = 0.0001, 
                m=para2[1], s=para[2], e=para[3], l=para[4], r=para[5])$root
delta <- function(eps) { 
  delta_fn(para2[1], para[2], para[3], para[4], para[5], eps) 
}
threshold_mu <- ((delta(eps)/para[5])* beta/(beta-1+eps))^(1/(1-eps))* (1-I*0.1)
# loop over epsilon
for (j in 1:3) {
  # Create adaptation matrix
  adaptation <- apply(x_c_t_mu, 1, FUN = function(x) {x < threshold_mu[j]})
  for (i in 1:ncol(adaptation)) {
    if(sum(adaptation[,i])>0) {
      adaptation[,i][min(which(adaptation[,i]==TRUE)):length(adaptation[,i])] <- TRUE
    }
  }
  # Create prob variable
  prob_mu <- list()
  for (i in 1:N) {
    prob_mu[[i]] <- rowSums(adaptation)[i]/ncol(adaptation)*100
  }
  # Assign prob variable to the corresponding epsilon
  assign(paste0("prob", j, "_mu"), prob_mu)
}

##### change in lambda #####
# threshold
beta <- uniroot(polynomial_fn, c(-10, 0), tol = 0.0001, 
                m=para[1], s=para[2], e=para[3], l=para2[4], r=para[5])$root
delta <- function(eps) { 
  delta_fn(para[1], para[2], para[3], para2[4], para[5], eps) 
}
threshold_lambda <- ((delta(eps)/para[5])* beta/(beta-1+eps))^(1/(1-eps))* (1-I*0.1)
# loop over epsilon
for (j in 1:3) {
  adaptation <- apply(x_c_t_lambda, 1, FUN = function(x) {x < threshold_lambda[j]})
  for (i in 1:ncol(adaptation)) {
    if(sum(adaptation[,i])>0) {
      adaptation[,i][min(which(adaptation[,i]==TRUE)):length(adaptation[,i])] <- TRUE
    }
  }
  prob_lambda <- list()
  for (i in 1:N) {
    prob_lambda[[i]] <- rowSums(adaptation)[i]/ncol(adaptation)*100
  }
  assign(paste0("prob", j, "_lambda"), prob_lambda)
}

##### change in eta #####
# threshold
beta <- uniroot(polynomial_fn, c(-10, 0), tol = 0.0001, m=para[1], 
                s=para[2], e=para2[3], l=para[4], r=para[5])$root
delta <- function(eps) { 
  delta_fn(para[1], para[2], para2[3], para[4], para[5], eps) 
}
threshold_eta <- ((delta(eps)/para[5])* beta/(beta-1+eps))^(1/(1-eps))* (1-I*0.1)
# loop over epsilon
for (j in 1:3) {
  # Create adaptation matrix
  adaptation <- apply(x_c_t_eta, 1, FUN = function(x) {x < threshold_eta[j]})
  for (i in 1:ncol(adaptation)) {
    if(sum(adaptation[,i])>0) {
      adaptation[,i][min(which(adaptation[,i]==TRUE)):length(adaptation[,i])] <- TRUE
    }
  }
  # Create prob variable
  prob_eta <- list()
  for (i in 1:N) {
    prob_eta[[i]] <- rowSums(adaptation)[i]/ncol(adaptation)*100
  }
  # Assign prob variable to the corresponding epsilon
  assign(paste0("prob", j, "_eta"), prob_eta)
}

#####change in sigma #####
# threshold
beta <- uniroot(polynomial_fn, c(-10, 0), tol = 0.0001, 
                m=para[1], s=para2[2], e=para[3], l=para[4], r=para[5])$root
delta <- function(eps) { 
  delta_fn(para[1], para2[2], para[3], para[4], para[5], eps) 
}
threshold_sigma <- ((delta(eps)/para[5])* beta/(beta-1+eps))^(1/(1-eps))* (1-I*0.1)
# loop over epsilon
for (j in 1:3) {
  # Create adaptation matrix
  adaptation <- apply(x_c_t_sigma, 1, FUN = function(x) {x < threshold_sigma[j]})
  for (i in 1:ncol(adaptation)) {
    if(sum(adaptation[,i])>0) {
      adaptation[,i][min(which(adaptation[,i]==TRUE)):length(adaptation[,i])] <- TRUE
    }
  }
  # Create prob variable
  prob_sigma <- list()
  for (i in 1:N) {
    prob_sigma[[i]] <- rowSums(adaptation)[i]/ncol(adaptation)*100
  }
  # Assign prob variable to the corresponding epsilon
  assign(paste0("prob", j, "_sigma"), prob_sigma)
}


### Graph ### 
t <- 1:N  

png("03_Programming/simulation/figure/sensitivity0.png", width = 5, height = 4, units = 'in', res = 500)
plot(t, prob1, type='l', lty=1, ylim=c(0, 100), xlab='Time', ylab='Probability')
lines(t, prob1_mu, lty = 2, col="blue")
lines(t, prob1_sigma, lty = 2, col="green")
lines(t, prob1_lambda, lty = 2, col="red")
lines(t, prob1_eta, lty = 2, col="orange")
legend("bottomright", cex = 0.8, lty=c(1, 2, 2, 2, 2), 
       title = c(expression(paste(epsilon, " = 0"))), 
       col=c("black", "blue", "green", "red", "orange"),
       legend = c(expression(paste("Baseline")), 
                  expression(paste(mu, "+", Delta, mu)),
                  expression(paste(sigma, "+", Delta, sigma)),
                  expression(paste(lambda, "+", Delta, lambda)),
                  expression(paste(eta, "+", Delta, eta))))
dev.off()


png("03_Programming/simulation/figure/sensitivity2.png", width = 5, height = 4, units = 'in', res = 500)
plot(t, prob3, type='l', lty=1, ylim=c(0, 100), xlab='Time', ylab='Probability')
lines(t, prob3_mu, lty = 2, col="blue")
lines(t, prob3_sigma, lty = 2, col="green")
lines(t, prob3_lambda, lty = 2, col="red")
lines(t, prob3_eta, lty = 2, col="orange")
legend("bottomright", cex = 0.8, lty=c(1, 2, 2, 2, 2), 
       title = c(expression(paste(epsilon, " = 2"))), 
       col=c("black", "blue", "green", "red", "orange"),
       legend = c(expression(paste("Baseline")), 
                  expression(paste(mu, "+", Delta, mu)),
                  expression(paste(sigma, "+", Delta, sigma)),
                  expression(paste(lambda, "+", Delta, lambda)),
                  expression(paste(eta, "+", Delta, eta))))
dev.off()

png("03_Programming/simulation/figure/sensitivity.png", width = 5, height = 4, units = 'in', res = 500)
plot(t, prob1, type='l', lty=1, ylim=c(0, 100), xlab='Time', ylab='Probability')
lines(t, prob3, lty = 2, col="black")
lines(t, prob3_mu, lty = 2, col="blue")
lines(t, prob3_sigma, lty = 2, col="green")
lines(t, prob3_lambda, lty = 2, col="red")
lines(t, prob3_eta, lty = 2, col="orange")
legend("bottomright", cex = 0.8, lty=c(1, 2, 2, 2, 2, 2), 
       col=c("black", "black", "blue", "green", "red", "orange"),
       legend = c(expression(paste("Baseline", " (", epsilon, " = 0)")), 
                  expression(paste("Baseline", " (", epsilon, " = 2)")), 
                  expression(paste(mu, "+", Delta, mu, " (", epsilon, " = 2)")),
                  expression(paste(sigma, "+", Delta, sigma, " (", epsilon, " = 2)")),
                  expression(paste(lambda, "+", Delta, lambda, " (", epsilon, " = 2)")),
                  expression(paste(eta, "+", Delta, eta, " (", epsilon, " = 2)"))))
dev.off()

############################################
### 3.1 Export to table
############################################
###  Changes in threshold compared to baseline
threshold2 <- c(threshold[1], threshold[2], threshold[3])
ct_mu <- (threshold_mu - threshold2 )/ threshold2 * 100
ct_sigma <- (threshold_sigma - threshold2 )/ threshold2 * 100
ct_lambda <- (threshold_lambda - threshold2 )/ threshold2 * 100
ct_eta <- (threshold_eta - threshold2 )/ threshold2 * 100

###  Changes in probability compared to baseline
for (i in 1:3) {
  # unlist
  assign(paste0("p_EPV", i), unlist(get(paste0("prob_EPV", i))))
  assign(paste0("p", i), unlist(get(paste0("prob", i))))
  assign(paste0("p", i, "_mu"), unlist(get(paste0("prob", i, "_mu"))))
  assign(paste0("p", i, "_sigma"), unlist(get(paste0("prob", i, "_sigma"))))
  assign(paste0("p", i, "_lambda"), unlist(get(paste0("prob", i, "_lambda"))))
  assign(paste0("p", i, "_eta"), unlist(get(paste0("prob", i, "_eta"))))
  # define changes in probability
  assign(paste0("cp", i, "_mu"), (get(paste0("p", i, "_mu")) - get(paste0("p", i)) ) / get(paste0("p", i)) * 100)
  assign(paste0("cp", i, "_sigma"), (get(paste0("p", i, "_sigma")) - get(paste0("p", i)) ) / get(paste0("p", i)) * 100)
  assign(paste0("cp", i, "_lambda"), (get(paste0("p", i, "_lambda")) - get(paste0("p", i)) ) / get(paste0("p", i)) * 100)
  assign(paste0("cp", i, "_eta"), (get(paste0("p", i, "_eta")) - get(paste0("p", i)) ) / get(paste0("p", i)) * 100)
}

###  Data frame for table
df_threshold <- data.frame(
  Epsilon = c(0, 1/2, 2),
  EPV = c(threshold_EPV[1], threshold_EPV[2], threshold_EPV[3]),
  ROA = c(threshold[1], threshold[2], threshold[3]),
  Mu = c(threshold_mu[1], threshold_mu[2], threshold_mu[3]),
  Mu_c = c(ct_mu[1], ct_mu[2], ct_mu[3]),
  Sigma = c(threshold_sigma[1], threshold_sigma[2], threshold_sigma[3]),
  Sigma_c = c(ct_sigma[1], ct_sigma[2], ct_sigma[3]),
  Lambda = c(threshold_lambda[1], threshold_lambda[2], threshold_lambda[3]),
  Lambda_c = c(ct_lambda[1], ct_lambda[2], ct_lambda[3]),
  Eta = c(threshold_eta[1], threshold_eta[2], threshold_eta[3]),
  Eta_c = c(ct_eta[1], ct_eta[2], ct_eta[3])
)

for (i in 1:3) {
  assign(paste0("df", i), data.frame(
    Year = c(5, 10, 15, 20, 25, 30),
    EPV = c(get(paste0("p_EPV", i))[6], get(paste0("p_EPV", i))[11], get(paste0("p_EPV", i))[16], 
            get(paste0("p_EPV", i))[21], get(paste0("p_EPV", i))[26], get(paste0("p_EPV", i))[31]),
    ROA = c(get(paste0("p", i))[6], get(paste0("p", i))[11], get(paste0("p", i))[16], 
            get(paste0("p", i))[21], get(paste0("p", i))[26], get(paste0("p", i))[31]),
    Mu = c(get(paste0("p", i, "_mu"))[6], get(paste0("p", i, "_mu"))[11], get(paste0("p", i, "_mu"))[16], 
            get(paste0("p", i, "_mu"))[21], get(paste0("p", i, "_mu"))[26], get(paste0("p", i, "_mu"))[31]),
    Mu_c = c(get(paste0("cp", i, "_mu"))[6], get(paste0("cp", i, "_mu"))[11], get(paste0("cp", i, "_mu"))[16], 
           get(paste0("cp", i, "_mu"))[21], get(paste0("cp", i, "_mu"))[26], get(paste0("cp", i, "_mu"))[31]),
    Sigma = c(get(paste0("p", i, "_sigma"))[6], get(paste0("p", i, "_sigma"))[11], get(paste0("p", i, "_sigma"))[16], 
            get(paste0("p", i, "_sigma"))[21], get(paste0("p", i, "_sigma"))[26], get(paste0("p", i, "_sigma"))[31]),
    Sigma_c = c(get(paste0("cp", i, "_sigma"))[6], get(paste0("cp", i, "_sigma"))[11], get(paste0("cp", i, "_sigma"))[16], 
              get(paste0("cp", i, "_sigma"))[21], get(paste0("cp", i, "_sigma"))[26], get(paste0("cp", i, "_sigma"))[31]),
    Lambda = c(get(paste0("p", i, "_lambda"))[6], get(paste0("p", i, "_lambda"))[11], get(paste0("p", i, "_lambda"))[16], 
              get(paste0("p", i, "_lambda"))[21], get(paste0("p", i, "_lambda"))[26], get(paste0("p", i, "_lambda"))[31]),
    Lambda_c = c(get(paste0("cp", i, "_lambda"))[6], get(paste0("cp", i, "_lambda"))[11], get(paste0("cp", i, "_lambda"))[16], 
               get(paste0("cp", i, "_lambda"))[21], get(paste0("cp", i, "_lambda"))[26], get(paste0("cp", i, "_lambda"))[31]),
    Eta = c(get(paste0("p", i, "_eta"))[6], get(paste0("p", i, "_eta"))[11], get(paste0("p", i, "_eta"))[16], 
               get(paste0("p", i, "_eta"))[21], get(paste0("p", i, "_eta"))[26], get(paste0("p", i, "_eta"))[31]),
    Eta_c = c(get(paste0("cp", i, "_eta"))[6], get(paste0("cp", i, "_eta"))[11], get(paste0("cp", i, "_eta"))[16], 
               get(paste0("cp", i, "_eta"))[21], get(paste0("cp", i, "_eta"))[26], get(paste0("cp", i, "_eta"))[31])
  )
  )
}
df_prob <- rbind(df1, df2, df3)

### Preparation
# Round the values to 2 decimal places
df_prob[, 2:11] <- round(df_prob[, 2:11], 2)
df_threshold[, 2:11] <- round(df_threshold[, 2:11], 3)
df_threshold[, c(5, 7, 9, 11)] <- round(df_threshold[, c(5, 7, 9, 11)], 2)

# Format the values
for (i in 2:11) {
  df_prob[, i] <- paste0(df_prob[, i], "%")
}
for (i in c(5, 7, 9, 11)) {
  df_prob[, i] <- paste0("(", df_prob[, i], ")")
  df_threshold[, i] <- paste0("(", df_threshold[, i], "%)")
}

### Export the data frame
# Create a new workbook
wb <- createWorkbook()
# Add a worksheet for each data frame
addWorksheet(wb, "Threshold")
addWorksheet(wb, "Probability")
# Write the data to the worksheets
writeData(wb, "Threshold", df_threshold)
writeData(wb, "Probability", df_prob)
# Save the workbook to a file
saveWorkbook(wb, "03_Programming/table/sensitivity.xlsx")



############################################
### 4. Different values of parameters
############################################
para <- c(0.03, 0.2, 0.25, 0.05, 0.2) # baseline (mu, sigma, eta, lambda, rho)
# range of parameters
cat("Max of drift rate is: ", para[5] - (para[2])^2 - para[3]* para[4]/ (1 - para[3]) )
cat("Max of volatility is: ", sqrt(para[5] - para[1] - para[3]* para[4]/ (1 - para[3])) )
cat("Max of frequency is: ", (1 - para[3])/ para[3]* (para[5] - para[1] - (para[2])^2) )
cat("Max of magnitude is: ", (para[5] - para[1] - (para[2])^2)/ (para[4] + para[5] - para[1] - (para[2])^2) )
# elasticity
para2 <- c(para[1] + 0.025/para[4]*para[1], para[2] + 0.025/para[4]*para[2], para[3] + 0.025/para[4]*para[3], para[4] + 0.025) 
para3 <- c(para[1] + 0.05/para[4]*para[1], para[2] + 0.05/para[4]*para[2], para[3] + 0.05/para[4]*para[3], para[4] + 0.05) 
para4 <- c(para[1] - 0.05/para[4]*para[1], para[2] - 0.05/para[4]*para[2], para[3] - 0.05/para[4]*para[3], para[4] - 0.05) 
# expected loss & elasticity
#para2 <- c(para[1] + para[3]*para[4], para[2] + para[2], para[3] + para[3], para[4] + para[4])


eps <- c(0, 0.5, 2) # list of epsilon
I <- 3

# Simulation parameters
n <- 10000                # Number of runs --> GO MUCH LARGER!
Tlength <- 30              # Time period
dt <- 1                    # Time step
t <- seq(0, Tlength, dt)  # steps
N <- length(t)            # Number of steps

# Simulation initialization
set.seed(1404)      # WK Birthday
variables <- c("mu", "sigma", "lambda", "eta")
x_ct <- matrix(0, n, N)
x_ct[, 1] <- x0


### 1. Relationship with drift rate ### 
drift <- seq(0, 0.05, 0.01) # drift rate range
for (i in 0:2) {
  assign(paste0("thr_list", i, "a"), list())
  assign(paste0("prob_list", i), list())
  assign(paste0("prob2_list", i), list())
  assign(paste0("prob5_mu", i), list())
  assign(paste0("prob10_mu", i), list())
} 
delta_list <- list()
hurdle_list <- list()
beta <- list()

for (x in 1:length(drift)){
  # threshold
  delta_list <- delta_fn(drift[x], para[2], para[3], para[4], para[5], eps)
  beta <- uniroot(polynomial_fn, c(-10, 0), tol = 0.0001, m = drift[x], 
                             s = para[2], e = para[3], l = para[4], r = para[5])$root
  hurdle_list <- (beta/ (beta - 1 + eps))^ (1/(1-eps))
  threshold_mu <- hurdle_list* (delta_list/para[5])^ (1/ (1 - eps))* (1 - I* 0.1)
  thr_list0a[[x]] <- hurdle_list[1]* (delta_list[1]/para[5])* (1 - I* 0.1)
  thr_list1a[[x]] <- hurdle_list[2]* (delta_list[2]/para[5])^2* (1 - I* 0.1)
  thr_list2a[[x]] <- hurdle_list[3]* (para[5]/delta_list[3])* (1 - I* 0.1)
  # simulation
  for (i in 1:n) {
    for (t in 1:(N-1)) {
      W <- rnorm(1, mean = 0, sd = sqrt(dt))
      J <- rpois(1, para[4] * dt)
      x_ct[i, t+1] <- x_ct[i, t] - drift[x] * x_ct[i, t] + W * para[2] * x_ct[i, t] - para[3] * J * x_ct[i, t]
    }
  }
  # loop over epsilon
  for (j in 1:3) {
    # Create adaptation matrix
    adaptation <- apply(x_ct, 1, FUN = function(x) {x < threshold_mu[j]})
    for (i in 1:ncol(adaptation)) {
      if(sum(adaptation[,i])>0) {
        adaptation[,i][min(which(adaptation[,i]==TRUE)):length(adaptation[,i])] <- TRUE
      }
    }
    # Create prob variable
    prob <- list()
    for (i in 1:N) {
      prob[[i]] <- rowSums(adaptation)[i]/ncol(adaptation)*100
    }
    # Assign prob variable to the corresponding epsilon
    assign(paste0("probability", j), prob)
  }
  prob_list0[[x]] <- probability1[6]
  prob_list1[[x]] <- probability2[6]
  prob_list2[[x]] <- probability3[6]
  prob2_list0[[x]] <- probability1[11]
  prob2_list1[[x]] <- probability2[11]
  prob2_list2[[x]] <- probability3[11]
}

# make them vector
prob5_mu0 <- unlist(lapply(prob_list0, `[[`, 1))
prob5_mu1 <- unlist(lapply(prob_list1, `[[`, 1))
prob5_mu2 <- unlist(lapply(prob_list2, `[[`, 1))
prob10_mu0 <- unlist(lapply(prob2_list0, `[[`, 1))
prob10_mu1 <- unlist(lapply(prob2_list1, `[[`, 1))
prob10_mu2 <- unlist(lapply(prob2_list2, `[[`, 1))



### 2. Relationship with volatility ### 
volatil <- seq(0, 0.5, 0.005) # volatility range
for (i in 0:2) {
  assign(paste0("thr_list", i, "b"), list())
  assign(paste0("prob_list", i), list())
  assign(paste0("prob2_list", i), list())
  assign(paste0("prob5_sigma", i), list())
  assign(paste0("prob10_sigma", i), list())
} 
delta_list <- list()
hurdle_list <- list()
beta <- list()

for (x in 1:length(volatil)){
  delta_list <- delta_fn(para[1], volatil[x], para[3], para[4], para[5], eps)
  beta <- uniroot(polynomial_fn, c(-10, 0), tol = 0.0001, m = para[1], 
                        s = volatil[x], e = para[3], l = para[4], r = para[5])$root
  hurdle_list <- (beta/ (beta - 1 + eps))^ (1/(1-eps))
  threshold_sigma <- hurdle_list* (delta_list/para[5])^ (1/ (1 - eps))* (1 - I* 0.1)
  thr_list0b[[x]] <- hurdle_list[1]* (delta_list[1]/para[5])* (1 - I* 0.1)
  thr_list1b[[x]] <- hurdle_list[2]* (delta_list[2]/para[5])^2* (1 - I* 0.1)
  thr_list2b[[x]] <- hurdle_list[3]* (para[5]/delta_list[3])* (1 - I* 0.1)
  # simulation
  for (i in 1:n) {
    for (t in 1:(N-1)) {
      W <- rnorm(1, mean = 0, sd = sqrt(dt))
      J <- rpois(1, para[4] * dt)
      x_ct[i, t+1] <- x_ct[i, t] - para[1] * x_ct[i, t] + W * volatil[x] * x_ct[i, t] - para[3] * J * x_ct[i, t]
    }
  }
  # loop over epsilon
  for (j in 1:3) {
    # Create adaptation matrix
    adaptation <- apply(x_ct, 1, FUN = function(x) {x < threshold_sigma[j]})
    for (i in 1:ncol(adaptation)) {
      if(sum(adaptation[,i])>0) {
        adaptation[,i][min(which(adaptation[,i]==TRUE)):length(adaptation[,i])] <- TRUE
      }
    }
    # Create prob variable
    prob <- list()
    for (i in 1:N) {
      prob[[i]] <- rowSums(adaptation)[i]/ncol(adaptation)*100
    }
    # Assign prob variable to the corresponding epsilon
    assign(paste0("probability", j), prob)
  }
  prob_list0[[x]] <- probability1[6]
  prob_list1[[x]] <- probability2[6]
  prob_list2[[x]] <- probability3[6]
  prob2_list0[[x]] <- probability1[11]
  prob2_list1[[x]] <- probability2[11]
  prob2_list2[[x]] <- probability3[11]
}

# make them vector
prob5_sigma0 <- unlist(lapply(prob_list0, `[[`, 1))
prob5_sigma1 <- unlist(lapply(prob_list1, `[[`, 1))
prob5_sigma2 <- unlist(lapply(prob_list2, `[[`, 1))
prob10_sigma0 <- unlist(lapply(prob2_list0, `[[`, 1))
prob10_sigma1 <- unlist(lapply(prob2_list1, `[[`, 1))
prob10_sigma2 <- unlist(lapply(prob2_list2, `[[`, 1))


### 3. Relationship with magnitude of weather extremes ### 
etalist <- seq(0.0, 0.5, 0.01) # magnitude of shock range
for (i in 0:2) {
  assign(paste0("thr_list", i, "c"), list())
  assign(paste0("prob_list", i), list())
  assign(paste0("prob2_list", i), list())
  assign(paste0("prob5_eta", i), list())
  assign(paste0("prob10_eta", i), list())
} 
delta_list <- list()
hurdle_list <- list()
beta <- list()

for (x in 1:length(etalist)){
  delta_list <- delta_fn(para[1], para[2], etalist[x], para[4], para[5], eps)
  beta <- uniroot(polynomial_fn, c(-10, 0), tol = 0.0001, m = para[1], 
                  s = para[2], e = etalist[x], l = para[4], r = para[5])$root
  hurdle_list <- (beta/ (beta - 1 + eps))^ (1/(1-eps))
  threshold_eta <- hurdle_list* (delta_list/para[5])^ (1/ (1 - eps))* (1 - I* 0.1)
  thr_list0c[[x]] <- hurdle_list[1]* (delta_list[1]/para[5])* (1 - I* 0.1)
  thr_list1c[[x]] <- hurdle_list[2]* (delta_list[2]/para[5])^2* (1 - I* 0.1)
  thr_list2c[[x]] <- hurdle_list[3]* (para[5]/delta_list[3])* (1 - I* 0.1)
  # simulation
  for (i in 1:n) {
    for (t in 1:(N-1)) {
      W <- rnorm(1, mean = 0, sd = sqrt(dt))
      J <- rpois(1, para[4] * dt)
      x_ct[i, t+1] <- x_ct[i, t] - para[1] * x_ct[i, t] + W * para[2] * x_ct[i, t] - etalist[x] * J * x_ct[i, t]
    }
  }
  # loop over epsilon
  for (j in 1:3) {
    # Create adaptation matrix
    adaptation <- apply(x_ct, 1, FUN = function(x) {x < threshold_eta[j]})
    for (i in 1:ncol(adaptation)) {
      if(sum(adaptation[,i])>0) {
        adaptation[,i][min(which(adaptation[,i]==TRUE)):length(adaptation[,i])] <- TRUE
      }
    }
    # Create prob variable
    prob <- list()
    for (i in 1:N) {
      prob[[i]] <- rowSums(adaptation)[i]/ncol(adaptation)*100
    }
    # Assign prob variable to the corresponding epsilon
    assign(paste0("probability", j), prob)
  }
  prob_list0[[x]] <- probability1[6]
  prob_list1[[x]] <- probability2[6]
  prob_list2[[x]] <- probability3[6]
  prob2_list0[[x]] <- probability1[11]
  prob2_list1[[x]] <- probability2[11]
  prob2_list2[[x]] <- probability3[11]
}

# make them vector
prob5_eta0 <- unlist(lapply(prob_list0, `[[`, 1))
prob5_eta1 <- unlist(lapply(prob_list1, `[[`, 1))
prob5_eta2 <- unlist(lapply(prob_list2, `[[`, 1))
prob10_eta0 <- unlist(lapply(prob2_list0, `[[`, 1))
prob10_eta1 <- unlist(lapply(prob2_list1, `[[`, 1))
prob10_eta2 <- unlist(lapply(prob2_list2, `[[`, 1))


### 4. Relationship with frequency of weather extreme ### 
lamlist <- seq(0, 0.2, 0.01) # magnitude of shock range
for (i in 0:2) {
  assign(paste0("thr_list", i, "d"), list())
  assign(paste0("prob_list", i), list())
  assign(paste0("prob2_list", i), list())
  assign(paste0("prob5_lambda", i), list())
  assign(paste0("prob10_lambda", i), list())
} 
delta_list <- list()
hurdle_list <- list()
beta <- list()

for (x in 1:length(lamlist)){
  delta_list <- delta_fn(para[1], para[2], para[3], lamlist[x], para[5], eps)
  beta <- uniroot(polynomial_fn, c(-10, 0), tol = 0.0001, m = para[1], 
                  s = para[2], e = para[3], l = lamlist[x], r = para[5])$root
  hurdle_list <- (beta/ (beta - 1 + eps))^ (1/(1-eps))
  threshold_lam <- hurdle_list* (delta_list/para[5])^ (1/ (1 - eps))* (1 - I* 0.1)
  thr_list0d[[x]] <- hurdle_list[1]* (delta_list[1]/para[5])* (1 - I* 0.1)
  thr_list1d[[x]] <- hurdle_list[2]* (delta_list[2]/para[5])^2* (1 - I* 0.1)
  thr_list2d[[x]] <- hurdle_list[3]* (para[5]/delta_list[3])* (1 - I* 0.1)
  # simulation
  for (i in 1:n) {
    for (t in 1:(N-1)) {
      W <- rnorm(1, mean = 0, sd = sqrt(dt))
      J <- rpois(1, lamlist[x] * dt)
      x_ct[i, t+1] <- x_ct[i, t] - para[1] * x_ct[i, t] + W * para[2] * x_ct[i, t] - para[3] * J * x_ct[i, t]
    }
  }
  # loop over epsilon
  for (j in 1:3) {
    # Create adaptation matrix
    adaptation <- apply(x_ct, 1, FUN = function(x) {x < threshold_lam[j]})
    for (i in 1:ncol(adaptation)) {
      if(sum(adaptation[,i])>0) {
        adaptation[,i][min(which(adaptation[,i]==TRUE)):length(adaptation[,i])] <- TRUE
      }
    }
    # Create prob variable
    prob <- list()
    for (i in 1:N) {
      prob[[i]] <- rowSums(adaptation)[i]/ncol(adaptation)*100
    }
    # Assign prob variable to the corresponding epsilon
    assign(paste0("probability", j), prob)
  }
  prob_list0[[x]] <- probability1[6]
  prob_list1[[x]] <- probability2[6]
  prob_list2[[x]] <- probability3[6]
  prob2_list0[[x]] <- probability1[11]
  prob2_list1[[x]] <- probability2[11]
  prob2_list2[[x]] <- probability3[11]
}

# make them vector
prob5_lambda0 <- unlist(lapply(prob_list0, `[[`, 1))
prob5_lambda1 <- unlist(lapply(prob_list1, `[[`, 1))
prob5_lambda2 <- unlist(lapply(prob_list2, `[[`, 1))
prob10_lambda0 <- unlist(lapply(prob2_list0, `[[`, 1))
prob10_lambda1 <- unlist(lapply(prob2_list1, `[[`, 1))
prob10_lambda2 <- unlist(lapply(prob2_list2, `[[`, 1))


############################################
### 4.1 Graph: threshold
############################################
setwd("C:/OwnCloud/A05/08_WP2_ROA/01_Irrigation_investment_numerical")
setwd("/Users/kodam1/ownCloud - wataru.kodama@uni-goettingen.de@owncloud.gwdg.de/A05/08_WP2_ROA/01_Irrigation_investment_numerical")
#par(mfrow = c(2, 2), cex.main = 0.8)
#par(mfrow=c(1,1), cex.main = 1)

# First plot
png("03_Programming/figure/threhold_mu.png", width = 5, height = 4, units = 'in', res = 500)
plot(drift, thr_list0a, type='l', lty=1, ylim=c(0.51,0.60), xlab='', ylab='')
lines(drift, thr_list1a, lty=2)
lines(drift, thr_list2a, lty=3)
#points(drift, thr_list0a, pch = 10, cex = 0.5)
legend("bottomright", cex = 0.8, lty=1:3,
       legend = c(expression(paste(epsilon," = 0")), 
                  expression(paste(epsilon," = 0.5 ")), 
                  expression(paste(epsilon," = 2"))))
dev.off()

# Second plot
png("03_Programming/figure/threhold_sigma.png", width = 5, height = 4, units = 'in', res = 500)
plot(volatil, thr_list0b, type='l', lty=1, ylim=c(0.32,0.7), xlab='', ylab='')
lines(volatil, thr_list1b, lty=2)
lines(volatil, thr_list2b, lty=3)
#points(volatil, thr_list0b, pch = 10, cex = 0.5)
legend("topright", cex = 0.8, lty=1:3,
       legend = c(expression(paste(epsilon," = 0")), 
                  expression(paste(epsilon," = 0.5 ")), 
                  expression(paste(epsilon," = 2"))))
dev.off()

# Third plot
png("03_Programming/figure/threhold_lambda.png", width = 5, height = 4, units = 'in', res = 500)
plot(lamlist, thr_list0d, type='l', lty=1, ylim=c(0.54,0.59), xlab='', ylab='')
lines(lamlist, thr_list1d, lty=2)
lines(lamlist, thr_list2d, lty=3)
#points(lamlist, thr_list0d, pch = 10, cex = 0.5)
legend("bottomright", cex = 0.8, lty=1:3,
       legend = c(expression(paste(epsilon," = 0")), 
                  expression(paste(epsilon," = 0.5 ")), 
                  expression(paste(epsilon," = 2"))))
dev.off()

# Fourth plot
png("03_Programming/figure/threhold_eta.png", width = 5, height = 4, units = 'in', res = 500)
plot(etalist, thr_list0c, type='l', lty=1, ylim=c(0.52,0.59), xlab='', ylab='')
lines(etalist, thr_list1c, lty=2)
lines(etalist, thr_list2c, lty=3)
#points(etalist, thr_list0c, pch = 10, cex = 0.5)
legend("topright", cex = 0.8, lty=1:3,
       legend = c(expression(paste(epsilon," = 0")), 
                  expression(paste(epsilon," = 0.5 ")), 
                  expression(paste(epsilon," = 2"))))
dev.off()


############################################
### 4.2 Graph: probability
############################################
setwd("C:/OwnCloud/A05/08_WP2_ROA/01_Irrigation_investment_numerical")
setwd("/Users/kodam1/ownCloud - wataru.kodama@uni-goettingen.de@owncloud.gwdg.de/A05/08_WP2_ROA/01_Irrigation_investment_numerical")

# First plot
png("03_Programming/figure/probability_mu.png", width = 5, height = 4, units = 'in', res = 500)
plot(drift, prob5_mu0, type='l', lty=1, ylim=c(15, 85), xlab='', ylab='')
#points(drift, prob5_mu0, pch = 10, cex = 0.5)
lines(drift, prob5_mu1, lty=2)
lines(drift, prob5_mu2, lty=3)
lines(drift, prob10_mu0, lty=1, col = "blue")
lines(drift, prob10_mu1, lty=2, col = "blue")
lines(drift, prob10_mu2, lty=3, col = "blue")
legend("bottomright", cex = 0.8, lty=1:3,
       legend = c(expression(paste(epsilon," = 0")), 
                  expression(paste(epsilon," = 0.5 ")), 
                  expression(paste(epsilon," = 2"))))
dev.off()

# Second plot
png("03_Programming/figure/probability_sigma.png", width = 5, height = 4, units = 'in', res = 500)
plot(volatil, prob5_sigma0, type='l', lty=1, ylim=c(15, 85), xlab='', ylab='')
#points(volatil, prob5_sigma0, pch = 10, cex = 0.5)
lines(volatil, prob5_sigma1, lty=2)
lines(volatil, prob5_sigma2, lty=3)
lines(volatil, prob10_sigma0, lty=1, col = "blue")
lines(volatil, prob10_sigma1, lty=2, col = "blue")
lines(volatil, prob10_sigma2, lty=3, col = "blue")
legend("bottomright", cex = 0.8, lty=1:3,
       legend = c(expression(paste(epsilon," = 0")), 
                  expression(paste(epsilon," = 0.5 ")), 
                  expression(paste(epsilon," = 2"))))
dev.off()

# Third plot
png("03_Programming/figure/probability_lambda.png", width = 5, height = 4, units = 'in', res = 500)
plot(lamlist, prob5_lambda0, type='l', lty=1, ylim=c(15, 85), xlab='', ylab='')
#points(lamlist, prob5_lambda0, pch = 10, cex = 0.5)
lines(lamlist, prob5_lambda1, lty=2)
lines(lamlist, prob5_lambda2, lty=3)
lines(lamlist, prob10_lambda0, lty=1, col = "blue")
lines(lamlist, prob10_lambda1, lty=2, col = "blue")
lines(lamlist, prob10_lambda2, lty=3, col = "blue")
legend("bottomright", cex = 0.8, lty=1:3,
       legend = c(expression(paste(epsilon," = 0")), 
                  expression(paste(epsilon," = 0.5 ")), 
                  expression(paste(epsilon," = 2"))))
dev.off()


# Fourth plot
png("03_Programming/figure/probability_eta.png", width = 5, height = 4, units = 'in', res = 500)
plot(etalist, prob5_eta0, type='l', lty=1, ylim=c(15, 85), xlab='', ylab='')
#points(etalist, prob5_eta0, pch = 10, cex = 0.5)
lines(etalist, prob5_eta1, lty=2)
lines(etalist, prob5_eta2, lty=3)
lines(etalist, prob10_eta0, lty=1, col = "blue")
lines(etalist, prob10_eta1, lty=2, col = "blue")
lines(etalist, prob10_eta2, lty=3, col = "blue")
legend("bottomright", cex = 0.8, lty=1:3,
       legend = c(expression(paste(epsilon," = 0")), 
                  expression(paste(epsilon," = 0.5 ")), 
                  expression(paste(epsilon," = 2"))))
dev.off()
