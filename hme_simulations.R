source("hme optim.R")
source("value_of_rule.R")
source("estimators_function.R")

library(pracma)
library(mvtnorm)
library(boot)

set.seed(1934)

# Create a large super population
largeN <- 10^6

# True parameters for the Xs
lx2 <- 1 ; lx3 <- 1; lx4 <- 1;
lx5 <- 1 ; lx6 <- 1; lx7 <- .001

# True parameters for Ya0|X=x
alpha <- c(0, -3, -.5, 5, -1.5, -2, 0)

# True parameters for Ya1|X=x
beta <- alpha * -1

# True parameters for S|X=x
gamma <- c(0, 0, 0, 0, 0, 0, 10)

# True parameters for r(X)
delta <- c(.05, -.5, .5, -.5, .5, .0, 0)

# True parameters for Psi(X)  
zeta <- delta 

#
nvar <- 6
O <- randortho(nvar, type = "orthonormal")
D_mat <- matrix(0, nrow = nvar, ncol = nvar)
diag(D_mat) <- seq(1, by=.2, length=nvar)
Sigma <- O %*% D_mat %*% t(O)

# Functions that takes sample size and true value of the parameters and outputs
# a simulated data in dataframe format

dgp <- function(n, alpha, beta, gamma, delta, zeta, Sigma) {
  
  # Generate Xs
  nvar <- 6
  X_temp <- rmvnorm(n, mean = rep(0, nrow(Sigma)), sigma = Sigma)
  X_temp[,1:2] <- as.numeric(X_temp[,1:2] < 0)
  X_temp[,3:5] <- exp(X_temp[,3:5])
  X <- cbind(1, X_temp)
  
  # Generate Ss
  rho <- expit(X %*% gamma)
  S <- rbinom(n, 1, rho)
  
  # Generate r(X)
  r <-  as.numeric(X %*% delta < 0)
  
  # Generate potential outcomes A_s0 and A_s1
  psi <-  expit(X %*% zeta)
  A_s0 <- rbinom(n, 1, psi)
  
  A_s1 <- r
  
  # Generate potential outcomes Y_a0 and Y_a0
  Pr_Ya0 <- expit(X %*% alpha)
  Pr_Ya1 <- expit(X %*% beta)
  
  Y_a0 <- rbinom(n, 1, Pr_Ya0)
  Y_a1 <- rbinom(n, 1, Pr_Ya1)
  
  # Generate potential outcomes Y_s0 and Y_s0
  Y_s0 <- A_s0*Y_a1 + (1-A_s0)*Y_a0
  Y_s1 <- A_s1*Y_a1 + (1-A_s1)*Y_a0
  
  # Generate observed outcome Y (under partial implementation of an ITR)
  Y <- S*Y_s1 + (1-S)*Y_s0
  
  # Generate observed outcome A (under partial implementation of an ITR)
  A <- S*A_s1 + (1-S)*A_s0
  
  # Create C
  C <- as.numeric(r==A)
  
  return(data.frame(X=X, rho=rho, S=S, Pr_Ya1=Pr_Ya1, Pr_Ya0=Pr_Ya0, Y_s1=Y_s1, Y_s0=Y_s0, Y=Y, r=r, A=A, psi=psi, C=C))
}

seed <- 1934
set.seed(seed)

### Simulation setup ###

n_sim <- 1000
boot_resample <- 999
prec <- 6
sample_sizes <- c(1000, 5000, 20000)

### Psi and r very different i.e., zeta = delta ### 

superpop_r_psi_different <- dgp(largeN, alpha, beta, gamma, delta, zeta=delta, Sigma)
TrueVals_r_psi_different <- apply(superpop_r_psi_different[, c("Y_s1", "Y_s0", "Y", "Pr_Ya0", "Pr_Ya1")], 2, mean) 
TrueVals_r_psi_different
TrueEstimands_r_psi_different <- data.frame(True_MIG=TrueVals_r_psi_different["Y_s1"]-TrueVals_r_psi_different["Y"],
           True_ARE=TrueVals_r_psi_different["Y_s1"]-TrueVals_r_psi_different["Y_s0"],
           True_AIE=TrueVals_r_psi_different["Y"]-TrueVals_r_psi_different["Y_s0"])

TrueEstimands_r_psi_different

estimates_r_psi_different <- list()
lb_ci_r_psi_different <- list()
ub_ci_r_psi_different <- list()

ss_it <- 0
for (sample_size in sample_sizes){
ss_it <- ss_it + 1

estimates_r_psi_different[[ss_it]] <- matrix(NA, nrow=n_sim, ncol=4)
lb_ci_r_psi_different[[ss_it]] <- matrix(NA, nrow=n_sim, ncol=4)
ub_ci_r_psi_different[[ss_it]] <- matrix(NA, nrow=n_sim, ncol=4)

colnames(estimates_r_psi_different[[ss_it]]) <- c("MIG", "ARE_MOE", "ARE_COMB", "AIE")
colnames(lb_ci_r_psi_different[[ss_it]]) <- c("MIG", "ARE_MOE", "ARE_COMB", "AIE")
colnames(ub_ci_r_psi_different[[ss_it]]) <- c("MIG", "ARE_MOE", "ARE_COMB", "AIE")

x <- 1:nrow(superpop_r_psi_different);
y <- seq(from = 1, to = nrow(superpop_r_psi_different), by = sample_size);
indices <- sapply(split(x, f = findInterval(x = x, vec = y)), c)

for (i in 1:n_sim) {
  gc() # Free memory
  res <- boot(superpop_r_psi_different[indices[, i], ],
              estimators_boot,
              R=boot_resample,
              epsilon = 10^-2, maxit=10, verbose = FALSE, seed=seed)
  
  estimates_r_psi_different[[ss_it]][i,] <- round(res$t0, prec)
  ci <- round(emp_boot_ci(t0 = res$t0, t = res$t), prec)
  lb_ci_r_psi_different[[ss_it]][i,] <- ci[1,]
  ub_ci_r_psi_different[[ss_it]][i,] <- ci[2,]
}
}

save.image(file="end_first_scenario.RData")

### Psi is random i.e., zeta = rep(0,7) ### 

superpop_psi_random <- dgp(largeN, alpha, beta, gamma, delta, zeta=rep(0,7), Sigma)
TrueVals_psi_random <- apply(superpop_psi_random[, c("Y_s1", "Y_s0", "Y", "Pr_Ya0", "Pr_Ya1")], 2, mean) 
TrueVals_psi_random
TrueEstimands_psi_random <- data.frame(True_MIG=TrueVals_psi_random["Y_s1"]-TrueVals_psi_random["Y"],
                                            True_ARE=TrueVals_psi_random["Y_s1"]-TrueVals_psi_random["Y_s0"],
                                            True_AIE=TrueVals_psi_random["Y"]-TrueVals_psi_random["Y_s0"])
TrueEstimands_psi_random

estimates_psi_random <- list()
lb_ci_psi_random <- list()
ub_ci_psi_random <- list()

ss_it <- 0
for (sample_size in sample_sizes){
  ss_it <- ss_it + 1
  
  estimates_psi_random[[ss_it]] <- matrix(NA, nrow=n_sim, ncol=4)
  lb_ci_psi_random[[ss_it]] <- matrix(NA, nrow=n_sim, ncol=4)
  ub_ci_psi_random[[ss_it]] <- matrix(NA, nrow=n_sim, ncol=4)
  
  colnames(estimates_psi_random[[ss_it]]) <- c("MIG", "ARE_MOE", "ARE_COMB", "AIE")
  colnames(lb_ci_psi_random[[ss_it]]) <- c("MIG", "ARE_MOE", "ARE_COMB", "AIE")
  colnames(ub_ci_psi_random[[ss_it]]) <- c("MIG", "ARE_MOE", "ARE_COMB", "AIE")
  
  x <- 1:nrow(superpop_psi_random);
  y <- seq(from = 1, to = nrow(superpop_psi_random), by = sample_size);
  indices <- sapply(split(x, f = findInterval(x = x, vec = y)), c)
  
  for (i in 1:n_sim) {
    gc() # Free memory
    res <- boot(superpop_psi_random[indices[, i], ],
                estimators_boot,
                R=boot_resample,
                epsilon = 10^-2, maxit=10, verbose = FALSE, seed=seed)
    
    estimates_psi_random[[ss_it]][i,] <- round(res$t0, prec)
    ci <- round(emp_boot_ci(t0 = res$t0, t = res$t), prec)
    lb_ci_psi_random[[ss_it]][i,] <- ci[1,]
    ub_ci_psi_random[[ss_it]][i,] <- ci[2,]
  }
}

save.image(file="end_second_scenario.RData")

### Psi is similar to r i.e., zeta = - delta ### 

superpop_psi_similar <- dgp(largeN, alpha, beta, gamma, delta, zeta=- delta, Sigma)
TrueVals_psi_similar <- apply(superpop_psi_similar[, c("Y_s1", "Y_s0", "Y", "Pr_Ya0", "Pr_Ya1")], 2, mean) 
TrueVals_psi_similar
TrueEstimands_psi_similar <- data.frame(True_MIG=TrueVals_psi_similar["Y_s1"]-TrueVals_psi_similar["Y"],
                                       True_ARE=TrueVals_psi_similar["Y_s1"]-TrueVals_psi_similar["Y_s0"],
                                       True_AIE=TrueVals_psi_similar["Y"]-TrueVals_psi_similar["Y_s0"])
TrueEstimands_psi_similar

estimates_psi_similar <- list()
lb_ci_psi_similar <- list()
ub_ci_psi_similar <- list()

ss_it <- 0
for (sample_size in sample_sizes){
  ss_it <- ss_it + 1
  
  estimates_psi_similar[[ss_it]] <- matrix(NA, nrow=n_sim, ncol=4)
  lb_ci_psi_similar[[ss_it]] <- matrix(NA, nrow=n_sim, ncol=4)
  ub_ci_psi_similar[[ss_it]] <- matrix(NA, nrow=n_sim, ncol=4)
  
  colnames(estimates_psi_similar[[ss_it]]) <- c("MIG", "ARE_MOE", "ARE_COMB", "AIE")
  colnames(lb_ci_psi_similar[[ss_it]]) <- c("MIG", "ARE_MOE", "ARE_COMB", "AIE")
  colnames(ub_ci_psi_similar[[ss_it]]) <- c("MIG", "ARE_MOE", "ARE_COMB", "AIE")
  
  x <- 1:nrow(superpop_psi_similar);
  y <- seq(from = 1, to = nrow(superpop_psi_similar), by = sample_size);
  indices <- sapply(split(x, f = findInterval(x = x, vec = y)), c)
  
  for (i in 1:n_sim) {
    gc() # Free memory
    res <- boot(superpop_psi_similar[indices[, i], ],
                estimators_boot,
                R=boot_resample,
                epsilon = 10^-2, maxit=10, verbose = FALSE, seed=seed)
    
    estimates_psi_similar[[ss_it]][i,] <- round(res$t0, prec)
    ci <- round(emp_boot_ci(t0 = res$t0, t = res$t), prec)
    lb_ci_psi_similar[[ss_it]][i,] <- ci[1,]
    ub_ci_psi_similar[[ss_it]][i,] <- ci[2,]
  }
}

save.image(file="end_third_scenario.RData")


