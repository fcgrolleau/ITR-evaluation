source("hme optim.R")

library(pracma)
library(mvtnorm)

set.seed(353)

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

# Functions that takes sample size and true value of the parameters and outputs
# a simulated data in dataframe format

dgp <- function(n, alpha, beta, gamma, delta, zeta) {
  # Generate Xs
  
  nvar <- 6
  O <- randortho(nvar, type = "orthonormal")
  D_mat <- matrix(0, nrow = nvar, ncol = nvar)
  diag(D_mat) <- seq(1, by=.2, length=nvar)
  Sigma <- O %*% D_mat %*% t(O)
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

superpop <- dgp(largeN, alpha, beta, gamma, delta, zeta)

TrueVals <- apply(superpop[, c("Y_s1", "Y_s0", "Y", "Pr_Ya0", "Pr_Ya1")], 2, mean) 
TrueVals







