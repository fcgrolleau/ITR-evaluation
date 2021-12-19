library(MASS)
library(mgcv)
expit <- function(x) 1 / (1+exp(-x))

set.seed(353)
# Sample size
n <- 10000

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
X <- cbind(1, rexp(n, lx2), rexp(n, lx3), rexp(n, lx4), rexp(n, lx5), rexp(n, lx6), rnorm(n, mean=0, sd = 1))

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


# Compute true value of E[Y_s1] and E[Y_s0]
largeN <- 10^6
TrueVals <- apply(dgp(largeN, alpha, beta, gamma, delta, zeta)[, c("Y_s1", "Y_s0", "Y", "Pr_Ya0", "Pr_Ya1")], 2, mean) 
TrueVals
dat <- dgp(10000, alpha, beta, gamma, delta, zeta)
table(dat$S)
gamma

hist(dat$rho)
with(dat, plot(X.7, rho))

### Test MOE No. 1
library(flexmix)
library(splines)

nodes <- 4

Model <- FLXMRglmfix(fixed = ~ 1, formula = ~ -1 + X.2 * X.3 * X.4 * X.5 * X.6, family = "binomial" )


concomitantModel <- FLXPmultinom(~ -1 + S)

Moe <- stepFlexmix(cbind(Y, 1 - Y) ~  -1,
                   model = Model, k=2,
                   concomitant = concomitantModel, data = dat, nrep = 3)


refit_Moe <- refit(Moe)
summary(refit_Moe)

exp_coef <- parameters(Moe, which="model")
gate_coef <- parameters(Moe, which="concomitant")
exp_coef
gate_coef

temp1 <- model.matrix(~ 1 + X.2 * X.3 * X.4 * X.5 * X.6, dat)
temp2 <- model.matrix(~ 1 + X.2 + X.3 + X.4 + X.5 + X.6, dat)

EY_s0s <- c(mean(expit(temp1 %*% exp_coef[,1] )), mean(expit(temp1 %*% exp_coef[,2] )) )
EY_s0s
TrueVals

