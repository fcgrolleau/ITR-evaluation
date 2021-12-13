library(MASS)
expit <- function(x) 1 / (1+exp(-x))

set.seed(657)
# Sample size
n <- 10000

# True parameters
lx1 <- 1 ; lx2 <- 2; lx3 <- 3;
lx4 <- 4 ; lx5 <- 3; lx6 <- 2

alpha <- c(0, -3, -.5, 5, -1.5, -2, 1)
beta <- alpha * c(1,1,1,3,1,1,2)
gamma <- c(-.1, 1, -1, 1, -1, 1, -1)*1

# Functions that takes sample size and true value of the parameters and outputs
# a simulated data in dataframe format

dgp <- function(n, Sigma, alpha, beta, gamma) {
# Generate Xs
X <- cbind(1, rexp(n, lx1), rexp(n, lx2), rexp(n, lx3), rexp(n, lx4), rexp(n, lx5), rexp(n, lx6))

# Generate Ss
rho <- expit(X %*% gamma)
S <- rbinom(n, 1, rho)

# Generate potential outcomes Y_s1 and Y_s0
Pr_Ys1 <- expit(X %*% alpha)
Pr_Ys0 <- expit(X %*% beta)

Y_s1 <- rbinom(n, 1, Pr_Ys1)
Y_s0 <- rbinom(n, 1, Pr_Ys0)

# Generate observed outcome Y (under partial implementation of an ITR)
Y <- S*Y_s1 + (1-S)*Y_s0

return(data.frame(X=X, rho=rho, S=S, Pr_Ys1=Pr_Ys1, Pr_Ys0=Pr_Ys0, Y_s1=Y_s1, Y_s0=Y_s0, Y=Y))
}

# Compute true value of E[Y_s1] and E[Y_s0]
largeN <- 10^6
TrueVals <- apply(dgp(largeN, Sigma, alpha, beta, gamma)[, c("Y_s1", "Y_s0", "Y")], 2, mean) 

dat <- dgp(n, Sigma, alpha, beta, gamma)
with(dat[1:1500,], plot(Pr_Ys1, Pr_Ys0, pch=16, col=rgb(1,0,0, alpha=.1), xlim=c(0,1), ylim=c(0,1)))
abline(lm(Pr_Ys0~Pr_Ys1, data=dat), lty=3)
abline(c(0,0), c(1,1))
table(dat$S)

### Test MOE No. 1
library(flexmix)

Model <- FLXMRglmfix(fixed = ~ -1, formula =  ~ -1 + X.2 + X.3 + X.4 + X.5 + X.6 + X.7, family = "binomial")
                         
concomitantModel <-FLXPmultinom(~ ~ 1 + X.2 + X.3 + X.4 + X.5 + X.6 + X.7)

Moe <- stepFlexmix(cbind(Y, 1 - Y) ~  1,
               k = 2, model = Model,
               concomitant = concomitantModel, data = dat, nrep = 3)

refit_Moe <- refit(Moe)
summary(refit_Moe)

exp_coef <- parameters(Moe, which="model")
gate_coef <- parameters(Moe, which="concomitant")
exp_coef
gate_coef

temp <- model.matrix(~ -1 + X.2 + X.3 + X.4 + X.5 + X.6 + X.7, dat)
EY_s0s <- c(mean(expit(temp %*% exp_coef[,1])), mean(expit(temp %*% exp_coef[,2])) )
EY_s0s
TrueVals

