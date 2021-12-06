library(boot)

mimic_si <- read.csv("/Users/francois/Desktop/github repos/ITR-evaluation/ITR-evaluation/mimic_si_preds.csv")

# Create concrdance variable
mimic_si$c <- with(mimic_si, a1 == r)

## Value of the rule APIWE E_Y^{s=1}

# Fit pronositc model
pr_mod <- glm(d60d ~ admission_age + a1*weight + a1*bun_k1 + a1*ph_k1 + a1*pot_k1 + a1*SOFA_24hours + a1*immunosuppressant, 
              data=mimic_si, family = "binomial")

# Get predictions under the recommended treatment
newdata <- mimic_si[, c("r", "admission_age", "weight", "SOFA_24hours", "immunosuppressant", "bun_k1", "ph_k1", "pot_k1" )]
colnames(newdata)[1] <- "a1"
mimic_si$pon_hat_d <- predict(pr_mod, newdata, type="response")


# Get predictions for the ITE
newdata <- mimic_si[, c("r", "admission_age", "weight", "SOFA_24hours", "immunosuppressant", "bun_k1", "ph_k1", "pot_k1" )]
newdata[,"a1"] <- 1
mimic_si$E_hat_Ya1_given_x <- predict(pr_mod, newdata, type="response")

newdata[,"a1"] <- 0
mimic_si$E_hat_Ya0_given_x <- predict(pr_mod, newdata, type="response")

mimic_si$ITE_hat <- with(mimic_si, E_hat_Ya1_given_x - E_hat_Ya0_given_x)

# Fit PS model
ps_mod <- glm(a1 ~ ph_k1 + bun_k1 + pot_k1, 
              data=mimic_si, family = "binomial")

# Get PS predictions
mimic_si$ps_hat <-  ps_mod$fitted.values

# Get predictions of compliance with the rule
mimic_si$ps_hat_d <-  mimic_si$ps_hat*mimic_si$a1 + (1-mimic_si$ps_hat)*(1-mimic_si$a1) 

# Compute AIPWE Vale of the rule r
EY_s1 <- mean( with(mimic_si, c*d60d / ps_hat_d - pon_hat_d*(c-ps_hat_d)/ps_hat_d ) )

## Compute EY
EY <- with(mimic_si, mean(d60d))
EY_s1 - EY

it <- 0 ; resamples <- 100;
res <- boot(mimic_si, are_aipw_new, R=resamples )
res
plot(res)
quantile(res$t, probs = c(.025, .975))

### Delta ITE 
res <- c()
sequence <- seq(0,.999, by=.1) 
for (i in sequence ) {
alpha <- i
mimic_si$rho <- with(mimic_si, (1- abs(r - ps_hat) )^(.5*log((alpha+1)/(1-alpha))) )
res <- c(res, mean( with(mimic_si, rho*(r - ps_hat ) * ITE_hat)) )
}

plot(res~sequence)

     