library(boot)

mimic_si <- read.csv("/Users/francois/Desktop/github repos/ITR-evaluation/ITR-evaluation/mimic_si_preds.csv")

# Create concrdance variable
mimic_si$c <- with(mimic_si, a1 == r)

## Value of the rule APIWE E_Y^{s=1}

# Fit pronositc model
pr_mod <- glm(d60d ~ a1 + admission_age + weight + SOFA_24hours + immunosuppressant + uo_k1 + bun_k1 + ph_k1 + pot_k1, 
              data=mimic_si, family = "binomial")

mimic_si$pon_hat <-  pr_mod$fitted.values


### PROBABILY WRONG HERE.. ###
# Get predictions under the recommended treatment
newdata <- mimic_si[, c("r", "admission_age", "weight", "SOFA_24hours", "immunosuppressant", "uo_k1", "bun_k1", "ph_k1", "pot_k1" )]
colnames(newdata)[1] <- "a1"
mimic_si$pon_hat_d <- predict(pr_mod, newdata, type="response")

# Fit PS model
ps_mod <- glm(a1 ~ pot_k1 + ph_k1 + bun_k1 + uo_k1, 
              data=mimic_si, family = "binomial")

# Get PS predictions
mimic_si$ps_hat <-  ps_mod$fitted.values

# Get predictions of compliance with the rule
mimic_si$ps_hat_d <-  mimic_si$ps_hat*mimic_si$r + (1-mimic_si$ps_hat)*(1-mimic_si$r) 

# Compute AIPWE Vale of the rule r_old
EY_s1 <- mean( with(mimic_si, c*d60d / ps_hat_d - pon_hat_d*(c-ps_hat_d)/ps_hat_d ) )

## Compute EY
EY <- with(mimic_si, mean(d60d))

it <- 0
res <- boot(mimic_si, are_aipw_new, R=500 )
boot.ci(res)

table(res$t>0)



