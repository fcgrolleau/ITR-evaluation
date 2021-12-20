value_of_r <- function(dat, r_var, ttt_var, out_y, 
                 pr_var_except_ttt,
                 f_pr_mod,
                 f_ps_mod ) {

  
## Value of the rule APIWE E_Y^{s=1}
##  ---- Inputs ----
## dat: data.frame
## r_var: string for rule variable
## ttt_var: string for treatment variable
## out_y: string for outcome variable Y
## pr_var_except_ttt: concatenation of strings: variables pronostic model not including treatment
## f_pr_mod: formula for the pronostic model always use -1 and add column name filled with ones if intercept is needed 
## f_ps_mod: formula for the propensity score model always use -1 and add column name filled with ones if intercept is needed
## 
  
  # Fit pronositc model
pr_mod <- glm(f_pr_mod, data=dat, family = "binomial")

dat$pon_hat <-  pr_mod$fitted.values

# Get predictions under the recommended treatment
newdata <- dat[, c(r_var, pr_var_except_ttt )]
colnames(newdata)[1] <- ttt_var
dat$pon_hat_d <- predict(pr_mod, newdata, type="response")

# Fit PS model
ps_mod <- glm(f_ps_mod, 
              data=dat, family = "binomial")

# Get PS predictions
dat$ps_hat <-  ps_mod$fitted.values

# Get predictions of compliance with the rule
dat$ps_hat_d <-  dat$ps_hat * dat[, r_var] + (1-dat$ps_hat)*(1-dat[, r_var]) 

# Compute AIPWE Vale of the rule r_old
dat$C <- as.numeric(dat[, r_var] == dat[, ttt_var]) # Get Concordance variable first
dat$Y <- dat[, out_y]
EY_s1 <- mean( with(dat, C*Y / ps_hat_d - pon_hat_d*(C-ps_hat_d)/ps_hat_d ) )

return(EY_s1)
}

# Example

# value_of_r(dat[10000:20000,], r_var = "r", ttt_var = "A", out_y = "Y", 
#           pr_var_except_ttt = c("X.2", "X.3", "X.4", "X.5", "X.6"),
#           f_pr_mod = formula(Y ~ -1 + A*X.2 + A*X.3 + A*X.4 + A*X.5 + A*X.6),
#           f_ps_mod = formula(A ~ 1 + X.1 + X.2 + X.3 + X.4 + X.5)) 
