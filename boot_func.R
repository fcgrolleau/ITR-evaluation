emp_boot_ci <-  function(t0, t, alpha=.05){
  # Empirical Bootstrap CI
  # as described in John Rice, Mathematical Statistics and Data Analysis, 3rd edition, p. 285.        
  return( as.numeric( c(2*t0 - quantile(t, probs = 1-alpha/2), 2*t0 - quantile(t, probs = alpha/2)) ) )
}

are_aipw_new <- function(d, i=1:nrow(d)) {
  z<-d[i,]
  
  it <<- it + 1  
  ## Value of the rule APIWE E_Y^{s=1}
  
  # Fit pronositc model
  pr_mod <- glm(d60d ~ a1*admission_age + a1*SOFA_24hours + weight + bun_k1 + a1*ph_k1 + a1*poly(pot_k1,3), 
                data=z, family = "binomial")
  
  z$pon_hat <-  pr_mod$fitted.values
  
  # Get predictions under the recommended treatment
  newdata <- z[, c("r", "admission_age", "weight", "SOFA_24hours", "immunosuppressant", "bun_k1", "ph_k1", "pot_k1" )]
  colnames(newdata)[1] <- "a1"
  z$pon_hat_d <- predict(pr_mod, newdata, type="response")
  
  # Fit PS model
  ps_mod <- glm(a1 ~ admission_age + weight + bun_k1 + ph_k1 + pot_k1 + SOFA_24hours + immunosuppressant, 
                data=z, family = "binomial")
  
  # Get PS predictions
  z$ps_hat <-  ps_mod$fitted.values
  
  # Get predictions of compliance with the rule
  z$ps_hat_d <-  z$ps_hat*z$r + (1-z$ps_hat)*(1-z$r) 
  
  # Compute AIPWE Vale of the rule r_old
  EY_s1 <- mean( with(z, c*d60d / ps_hat_d - pon_hat_d*(c-ps_hat_d)/ps_hat_d ) )
  
  ## Compute EY
  EY <- with(z, mean(d60d))

  if ( (100*it/resamples) %% 10 == 0) {
    print(paste0('Iteration ', it, ': ', 100*it/resamples, '%')) }
  
  return(EY_s1 - EY)
}


cb_new <- function(d, i=1:nrow(d)) {
  z<-d[i,]
  
  it_boot <<- it_boot + 1 
  
  Deltas <- c()
  Alpha_seq <- seq(0,.999, by=.1) 
  mean_imp <- c()
  
  
  # cognitive bias scenario

  z$rho <- with(z, (1- abs(r - ps_hat) )^(.5*log((alpha_cb+1)/(1-alpha_cb))) )
  Delta <-  mean( with(z, rho*(r - ps_hat ) * ITE_hat))
  
  if ( (100*it_boot/resamples) %% 10 == 0) {
    print(paste0('Iteration alpha ', it, ' Iteration bootstrap ', it_boot, ': ', 100*it_boot/resamples, '%')) }
  return (Delta)
}

cl_new <- function(d, i=1:nrow(d)) {
  z<-d[i,]
  
  it_boot <<- it_boot + 1 
  
  Deltas <- c()
  Alpha_seq <- seq(0,.999, by=.1) 
  mean_imp <- c()
  
  # confidence level scenario
  
  z$rho <- with(z, (iard - qnorm(1 - alpha_cl/2))*(iard + qnorm(1 - alpha_cl/2)) > 0 )
  Delta <-  mean( with(z, rho*(r - ps_hat ) * ITE_hat)) 
  
  if ( (100*it_boot/resamples) %% 10 == 0) {
    print(paste0('Iteration alpha ', it, ' Iteration bootstrap ', it_boot, ': ', 100*it_boot/resamples, '%')) }
  return (Delta)
} 

##
ipl_itr_boot <- function(d, i=1:nrow(d)) {
  z<-d[i,]
  
  it <<- it + 1 

  #create a variable r_old for an old impolemented rule (SOFA>10)
  z$r_old <- as.numeric(z$SOFA_24hours >11 )
  
  #create a variable c for concordance between r and the delivered treatment 
  z$c <- z$r_old == z$a1
  
  ## Value of the rule APIWE E_Y^{s=1}
  
  # Fit pronositc model
  pr_mod <- glm(d60d ~ a1*admission_age + a1*SOFA_24hours + weight + bun_k1 + a1*ph_k1 + a1*poly(pot_k1,3), 
                data=z, family = "binomial")
  
  z$pon_hat <-  pr_mod$fitted.values
  
  # Get predictions under the recomended treatment
  newdata <- z[, c("r_old", "admission_age", "weight", "SOFA_24hours", "immunosuppressant", "uo_k1", "bun_k1", "ph_k1", "pot_k1" )]
  colnames(newdata)[1] <- "a1"
  z$pon_hat_d <- predict(pr_mod, newdata, type="response")
  
  # Fit PS model
  ps_mod <- glm(a1 ~ admission_age + weight + bun_k1 + ph_k1 + pot_k1 + SOFA_24hours + immunosuppressant, 
                data=z, family = "binomial")
  
  # Get PS predictions
  z$ps_hat <-  ps_mod$fitted.values
  
  # Get predictions of compliance with the rule
  z$ps_hat_d <-  z$ps_hat*z$r_old + (1-z$ps_hat)*(1-z$r_old) 
  
  # Compute AIPWE Vale of the rule r_old
  EY_s1 <- mean( with(z, c*d60d / ps_hat_d - pon_hat_d*(c-ps_hat_d)/ps_hat_d ) )
  
  
  ## E_Y^{s=0}
  Model <- FLXMRglmfix(fixed = ~ 1, formula =  ~ -1 + admission_age +  SOFA_24hours + weight + bun_k1 + ph_k1 + pot_k1, family = "binomial")
  concomitantModel <- FLXPmultinom(~ ~ 1 + admission_age + bun_k1 + ph_k1 + pot_k1)
  
  Moe <- stepFlexmix(cbind(d60d, 1 - d60d) ~ 1,
                     k = 2, model = Model,
                     concomitant = concomitantModel, data = z, nrep = 1)

  exp_coef <- parameters(Moe, which="model")
  gate_coef <- parameters(Moe, which="concomitant")
  exp_coef
  gate_coef
  
  
  temp <- model.matrix(~ 1 + admission_age +  SOFA_24hours + weight + bun_k1 + ph_k1 + pot_k1, z)
  
  expit <- function(x) 1 / (1+exp(-x))
  
  EY_s0s <- c(mean(expit(temp %*% exp_coef[,1])), mean( expit(temp %*% exp_coef[,2])) )
  
  EY_s0 <- EY_s0s[which.max(abs(EY_s1-EY_s0s))]
  EY_s1_moe <- EY_s0s[which.min(abs(EY_s1-EY_s0s))]
  
  EY <- mean(z$d60d)
  
  # MIG
  MIG <- EY_s1 - EY
  
  # ARE
  ARE_moe <- EY_s1  - EY_s0
  ARE_comb <- EY_s1_moe  - EY_s0
  
  # AIE
  AIE <- EY  - EY_s0

if ( (100*it/resamples) %% 10 == 0) {
  print(paste0('Iteration bootstrap ', it, ': ', 100*it/resamples, '%')) }
return ( c(ARE_moe, ARE_comb, AIE, MIG) )
}

