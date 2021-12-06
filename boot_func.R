are_aipw_new <- function(d, i=1:nrow(d)) {
  z<-d[i,]
  
  it <<- it + 1  
  ## Value of the rule APIWE E_Y^{s=1}
  
  # Fit pronositc model
  pr_mod <- glm(d60d ~ admission_age + a1*weight + a1*bun_k1 + a1*ph_k1 + a1*pot_k1 + a1*SOFA_24hours + a1*immunosuppressant, 
                data=z, family = "binomial")
  
  z$pon_hat <-  pr_mod$fitted.values
  
  # Get predictions under the recommended treatment
  newdata <- z[, c("r", "admission_age", "weight", "SOFA_24hours", "immunosuppressant", "bun_k1", "ph_k1", "pot_k1" )]
  colnames(newdata)[1] <- "a1"
  z$pon_hat_d <- predict(pr_mod, newdata, type="response")
  
  # Fit PS model
  ps_mod <- glm(a1 ~ ph_k1 + bun_k1 + pot_k1, 
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

