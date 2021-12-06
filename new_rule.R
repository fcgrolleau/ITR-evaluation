library(boot)
library(gplots)

mimic_si <-read.csv("/Users/francois/Desktop/github repos/ITR-evaluation/ITR-evaluation/mimic_si_preds.csv")

# Create concrdance variable hfd
mimic_si$c <- with(mimic_si, a1 == r)

## Value of the rule APIWE E_Y^{s=1}

# Fit pronositc model
pr_mod <- glm(d60d ~ a1*admission_age + a1*SOFA_24hours + weight + bun_k1 + a1*ph_k1 + a1*poly(pot_k1,3), 
              data=mimic_si, family = "binomial")


# Get predictions under the recommended treatment
newdata <- mimic_si[, c("r", "admission_age", "weight", "SOFA_24hours", "immunosuppressant", "bun_k1", "ph_k1", "pot_k1" )]
colnames(newdata)[1] <- "a1"
mimic_si$pon_hat_d <- predict(pr_mod, newdata, type="response")


# Get predictions for the ITE
newdata <- mimic_si[, c("a1", "admission_age", "weight", "SOFA_24hours", "immunosuppressant", "bun_k1", "ph_k1", "pot_k1" )]
newdata[,"a1"] <- 1
mimic_si$E_hat_Ya1_given_x <- predict(pr_mod, newdata, type="response")

newdata[,"a1"] <- 0
mimic_si$E_hat_Ya0_given_x <- predict(pr_mod, newdata, type="response")

mimic_si$ITE_hat <- with(mimic_si, E_hat_Ya1_given_x - E_hat_Ya0_given_x)

# Fit PS model
ps_mod <- glm(a1 ~ admission_age + weight + bun_k1 + ph_k1 + pot_k1 + SOFA_24hours + immunosuppressant, 
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
#res <- boot(mimic_si, are_aipw_new, R=resamples )
res
plot(res)
quantile(res$t, probs = c(.025, .975))

### Delta ITE 
mean( with(mimic_si, (r - ps_hat ) * ITE_hat))

Deltas <- c()
Alpha_seq <- abs( seq(0,1, by=.1)-10^(-10) ) 
mean_imp <- c()
resamples <- 100
boot_ci <- list()
it <- 1
it_boot <- 0
resamples <- 100

# cognitive bias scenario
for (i in Alpha_seq ) {
  alpha_cb <- i
  mimic_si$rho <- with(mimic_si, (1- abs(r - ps_hat) )^(.5*log((alpha_cb+1)/(1-alpha_cb))) )
  Deltas <- c(Deltas, mean( with(mimic_si, rho*(r - ps_hat ) * ITE_hat)) )
  res <- boot(mimic_si, cb_new, R=resamples )
  it_boot <- 0
  boot_ci[[it]] <- quantile(res$t, probs = c(.025, .975))
  mean_imp <- c(mean_imp, mean(mimic_si$rho))
  it <- it+1
}
boot_ci <- sapply(boot_ci, c)

par(mfrow=c(1,2))
plot(NULL, xlim=c(0,1), ylim=c(min(boot_ci[1,]), max(boot_ci[2,])), bty="n", las=1, xlab="Cognitive Biais Parameter Alpha", ylab="Deltas")
abline(h = 0, lty=1)
plotCI(Alpha_seq, Deltas, ui=boot_ci[2,], li=boot_ci[1,], pch=18, gap=0, sfrac=0.003, col="#00a1d5", barcol="black", xlab="", ylab="Deltas", bty="n", add=TRUE)

plot(NULL, xlim=c(0,1), ylim=c(min(boot_ci[1,]), max(boot_ci[2,])), bty="n", yaxt="n", las=1, xlab="Proportion of Patients Implementing The Rule", ylab="")
abline(h = 0, lty=1)
plotCI(mean_imp, Deltas, ui=boot_ci[2,], li=boot_ci[1,], pch=18, gap=0, sfrac=0.003, col="#00a1d5", barcol="black", xlab="", ylab="",
       xlim=c(0, 1), bty="n", yaxt="n", add=TRUE)

# confidence level scenario
Alpha_seq <- seq(0,1, by=.1) 
Deltas <- c()
mean_imp <- c()
boot_ci <- list()
it <- 1
it_boot <- 0
resamples <- 50

for (i in Alpha_seq ) {
  alpha_cl <- i
  mimic_si$rho <- with(mimic_si, (iard - qnorm(1 - alpha_cl/2))*(iard + qnorm(1 - alpha_cl/2)) > 0 )
  Deltas <- c(Deltas, mean( with(mimic_si, rho*(r - ps_hat ) * ITE_hat)) )
  res <- boot(mimic_si, cl_new, R=resamples )
  it_boot <- 0
  boot_ci[[it]] <- quantile(res$t, probs = c(.025, .975))
  mean_imp <- c(mean_imp, mean(mimic_si$rho))
  it <- it+1
}     
boot_ci <- sapply(boot_ci, c)

par(mfrow=c(1,2))
plot(Alpha_seq, Deltas)
plot(mean_imp, Deltas)

temp_df <- data.frame(Alpha_cl=sequence, Delta=Deltas)
temp_df[temp_df$Delta<0, ][1,]
