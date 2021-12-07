library(boot)
library(gplots)
library(latex2exp)

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

Deltas_cb <- c()
Alpha_seq_cb <- seq(0, 1-10^-100, length=10)
mean_imp_cb <- c()
resamples <- 100
boot_ci_cb <- list()
it <- 1
it_boot <- 0
resamples <- 100

# cognitive bias scenario
for (i in Alpha_seq_cb ) {
  alpha_cb <- i
  mimic_si$rho <- with(mimic_si, (1- abs(r - ps_hat) )^(.5*log((alpha_cb+1)/(1-alpha_cb))) )
  Deltas_cb <- c(Deltas_cb, mean( with(mimic_si, rho*(r - ps_hat ) * ITE_hat)) )
  res <- boot(mimic_si, cb_new, R=resamples )
  it_boot <- 0
  boot_ci_cb[[it]] <- quantile(res$t, probs = c(.025, .975))
  mean_imp_cb <- c(mean_imp_cb, mean(mimic_si$rho))
  it <- it+1
}
boot_ci_cb <- sapply(boot_ci_cb, c)

# confidence level scenario
Alpha_seq_cl <- c(seq(0, .93, length.out = 3), seq(.94, 1, by=.01) ) 
Deltas_cl <- c()
mean_imp_cl <- c()
boot_ci_cl <- list()
it <- 1
it_boot <- 0
#resamples <- 50

for (i in Alpha_seq_cl ) {
  alpha_cl <- i
  mimic_si$rho <- with(mimic_si, (iard - qnorm(1 - alpha_cl/2))*(iard + qnorm(1 - alpha_cl/2)) > 0 )
  Deltas_cl <- c(Deltas_cl, mean( with(mimic_si, rho*(r - ps_hat ) * ITE_hat)) )
  res <- boot(mimic_si, cl_new, R=resamples )
  it_boot <- 0
  boot_ci_cl[[it]] <- quantile(res$t, probs = c(.025, .975))
  mean_imp_cl <- c(mean_imp_cl, mean(mimic_si$rho))
  it <- it+1
}     
boot_ci_cl <- sapply(boot_ci_cl, c)

ci_width <- 0.004
diamond_size <- 2
lab_size <- 2
xlab_pos <- 3.7
ylab_pos <- 3.4
wl <- 8
dev.new(width=wl, height=wl, pointsize=7, noRStudioGD = TRUE)
par(mfcol=c(2,2), mar = c(5.5, 5.9, 2, 2))
plot(NULL, xlim=c(0,1), ylim=round(c(min(c(boot_ci_cl, boot_ci_cb)), max(c(boot_ci_cl, boot_ci_cb))), 3), bty="n", las=1,
  xlab="", ylab="")
title(xlab=TeX('$\\alpha$'), line=xlab_pos, cex.lab=lab_size)
title(ylab=TeX('$\\widehat{\\Delta}_{ITE}(r,\\rho_{cb,\\alpha})$'), line=ylab_pos, cex.lab=lab_size)
abline(h = 0, lty=1, lwd = 2)
plotCI(Alpha_seq_cb, Deltas_cb, ui=boot_ci_cb[2,], li=boot_ci_cb[1,], pch=18, gap=0, cex=diamond_size, sfrac=ci_width, col="#00a1d5ff", barcol="black", add=TRUE)
plot(NULL, xlim=c(0,1), ylim=round(c(min(c(boot_ci_cl, boot_ci_cb)), max(c(boot_ci_cl, boot_ci_cb))), 3), bty="n", las=1,
     xlab="", ylab="")
title(xlab=TeX('$\\alpha$'), line=xlab_pos, cex.lab=lab_size)
title(ylab=TeX('$\\widehat{\\Delta}_{ITE}(r,\\rho_{cl,\\alpha})$'), line=ylab_pos, cex.lab=lab_size)
abline(h = 0, lty=1, lwd = 2)
plotCI(Alpha_seq_cl, Deltas_cl, ui=boot_ci_cl[2,], li=boot_ci_cl[1,], pch=18, gap=0, cex=diamond_size, sfrac=ci_width, col="#df8f44ff", barcol="black", add=TRUE)


plot(NULL, xlim=c(0,1), ylim=round(c(min(c(boot_ci_cl, boot_ci_cb)), max(c(boot_ci_cl, boot_ci_cb))), 3), bty="n", las=1,
     xlab="", ylab="")
title(xlab=TeX(r'($\widehat{\mathbf{E}}\{\rho_{cb,\alpha}(X)\}$)'), line=xlab_pos, cex.lab=lab_size)
title(ylab=TeX('$\\widehat{\\Delta}_{ITE}(r,\\rho_{cb,\\alpha})$'), line=ylab_pos, cex.lab=lab_size)
abline(h = 0, lty=1, lwd = 2)
plotCI(mean_imp_cb, Deltas_cb, ui=boot_ci_cb[2,], li=boot_ci_cb[1,], pch=18, gap=0, cex=diamond_size, sfrac=ci_width, col="#00a1d5ff", barcol="black", add=TRUE)
plot(NULL, xlim=c(0,1), ylim=round(c(min(c(boot_ci_cl, boot_ci_cb)), max(c(boot_ci_cl, boot_ci_cb))), 3), bty="n", las=1,
     xlab="", ylab="")
title(xlab=TeX(r'($\widehat{\mathbf{E}}\{\rho_{cl,\alpha}(X)\}$)'), line=xlab_pos, cex.lab=lab_size)
title(ylab=TeX('$\\widehat{\\Delta}_{ITE}(r,\\rho_{cl,\\alpha})$'), line=ylab_pos, cex.lab=lab_size)
abline(h = 0, lty=1, lwd = 2)
plotCI(mean_imp_cl, Deltas_cl, ui=boot_ci_cl[2,], li=boot_ci_cl[1,], pch=18, gap=0, cex=diamond_size, sfrac=ci_width, col="#df8f44ff", barcol="black", add=TRUE)



