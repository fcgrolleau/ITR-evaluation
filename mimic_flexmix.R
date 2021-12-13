library(flexmix)
library(boot)
library(latex2exp)

set.seed(4321)

mimic_si <- read.csv("~/Desktop/github repos/Emulated-ITR/Python/data/mimic_si_preds.csv")

#create a variable r_old for an old implemented rule (SOFA>10)
mimic_si$r_old <- as.numeric(mimic_si$SOFA_24hours >11 )

#create a variable c for concordance between r and the delivered treatment 
mimic_si$c <- mimic_si$r_old == mimic_si$a1

## Value of the rule APIWE E_Y^{s=1}

# Fit pronositc model
pr_mod <- glm(d60d ~ a1*admission_age + a1*SOFA_24hours + weight + bun_k1 + a1*ph_k1 + a1*poly(pot_k1,3), 
              data=mimic_si, family = "binomial")

mimic_si$pon_hat <-  pr_mod$fitted.values

# Get predictions under the recomended treatment
newdata <- mimic_si[, c("r_old", "admission_age", "weight", "SOFA_24hours", "immunosuppressant", "uo_k1", "bun_k1", "ph_k1", "pot_k1" )]
colnames(newdata)[1] <- "a1"
mimic_si$pon_hat_d <- predict(pr_mod, newdata, type="response")

# Fit PS model
ps_mod <- glm(a1 ~ admission_age + weight + bun_k1 + ph_k1 + pot_k1 + SOFA_24hours + immunosuppressant, 
              data=mimic_si, family = "binomial")

# Get PS predictions
mimic_si$ps_hat <-  ps_mod$fitted.values

# Get predictions of compliance with the rule
mimic_si$ps_hat_d <-  mimic_si$ps_hat*mimic_si$r_old + (1-mimic_si$ps_hat)*(1-mimic_si$r_old) 

# Compute AIPWE Vale of the rule r_old
EY_s1 <- mean( with(mimic_si, c*d60d / ps_hat_d - pon_hat_d*(c-ps_hat_d)/ps_hat_d ) )


## E_Y^{s=0}
Model <- FLXMRglmfix(fixed = ~ 1, formula =  ~ -1 + admission_age +  SOFA_24hours + weight + bun_k1 + ph_k1 + pot_k1, family = "binomial")
concomitantModel <- FLXPmultinom(~ ~ 1 + admission_age + bun_k1 + ph_k1 + pot_k1 + c)

Moe <- stepFlexmix(cbind(d60d, 1 - d60d) ~ 1,
                   k = 2, model = Model,
                   concomitant = concomitantModel, data = mimic_si, nrep = 3)

refit_Moe <- refit(Moe)
summary(refit_Moe)

exp_coef <- parameters(Moe, which="model")
gate_coef <- parameters(Moe, which="concomitant")
exp_coef
gate_coef


temp <- model.matrix(~ 1 + admission_age +  SOFA_24hours + weight + bun_k1 + ph_k1 + pot_k1, mimic_si)

expit <- function(x) 1 / (1+exp(-x))
  
EY_s0s <- c(mean(expit(temp %*% exp_coef[,1])), mean( expit(temp %*% exp_coef[,2])) )

EY_s0 <- EY_s0s[which.max(abs(EY_s1-EY_s0s))]
EY_s1_moe <- EY_s0s[which.min(abs(EY_s1-EY_s0s))]

EY <- mean(mimic_si$d60d)

# MIG
MIG <- EY_s1 - EY

# ARE
ARE_moe <- EY_s1  - EY_s0
ARE_comb <- EY_s1_moe  - EY_s0

# AIE
AIE <- EY  - EY_s0

# Recall the comparison
EY_s0s
EY_s0
EY_s1
EY

resamples <- 1000; it <- 0
res <- boot(mimic_si, ipl_itr_boot, R=resamples)

estimates <- c(ARE_moe, ARE_comb, AIE, MIG)
est_prep <- lapply(1:length(estimates), function(x) list(estimates[x], res$t[,x]))
estimates_ci <- sapply(est_prep, function(x) emp_boot_ci(x[[1]], x[[2]]))

##### Make a nice plot
dev.new(width=6,height=6,pointsize=9)
par(mar=c(3,10,2,2)+.1, mgp=c(1.75,0.5,0), tcl=-.4)
xl <- c(-.2, .5)
yl <- c(0.5,4.5)
#
plot(xl,yl,axes=F,type="n",xlab="",ylab="", xaxs="i")
segments(0,par('usr')[3],0,11)
segments(estimates_ci[1,],length(estimates):1,estimates_ci[2,],length(estimates):1)
points(estimates,length(estimates):1, pch=15, col="#00A1D5FF", cex=2)
axis(1, at=c(seq(from=xl[1], to=xl[2], by=.1),0), labels=c(seq(from=xl[1], to=xl[2], by=.1),0), cex.axis=.8, lwd=0.75, lwd.tick=0.75)
mtext("Absolute Risk Difference", side=1, line=1.75, cex=1)
#
arrows(c(-.01,.01), 4.5, 25*c(-.01,.01), 4.5, length = 0.07, xpd=TRUE)
mtext("Favors ITR implementation",side=3, line=-.6, at=-.01, font=1, las=1, adj=1, cex=.8)
mtext("Favors no ITR implementation",side=3, line=-.6, at=.01, font=1, las=1, adj=0, cex=.8)
mtext("Estimates",side=1, line=0, at=-.45, font=2, las=1, adj=0)
mtext(c(TeX(r'($\widehat{\Delta}_{MOE}(r)$)'),
            TeX(r'($\widehat{\Delta}_{COMB}(r)$)'),
                TeX(r'($\widehat{\Lambda}(r, \rho)$)'),
                    TeX(r'($\widehat{\Gamma}(r, \rho)$)')
                         ) ,side=2, line=7, at=length(estimates):1, las=1, adj=0)

save.image("implemented_mimic_anlysis.RData")

