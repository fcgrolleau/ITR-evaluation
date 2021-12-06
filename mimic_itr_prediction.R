mimic_si <- read.csv("~/Desktop/github repos/Emulated-ITR/Python/data/mimic_si.csv")

library(survival)
library(splines)

load(url("https://github.com/fcgrolleau/Functions-and-explanatory-analysis/raw/master/testenvironment.RData"))

ite_se <- function(pot, sofa, weight, is, bun, ph) {
  
  expit <- function(x){1/(1+exp(-x))}
  coefs<-finalmodel$estimate

  ### compute predicted probablity of RRT initiation
x_vect<-as.numeric(c(1, pot, sofa, weight, is, bun, ph))
pred<- expit(x_vect %*% coefs )
  #trunc_pred<-ifelse(pred>.55, .55, pred)
  if(pred>.55) {trunc_pred<-.55
    } else { if(pred<0.02468) {trunc_pred<-0.02468
          } else{trunc_pred<-pred} }
  
  
  ### compute other indivudal event rates and absiolute risk difference
  
  ## ER plots  
  dr_early.survfitobj<-survfit(fit.int.er, newdata=data.frame(bras="STRATEGIE PRECOCE", probs.c=trunc_pred-.25, probs.int.t=0))
  dr_late.survfitobj<-survfit(fit.int.er, newdata=data.frame(bras="STRATEGIE D ATTENTE", probs.c=trunc_pred-.25, probs.int.t=trunc_pred-.25))
  
  dr_early_ci <- c(1-tail(dr_early.survfitobj$surv,1),
                   1-tail(dr_early.survfitobj$upper,1),
                   1-tail(dr_early.survfitobj$lower,1))
  
  dr_late_ci<-c(1-tail(dr_late.survfitobj$surv,1), 
                1-tail(dr_late.survfitobj$upper,1),
                1-tail(dr_late.survfitobj$lower,1))
  
  ## ARD plots
  #closest position
  closest<-which.min(abs(rep(trunc_pred,106)-seq(0.02468, .55, by=.005)))
  
  
  pred.ard<-c(dr_early_ci[1]-dr_late_ci[1],
              dr_early_ci[1]-dr_late_ci[1]-qnorm(.975)*se_s_pooled[closest],
              dr_early_ci[1]-dr_late_ci[1]+qnorm(.975)*se_s_pooled[closest])

i_preds <- c(trunc_pred, dr_early_ci[1]-dr_late_ci[1], se_s_pooled[closest])

#print(c(pot, sofa, weight, is, bun, ph, i_preds)) # for debugging
return (i_preds)
}

# Compute predictions
external_preds <- apply(mimic_si, 1, function(x)  ite_se(pot = as.numeric(x["pot_k1"]),
                                             sofa = as.numeric(x["SOFA_24hours"]),
                                             weight = as.numeric(x["weight"]),
                                             is = as.logical(x["immunosuppressant"]),
                                             bun = as.numeric(x["bun_k1"]),
                                             ph = as.numeric(x["ph_k1"]) ) ) 

external_preds <- as.data.frame(t(external_preds))
colnames(external_preds) <-c("pred", "iard", "ise")

# Store predictions in a new data frame
mimic_si_preds <- cbind(mimic_si, external_preds )

#Compute recomendations from predictions
mimic_si_preds$r <- ifelse(mimic_si_preds$iard<0, 1, 0)

#Compute 60-day mortality
mimic_si_preds$d60d <- with(mimic_si_preds, difftime(dod, charttime.kdigo3, units = "days") )<=60 
mimic_si_preds$d60d[is.na(mimic_si_preds$d60d)] <- 0

# Export
#write.csv(x = mimic_si_preds, "mimic_si_preds.csv", row.names=FALSE)

