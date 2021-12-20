library(latex2exp)

expit <- function(x) 1 / (1+exp(-x))

###
hme <- function(dat, 
         expert_vars = c("X.2", "X.3", "X.4", "X.5", "X.6"), 
         top_gate_vars = c("X.7"), 
         o_g_z_gate_vars = c("X.1", "X.2", "X.3", "X.4", "X.5"), 
         o_g_o_gate_vars = c("X.1", "X.2", "X.3", "X.4", "X.5"),
         epsilon=10^-3, maxit=40, seed=123, verbose=TRUE) {

set.seed(seed)  
sd <- 100000000

Ys0_hat <- c()
Ys1_hat <- c()

reached_maxit <- FALSE
converged <- FALSE

while(reached_maxit==FALSE & converged==FALSE) {
try({
  
  z_g_z_expert_param <- rnorm(length(expert_vars), sd = sd) 
  o_g_z_expert_param <- rnorm(length(expert_vars), sd = sd) 
  
  z_g_o_expert_param <- rnorm(length(expert_vars), sd = sd) 
  o_g_o_expert_param <- rnorm(length(expert_vars), sd = sd) 
  
  top_gate_param <- rnorm(length(top_gate_vars), sd = sd)
  
  o_g_z_gate_param <- rnorm(length(o_g_z_gate_vars), sd = sd) 
  o_g_o_gate_param <- rnorm(length(o_g_o_gate_vars), sd = sd)
  ###
  
  for(i in 1:maxit) {
    ### h_i's
  
    g_o <- expit(as.matrix(dat[top_gate_vars]) %*% top_gate_param)
    g_z <- 1 - g_o
    
    g_o_g_o <- expit(as.matrix(dat[o_g_o_gate_vars]) %*% o_g_o_gate_param)
    g_z_g_o <- 1 - g_o_g_o
    
    g_o_g_z <- expit(as.matrix(dat[o_g_z_gate_vars]) %*% o_g_z_gate_param)
    g_z_g_z <- 1 - g_o_g_z
    
    p_o_o <- expit(as.matrix(dat[expert_vars]) %*% o_g_o_expert_param)
    p_o_z <- expit(as.matrix(dat[expert_vars]) %*% z_g_o_expert_param)
    
    p_z_o <- expit(as.matrix(dat[expert_vars]) %*% o_g_z_expert_param)
    p_z_z <- expit(as.matrix(dat[expert_vars]) %*% z_g_z_expert_param)
    
    numerator <- g_o * (g_o_g_o * p_o_o + g_z_g_o * p_o_z )
    denominator <- numerator + g_z * (g_o_g_z * p_z_o + g_z_g_z * p_z_z )
    dat$h_o <- numerator / denominator
    dat$h_z <- 1 - dat$h_o 
  
    ### h_j_g_i's
    numerator <- g_z_g_z * p_z_z
    denominator <- numerator + g_o_g_z * p_z_o
    dat$h_z_g_z <- numerator / denominator
    dat$h_o_g_z <-  1 - dat$h_z_g_z
    
    numerator <- g_z_g_o * p_o_z
    denominator <- numerator + g_o_g_o * p_o_o
    dat$h_z_g_o <- numerator / denominator
    dat$h_o_g_o <- 1 - dat$h_z_g_o
  
  ### M STEP
  
  # step 2
  f_experts <- formula(paste0("Y ~ -1 +", paste(expert_vars, collapse=" + ")))
  
  z_g_z_expert_mod <- glm(f_experts, weights = h_z * h_z_g_z, family = "quasibinomial", data=dat, method = glm.fit)
  o_g_z_expert_mod <- glm(f_experts, weights = h_z * h_o_g_z, family = "quasibinomial", data=dat, method = glm.fit)
  
  z_g_o_expert_mod <- glm(f_experts, weights = h_o * h_z_g_o, family = "quasibinomial", data=dat, method = glm.fit)
  o_g_o_expert_mod <- glm(f_experts, weights = h_o * h_o_g_o, family = "quasibinomial", data=dat, method = glm.fit)
  
  # step 3
  f_topgate <- formula(paste0("h_o ~ -1 +", paste(top_gate_vars, collapse=" + ")))
  top_gate_mod <- glm(f_topgate, family = "quasibinomial", data=dat, method = glm.fit) # Not fit for debugging
  
  # step 4
  f_o_g_z_gate <- formula(paste0("h_o_g_z ~ -1 +", paste(o_g_z_gate_vars, collapse=" + ")))
  o_g_z_gate_mod <- glm(f_o_g_z_gate, weights = h_z, family = "quasibinomial", data=dat, method = glm.fit)
  
  f_o_g_o_gate <- formula(paste0("h_o_g_o ~ -1 +", paste(o_g_o_gate_vars, collapse=" + ")))  # only if this need be updated
  o_g_o_gate_mod <- glm(f_o_g_o_gate, family = "quasibinomial", data=dat, method = glm.fit)  # only if this need be updated
  
  
  # step 5 update parameters
  
  #new_thetas <- c(coef(z_g_z_expert_mod), coef(o_g_z_expert_mod), coef(z_g_o_expert_mod), coef(o_g_o_expert_mod), coef(o_g_z_gate_mod), coef(top_gate_mod), coef(o_g_o_gate_mod))
  #old_thetas <- c(z_g_z_expert_param, o_g_z_expert_param, z_g_o_expert_param, o_g_o_expert_param, o_g_z_gate_param, top_gate_param, o_g_o_gate_param)
  
  Ys1_hat <- c(Ys1_hat, mean(g_o_g_o * p_o_o + g_z_g_o * p_o_z))
  Ys0_hat <- c(Ys0_hat, mean(g_o_g_z * p_z_o + g_z_g_z * p_z_z))
  
  conv <- norm(diff(tail(cbind(Ys1_hat, Ys0_hat),2)), type = "2")
  if (conv > epsilon) {
  
  z_g_z_expert_param <- coef(z_g_z_expert_mod)
  o_g_z_expert_param <- coef(o_g_z_expert_mod)
  
  z_g_o_expert_param <- coef(z_g_o_expert_mod)
  o_g_o_expert_param <- coef(o_g_o_expert_mod)
  
  top_gate_param <- coef(top_gate_mod) # No update for debugging
  
  o_g_z_gate_param <- coef(o_g_z_gate_mod)
  o_g_o_gate_param <- coef(o_g_o_gate_mod) # only if this need be updated
  } else { converged <- TRUE 
          break }
  
  ### Debugging
  if(verbose==TRUE) { print(paste(i, format(round(conv,2), nsmall = 2) )) }
  if(i==maxit) { reached_maxit <- TRUE  }
  }

}, silent=TRUE)

}
if(reached_maxit) {print("Reached maximum number of iterations")}
if(converged) {print("Algorithm converged")}
return(data.frame(Yk0_hat = Ys0_hat, Yk1_hat = Ys1_hat))
}


plot.hme <- function(res){
    par(mar = c(5, 5, 2, 2))
    matplot(as.matrix(res), type = c("b"), pch=16, cex=.5, col = c("#DF8F44FF", "#00A1D5FF"), las=1, bty="n",
            ylim=c(0,1), xlim=c(0, nrow(res)), 
            xlab = "EM iteration", ylab="Estimates", 
            main="HME convergence")
    legend("topright", c(TeX(r'($\widehat{E}\[Y^{k=1}\]$)'), TeX(r'($\widehat{E}\[Y^{k=0}\]$)')),  pch = 16, col = c("#00A1D5FF", "#DF8F44FF"), bty = "n")
}

## Example
#res <- hme(superpop[10000:20000,], 
#           expert_vars = c("X.2", "X.3", "X.4", "X.5", "X.6"), 
#           top_gate_vars = c("X.7"), 
#           o_g_z_gate_vars = c("X.1", "X.2", "X.3", "X.4", "X.5"), 
#           o_g_o_gate_vars = c("X.1", "X.2", "X.3", "X.4", "X.5"),
#           epsilon=10^-3, maxit=40, seed=66, verbose=TRUE)

#plot.hme(res)
#abline(h=TrueVals["Y_s0"], lty = 3)
#abline(h=TrueVals["Y_s1"], lty = 3)

#data.frame(Hats=as.numeric(tail(res,1)), True=c(TrueVals["Y_s0"], TrueVals["Y_s1"]) )


