estimators <- function(dat, 
         expert_vars = c("X.2", "X.3", "X.4", "X.5", "X.6"),
         top_gate_vars = c("X.7"), 
         o_g_z_gate_vars = c("X.1", "X.2", "X.3", "X.4", "X.5"), 
         o_g_o_gate_vars = c("X.1", "X.2", "X.3", "X.4", "X.5"),
         epsilon=10^-2, maxit=40, seed=999, verbose=TRUE,
         
         r_var = "r", ttt_var = "A", out_y = "Y", 
         pr_var_except_ttt = c("X.2", "X.3", "X.4", "X.5", "X.6"),
         f_pr_mod = formula(Y ~ -1 + A*X.2 + A*X.3 + A*X.4 + A*X.5 + A*X.6),
         f_ps_mod = formula(A ~ 1 + X.1 + X.2 + X.3 + X.4 + X.5)
         
         ) {

res <- hme(dat, 
             expert_vars = expert_vars, 
             top_gate_vars = top_gate_vars, 
             o_g_z_gate_vars = o_g_z_gate_vars, 
             o_g_o_gate_vars = o_g_o_gate_vars,
             epsilon=epsilon, maxit=maxit, seed=seed, verbose=verbose)

EY_s0s <- as.numeric(tail(res, 1))

V_AIPW_r <- value_of_r(dat, r_var = r_var, ttt_var = ttt_var, out_y = out_y, 
           pr_var_except_ttt = pr_var_except_ttt,
           f_pr_mod = f_pr_mod,
           f_ps_mod = f_ps_mod,
           seed=seed) 

EY_s0_MOE <- EY_s0s[which.max(abs(V_AIPW_r-EY_s0s))]
EY_s1_MOE <- EY_s0s[which.min(abs(V_AIPW_r-EY_s0s))]

Y <- mean(dat[,out_y])

MIG <- V_AIPW_r - Y
ARE_MOE <- EY_s1_MOE - EY_s0_MOE
ARE_COMB <- V_AIPW_r - EY_s0_MOE
AIE <- Y - EY_s0_MOE

return (cbind(MIG=MIG, ARE_MOE=ARE_MOE, ARE_COMB=ARE_COMB, AIE=AIE) )
}

estimators_boot <- function(d, i=1:nrow(d),
                            expert_vars = c("X.2", "X.3", "X.4", "X.5", "X.6"),
                            top_gate_vars = c("X.7"), 
                            o_g_z_gate_vars = c("X.1", "X.2", "X.3", "X.4", "X.5"), 
                            o_g_o_gate_vars = c("X.1", "X.2", "X.3", "X.4", "X.5"),
                            epsilon=10^-2, maxit=40, seed=999, verbose=TRUE,
                            
                            r_var = "r", ttt_var = "A", out_y = "Y", 
                            pr_var_except_ttt = c("X.2", "X.3", "X.4", "X.5", "X.6"),
                            f_pr_mod = formula(Y ~ -1 + A*X.2 + A*X.3 + A*X.4 + A*X.5 + A*X.6),
                            f_ps_mod = formula(A ~ 1 + X.1 + X.2 + X.3 + X.4 + X.5)) {
  z<-d[i,]
  return(estimators(z,
                    expert_vars = expert_vars,
                    top_gate_vars = top_gate_vars, 
                    o_g_z_gate_vars = o_g_z_gate_vars, 
                    o_g_o_gate_vars = o_g_o_gate_vars,
                    epsilon=epsilon, maxit=maxit, seed=seed, verbose=verbose,
                    
                    r_var = r_var, ttt_var = ttt_var, out_y = out_y, 
                    pr_var_except_ttt = pr_var_except_ttt,
                    f_pr_mod = f_pr_mod,
                    f_ps_mod = f_ps_mod))
}

emp_boot_ci <-  function(t0, t, alpha=.05){
  # Empirical Bootstrap CI
  # as described in John Rice, Mathematical Statistics and Data Analysis, 3rd edition, p. 285.
  temp <- rbind(
    2*t0 - apply(t , 2, function(x) quantile(x, probs = 1-alpha/2) ),
    2*t0 - apply(t , 2, function(x) quantile(x, probs = alpha/2) )
  )
  
  return( temp ) 
}

t0 <- res$t0
t <- res$t


  
# Example
#res1 <- estimators(dat)
#res1
#TrueEstimands


#library(boot)
#res<-boot(dat, estimators_boot, R=100, verbose=FALSE)
#res
#plot(res, index=2)
#plot(res, index = 2)

