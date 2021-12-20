res <- hme(superpop[1000:2000,], 
             expert_vars = c("X.2", "X.3", "X.4", "X.5", "X.6"), 
             top_gate_vars = c("X.7"), 
             o_g_z_gate_vars = c("X.1", "X.2", "X.3", "X.4", "X.5"), 
             o_g_o_gate_vars = c("X.1", "X.2", "X.3", "X.4", "X.5"),
             epsilon=10^-2, maxit=40, seed=999, verbose=FALSE)
EY_s0s <- tail(res, 1)

V_AIPW_r <- value_of_r(superpop[1000:2000,], r_var = "r", ttt_var = "A", out_y = "Y", 
           pr_var_except_ttt = c("X.2", "X.3", "X.4", "X.5", "X.6"),
           f_pr_mod = formula(Y ~ -1 + A*X.2 + A*X.3 + A*X.4 + A*X.5 + A*X.6),
           f_ps_mod = formula(A ~ 1 + X.1 + X.2 + X.3 + X.4 + X.5)) 

EY_s0 <- EY_s0s[which.max(abs(V_AIPW_r-EY_s0s))]
