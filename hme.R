Ys0_hat <- c()
Ys1_hat <- c()

dat$r_sign <- with(dat, ifelse(r==0, -1, r))
dat$S_sign <- with(dat, ifelse(S==0, -1, S))

expit <- function(x) 1 / (1+exp(-x))

###
expert_vars <- c("X.2", "X.3", "X.4", "X.5", "X.6")
o_g_z_gate_vars <- c("X.1", "X.2", "X.3", "X.4", "X.5")
o_g_o_gate_vars <- c("X.1", "X.2", "X.3", "X.4", "X.5")
top_gate_vars <- c("X.7")

sd <- 1000000
  
z_g_z_expert_param <- rnorm(length(expert_vars), sd = sd) # c(-3, -.5, 5, -1.5, -2)  + rnorm(length(expert_vars), sd = sd)
o_g_z_expert_param <- rnorm(length(expert_vars), sd = sd) # c(-3, -.5, 5, -1.5, -2) *-1 + rnorm(length(expert_vars), sd = sd)

z_g_o_expert_param <- rnorm(length(expert_vars), sd = sd) # c(-3, -.5, 5, -1.5, -2) + rnorm(length(expert_vars), sd = sd)
o_g_o_expert_param <- rnorm(length(expert_vars), sd = sd) # c(-3, -.5, 5, -1.5, -2) *-1 + rnorm(length(expert_vars), sd = sd)

top_gate_param <- rnorm(length(top_gate_vars), sd = sd)

o_g_z_gate_param <- rnorm(length(o_g_z_gate_vars), sd = sd) #c(.05, -.5, .5, -.5, .5) + rnorm(length(o_g_z_gate_vars), sd = sd)
o_g_o_gate_param <- rnorm(length(o_g_o_gate_vars), sd = sd)
###

epsilon <- .05
maxit <- 100
for(i in 1:maxit) {
# h_i's

dat$h_o <- expit(as.matrix(dat[top_gate_vars]) %*% top_gate_param) *
  (
    expit(as.matrix(dat[o_g_o_gate_vars]) %*% o_g_o_gate_param) *
      expit(as.matrix(dat[expert_vars]) %*% o_g_o_expert_param)
    
    + 
      
      ( 1 - expit(as.matrix(dat[o_g_o_gate_vars]) %*% o_g_o_gate_param) ) *
      expit(as.matrix(dat[expert_vars]) %*% z_g_o_expert_param)
  )  / (
( 1 - expit(as.matrix(dat[top_gate_vars]) %*% top_gate_param) ) *
  (
expit(as.matrix(dat[o_g_z_gate_vars]) %*% o_g_z_gate_param) *
expit(as.matrix(dat[expert_vars]) %*% o_g_z_expert_param)
 
+ 

( 1 - expit(as.matrix(dat[o_g_z_gate_vars]) %*% o_g_z_gate_param) ) *
expit(as.matrix(dat[expert_vars]) %*% z_g_z_expert_param)
  ) +

 expit(as.matrix(dat[top_gate_vars]) %*% top_gate_param) *
   (
    expit(as.matrix(dat[o_g_o_gate_vars]) %*% o_g_o_gate_param) *
      expit(as.matrix(dat[expert_vars]) %*% o_g_o_expert_param)
    
    + 
      
      ( 1 - expit(as.matrix(dat[o_g_o_gate_vars]) %*% o_g_o_gate_param) ) *
      expit(as.matrix(dat[expert_vars]) %*% z_g_o_expert_param)
  ) )

dat$h_z <- 1 - dat$h_o

### h_j_g_i's

dat$h_z_g_z <- 

( 1 - expit(as.matrix(dat[o_g_z_gate_vars]) %*% o_g_z_gate_param) ) *
  expit(as.matrix(dat[expert_vars]) %*% z_g_z_expert_param) / (
    
  ( 1 - expit(as.matrix(dat[o_g_z_gate_vars]) %*% o_g_z_gate_param) ) *
    expit(as.matrix(dat[expert_vars]) %*% z_g_z_expert_param)
  
  + 
    
    expit(as.matrix(dat[o_g_z_gate_vars]) %*% o_g_z_gate_param) *
    expit(as.matrix(dat[expert_vars]) %*% o_g_z_expert_param)
)

dat$h_o_g_z <-  1 - dat$h_z_g_z

dat$h_z_g_o <- 
  
  ( 1 - expit(as.matrix(dat[o_g_o_gate_vars]) %*% o_g_o_gate_param) ) *
  expit(as.matrix(dat[expert_vars]) %*% z_g_o_expert_param) / (
    
    ( 1 - expit(as.matrix(dat[o_g_o_gate_vars]) %*% o_g_o_gate_param) ) *
      expit(as.matrix(dat[expert_vars]) %*% z_g_o_expert_param)
    
    + 
      
      expit(as.matrix(dat[o_g_o_gate_vars]) %*% o_g_o_gate_param) *
      expit(as.matrix(dat[expert_vars]) %*% o_g_o_expert_param)
  )

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
top_gate_mod <- glm(f_topgate, family = "quasibinomial", data = dat, method = glm.fit) # Not fit for debugging

# step 4
f_o_g_z_gate <- formula(paste0("h_o_g_z ~ -1 +", paste(o_g_z_gate_vars, collapse=" + ")))
o_g_z_gate_mod <- glm(f_o_g_z_gate, weights = h_z, family = "quasibinomial", data = dat)

#f_o_g_o_gate <- formula(paste0("h_o_g_o ~ -1 +", paste(o_g_o_gate_vars, collapse=" + ")))  # only if this need be updated
#o_g_o_gate_mod <- lm(f_o_g_o_gate, weights = h_o, data = dat)  # only if this need be updated


# step 5 update parameters

new_thetas <- c(coef(z_g_z_expert_mod), coef(o_g_z_expert_mod), coef(z_g_o_expert_mod), coef(o_g_o_expert_mod), coef(o_g_z_gate_mod))# coef(top_gate_mod), coef(o_g_o_gate_mod))
old_thetas <- c(z_g_z_expert_param, o_g_z_expert_param, z_g_o_expert_param, o_g_o_expert_param, o_g_z_gate_param)# top_gate_param, o_g_o_gate_param)

conv <- norm(new_thetas-old_thetas, type = "2")
if (conv > epsilon) {

z_g_z_expert_param <- coef(z_g_z_expert_mod)
o_g_z_expert_param <- coef(o_g_z_expert_mod)

z_g_o_expert_param <- coef(z_g_o_expert_mod)
o_g_o_expert_param <- coef(o_g_o_expert_mod)

top_gate_param <- coef(top_gate_mod) # No update for debugging

o_g_z_gate_param <- coef(o_g_z_gate_mod)
#o_g_o_gate_param <- coef(o_g_o_gate_mod) # only if this need be updated
} else break

### Debugging
print(paste(i, format(round(conv,2), nsmall = 2) ))

Ys1_hat <- c(Ys1_hat, 
mean   (
  expit(as.matrix(dat[o_g_o_gate_vars]) %*% o_g_o_gate_param) *
    expit(as.matrix(dat[expert_vars]) %*% o_g_o_expert_param)
  
  + 
    
    ( 1 - expit(as.matrix(dat[o_g_o_gate_vars]) %*% o_g_o_gate_param) ) *
    expit(as.matrix(dat[expert_vars]) %*% z_g_o_expert_param)
)
)

Ys0_hat <- c(Ys0_hat,
mean(  
  expit(as.matrix(dat[o_g_z_gate_vars]) %*% o_g_z_gate_param) *
    expit(as.matrix(dat[expert_vars]) %*% o_g_z_expert_param)
  
  + 
    
    ( 1 - expit(as.matrix(dat[o_g_z_gate_vars]) %*% o_g_z_gate_param) ) *
    expit(as.matrix(dat[expert_vars]) %*% z_g_z_expert_param)
)
)
}

matplot(as.matrix(data.frame(Ys0_hat = Ys0_hat, Ys1_hat = Ys1_hat)), type = c("b"),pch=1, col = 1:2)
data.frame(Hats=c(tail(Ys0_hat,1), tail(Ys1_hat,1)), True=c(TrueVals["Y_s0"], TrueVals["Y_s1"]) )

