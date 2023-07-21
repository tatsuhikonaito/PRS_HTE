
library(grf)

# Specify a random seed.
set.seed()
character_list = c("PRS", "x_int", paste0("x_", seq(1, 10)))

args <- commandArgs(trailingOnly = T)
n <- as.integer(args[1])
f_disease <- as.numeric(args[2])
f_env <- as.numeric(args[3])
a_PRS <- as.numeric(args[4])
a_env <- as.numeric(args[5])
a_env_PRS <- as.numeric(args[6])
a_int <- as.numeric(args[7])
a_int_env_x <- as.numeric(args[8])
a_x <- as.numeric(args[9])

run_grf <- function(data, character_list){
  # Prepare dataset for causal forest.
  # We assume a dataframe with a disease, exposure, and various characteristics including PRS in each column.
  X <- data[, character_list]
  if (length(character_list)==1){
    X = matrix(X)
  }
  W <- data$env
  Y <- data$disease
  # Prepare for data.splitting.
  # Assign a fold nber to each observation.
  n.folds <- 10
  n <- nrow(data)
  folds <- sample(c(1:n.folds), replace=TRUE, size=n)
  # Apply causal forest to estimate ITEs of the exposure on the disease. 
  # Randomized settings with fixed and known probabilities (here: 0.5).
  forest <- causal_forest(X, Y, W, 
                          clusters=folds,
                          tune.parameters="all")
  return(forest)
}

# Detemine an intercept based on the gradient descent method
determine_intercept <- function(f, array_covar, n_iter=50, e=1){
  intercept <- 0
  res_test <- c()
  for (i in 1:n_iter){
    exp_denom_prob <- intercept+array_covar
    f_tmp <- 1/(1+exp(-exp_denom_prob))
    Y <- rbinom(n, 1, f_tmp)
    diff <- mean(Y) - f
    print(mean(Y))
    print(intercept)
    intercept <- intercept - e*diff
    res_test <- rbind(res_test, c(diff, intercept))
  }
  return (res_test[which.min(abs(res_test[,1])), 2])
}

simulation <- function(n, f_disease, f_env, a_PRS, a_env, a_env_PRS, a_int, a_int_env_x, int_PRS_env_x, a_x, n_variable=10){
  prefix = paste(paste0("n", n), 
                 paste0("f_disease", f_disease), 
                 paste0("f_env", f_env), 
                 paste0("a_PRS", a_PRS), 
                 paste0("a_env", a_env), 
                 paste0("a_env_PRS", a_env_PRS), 
                 paste0("a_int", a_int), 
                 paste0("a_int_env_x", a_int_env_x), 
                 paste0("a_x", a_x), 
                 sep = ".")

  # Sample size
  n_case <- as.integer(n*f_disease)
  n_control <- n - n_case
  n1 <- as.integer(n/2)  # The number of samples whose x is assigned 1
  
  # Dataframe to store results
  df <- data.frame(matrix(rep(NaN, 6*n), ncol=6))
  colnames(df) <- c("disease", "env", "PRS", "ITE", "ATE_strat_10", "TE_interact")

  # Generate PRS values based on standardized normal distribution
  df$PRS <- rnorm(n, 0, 1)

  # Assign an environmental factor
  if (log(a_env_PRS)==0){
    df$env <- rbinom(n, 1, f_env)
  } else {
    array_covar <- log(a_env_PRS)*df$PRS
    for (i in 1:5){
      df[, paste0("xenv_", i)] <- rnorm(n, 0, 1)
    }
    for (i in 1:5){
      array_covar <- array_covar + 1.5*df[, paste0("xenv_", i)]
    }
    intercept_env <- determine_intercept(f_env, array_covar)
    exp_denom_prob_PRS <- intercept_env+array_covar
    f_env_PRS <- 1/(1+exp(-exp_denom_prob_PRS))
    df$env <- rbinom(n, 1, f_env_PRS)
  }

  # Generate other variables
  for (i in 1:n_variable){
    df[, paste0("x_", i)] <- rnorm(n, 0, 1)
  }
  df$x_int <- c(rep(0, n1), rep(1, n-n1))  # This varible can generate two models with different environmental factor effects

  # Generate disease status
  array_covar <- log(a_PRS)*df$PRS+log(a_env)*df$env+log(a_int)*df$PRS*df$env+log(a_int_env_x)*df$env*df$x_int
  for (i in 1:n_variable){
    array_covar <- array_covar + log(a_x)*df[, paste0("x_", i)]
  }
  intercept <- determine_intercept(f_disease, array_covar)
  exp_denom_prob <- intercept+array_covar
  prob <- 1/(1+exp(-exp_denom_prob))
  df$disease <- rbinom(n, 1, prob)


  ## Predict ITEs using causal forest
  forest <- run_grf(df, character_list)
  predictions <- predict(forest)
  df$ITE <- predictions$predictions
  res_calib = test_calibration(forest)

  # Generate expected treatment effects
  n_iter <- 10000
  n_point <- 1000
  data_exp <- data.frame(matrix(rep(NaN, 2*n_point), ncol=2))
  colnames(data_exp) <- c("x", "y")
  data_exp$x <- seq(-4, 4, length=n_point)
  for (i in 1:dim(data_exp)[1]){
    x <- data_exp$x[i]
    rd <- 0
    for (w in c(0, 1)){
      exp_denom_prob <- intercept+log(a_PRS)*x+log(a_env)*w+log(a_int)*x*w
      prob <- mean(sapply(seq(n_iter), 
                         function(x){1/(1+exp(-(exp_denom_prob+sum(log(a_x)*rnorm(n_variable, 0, 1))+log(a_int_env_x)*w*rbinom(1, 1, 0.5))))}))
      if(w==0){
        rd <- rd - prob
      } else if(w==1){
        rd <- rd + prob
      }
    }
    data_exp[i, "y"] <- rd
  }

  # Generate expected treatment effects in stratified groups
  if (a_int_env_x!=1){
    data_exp_strat <- data.frame(matrix(rep(NaN, 3*n_point), ncol=3))
    colnames(data_exp_strat) <- c("x", "y0", "y1")
    data_exp_strat$x <- seq(-4, 4, length=n_point)
    for (j in c(0, 1)){
      for (i in 1:dim(data_exp_strat)[1]){
        x <- data_exp_strat$x[i]
        rd <- 0
        for (w in c(0, 1)){
          exp_denom_prob <- intercept+log(a_PRS)*x+log(a_env)*w+log(a_int)*x*w+log(a_int_env_x)*w*j
          prob <- mean(sapply(seq(n_iter), 
                             function(x){1/(1+exp(-(exp_denom_prob+sum(log(a_x)*rnorm(n_variable, 0, 1)))))}))
          if(w==0){
            rd <- rd - prob
          } else if(w==1){
            rd <- rd + prob
          }
        }
        data_exp_strat[i, paste0("y", j)] = rd
      }
    }
  }

  ## Evaluate ATEs in groups stratified by PRS values
  df <- df[order(df$PRS),]
  n_split <- 10
  for (split in 1:n_split){
    df_split <- df[((n%/%n_split)*(split-1)+1):(min(c(n%/%n_split)*split, n)), ]
    model <- lm(as.formula(paste("disease", paste(c("PRS", "env", "x_int", paste0("x_", seq(n_variable))), collapse="+"), sep="~")), 
               data=df_split)
    ATE_strat <- summary(model)$coef["env", "Estimate"]
    df[((n%/%n_split)*(split-1)+1):(min(c(n%/%n_split)*split, n)), paste0("ATE_strat_", n_split)] <- ATE_strat
  }

  # Evaluate TEs based on additive interaction analysis
  model <- lm(as.formula(paste("disease", paste(c("PRS", "env", "int_PRS_env", "x_int", paste0("x_", seq(n_variable))), collapse="+"), sep="~")), 
             data=df)
  res_int <- summary(model)$coef[, c("Estimate", "Std. Error", "Pr(>|t|)")]
  colnames(res_int) <- c("beta", "SE", "P")
  beta_interact <- res_int["env", "beta"]
  slope_interact <- res_int["int_PRS_env", "beta"]
  df$TE_interact <- beta_interact + slope_interact * df$PRS


  # Save simutation results.
  write.table(df, 
              paste("simulation_results", prefix, "txt", sep=".")),
              quote=F, row.names = F, col.names =T)
  write.table(res_calib, 
              paste("res_calib", prefix, "txt", sep=".")),
              quote=F, row.names = F, col.names =T)
  write.table(data_exp, 
              paste("exp", prefix, "txt", sep=".")),
              quote=F, row.names = F, col.names =T)
  if (a_int_env_x!=1){
    write.table(data_exp_strat, 
                paste("exp_strat", prefix, "txt", sep=".")),
                quote=F, row.names = F, col.names =T)
  }
}

simulation(n, f_disease, f_env, a_PRS, a_env, a_env_PRS, a_int, a_int_env_x, a_x)
