library(grf)

set.seed()

# Prepare dataset for causal forest.
# We assume a dataframe with a disease, exposure, and various characteristics including PRS in each column.
X <- data[, character_list]
W <- data$exposure
Y <- data$disease

# Prepare for data.splitting.
# Assign a fold number to each observation.
num.folds <- 10
n<-nrow(data)
folds <- sample(c(1:num.folds), replace=TRUE, size=n)

# Apply causal forest to estimate ITEs of the exposure on the disease. 
# Randomized settings with fixed and known probabilities (here: 0.5).
forest <- causal_forest(X, Y, W, 
                        clusters=folds,
                        tune.parameters="all")
# Predict the ITEs.
predictions <- predict(forest)
tau.hat <- predictions$predictions

# Plot a partial dependence plot of PRSs on ITEs.
plot(X[, "PRS"], tau.hat, 
     xlab = "PRS", ylab = "ITE")

# Assess the slope of the calibration line between predicted ARR and observed ARR.
test_calibration(forest) 
