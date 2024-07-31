# Define the Lagrangian function (objective)
lagrangian <- function(lambda, y, y_tr) {
  T <- length(y)  # Number of time steps
  lagrangian_value <- 0
  
  for (t in 1:T) {
    lagrangian_value <- lagrangian_value + (y[t] - y_tr[t])^2
  }
  
  for (t in 2:(T - 1)) {
    lagrangian_value <- lagrangian_value + ((y_tr[t + 1] - y_tr[t]) - (y_tr[t] - y_tr[t - 1]))^2
  }
  
  V <- sum(((diff(y_tr) - diff(diff(y_tr)))^2)) / sum((y - y_tr)^2)
  
  # Calculate the Lagrangian equation
  lagrangian_value <- (1 - lambda * V) * lagrangian_value + lambda * sum(((diff(y_tr) - diff(diff(y_tr)))^2))
  
  return(lagrangian_value)
}

# Define your time series data y and y_tr
y <- c(1, 2, 3, 4, 5)
y_tr <- c(0.5, 1.5, 2.5, 3.5, 4.5)

# Calculate V
V <- sum(((diff(y_tr) - diff(diff(y_tr)))^2)) / sum((y - y_tr)^2)

# Set a range of lambda values from 0 to 1/V
lambda_range <- seq(0, 1/V, length.out = 100)

# Initialize a vector to store Lagrangian values
lagrangian_values <- numeric(length(lambda_range))

# Calculate Lagrangian values for each lambda in the range
for (i in 1:length(lambda_range)) {
  lagrangian_values[i] <- lagrangian(lambda_range[i], y = y, y_tr = y_tr)
}

# Find the lambda that minimizes the Lagrangian value
optimal_lambda_index <- which.min(lagrangian_values)
optimal_lambda <- lambda_range[optimal_lambda_index]
minimized_lagrangian_value <- lagrangian_values[optimal_lambda_index]

# Print the results
cat("Optimal lambda:", optimal_lambda, "\n")
cat("Minimized Lagrangian value:", minimized_lagrangian_value, "\n")
