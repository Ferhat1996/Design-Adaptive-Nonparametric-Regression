# install.packages("stats")
library(stats)

# Define the Gaussian kernel function
K <- function(x) {
  return((1/sqrt(2*pi))*exp(-x^2/2))
}

# Define the local linear estimator
LL <- function(x, X, Y, h) {
  W <- K((x - X) / h)
  b <- lm(Y ~ X, weights = W)
  return(predict(b, newdata = list(X = x)))
}

# Define the Nadaraya-Watson estimator
NW <- function(x, X, Y, h) {
  W <- K((x - X) / h)
  return(sum(W * Y) / sum(W))
}

# Define the true regression function
M <- function(x) {
  return(sin(2.5 * x))
}

# Define the AMISE function
AMISE <- function(h, n, x) {
  return(0.006 * cos(2.5 * x)^2 * h^(4/5) + 0.08 / (n * h))
}

# Define the golden section search method
golden_section_search <- function(f, a, b, tol = 1e-5) {
  phi <- (1 + sqrt(5)) / 2
  c <- b - (b - a) / phi
  d <- a + (b - a) / phi
  while (abs(b - a) > tol) {
    if (f(c) < f(d)) {
      b <- d
    } else {
      a <- c
    }
    c <- b - (b - a) / phi
    d <- a + (b - a) / phi
  }
  return((a + b) / 2)
}

# Set parameters for different sample sizes
sample_sizes <- c(100, 200, 400)

# Initialize an empty data frame to store results
results_table <- data.frame(Sample_Size = numeric(),
                            LL_MISE = numeric(),
                            NW_MISE = numeric(),
                            Efficiency = numeric())

for (n in sample_sizes) {
  # Find the optimal bandwidth by minimizing the AMISE at x = 0
  h <- golden_section_search(function(h) AMISE(h, n, 0), 0.01, 0.5)
  R <- 300
  
  # Initialize MISEs
  LL_MISE <- 0
  NW_MISE <- 0
  
  # Run simulations
  for (r in 1:R) {
    
    set.seed(r)
    
    # Generate data
    s <- rbinom(n, 1, 0.5)
    X <- s * rnorm(n, mean = -1, sd = 1) + (1 - s) * rnorm(n, mean = 1.75, sd = 0.25)
    Y <- sin(2.5 * X) + 0.4 * rnorm(n)
    
    # Set grid of points
    x_grid <- seq(-2, 2, 0.1)
    
    # Initialize ISEs
    LL_ISE <- 0
    NW_ISE <- 0
    
    # Apply estimators
    for (x in x_grid) {
      LL_hat <- LL(x, X, Y, h)
      NW_hat <- NW(x, X, Y, h)
      M_true <- M(x)
      LL_ISE <- LL_ISE + (LL_hat - M_true)^2
      NW_ISE <- NW_ISE + (NW_hat - M_true)^2
    }
    
    # Update MISEs
    LL_MISE <- LL_MISE + LL_ISE / length(x_grid)
    NW_MISE <- NW_MISE + NW_ISE / length(x_grid)
  }
  
  # Calculate final MISEs
  LL_MISE <- LL_MISE / R
  NW_MISE <- NW_MISE / R
  
  # Calculate efficiency
  Efficiency <- sqrt(LL_MISE / NW_MISE)
  
  # Add results to the table
  results_table <- rbind(results_table, c(n, LL_MISE, NW_MISE, Efficiency))
}

# Add column headers
colnames(results_table) <- c("Sample_Size", "LL_MISE", "NW_MISE", "Efficiency")

# Print the results table
print(results_table)
