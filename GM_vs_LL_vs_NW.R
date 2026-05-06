########################################################################
## A comparison of Dirichlet-kernel regression methods on the simplex ##
########################################################################

## Written by F. Ouimet (April 2026)

require("cubature")             # for integrals
require("doFuture")             # for parallel execution with foreach
require("fs")                   # for filesystem operations (dependency)
require("future.batchtools")    # for batchtools integration with future
require("ggplot2")              # for plotting
require("LaplacesDemon")        # for the Dirichlet distribution
require("parallel")             # for parallel execution of calculations
require("tidyverse")            # for data manipulation and visualization
require("writexl")              # to export output to Excel files

##############################
## Parallelization on cores ##
##############################

# Define the list of libraries to load on each cluster node

libraries_to_load <- c(
  "cubature",
  "doFuture",
  "fs",
  "future.batchtools",
  "ggplot2",
  "LaplacesDemon",
  "parallel",
  "tidyverse",
  "writexl"
)

# Define the list of variables/functions to export to the worker nodes

vars_to_export <- c(
  "BB",
  "JJ",
  "KK",
  "LOOCV",
  "ISE_MC",
  "MCsim",
  "MM",
  "RR",
  "adaptIntegrate",
  "b_opt_grid",
  "d",
  "hat_m",
  "m",
  "mesh",
  "mm",
  "n_to_k",
  "tol1"
)

# Sets up a parallel cluster, loads necessary libraries, and exports required variables globally

setup_parallel_cluster <- function() {
  num_cores <<- detectCores() - 1
  cl <- makeCluster(num_cores) # Create the cluster
  
  # Export the list of libraries to the worker nodes
  clusterExport(cl, varlist = "libraries_to_load")
  
  # Load necessary libraries on each cluster node
  invisible(clusterEvalQ(cl, {
    lapply(libraries_to_load, library, character.only = TRUE)
  }))
  
  # Export all necessary objects, functions, and parameters to the worker nodes
  clusterExport(cl, varlist = vars_to_export)
  
  return(cl) # Return the cluster object
}

# Initialize all variables in the list as NULL except vars_to_export and setup_parallel_cluster

invisible(
  lapply(
    vars_to_export[!(vars_to_export %in% c("vars_to_export", "setup_parallel_cluster"))],
    function(x) assign(x, NULL, envir = .GlobalEnv)
  )
)

################
# Set the path #
################

# path <- file.path(
#   "C://Users//fred1//Desktop//Github_DirichletKernelRegression",
#   fsep = .Platform$file.sep
# )
path <- getwd()
# setwd(path)

##############
# Parameters #
##############

d <- 2 # dimension of simplex
MCsim <- 10 ^ 3 # number of uniforms sampled for integral MC estimates

cores_per_node <- 64 # number of cores for each node in the super-computer
BB <- seq(0.01, 1, length.out = cores_per_node) # bandwidths for LOOCV graphs

MM <- list("GM", "NW", "LL") # list of Dirichlet-kernel methods
KK <- c(7, 10, 14) # indices for the mesh
JJ <- 1:6 # target regression function indices
RR <- 1:100 # replication indices

tol1 <- 1e-3

##############################
## Parallelization on nodes ##
##############################

resources_list <- list(
  cpus_per_task = cores_per_node,
  mem = "180G",
  walltime = "24:00:00",
  nodes = 1
  # Omit 'partition' to let SLURM choose
)

###########################
## Mesh of design points ##
###########################

# Mesh of the 2-dim simplex

mesh <- function(k) {
  w <- (k - 1 / sqrt(2)) / (k - 1)
  res <- list()
  for (i in 1:k) {
    for (j in i:k) {
      res <- append(res, list((c(w * (i - 1) + 0.5, w * (k - j) + 0.5) / (k + 1))))
    }
  }
  return(res)
}

# Number of points on each side of the mesh

n_to_k <- function(n) {
  (-1 + sqrt(1 + 8 * n)) / 2
}

#################################
## Target regression functions ##
#################################

# for one design point x

m <- function(j, x) { # x is a d-dim vector on the simplex
  if (j == 1) {
    # Case when j = 1
    res <- log(1 + x[1] + x[2])
  } else if (j == 2) {
    # Case when j = 2
    res <- sin(x[1]) + cos(x[2])
  } else if (j == 3) {
    # Case when j = 3
    res <- sqrt(x[1]) + sqrt(x[2])
  } else if (j == 4) {
    # Case when j = 4
    res <- x[1] * (1 + x[2])
  } else if (j == 5) {
    # Case when j = 5
    res <- (x[1] + 0.25) ^ 2 + (x[2] + 0.75) ^ 2
  } else if (j == 6) {
    # Case when j = 6
    res <- (1 + x[1]) * exp(x[2])
  } else {
    # Default case if j is not in 1:6
    warning("Invalid value of j. Should be between 1 and 6.")
    res <- NULL
  }
  return(res)
}

# for a mesh of design points x_1, ..., x_n

mm <- function(j, xx) { # xx is a list of d-dim vectors on the simplex
  res <- list()
  n <- length(xx)
  for (i in 1:n) {
    res <- append(res, list(m(j, xx[[i]])))
  }
  return(res)
}

################
## Estimators ##
################

hat_m <- function(xx, b, s, j, method, y = NULL) {
  # xx is a list of d-dim vectors on the simplex, s is a d-dim vector on the simplex
  n <- length(xx)
  d <- length(xx[[1]])
  
  if (is.null(y)) {
    y <- as.numeric(mm(j, xx)) # without random noise (this is not observed)
    y <- y + sqrt(0.1 * IQR(y)) * rnorm(n, 0, 1) # with random noise (this is observed)
  }
  
  u <- s / b + rep(1, d)
  v <- (1 - sum(s)) / b + 1
  
  switch(method,
         "GM" = { # GM is implemented specifically for the fixed mesh here, the general definition is different
           # GM method (Gasser-Muller)
           k <- n_to_k(n)
           k_round <- round(k)
           if (abs(k - k_round) < sqrt(.Machine$double.eps) && k_round * (k_round + 1) / 2 == n) {
             w <- (k_round - 1 / sqrt(2)) / (k_round - 1)
             half_width <- w / (k_round + 1) / 2
             # k indices for the closest points to the line y(x) = 1 - x
             cpi <- order(sapply(xx, function(x) abs(1 - x[1] - x[2]) / sqrt(2)))[1:k_round]
             integrand <- function(x) {
               # Check if x is inside the simplex (sum(x) <= 1)
               if (sum(x) >= 1) {
                 return(0)  # Return 0 if x is outside the simplex
               } else {
                 if (abs(1 - x[1] - x[2]) / sqrt(2) >= (1 / (2 * (k_round + 1)))) {  # if x is regular
                   for (i in 1:n) {
                     if (abs(x[1] - xx[[i]][1]) <= half_width && abs(x[2] - xx[[i]][2]) <= half_width) {
                       return(y[i] * LaplacesDemon::ddirichlet(c(x, 1 - sum(x)), c(u, v), log = FALSE))
                     }
                   }
                 } else { # if x is not regular
                   # Find the index i for which x is closest to xx[[i]] among all j in cpi
                   closest_i <- cpi[which.min(sapply(cpi, function(j) dist(rbind(x, xx[[j]]))))]
                   
                   # Return the result for the closest index
                   return(y[closest_i] * LaplacesDemon::ddirichlet(c(x, 1 - sum(x)), c(u, v), log = FALSE))
                 }
               }
             }
           } else {
             integrand <- function(x) {
               if (sum(x) >= 1) {
                 return(0)
               } else {
                 closest_i <- which.min(vapply(xx, function(z) sum((x - z) ^ 2), numeric(1)))
                 return(y[closest_i] * LaplacesDemon::ddirichlet(c(x, 1 - sum(x)), c(u, v), log = FALSE))
               }
             }
           }
           return(adaptIntegrate(integrand, lowerLimit = c(0, 0), upperLimit = c(1, 1), tol = tol1)$integral)
         },
         
         "LL" = {
           # LL method (Local Linear)
           kernel_vec <- rep(NA, n)
           for (i in 1:n) {
             kernel_vec[i] <- LaplacesDemon::ddirichlet(c(xx[[i]], 1 - sum(xx[[i]])), c(u, v), log = FALSE)
           }
           design_mat <- matrix(1, nrow = n, ncol = d + 1)
           for (i in 1:n) {
             design_mat[i, -1] <- xx[[i]] - s
           }
           W <- diag(kernel_vec)
           return(solve(t(design_mat) %*% W %*% design_mat, t(design_mat) %*% W %*% y)[1])
         },
         
         "NW" = {
           # NW method (Nadaraya-Watson)
           kernel_vec <- rep(NA, n)
           for (i in 1:n) {
             kernel_vec[i] <- LaplacesDemon::ddirichlet(c(xx[[i]], 1 - sum(xx[[i]])), c(u, v), log = FALSE)
           }
           return(sum(y * kernel_vec) / sum(kernel_vec))
         },
         
         stop("Invalid method. Choose either 'GM', 'LL', or 'NW'.")
  )
}

#########################################################
## Leave-One-Out Cross-Validation (LOOCV) exact version##
#########################################################

LOOCV <- function(xx, b, j, method, y = NULL) {
  n <- length(xx)
  d <- length(xx[[1]])
  
  if (is.null(y)) {
    y <- as.numeric(mm(j, xx))
    y <- y + sqrt(0.1 * IQR(y)) * rnorm(length(y), 0, 1)
  }
  
  switch(method,
         "GM" = {
           rss <- 0
           for (i in 1:n) {
             rss <- rss + (y[i] - hat_m(xx[-i], b, xx[[i]], j, method, y[-i])) ^ 2
           }
           return(rss / n)
         },
         
         "LL" = {
           kernel_mat <- matrix(NA, nrow = n, ncol = n)
           for (i in 1:n) {
             s <- xx[[i]]
             u <- s / b + rep(1, d)
             v <- (1 - sum(s)) / b + 1
             for (l in 1:n) {
               kernel_mat[i, l] <- LaplacesDemon::ddirichlet(c(xx[[l]], 1 - sum(xx[[l]])), c(u, v), log = FALSE)
             }
           }
           y_hat <- rep(NA, n)
           s_diag <- rep(NA, n)
           for (i in 1:n) {
             design_mat <- matrix(1, nrow = n, ncol = d + 1)
             for (l in 1:n) {
               design_mat[l, -1] <- xx[[l]] - xx[[i]]
             }
             W <- diag(kernel_mat[i, ])
             A_inv <- solve(t(design_mat) %*% W %*% design_mat)
             beta_hat <- A_inv %*% (t(design_mat) %*% W %*% y)
             y_hat[i] <- beta_hat[1]
             s_diag[i] <- kernel_mat[i, i] * A_inv[1, 1]
           }
           return(mean(((y - y_hat) / (1 - s_diag)) ^ 2))
         },
         
         "NW" = {
           kernel_mat <- matrix(NA, nrow = n, ncol = n)
           for (i in 1:n) {
             s <- xx[[i]]
             u <- s / b + rep(1, d)
             v <- (1 - sum(s)) / b + 1
             for (l in 1:n) {
               kernel_mat[i, l] <- LaplacesDemon::ddirichlet(c(xx[[l]], 1 - sum(xx[[l]])), c(u, v), log = FALSE)
             }
           }
           numerator <- as.numeric(kernel_mat %*% y)
           denominator <- rowSums(kernel_mat)
           return(mean(((y * denominator - numerator) / (denominator - diag(kernel_mat))) ^ 2))
         },
         
         stop("Invalid method. Choose either 'GM', 'LL', or 'NW'.")
  )
}

######################################
## Optimal Bandwidth (grid version) ##
######################################

# Function to find the optimal bandwidth using a grid search for LOOCV
b_opt_grid <- function(xx, j, method, return_LOOCV = FALSE, y = NULL) {
  if (is.null(y)) {
    y <- as.numeric(mm(j, xx))
    y <- y + sqrt(0.1 * IQR(y)) * rnorm(length(y), 0, 1)
  }
  
  # Determine the number of cores and set the grid size
  num_cores <- detectCores() - 1
  grid_size <- 1 * num_cores
  
  # Generate grid points between the lower and upper bounds of BB
  b_grid <- seq(min(BB), max(BB), length.out = grid_size)
  
  # Initialize parallel cluster
  cl <- setup_parallel_cluster()
  
  # Parallelize computation for all b values on the grid
  LOOCV_values <- parSapply(cl, b_grid, function(b, xx, j, method, y) {
    LOOCV(xx, b, j, method, y)
  }, xx = xx, j = j, method = method, y = y)
  
  # Stop the cluster after computation
  stopCluster(cl)
  
  # Find the index of the b that minimizes LOOCV_values
  min_index <- which.min(LOOCV_values)
  
  # Get the optimal b value and the corresponding LOOCV value
  b_opt_value <- b_grid[min_index]
  min_LOOCV_value <- LOOCV_values[min_index]
  
  # Return the desired value(s) based on the return_LOOCV argument
  if (return_LOOCV) {
    return(min_LOOCV_value)
  } else {
    return(b_opt_value)
  }
}

###############################################
## Integrated Squared Error (ISE_MC)         ##
###############################################

ISE_MC <- function(xx, j, method, y = NULL) {
  # 1. Generate the response variable y if not provided
  if (is.null(y)) {
    y <- as.numeric(mm(j, xx))
    y <- y + sqrt(0.1 * IQR(y)) * rnorm(length(y), 0, 1)
  }
  
  # 2. Find the optimal bandwidth (minimizing LOOCV on the design points)
  b_hat <- b_opt_grid(xx, j, method, return_LOOCV = FALSE, y = y)
  
  # 3. Generate 1,000 uniform points on the simplex
  d <- length(xx[[1]])
  # rdirichlet generates rows of (x1, x2, x3) that sum to 1. We just need the first d components.
  U_matrix <- LaplacesDemon::rdirichlet(MCsim, rep(1, d + 1)) 
  
  # 4. Calculate the squared errors at these uniform points
  sq_errors <- numeric(MCsim)
  for (i in 1:MCsim) {
    u <- U_matrix[i, 1:d]
    m_hat_val <- hat_m(xx, b_hat, s = u, j = j, method = method, y = y)
    m_true_val <- m(j, u)
    sq_errors[i] <- (m_hat_val - m_true_val)^2
  }
  
  # 5. Return the MC estimate of the ISE
  # The factor 2 in the denominator is the normalization constant for the uniform distribution on the 2D simplex
  return(mean(sq_errors) / 2)
}

###############
## Main code ##
###############

.libPaths("~/R/library")

# Disable the check for random number generation misuse in doFuture
options(doFuture.rng.onMisuse = "ignore")

# Register the doFuture parallel backend
registerDoFuture()

# Tweak the batchtools_slurm with the custom template and resources
myslurm <- tweak(
  batchtools_slurm,
  template = "batchtools.slurm.tmpl",
  resources = resources_list
)

# Set the plan for future
plan(list(myslurm, multisession))

raw_results <- data.frame(
  j = integer(),
  k = integer(),
  method = character(),
  ISE_MC = numeric(),
  stringsAsFactors = FALSE
)

# Capture the start time
start_time <- Sys.time()

# Parallel loop over the replications (RR), each node processes one set of RR values
res <- foreach(r = RR, .combine = "rbind", 
               .export = vars_to_export,
               .packages = libraries_to_load) %dopar% {
                 # Set a unique seed for each node (replication)
                 set.seed(r)
                 
                 # Set library paths within each worker node
                 .libPaths("~/R/library")
                 
                 local_raw_results <- data.frame(
                   j = integer(),
                   k = integer(),
                   method = character(),
                   ISE_MC = numeric(),
                   stringsAsFactors = FALSE
                 )
                 
                 # Loop over combinations of j, k, and method within each worker
                 for (j in JJ) {
                   for (k in KK) {
                     # Generate the mesh of design points once for each k
                     xx <- mesh(k)
                     y <- as.numeric(mm(j, xx))
                     y <- y + sqrt(0.1 * IQR(y)) * rnorm(length(y), 0, 1)
                     
                     for (method in MM) {
                       # Compute LSCV_MC for the current combination
                       ISE_MC_value <- ISE_MC(xx, j, method, y)
                       
                       # Store the raw results for this replication
                       local_raw_results <- rbind(
                         local_raw_results,
                         data.frame(
                           j = j,
                           k = k,
                           method = method,
                           ISE_MC = ISE_MC_value,
                           stringsAsFactors = FALSE
                         )
                       )
                     }
                   }
                 }
                 
                 # Return the raw results for this replication
                 return(local_raw_results)
               }

# Combine results from all nodes
raw_results <- res

# Stop parallel execution
plan(sequential)

# Calculate the duration in minutes
elapsed_time_minutes <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
print(paste("Elapsed time:", round(elapsed_time_minutes, 2), "minutes"))

# Save the raw results to an Excel file in the specified path
raw_output_file <- file.path(path, "raw_ISE_MC_results.csv")
write.csv(raw_results, raw_output_file, row.names = FALSE)

print("Raw results saved to raw_ISE_MC_results.csv")

#########################
## Process the results ##
#########################

# Create a data frame to store the summary results
summary_results <- data.frame(
  j = integer(),
  k = integer(),
  method = character(),
  mean_ISE_MC = numeric(),
  sd_ISE_MC = numeric(),
  median_ISE_MC = numeric(),
  IQR_ISE_MC = numeric(),
  stringsAsFactors = FALSE
)

# Loop through the results to compute the summary statistics
for (j in JJ) {
  for (k in KK) {
    for (method in MM) {
      # Filter the raw results by j, k, and method
      filtered_results <- raw_results %>%
        filter(j == !!j, k == !!k, method == !!method)
      
      ISE_values <- filtered_results$ISE_MC
      mean_ISE_MC <- mean(ISE_values)
      sd_ISE_MC <- sd(ISE_values)
      median_ISE_MC <- median(ISE_values)
      IQR_ISE_MC <- IQR(ISE_values)
      
      # Store the summary results
      summary_results <- rbind(
        summary_results,
        data.frame(
          j = j,
          k = k,
          method = method,
          mean_ISE_MC = mean_ISE_MC,
          sd_ISE_MC = sd_ISE_MC,
          median_ISE_MC = median_ISE_MC,
          IQR_ISE_MC = IQR_ISE_MC,
          stringsAsFactors = FALSE
        )
      )
    }
  }
}

# Save the summary results to an Excel file in the specified path
summary_output_file <- file.path(path, "ISE_MC_results.csv")
write.csv(summary_results, summary_output_file, row.names = FALSE)

print("Summary results saved to ISE_MC_results.csv")
