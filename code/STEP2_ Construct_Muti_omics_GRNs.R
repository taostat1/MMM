rm(list=ls())  # Clear current environment variables
renv::install("rTRNG")
renv::install("RcppArmadillo")
library("rTRNG")
library("RcppArmadillo")

Rcpp::sourceCpp("code/function_real.cpp")  # Load C++ helper functions

# Configuration: Define MR dataset parameters 
# Core parameters for MR dataset 
ds_name <- "MR"               # Dataset identifier
data_file <- "output/MRdata.RData"   
ds_seed <- 20250102           # Random seed for MR
K <- 4                        # Number of cell types for MR 

# Shared parameters 
num_iter <- 10000             # Total MCMC iterations
num_save <- 5000              # Number of iterations to retain (post-burn-in, i.e., keep last 5000)
protein_file <- "data/protein_matrix.csv"  # Path to protein interaction matrix file
ssp_v0 <- 0.02                # Spike-slab prior parameters (unchanged from original)
ssp_v1 <- 1
ssp_l <- 1


##############################################################################
# Step 1: Load and preprocess protein interaction matrix (retains original logic)
##############################################################################
# Read protein interaction matrix
protein_df <- read.csv(file = protein_file, header = FALSE)
protein_mat <- as.matrix(protein_df)
# Ensure matrix is numeric (binary 0/1: 0=no interaction, 1=interaction)
protein_mat <- as.numeric(protein_mat)
# Restore matrix dimensions (prevent dimension loss after vector conversion)
dim(protein_mat) <- c(nrow(protein_df), ncol(protein_df))

# Enforce binary format for protein matrix (prevent interference from non-0/1 values)
protein_mat[protein_mat != 0] <- 1  # Treat non-0 values as "interaction" (1)
protein_mat[protein_mat == 0] <- 0  # Treat 0 values as "no interaction"
cat(paste0("Protein interaction matrix dimensions: ", nrow(protein_mat), "×", ncol(protein_mat), "\n"))


##############################################################################
# Step 2: Process MR dataset
##############################################################################
cat(paste0("\n===== Processing ", ds_name, " Data =====\n"))

# --------------------------
# 1. Load and prepare MR data
# --------------------------
load(data_file)  # Load MR data (assumes gene expression matrix 'X' is included)
X <- as.matrix(X)  # Ensure X is in matrix format
G <- dim(X)[1]     # Number of genes (rows: genes)
N <- dim(X)[2]     # Number of cells (columns: cells)
cat(paste0(ds_name, " dataset dimensions: ", G, " genes × ", N, " cells\n"))

# Validate compatibility between protein matrix and MR gene count (critical check to avoid dimension mismatch)
if (nrow(protein_mat) != G || ncol(protein_mat) != G) {
  stop(paste0(
    "Protein matrix dimension mismatch with ", ds_name, "!\n",
    "Protein matrix: ", nrow(protein_mat), "×", ncol(protein_mat), 
    "; ", ds_name, " genes: ", G
  ))
}
cat(paste0(ds_name, " gene count matches protein matrix dimensions. Proceeding...\n"))


# --------------------------
# 2. Initialize cell type labels (k-means clustering)
# --------------------------
set.seed(ds_seed)  # Use MR-specific seed to ensure reproducible clustering results
# Perform k-means on log-transformed expression matrix (transpose to cells×genes for k-means input format)
group_t <- kmeans(t(log(X + 1)), centers = K, nstart = 50)$cluster
# Print initial cell type counts (check if clustering distribution is reasonable)
cat(paste0(ds_name, " initial cell type counts:\n"))
print(table(group_t))


# --------------------------
# 3. Initialize theta_t (zero imputation and log transformation for expression matrix)
# --------------------------
# Matrix marking zeros in expression data (for subsequent imputation)
ind_zero <- (X == 0)
theta_t <- X  # Initialize theta_t with original expression matrix

# Perform zero imputation by gene and cell type (fill zeros with 5th percentile of non-zero values in the cell type)
for (g in 1:G) {  # Iterate over each gene
  ind_zero_g <- ind_zero[g, ]  # Zero positions for the current gene
  X_g <- X[g, ]                # Original expression values for the current gene
  
  for (k in 1:K) {  # Iterate over each cell type
    # Extract expression values of the current gene in cell type k
    X_gk <- X_g[group_t == k]
    # Zero positions of the current gene in cell type k
    ind_zero_gk <- ind_zero_g[group_t == k]
    
    # Imputation logic: If non-zero values exist, fill zeros with 5th percentile; if all values are zero, retain 0
    if (sum(!ind_zero_gk) > 0) {
      theta_gk <- X_gk
      theta_gk[ind_zero_gk] <- quantile(X_gk[!ind_zero_gk], na.rm = TRUE, probs = 0.05)
      theta_t[g, group_t == k] <- theta_gk
    } else {
      theta_t[g, group_t == k] <- 0
    }
  }
}

# Log transformation (compress value range to fit normal distribution assumptions)
theta_t <- log(theta_t)
# Handle potential NA values (e.g., -Inf from log(0), replace with 0)
theta_t[is.na(theta_t) | theta_t == -Inf] <- 0
cat(paste0(ds_name, " theta_t initialization completed (zero imputation + log transformation)\n"))


# --------------------------
# 4. Initialize statistical parameters (mean, covariance, precision matrices)
# --------------------------
# mu_t: Gene×cell type mean matrix (average gene expression per cell type)
mu_t <- matrix(0, nrow = G, ncol = K)
# cov_t: Gene×gene×cell type covariance matrix (diagonal matrix, retains only gene-specific variance)
cov_t <- array(diag(1, G), dim = c(G, G, K))
# invcov_t: Inverse of covariance matrix (precision matrix, used for subsequent MCMC sampling)
invcov_t <- array(diag(1, G), dim = c(G, G, K))

# Calculate mean and covariance by cell type
for (k in 1:K) {
  # Extract theta_t for cell type k (genes×cells)
  theta_k <- theta_t[, group_t == k]
  # Gene mean for cell type k
  mu_t[, k] <- rowMeans(theta_k)
  
  # Calculate gene-specific variance (diagonal covariance, ignore inter-gene covariance for simplicity)
  # Formula: Sum of squared deviations from mean / (number of cells - 1)
  gene_var <- diag((theta_k - mu_t[, k]) %*% t(theta_k - mu_t[, k]) / (sum(group_t == k) - 1))
  # Avoid zero variance (replace with 1 to prevent errors in inverse matrix calculation)
  gene_var[gene_var == 0] <- 1
  # Assign to covariance and precision matrices
  diag(cov_t[,, k]) <- gene_var
  diag(invcov_t[,, k]) <- 1.0 / gene_var  # Precision = 1/variance
}
cat(paste0(ds_name, " Mean/covariance/precision matrix initialization completed\n"))


# --------------------------
# 5. Initialize edge matrix (edge_t): Based on protein interaction matrix
# --------------------------
# edge_t: Gene×gene×cell type edge matrix (0=no connection, 1=connection)
edge_t <- array(0, dim = c(G, G, K))
# Initial edges for each cell type are based on the protein interaction matrix (updated later in MCMC)
for (k in 1:K) {
  edge_t[,, k] <- protein_mat
}
cat(paste0(ds_name, " Edge matrix initialization completed (based on protein interaction matrix)\n"))


# --------------------------
# 6. Final parameter initialization (for MCMC sampling)
# --------------------------
# Transpose theta_t to cells×genes (matches MCMC input format)
theta_t <- t(theta_t)
# Transpose zero-marking matrix (matches theta_t dimensions)
ind_zero <- t(ind_zero)

# gam: Prior weights for cell types (initialized to 1 for no bias)
gam <- rep(1, K)
# pi_t: Cell type probabilities (updated based on initial clustering results)
pi_t <- update_pi_R(group_t, gam, K)  # Ensure the update_pi_R function is pre-defined

# Spike-slab prior parameters (gene-level priors to control edge existence probability)
lambda0_t <- rnorm(G, mean = 3, sd = 0.2)    # Prior parameter for non-edges (spike)
lambda1_t <- rnorm(G, mean = -2, sd = 0.2)   # Prior parameter for edges (slab)
ssp_xi <- 2 / (G - 1)                        # Adjustment term based on number of genes
cat(paste0(ds_name, " MCMC sampling parameter initialization completed\n"))


# --------------------------
# 7. Run MCMC sampling (core step, retains original parameter settings)
# --------------------------
cat(paste0("Starting ", ds_name, " MCMC sampling...\n"))
start_time <- Sys.time()

# Call MCMC_full function (ensure this function is pre-defined with consistent parameters)
Result <- MCMC_full(
  num_iter = num_iter,          # Total number of iterations
  num_save = num_save,          # Number of iterations to retain
  theta_t = theta_t,            # Imputed expression matrix (cells×genes)
  ind_zero = ind_zero,          # Zero-marking matrix
  mu_t = mu_t,                  # Gene mean matrix
  invcov_t = invcov_t,          # Precision matrix
  cov_t = cov_t,                # Covariance matrix
  edge_t = edge_t,              # Initial edge matrix
  group_t = group_t,            # Initial cell type labels
  lambda0_t = lambda0_t,        # Spike-slab prior lambda0
  lambda1_t = lambda1_t,        # Spike-slab prior lambda1
  pi_t = pi_t,                  # Cell type probabilities
  gam = gam,                    # Prior weights for cell types
  G = G,                        # Number of genes
  N = N,                        # Number of cells
  K = K,                        # Number of cell types
  ssp_v0 = ssp_v0,              # Spike-slab prior v0
  ssp_v1 = ssp_v1,              # Spike-slab prior v1
  ssp_l = ssp_l,                # Spike-slab prior l
  ssp_xi,              # Spike-slab prior xi
  epsilon_theta = 0.2,          # Step size for theta updates
  num_step_theta = 20,          # Number of steps for theta updates
  eta_mu = 0,                   # Prior mean for mu
  tau_sq_mu = 1,                # Prior variance for mu
  lam0_0 = 3,                   # Prior mean for lambda0
  lam1_0 = -2,                  # Prior mean for lambda1
  sigma2_lam0 = 0.2,            # Prior variance for lambda0
  sigma2_lam1 = 0.2,            # Prior variance for lambda1
  epsilon_lam = 0.01,           # Step size for lambda updates
  num_step_lam = 10             # Number of steps for lambda updates
)

# Calculate and print MCMC runtime
end_time <- Sys.time()
runtime <- difftime(end_time, start_time, units = "mins")
cat(paste0(ds_name, " MCMC runtime: ", round(runtime, 2), " minutes\n"))


# --------------------------
# 8. Extract and save MR results
# --------------------------
# Extract posterior estimates from MCMC results (mean of retained iterations)
lam0_est <- rowMeans(Result$lam0_iter)    # Posterior estimate of lambda0
lam1_est <- rowMeans(Result$lam1_iter)    # Posterior estimate of lambda1
mu_est <- Result$mu_est                   # Posterior estimate of gene means
invcov_est <- Result$invcov_est           # Posterior estimate of precision matrix
edge_est <- Result$edge_est               # Posterior estimate of edge matrix (probability values)
group_est <- Result$group_est             # Posterior estimate of cell type labels
Pi_est <- Result$Pi_est                   # Posterior estimate of cell type probabilities
theta_est <- Result$theta_est             # Posterior estimate of expression matrix

# Save results (filename fixed for MR)
result_file <- paste0("output/",paste0(ds_name, "Data_result_with_protein.RData"))
save(
  list = c("X", "invcov_est", "edge_est", "group_est", "lam0_est", 
           "lam1_est", "mu_est", "Pi_est", "theta_est", "protein_mat"),
  file = result_file
)
cat(paste0(ds_name, " results saved to: ", result_file, "\n"))




