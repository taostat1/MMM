rm(list=ls())  # Clear current environment variables
K <- 4 # Change according to needs
library("Matrix")
library("rTRNG")
library(igraph)
Rcpp::sourceCpp("code/function_real.cpp")
# Load core data 
load("output/goal_gene_mt.RData")
load("output/Matched_goal_expr_PR.RData")

#########################Constructing PR GRNs######################################################
# Shared parameters (used for both datasets)
num_iter <- 10000       # Total MCMC iterations
num_save <- 5000        # Iterations to retain (post-burn-in)
protein_file <- "data/protein_matrix.csv"  # Protein interaction matrix
ssp_v0 <- 0.02          # Spike-slab prior parameters
ssp_v1 <- 1
ssp_l <- 1


##############################################################################
# Step 1: Load and preprocess protein interaction matrix
##############################################################################
# Read protein matrix
protein_df <- read.csv(file = protein_file, header = FALSE)
protein_mat <- as.matrix(protein_df)
protein_mat <- as.numeric(protein_mat)  # Ensure numeric values (0/1)
dim(protein_mat) <- c(nrow(protein_df), ncol(protein_df))  # Restore dimensions

# Verify protein matrix is binary
protein_mat[protein_mat != 0] <- 1
protein_mat[protein_mat == 0] <- 0


# --------------------------
# 2. Process each PR-specific group (PR1 to PR4) in a loop.
# --------------------------
pr_network_results <- list()
pr_classes <- names(matched_goal_expr_PR)  
pr_class <- "1"

K_pr <- 1
for (pr_class in pr_classes) {
  cat(paste0("\n===== Start constructing PR GRNs", pr_class, " =====\n"))
  
  # --------------------------
  # Step 1: Extract expression data from the current PR group and preprocess it.
  # --------------------------
  X <- t(as.matrix(matched_goal_expr_PR[[pr_class]]))  # 核心数据：PR细胞的目标基因表达
  G <- dim(X)[1]     # Number of genes (rows: genes)
  N <- dim(X)[2]     # Number of cells (columns: cells)
  
  # Validate compatibility between protein matrix and MR gene count (critical check to avoid dimension mismatch)
  if (nrow(protein_mat) != G || ncol(protein_mat) != G) {
    stop(paste0(
      "Protein matrix dimension mismatch with ", pr_class, "!\n",
      "Protein matrix: ", nrow(protein_mat), "×", ncol(protein_mat), 
      "; ", pr_class, " genes: ", G
    ))
  }
  cat(paste0(pr_class, "class gene count matches protein matrix dimensions. Proceeding...\n"))
  
  
  # --------------------------
  # 2. Initialize cell type labels (k-means clustering)
  # --------------------------
  set.seed(20250102)  # Use MR-specific seed to ensure reproducible clustering results
  # Perform k-means on log-transformed expression matrix (transpose to cells×genes for k-means input format)
  group_t <- kmeans(t(log(X + 1)), centers = K_pr, nstart = 50)$cluster
  # Print initial cell type counts (check if clustering distribution is reasonable)
  cat(paste0(pr_class, " initial cell type counts:\n"))
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
    
    for (k in 1:K_pr) {  # Iterate over each cell type
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
  cat(paste0(pr_class, " theta_t initialization completed (zero imputation + log transformation)\n"))
  
  
  # --------------------------
  # 4. Initialize statistical parameters (mean, covariance, precision matrices)
  # --------------------------
  # mu_t: Gene×cell type mean matrix (average gene expression per cell type)
  mu_t <- matrix(0, nrow = G, ncol = K_pr)
  # cov_t: Gene×gene×cell type covariance matrix (diagonal matrix, retains only gene-specific variance)
  cov_t <- array(diag(1, G), dim = c(G, G, K_pr))
  # invcov_t: Inverse of covariance matrix (precision matrix, used for subsequent MCMC sampling)
  invcov_t <- array(diag(1, G), dim = c(G, G, K_pr))
  
  # Calculate mean and covariance by cell type
  for (k in 1:K_pr) {
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
  cat(paste0(pr_class, " Mean/covariance/precision matrix initialization completed\n"))
  
  
  # --------------------------
  # 5. Initialize edge matrix (edge_t): Based on protein interaction matrix
  # --------------------------
  # edge_t: Gene×gene×cell type edge matrix (0=no connection, 1=connection)
  edge_t <- array(0, dim = c(G, G, K_pr))
  # Initial edges for each cell type are based on the protein interaction matrix (updated later in MCMC)
  for (k in 1:K_pr) {
    edge_t[,, k] <- protein_mat
  }
  cat(paste0(pr_class, " Edge matrix initialization completed (based on protein interaction matrix)\n"))
  
  
  # --------------------------
  # 6. Final parameter initialization (for MCMC sampling)
  # --------------------------
  # Transpose theta_t to cells×genes (matches MCMC input format)
  theta_t <- t(theta_t)
  # Transpose zero-marking matrix (matches theta_t dimensions)
  ind_zero <- t(ind_zero)
  
  # gam: Prior weights for cell types (initialized to 1 for no bias)
  gam <- rep(1, K_pr)
  # pi_t: Cell type probabilities (updated based on initial clustering results)
  pi_t <- update_pi_R(group_t, gam, K_pr)  # Ensure the update_pi_R function is pre-defined
  
  # Spike-slab prior parameters (gene-level priors to control edge existence probability)
  lambda0_t <- rnorm(G, mean = 3, sd = 0.2)    # Prior parameter for non-edges (spike)
  lambda1_t <- rnorm(G, mean = -2, sd = 0.2)   # Prior parameter for edges (slab)
  ssp_xi <- 2 / (G - 1)                        # Adjustment term based on number of genes
  cat(paste0(pr_class, " MCMC sampling parameter initialization completed\n"))
  
  
  # --------------------------
  # 7. Run MCMC sampling (core step, retains original parameter settings)
  # --------------------------
  cat(paste0("Starting ", pr_class, " MCMC sampling...\n"))
  start_time <- Sys.time()
  
  # Call MCMC_full function 
  Result <- MCMC_full(num_iter, num_save, theta_t, ind_zero, mu_t, invcov_t, cov_t, edge_t = edge_t,
                      group_t, lambda0_t, lambda1_t, pi_t, gam, G, N, K_pr,
                      ssp_v0, ssp_v1, ssp_l, ssp_xi,
                      epsilon_theta = 0.2, num_step_theta = 20, eta_mu = 0, tau_sq_mu = 1,
                      lam0_0 = 3, lam1_0 = -2, sigma2_lam0 = 0.2, sigma2_lam1 = 0.2,
                      epsilon_lam = 0.01, num_step_lam = 10)
  
  # Calculate and print MCMC runtime
  end_time <- Sys.time()
  runtime <- difftime(end_time, start_time, units = "mins")
  cat(paste0(pr_class, "class MCMC runtime: ", round(runtime, 2), " minutes\n"))
  
  
  # --------------------------
  # 8. Extract and save PR results
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
  
  # Save results 
  pr_network_results[[pr_class]] <- edge_est
  
  # --------------------------
  # 9. Generate and save PR cell type-specific network plots
  # --------------------------
  # Binarize edge matrix: Treat posterior probability ≥0.5 as "edge exists" (1), else 0
  edge_est_mat <- array(0, dim = c(G, G, dim(edge_est)[3]))
  edge_est_mat[edge_est >= 0.5] <- 1
  
  # Prepare edge matrices with gene names for each cell type
  gene_names <- rownames(X)  # Get gene names from original expression matrix
  edge_full <- list()        # Store edge matrix for each cell type
  for (k in 1:K_pr) {
    edge_full[[k]] <- edge_est_mat[,, k]
    rownames(edge_full[[k]]) <- gene_names
    colnames(edge_full[[k]]) <- gene_names
  }
  
  # Define network plotting function
  plot_mr_network <- function(edge_data, num) {
    # Construct undirected graph from edge matrix
    G_graph <- graph_from_adjacency_matrix(
      edge_data,
      mode = "undirected",  # Protein interaction networks are undirected
      add.rownames = "gene" # Use row names (gene names) as node labels
    )
    
    # Set plot layout (circular layout for global connection visualization)
    plot_layout <- layout.circle(G_graph)
    node_labels <- V(G_graph)$gene  # Node labels (gene names)
    node_color <- rep("#EC7357", nrow(edge_data))  # Node color (orange, consistent with original settings)
    
    # Save as PDF (MR-specific filename)
    pdf_file <- paste0("output/figure/",paste0(pr_class, "pr_edge", num, ".pdf"))
    pdf(file = pdf_file, width = 14, height = 12)  # Set plot dimensions
    
    # Plot network
    plot(
      G_graph,
      layout = plot_layout,
      vertex.label = node_labels,    # Show gene names
      vertex.label.color = node_color, # Gene name color
      vertex.label.cex = 1.1,        # Gene name font size
      vertex.size = 8,               # Node size
      vertex.color = "white",        # Node fill color
      vertex.frame.color = "#030303",# Node border color (black)
      edge.color = "#030303",        # Edge color (black)
      edge.width = 2,                # Edge width
      margin = c(1, 1, 1, 1)         # Plot margins
    )
    
    dev.off()  # Close PDF device
    cat(paste0(pr_class, " network ", num, " saved to: ", pdf_file, "\n"))
  }
  
  # Plot and save network graphs for each MR cell type
  for (k in 1:K_pr) {
    # Adjust edge matrix order 
    plot_data <- edge_full[[k]][c(1, nrow(edge_full[[k]]):2), c(1, nrow(edge_full[[k]]):2)]
    # Call plotting function
    plot_mr_network(edge_data = plot_data, num = k)
  }
}


# --------------------------
# 3. Save all PR group network results 
# --------------------------
save(
  pr_network_results,
  file = "output/PR_4Groups_Network_Results.RData"
)

#########################Constructing Specific GRNs######################################################
# Load grouped data
load("output/MRData_result_with_protein.RData")
mr_edge_est <-edge_est

i <- 1
edge_diff_list <- list()
for (i in pr_classes) {
  mr_edge <- mr_edge_est[,, as.numeric(i)]            
  pr_edge <- pr_network_results[[i]][,,1]  
  
  # Calculate the difference (MR - PR)
  edge_diff <- mr_edge - pr_edge
  
  # Store the difference data containing pairing information
  edge_diff_list[[i]] <- list(
    pair_id = i,
    diff_matrix = edge_diff
  )
}

cat(sprintf("\n Calculated %d pairs of MR-PR adjacency matrix differences\n", length(edge_diff_list)))


# --------------------------
# Begin drawing the threshold network diagram
# --------------------------


diff_data_idx <- 1
library(tidyverse)   
library(igraph) 

plot_network <- function(edge_data, prefix, num) {
  edge_list <- data.frame(
    from = character(),
    to = character(),
    lty = integer(),
    value = numeric(),  
    threshold_result = character(),  
    stringsAsFactors = FALSE
  )
  
  # Iterate through the matrix to filter edges meeting the threshold (using the upper triangle to avoid duplicate edges in undirected graphs)
  for (i in 1:length(goal_gene)) {
    for (j in i:length(goal_gene)) {
      val <- edge_data[i, j]
      if (val > 0.5) {
        edge_list <- rbind(edge_list, data.frame(
          from = goal_gene[i], 
          to = goal_gene[j], 
          lty = 1,  # solid line
          value = val,  
          threshold_result = "1" # >0.5
        ))
      } else if (val < -0.5) {
        edge_list <- rbind(edge_list, data.frame(
          from = goal_gene[i], 
          to = goal_gene[j], 
          lty = 2,  # dotted line
          value = val,  
          threshold_result = "-1"  # <-0.5
        ))
      }
    }
  }
  
  if (nrow(edge_list) == 0) {
    warning(paste0("Edges with values > 0.5 or < -0.5 are absent; skipping plotting and saving for the ", num, "th network"))
    return(NULL)
  }
  
  # Construct an undirected network object (using goal_gene as nodes)
  G_graph <- graph_from_data_frame(
    d = edge_list,
    directed = FALSE,
    vertices = data.frame(gene = goal_gene)
  )
  
  # Set circular layout
  plot_layout <- layout.circle(G_graph)
  
  # Define output file paths
  pdf_file <- paste0("output/figure/",paste0(prefix, "_edge", num, ".pdf"))
  csv_file <- paste0("output/",paste0(prefix, "_edge", num, ".csv"))  # CSV file path
  
  # Save as PDF and plot (plotting parameters remain unchanged)
  pdf(file = pdf_file, width = 14, height = 12)
  plot(
    G_graph,
    layout = plot_layout,
    vertex.label = V(G_graph)$gene,    # Node labels (gene names)
    vertex.label.color = "#EC7357",    # Node label color (orange)
    vertex.label.cex = 1.1,            # Label size
    vertex.size = 8,                   # Node size
    vertex.color = "white",            # Node fill color
    vertex.frame.color = "#030303",    # Node border color (black)
    edge.lty = E(G_graph)$lty,         # Line type: 1=solid (>0.5), 2=dotted (<-0.5)
    edge.color = "#030303",            # Edges uniformly black
    edge.width = 2,                    # Edge width
    margin = c(1, 1, 1, 1)             # Margins
  )
  dev.off()
  
  # Save edge data as CSV file (including the new threshold_result column)
  write.csv(edge_list, file = csv_file, row.names = FALSE)
  
  cat(paste0("The ", num, "th network has been saved to: ", pdf_file, "\n"))
  cat(paste0("Edge data for the ", num, "th network has been saved to: ", csv_file, "\n"))
}


# Output file prefix (to distinguish different analyses)
output_prefix <- "Specific"

# Loop from 1 to K to plot and save data
for (num in 1:K) {
  # Extract the difference matrix for the num-th class
  current_diff_matrix <- edge_diff_list[[num]]$diff_matrix
  
  # Call the plotting function (saves CSV simultaneously)
  plot_network(
    edge_data = current_diff_matrix,
    prefix = output_prefix,
    num = num
  )
}

plot_network_separate <- function(edge_data, prefix, num) {
  # Initialize total edge list (retain original logic)
  edge_list <- data.frame(
    from = character(),
    to = character(),
    lty = integer(),
    value = numeric(),
    threshold_result = character(),
    stringsAsFactors = FALSE
  )
  
  # Iterate through the matrix to filter edges meeting the threshold (upper triangle to avoid duplicates)
  for (i in 1:length(goal_gene)) {
    for (j in i:length(goal_gene)) {
      val <- edge_data[i, j]
      if (val > 0.5) {
        edge_list <- rbind(edge_list, data.frame(
          from = goal_gene[i], 
          to = goal_gene[j], 
          lty = 1,  # solid line (positive correlation)
          value = val,
          threshold_result = "1"
        ))
      } else if (val < -0.5) {
        edge_list <- rbind(edge_list, data.frame(
          from = goal_gene[i], 
          to = goal_gene[j], 
          lty = 2,  # dashed line (negative correlation)
          value = val,
          threshold_result = "-1"
        ))
      }
    }
  }
  
  # Separate solid edges and dashed edges
  solid_edges <- edge_list[edge_list$lty == 1, ]  # solid lines (positive correlation)
  dashed_edges <- edge_list[edge_list$lty == 2, ]  # dashed lines (negative correlation)
  
  # Define plotting and saving function (internal subfunction to avoid code duplication)
  plot_and_save <- function(edges, type_suffix, line_type) {
    if (nrow(edges) == 0) {
      warning(paste0("No ", type_suffix, " edges (outside threshold); skipping plotting for the ", num, "th ", type_suffix, " network"))
      return(NULL)
    }
    
    # Construct network object
    G_graph <- graph_from_data_frame(
      d = edges,
      directed = FALSE,
      vertices = data.frame(gene = goal_gene)
    )
    
    # Circular layout
    plot_layout <- layout.circle(G_graph)
    
    # File names (with suffixes for distinction)
    pdf_file <- paste0("output/figure/",paste0(prefix, "_", type_suffix, "_edge", num, ".pdf"))
    csv_file <- paste0("output/",paste0(prefix, "_", type_suffix, "_edge", num, ".csv"))
    
    # Plot and save as PDF
    pdf(file = pdf_file, width = 14, height = 12)
    plot(
      G_graph,
      layout = plot_layout,
      vertex.label = V(G_graph)$gene,
      vertex.label.color = "#EC7357",
      vertex.label.cex = 1.1,
      vertex.size = 8,
      vertex.color = "white",
      vertex.frame.color = "#030303",
      edge.lty = line_type,  # Fixed line type for current edge category
      edge.color = "#030303",
      edge.width = 2,
      margin = c(1, 1, 1, 1)
    )
    dev.off()
    
    # Save as CSV
    write.csv(edges, file = csv_file, row.names = FALSE)
    
    cat(paste0("The ", num, "th ", type_suffix, " network has been saved to: ", pdf_file, "\n"))
    cat(paste0("Edge data for the ", num, "th ", type_suffix, " network has been saved to: ", csv_file, "\n"))
  }
  
  # Plot and save solid and dashed edge networks separately
  plot_and_save(solid_edges, type_suffix = "solid", line_type = 1)  # solid lines (positive correlation)
  plot_and_save(dashed_edges, type_suffix = "dashed", line_type = 2)  # dashed lines (negative correlation)
}


# Output file prefix
output_prefix <- "Specific"

# Loop to call the new function, outputting solid and dashed edge networks separately
for (num in 1:K) {
  current_diff_matrix <- edge_diff_list[[num]]$diff_matrix
  plot_network_separate(
    edge_data = current_diff_matrix,
    prefix = output_prefix,
    num = num
  )
}
