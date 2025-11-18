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
protein_file <- "protein_matrix.csv"  # Protein interaction matrix
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
# 2. 循环处理每个PR专属群（PR1~PR4）
# --------------------------
pr_network_results <- list()
# PR群名称（对应MR1~MR4的专属PR群）
pr_classes <- names(matched_goal_expr_PR)  # 如c("1","2","3","4")
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
  result_file <- paste0(pr_class, "Pr_Data_result_with_protein.RData")
  save(
    list = c("X", "invcov_est", "edge_est", "group_est", "lam0_est", 
             "lam1_est", "mu_est", "Pi_est", "theta_est", "protein_mat"),
    file = result_file
  )
  cat(paste0(pr_class, " results saved to: ", result_file, "\n"))
  
  
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
    pdf_file <- paste0(pr_class, "pr_edge", num, ".pdf")
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
# 3. 保存所有PR群的网络结果（便于后续与MR对比）
# --------------------------
save(
  pr_network_results,
  file = "PR_4Groups_Network_Results.RData"
)

# Load grouped data
load("MRData_result_with_protein.RData")
mr_edge_est <-edge_est

i <- 1
edge_diff_list <- list()
for (i in pr_classes) {
  mr_edge <- mr_edge_est[,, as.numeric(i)]             # MR对应子类的邻接矩阵
  pr_edge <- pr_network_results[[i]][,,1]  # PR对应子类的邻接矩阵
  
  # 计算差值（MR - PR）
  edge_diff <- mr_edge - pr_edge
  
  # 存储包含配对信息的差值数据
  edge_diff_list[[i]] <- list(
    pair_id = i,
    diff_matrix = edge_diff
  )
}

cat(sprintf("\n✅ 共计算%d对MR-PR邻接矩阵差值\n", length(edge_diff_list)))


# --------------------------
# 批量绘制所有配对的阈值化网络
# --------------------------
cat("\n=== 开始绘制阈值化网络图 ===\n")

diff_data_idx<- 1
library(tidyverse)   # 如果你用了 dplyr 的 %>% 等函数
library(igraph) 
plot_network <- function(edge_data, prefix, num) {
  # 初始化边列表（新增threshold_result列记录阈值处理结果）
  edge_list <- data.frame(
    from = character(),
    to = character(),
    lty = integer(),
    value = numeric(),  # 原始值
    threshold_result = character(),  # 阈值处理结果（新增）
    stringsAsFactors = FALSE
  )
  
  # 遍历矩阵筛选符合阈值的边（上三角避免无向图重复边）
  for (i in 1:length(goal_gene)) {
    for (j in i:length(goal_gene)) {
      val <- edge_data[i, j]
      if (val > 0.5) {
        edge_list <- rbind(edge_list, data.frame(
          from = goal_gene[i], 
          to = goal_gene[j], 
          lty = 1,  # 实线
          value = val,  # 原始值
          threshold_result = "1"  # 阈值处理结果：正相关且超过0.5
        ))
      } else if (val < -0.5) {
        edge_list <- rbind(edge_list, data.frame(
          from = goal_gene[i], 
          to = goal_gene[j], 
          lty = 2,  # 虚线
          value = val,  # 原始值
          threshold_result = "-1"  # 阈值处理结果：负相关且小于-0.5
        ))
      }
    }
  }
  
  # 容错：无符合条件的边则跳过
  if (nrow(edge_list) == 0) {
    warning(paste0("无差值>0.5或<-0.5的边，跳过第", num, "个网络绘制和保存"))
    return(NULL)
  }
  
  # 构建无向网络对象（使用goal_gene作为节点）
  G_graph <- graph_from_data_frame(
    d = edge_list,
    directed = FALSE,
    vertices = data.frame(gene = goal_gene)
  )
  
  # 设置圆形布局
  plot_layout <- layout.circle(G_graph)
  
  # 定义输出文件路径
  pdf_file <- paste0(prefix, "_edge", num, ".pdf")
  csv_file <- paste0(prefix, "_edge", num, ".csv")  # CSV文件路径
  
  # 保存为PDF并绘图（绘图部分不变）
  pdf(file = pdf_file, width = 14, height = 12)
  plot(
    G_graph,
    layout = plot_layout,
    vertex.label = V(G_graph)$gene,    # 节点label（基因名）
    vertex.label.color = "#EC7357",    # 节点label颜色（橙色）
    vertex.label.cex = 1.1,            # label大小
    vertex.size = 8,                   # 节点大小
    vertex.color = "white",            # 节点填充色
    vertex.frame.color = "#030303",    # 节点边框色（黑色）
    edge.lty = E(G_graph)$lty,         # 线型：1=实线（>0.5），2=虚线（<-0.5）
    edge.color = "#030303",            # 边统一为黑色
    edge.width = 2,                    # 边宽度
    margin = c(1, 1, 1, 1)             # 边距
  )
  dev.off()
  
  # 保存边数据为CSV文件（包含新增的threshold_result列）
  write.csv(edge_list, file = csv_file, row.names = FALSE)
  
  cat(paste0("第", num, "个网络已保存至: ", pdf_file, "\n"))
  cat(paste0("第", num, "个网络的边数据已保存至: ", csv_file, "\n"))
}


# 输出文件前缀（用于区分不同分析）
output_prefix <- "Specific"

# 从1到K循环绘图并保存数据
for (num in 1:K) {
  # 提取第num类的差值矩阵
  current_diff_matrix <- edge_diff_list[[num]]$diff_matrix
  
  # 调用绘图函数（同时保存CSV）
  plot_network(
    edge_data = current_diff_matrix,
    prefix = output_prefix,
    num = num
  )
}

plot_network_separate <- function(edge_data, prefix, num) {
  # 初始化总边列表（保留原始逻辑）
  edge_list <- data.frame(
    from = character(),
    to = character(),
    lty = integer(),
    value = numeric(),
    threshold_result = character(),
    stringsAsFactors = FALSE
  )
  
  # 遍历矩阵筛选符合阈值的边（上三角避免重复）
  for (i in 1:length(goal_gene)) {
    for (j in i:length(goal_gene)) {
      val <- edge_data[i, j]
      if (val > 0.5) {
        edge_list <- rbind(edge_list, data.frame(
          from = goal_gene[i], 
          to = goal_gene[j], 
          lty = 1,  # 实线（正相关）
          value = val,
          threshold_result = "1"
        ))
      } else if (val < -0.5) {
        edge_list <- rbind(edge_list, data.frame(
          from = goal_gene[i], 
          to = goal_gene[j], 
          lty = 2,  # 虚线（负相关）
          value = val,
          threshold_result = "-1"
        ))
      }
    }
  }
  
  # 分离实线边和虚线边
  solid_edges <- edge_list[edge_list$lty == 1, ]  # 实线（正相关）
  dashed_edges <- edge_list[edge_list$lty == 2, ]  # 虚线（负相关）
  
  # 定义绘图和保存函数（内部子函数，避免重复代码）
  plot_and_save <- function(edges, type_suffix, line_type) {
    if (nrow(edges) == 0) {
      warning(paste0("无", type_suffix, "边（阈值外），跳过第", num, "个", type_suffix, "网络绘制"))
      return(NULL)
    }
    
    # 构建网络对象
    G_graph <- graph_from_data_frame(
      d = edges,
      directed = FALSE,
      vertices = data.frame(gene = goal_gene)
    )
    
    # 圆形布局
    plot_layout <- layout.circle(G_graph)
    
    # 文件名（加后缀区分）
    pdf_file <- paste0(prefix, "_", type_suffix, "_edge", num, ".pdf")
    csv_file <- paste0(prefix, "_", type_suffix, "_edge", num, ".csv")
    
    # 绘图并保存PDF
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
      edge.lty = line_type,  # 固定当前类型的线型
      edge.color = "#030303",
      edge.width = 2,
      margin = c(1, 1, 1, 1)
    )
    dev.off()
    
    # 保存CSV
    write.csv(edges, file = csv_file, row.names = FALSE)
    
    cat(paste0("第", num, "个", type_suffix, "网络已保存至: ", pdf_file, "\n"))
    cat(paste0("第", num, "个", type_suffix, "边数据已保存至: ", csv_file, "\n"))
  }
  
  # 分别绘制实线和虚线网络
  plot_and_save(solid_edges, type_suffix = "solid", line_type = 1)  # 实线（正相关）
  plot_and_save(dashed_edges, type_suffix = "dashed", line_type = 2)  # 虚线（负相关）
}


# 输出文件前缀
output_prefix <- "Specific"

# 循环调用新函数，分别输出实线和虚线网络
for (num in 1:K) {
  current_diff_matrix <- edge_diff_list[[num]]$diff_matrix
  plot_network_separate(
    edge_data = current_diff_matrix,
    prefix = output_prefix,
    num = num
  )
}