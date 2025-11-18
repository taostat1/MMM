rm(list=ls())  # Clear current environment variables
setwd("D:/MMM-Bio")# Set working directory
K <- 4 # Change according to needs
library("Matrix")
library("rTRNG")
library(igraph)
Rcpp::sourceCpp("D:/MMM-Bio/function_real.cpp")
# Load core data 
load("goal_gene_mt.RData")
# Goal genes removed, only non-goal genes remain
rows_to_delete <- rownames(mt) %in% goal_gene
pro_MR_mt<-t(mt[!rows_to_delete,cols_to_MR])
pro_PR_mt<-t(mt[!rows_to_delete,cols_to_PR])
pro_cell_label<-cell$Labels[cols_to_MR] # Types of cells
pro_cell_complete<-cell$scater_qc.all.total_features_by_counts[c(cols_to_MR,cols_to_PR)]# Cellular transcriptome integrity
pro_cell_expression<-cell$scater_qc.all.total_counts[c(cols_to_MR,cols_to_PR)]# Total cellular expression levels
pro_cell_reliability<-cell$scater_qc.endogenous.pct_counts[c(cols_to_MR,cols_to_PR)]# Cellular autogenous reliability
pro_cell_apoptosis<-cell$scater_qc.feature_control_MT.pct_counts[c(cols_to_MR,cols_to_PR)]# Cell Apoptosis / Stress
pro_cell_homogeneity<-cell$scater_qc.all.pct_counts_in_top_50_features[c(cols_to_MR,cols_to_PR)]# Cellular expression homogeneity
pro_cell_individual<-cell$individual[c(cols_to_MR,cols_to_PR)] # Individual patients from whom the cells were derived
table(pro_cell_label)

goal_gene_expr <- mt[goal_gene, ]  # Expression matrix of target genes (1 row x n columns, n = total number of cells)
all_cells <- c(cols_to_MR, cols_to_PR)  # Arrange all cells in order MR → PR
goal_gene_expr <- t(goal_gene_expr[, all_cells, drop = FALSE])  # After conversion: 10957 rows x 27 columns
# Convert to a numeric matrix
goal_gene_expr <- as.numeric(goal_gene_expr)  
goal_gene_expr <- matrix(goal_gene_expr, nrow = length(all_cells), ncol = length(goal_gene))  # 恢复矩阵结构
colnames(goal_gene_expr) <- goal_gene 
rownames(goal_gene_expr) <- all_cells  

# 1. Constructing the propensity matching score dataset
# Add 0 and 1 to discriminate between processing and control groups
# Create a column vector with all 1
new_column_MR <- rep(1, nrow(pro_MR_mt))
# Add the new column vector to the matrix
pro_MR_mt <- cbind(pro_MR_mt, new_column_MR)
# Set the column name of the new column to "control"
colnames(pro_MR_mt)[ncol(pro_MR_mt)] <- "control"

# Create a column vector with all 0
new_column_PR <- rep(0, nrow(pro_PR_mt))
pro_PR_mt <- cbind(pro_PR_mt, new_column_PR)
colnames(pro_PR_mt)[ncol(pro_PR_mt)] <- "control"

# Combine the two matrices
pro_mt<-rbind(pro_MR_mt,pro_PR_mt)

# 2. Make a match score
# Convert Matrix objects to dataframes
pro_mt_data <- as.data.frame(as.matrix(pro_mt),row.names = T)

control_vec <- pro_mt[, "control"]
ps_data <- data.frame(
  control = as.numeric(control_vec),  # Processing group labels (1=MR, 0=PR)
  complete = as.numeric(pro_cell_complete),  # Cellular transcriptome integrity
  expression = as.numeric(pro_cell_expression),  # Total cellular expression levels
  reliability = as.numeric(pro_cell_reliability),  # Cellular autogenous reliability
  apoptosis = as.numeric(pro_cell_apoptosis),  # Cell Apoptosis / Stress
  homogeneity = as.numeric(pro_cell_homogeneity),  # Cellular expression homogeneity
  individual = as.factor(pro_cell_individual) # Individual patients from whom the cells were derived
)

# Calculation of the propensity score
ps_model <- glm(control ~ ., 
                data = ps_data, 
                family = binomial(link = "logit"))
fitted <- ps_model$fitted.values 

# Load grouped data
load("MRData_result_with_protein.RData")
MR_group<-group_est

# Construct the score dataset
pro_group<-cbind(MR_group,fitted[cols_to_MR])
cell_group<- cbind(MR_group,pro_cell_label)
head(cell_group)
class_count <- table(pro_group[, 1])
print(class_count)

library(dplyr)
mr_fitted_data <- pro_group[pro_group[, 1] %in% 1:K, ]  # 仅保留MR细胞（1~4类）

# 2. 按MR类别分组，计算每个类别的fitted均值
mr_fitted_stats <- mr_fitted_data %>%
  as.data.frame() %>%
  group_by(group = V1) %>%  # V1是pro_group的第1列（分组标签）
  summarise(
    fitted_mean = mean(fitted),  # 该MR类别的fitted均值
    dist_to_0.5 = abs(fitted_mean - 0.5),  # 与0.5的绝对距离
    .groups = "drop"
  ) %>%
  arrange(desc(dist_to_0.5))  # 按距离从大到小排序（优先处理难匹配的）
# 若想按距离从小到大排序，改为：arrange(dist_to_0.5)

# 查看计算结果
print("每个MR类别的fitted均值与0.5的距离（按距离排序）：")
print(mr_fitted_stats)
sorted_mr_classes <- as.character(mr_fitted_stats$group)  # 转换为字符型

library(MatchIt)
library(cobalt)

all_group <- tissue_group

mr_classes <- 1:K  # K=4，MR的4个类别标签
pr_pool_all <- which(all_group %in% (K+1):(2*K))  # PR细胞的索引（初始为所有PR细胞）
pr_used <- c()  # 记录已匹配的PR细胞索引

# 3准备存储每个MR类别的匹配结果
match_results <- list()  # 存储每个MR类别的匹配对象
pr_matched_list <- list()  # 存储每个MR类别匹配到的PR细胞索引


# 循环处理每个MR类别
for (class in sorted_mr_classes) {
  cat("开始匹配MR类别", class, "...\n")
  
  # 1. 提取当前MR类别的细胞和未被使用的PR细胞
  mr_idx <- which(all_group == class)  # 当前MR类别的细胞索引
  pr_available_idx <- setdiff(pr_pool_all, pr_used)  # 未被匹配的PR细胞索引
  
  # 检查PR可用细胞数量是否充足（至少与当前MR类别细胞数相等）
  if (length(pr_available_idx) < length(mr_idx)) {
    warning("MR类别", class, "的可用PR细胞不足，可能无法完成1:1匹配！")
  }
  
  # 2. 构建当前匹配的数据集（仅包含当前MR类别和可用PR细胞）
  # 提取对应细胞的协变量（control=1为MR，0为PR）
  current_data <- ps_data[c(mr_idx, pr_available_idx), ]
  current_data$control <- ifelse(seq_along(c(mr_idx, pr_available_idx)) <= length(mr_idx), 1, 0)
  
  # 3. 计算倾向得分
  ps_model_current <- glm(control ~ complete + expression + reliability + 
                            apoptosis + homogeneity + individual,
                          data = current_data, 
                          family = binomial(link = "logit"))
  current_data$fitted <- ps_model_current$fitted.values
  
  # 4. 1:1卡尺匹配（卡尺=0.2，确保平衡）
  # 注意：需要将数据转换为matchit可识别的格式（处理变量为二元）
  match_obj <- matchit(
    control ~ complete + expression + reliability + apoptosis + homogeneity + individual,
    data = current_data,
    method = "nearest",  # 最近邻匹配
    ratio = 1,           # 1:1匹配
    caliper = 0.2,       # 卡尺阈值（SMD<0.2）
    replace = FALSE      # 不重复使用PR细胞
  )
  
  # 5. 提取匹配结果，记录匹配的PR细胞索引
  pr_abs_idx <- which(current_data$control == 0 & match_obj$weights > 0)
  
  # 转换为pr_available_idx中的相对索引（减去MR细胞的数量）
  # 因为current_data中前length(mr_idx)行是MR，之后才是PR
  mr_count <- length(mr_idx)  # MR细胞数量（如973）
  matched_pr_current <- pr_abs_idx - mr_count  # 相对索引（1~5019）
  
  # 提取原始PR细胞索引（此时无越界，避免NA）
  matched_pr_original <- pr_available_idx[matched_pr_current]
  
  # 过滤可能的NA（理论上此时已无NA，保险起见保留）
  matched_pr_original_clean <- na.omit(matched_pr_original)
  cat("MR类别", class, "：匹配到", length(matched_pr_original), "个PR细胞，过滤NA后剩余", length(matched_pr_original_clean), "个\n")
  # 6. 更新已使用的PR细胞池（用过滤后的索引）
  pr_used <- c(pr_used, matched_pr_original_clean)
  
  # 7. 验证匹配平衡
  bal <- bal.tab(match_obj,m.threshold = 0.2, var.names = c("complete", "expression", "reliability", 
                                          "apoptosis", "homogeneity", "individual"))
  cat("MR类别", class, "匹配后的SMD:\n")
  print(bal$Balance)
  
  # 8. 保存结果（用过滤后的索引）
  match_results[[as.character(class)]] <- match_obj
  pr_matched_list[[as.character(class)]] <- matched_pr_original_clean  # 存入清理后的结果
}

group_pairs <- list()
for (class in mr_classes) {
  mr_idx <- which(all_group == class)  # 当前MR类别的细胞索引
  pr_idx <- pr_matched_list[[as.character(class)]]  # 匹配的PR细胞索引
  group_pairs[[as.character(class)]] <- list(MR = mr_idx, PR = pr_idx)
  cat("MR类别", class, "匹配成功：", length(mr_idx), "个MR细胞 →", length(pr_idx), "个PR细胞\n")
}

matched_goal_expr_PR <- list()
for (class in mr_classes) {
  idx <- group_pairs[[class]]$PR# MR对应的PR
  matched_goal_expr_PR[[as.character(class)]] <- goal_gene_expr[idx, , drop = FALSE]
}

#########################构造PR网络######################################################
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
  cat(paste0("\n===== 开始处理 PR专属群", pr_class, " =====\n"))
  
  # --------------------------
  # 步骤1：提取当前PR群的表达数据并预处理
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


#############################################################################################
# Load the dplyr library
library(dplyr)
cell_group_df <- as.data.frame(cell_group, stringsAsFactors = FALSE)
# 命名列（第一列是label，第二列是细胞类型）
colnames(cell_group_df) <- c("type_label", "pro_cell_label")
# 去除可能的引号（若存在）
cell_group_df <- cell_group_df %>%
  mutate(
    type_label = gsub("\"", "", type_label),  # 清理label中的引号
    pro_cell_label = gsub("\"", "", pro_cell_label)  # 清理细胞类型中的引号
  )

# 2. 按细胞类型分组，计算每种label的数量和比例
label_proportion <- cell_group_df %>%
  group_by(pro_cell_label, type_label) %>%  # 先按细胞类型，再按label分组
  summarise(
    数量 = n(),  # 每组内的数量
    .groups = "drop"  # 取消分组状态
  ) %>%
  group_by(pro_cell_label) %>%  # 按细胞类型单独分组，计算总数量
  mutate(
    总细胞数 = sum(数量),  # 该细胞类型的总细胞数
    比例 = round(数量 / 总细胞数 , 4)  # 计算比例（保留2位小数）
  ) %>%
  arrange(pro_cell_label, desc(比例))  # 按细胞类型排序，同类型内按比例降序

# 3. 输出结果
print("每类细胞中各种label的比例统计：")
print(label_proportion, n = Inf)  # n=Inf确保显示所有结果
# write.csv(label_proportion,"label_proportion.csv")
