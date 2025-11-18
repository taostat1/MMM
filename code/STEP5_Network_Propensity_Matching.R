rm(list=ls())  # Clear current environment variables
K <- 4 # Change according to needs
library("Matrix")
# Load core data 
load("output/goal_gene_mt.RData")
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
goal_gene_expr <- matrix(goal_gene_expr, nrow = length(all_cells), ncol = length(goal_gene))  
colnames(goal_gene_expr) <- goal_gene 
rownames(goal_gene_expr) <- all_cells  

# --------------------------
# 1. Constructing the propensity matching score dataset
# --------------------------
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

# --------------------------
# 2. Make a match score
# --------------------------
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
fitted <- ps_model$fitted.values[ps_model$model$control=="1"] # MR_fitted
# Load grouped data
load("output/MRData_result_with_protein.RData")
MR_group<-group_est

# Construct the score dataset
pro_group<-cbind(MR_group,fitted)
cell_group<- cbind(MR_group,pro_cell_label)
head(cell_group)
class_count <- table(pro_group[, 1])
print(class_count)


# --------------------------
# 3. Assessing Matching Difficulty
# --------------------------
# Group by MR category and calculate the fitted mean for each category
library(dplyr)
mr_fitted_stats <- pro_group %>%
  as.data.frame() %>%  
  group_by(group = V1) %>%  
  summarise(
    fitted_mean = mean(fitted, na.rm = TRUE),  
    dist_to_MS = abs(fitted_mean -mean(ps_model$fitted.values)),  # The absolute distance between the fitted mean of this MR category and Mean Score of Total Sample
    .groups = "drop"
  ) %>%
  arrange(desc(dist_to_MS)) 
# Sort by distance in descending order (prioritize difficult-to-match items)

print("Distance between the fitted mean of each MR category and Mean Score of Total Sample (sorted by distance):")
print(mr_fitted_stats)
sorted_mr_classes <- as.character(mr_fitted_stats$group)  

# --------------------------
# 4. Match corresponding cells
# --------------------------

library(MatchIt)
library(cobalt)

mr_classes <- 1:K  
pr_pool_all <- which(ps_model$model$control=="0")  # Index of PR Cells
pr_used <- c()  # Record matched PR cell indices

# Record matched PR cell indices
match_results <- list()  
pr_matched_list <- list()  


# Iterate through each MR category
for (class in sorted_mr_classes) {
  cat("Start Matching", class, "...\n")
  
  # 1. Extract cells of the current MR class and unused PR cells
  mr_idx <- which(MR_group == class)  # Indices of cells in the current MR class
  pr_available_idx <- setdiff(pr_pool_all, pr_used)  # Indices of unmatched PR cells
  
  # Check if the number of available PR cells is sufficient
  if (length(pr_available_idx) < length(mr_idx)) {
    warning("Insufficient available PR cells for MR class ", class, " - 1:1 matching may not be completed!")
  }
  
  # 2. Construct dataset for current matching (only includes current MR class and available PR cells)
  # Extract covariates for corresponding cells (control=1 for MR, 0 for PR)
  current_data <- ps_data[c(mr_idx, pr_available_idx), ]
  current_data$control <- ifelse(seq_along(c(mr_idx, pr_available_idx)) <= length(mr_idx), 1, 0)
  
  # 3. Calculate propensity scores
  ps_model_current <- glm(control ~ complete + expression + reliability + 
                            apoptosis + homogeneity + individual,
                          data = current_data, 
                          family = binomial(link = "logit"))
  current_data$fitted <- ps_model_current$fitted.values
  
  # 4. 1:1 caliper matching (caliper=0.2 to ensure balance)
  # Note: Data needs to be converted to a format recognizable by matchit (treatment variable is binary)
  match_obj <- matchit(
    control ~ complete + expression + reliability + apoptosis + homogeneity + individual,
    data = current_data,
    method = "nearest",  # Nearest neighbor matching
    ratio = 1,           # 1:1 matching
    caliper = 0.2,       # Caliper threshold (SMD < 0.2)
    replace = FALSE      # Do not reuse PR cells
  )
  
  # 5. Extract matching results and record indices of matched PR cells
  pr_abs_idx <- which(current_data$control == 0 & match_obj$weights > 0)
  
  # Convert to relative indices within pr_available_idx
  mr_count <- length(mr_idx)  # Number of MR cells
  matched_pr_current <- pr_abs_idx - mr_count  # Relative indices (1~5019)
  
  # Extract original PR cell indices
  matched_pr_original <- pr_available_idx[matched_pr_current]
  
  # Filter out potential NA values
  matched_pr_original_clean <- na.omit(matched_pr_original)
  cat("MR class", class, ": Matched", length(matched_pr_original), "PR cells,", length(matched_pr_original_clean), "remaining after NA filtering\n")
  
  # 6. Update the pool of used PR cells (using filtered indices)
  pr_used <- c(pr_used, matched_pr_original_clean)
  
  # 7. Verify matching balance
  bal <- bal.tab(match_obj, m.threshold = 0.2, var.names = c("complete", "expression", "reliability", 
                                                             "apoptosis", "homogeneity", "individual"))
  cat("SMD after matching for MR class", class, ":\n")
  print(bal$Balance)
  
  # 8. Save results (using filtered indices)
  match_results[[as.character(class)]] <- match_obj
  pr_matched_list[[as.character(class)]] <- matched_pr_original_clean  # Store cleaned results
}

# Match corresponding cells
group_pairs <- list()
for (class in mr_classes) {
  mr_idx <- which(MR_group == class)  
  pr_idx <- pr_matched_list[[as.character(class)]]  
  group_pairs[[as.character(class)]] <- list(MR = mr_idx, PR = pr_idx)
  cat("MR Class", class, "Sucessfuily mactch：", length(mr_idx), "MR cells →", length(pr_idx), "PR cells\n")
}

matched_goal_expr_PR <- list()
for (class in mr_classes) {
  idx <- group_pairs[[class]]$PR # MR Corresponding PR
  matched_goal_expr_PR[[as.character(class)]] <- goal_gene_expr[idx, , drop = FALSE]
}

# Save results 
result_file <- paste0("output/", "Matched_goal_expr_PR.RData")
save("matched_goal_expr_PR", file = result_file)




