rm(list=ls())  # Clear current environment variables
result_file <- "output/MRData_result_with_protein.RData"
load(result_file)
# Load required library
library(igraph)

# --------------------------
# 1. Generate and save protain network plots
# --------------------------
protein_df <- read.csv(file = "data/protein_matrix.csv", header = FALSE)
protein_mat <- as.matrix(protein_df)
# Ensure matrix is numeric (binary 0/1: 0=no interaction, 1=interaction)
protein_mat <- as.numeric(protein_mat)
dim(protein_mat) <- c(nrow(protein_df), ncol(protein_df))

# Enforce binary format for protein matrix (prevent interference from non-0/1 values)
protein_mat[protein_mat != 0] <- 1  # Treat non-0 values as "interaction" (1)
protein_mat[protein_mat == 0] <- 0  # Treat 0 values as "no interaction"
gene_names <- rownames(X)  # Get gene names from original expression matrix
rownames(protein_mat) <- gene_names
colnames(protein_mat) <- gene_names

# Define network plotting function
plot_network <- function(edge_data, ds_name, num) {
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
  pdf_file <- paste0("output/figure/",paste0(ds_name, "_edge", num, ".pdf"))
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
  cat(paste0(ds_name, " network ", num, " saved to: ", pdf_file, "\n"))
}

plot_network(
  edge_data= protein_mat,
  "Protein",
 1
)

# --------------------------
# 2. Generate and save MR cell type-specific network plots
# --------------------------
# Binarize edge matrix: Treat posterior probability ≥0.5 as "edge exists" (1), else 0
X <- as.matrix(X)  # Ensure X is in matrix format
G <- dim(X)[1]     # Number of genes (rows: genes)
N <- dim(X)[2]     # Number of cells (columns: cells)
K <- dim(edge_est)[3] # Number of kinds
edge_est_mat <- array(0, dim = c(G, G, dim(edge_est)[3]))
edge_est_mat[edge_est >= 0.5] <- 1

# Prepare edge matrices with gene names for each cell type
edge_full <- list()        # Store edge matrix for each cell type
for (k in 1:K) {
  edge_full[[k]] <- edge_est_mat[,, k]
  rownames(edge_full[[k]]) <- gene_names
  colnames(edge_full[[k]]) <- gene_names
}

# Plot and save network graphs for each MR cell type
for (k in 1:K) {
  # Adjust edge matrix order 
  plot_data <- edge_full[[k]][c(1, nrow(edge_full[[k]]):2), c(1, nrow(edge_full[[k]]):2)]
  # Call plotting function
  plot_network(edge_data = plot_data, ds_name= "MR", num = k)
}

# --------------------------
# 3. Generate and save MR heatmap
# --------------------------

# Start Timing
library(ggplot2)
library(gplots)   
library(RColorBrewer)  
start_time <- Sys.time()

# Convert expression matrix to cell×gene format
expr_matrix <- as.data.frame(t(X))  
gene_names <- rownames(X)           
cell_names <- colnames(X)           
rownames(expr_matrix) <- cell_names 
colnames(expr_matrix) <- gene_names 

# Retrieve Cell Type Grouping
cell_kinds <- as.factor(group_est)
cat("MR Data Cell Type Distribution：\n")
print(table(cell_kinds))

# Sort cells by cell type (ensure same types cluster)
sorted_idx <- order(cell_kinds)
expr_sorted <- expr_matrix[sorted_idx, ]  
cell_kinds_sorted <- cell_kinds[sorted_idx]  

# Log2 transformation of expression values
expr_log <- log2(expr_sorted + 1)

# Define color palette (blue-white-red gradient)
color_palette <- colorRampPalette(brewer.pal(9, "RdBu"))(100)
color_palette <- rev(color_palette)  # High expression: red; Low expression: blue

# Define cell type colors
type_colors <- c(
  "#FF6B6B",  # Cell Type 1: Pink
  "#4ECDC4",  # Cell Type 2: Blue-Green
  "#45B7D1",  # Cell Type 3: Sky Blue
  "#FFA07A"   # Cell Type 4: Light Orange
)
names(type_colors) <- levels(cell_kinds_sorted)
row_side_colors <- type_colors[as.character(cell_kinds_sorted)]  # Colors for row sidebar

# Plot and save heatmap
pdf("output/figure/MR_Expression_Heatmap.pdf", width = 12, height = 8)

# Generate heatmap (no cell barcodes displayed)
heatmap.2(
  as.matrix(expr_log),        
  col = color_palette,        
  scale = "row",              
  RowSideColors = row_side_colors,  
  Rowv = FALSE,                    
  Colv = TRUE,                     
  dendrogram = "column",           
  key = TRUE,                      
  key.title = "log2(Expression + 1)",  
  density.info = "none",          
  trace = "none",                 
  margins = c(8, 4), 
  xlab = "Genes",                 
  ylab = "Cells",                  
  axisRow = FALSE,    
  cexCol = 0.7,                    
  srtCol = 45,                     
  adjCol = c(1, 0.5),
  labRow = FALSE     
)

mtext("Cells", side = 2, line = 3, cex = 1)  
mtext("Genes", side = 1, line = 5, cex = 1)

# Add cell type legend
legend(
  "bottomleft",
  legend = levels(cell_kinds_sorted),
  fill = type_colors,
  border = FALSE,
  bty = "n",  
  cex = 0.8,
  title = "Cell Types"
)

dev.off() 

# End Timer and calculate runtime
end_time <- Sys.time()
execution_time <- end_time - start_time
print(execution_time)
# Time difference of 10.56918 mins





