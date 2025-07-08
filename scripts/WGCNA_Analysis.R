# Setup and Data Preparation
# Set working directory (change as appropriate)
setwd("your/working/directory")

# Load required packages
library(WGCNA)
library(rtracklayer)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(pheatmap)
library(DESeq2)  # For VST normalization

# Avoid factor conversion issues
options(stringsAsFactors = FALSE)

# Increase memory limit for large datasets
options(future.globals.maxSize = 250000 * 1024^2)

# Step 1: Prepare gene length data from GFF annotation
my_tags <- c("Name", "Dbxref", "gene")
my_columns <- c("seqid", "start", "end", "strand", "type")
my_filter <- list(type = "gene")

dat <- readGFF("GCF_003957565.2_bTaeGut1.4.pri_genomic.gff", 
               tags = my_tags, 
               columns = my_columns, 
               filter = my_filter)

dat$gene_name <- dat$gene
dat$interval_start <- dat$start
dat$interval_stop <- dat$end

dat_1 <- as.data.frame(dat)
all_data_from_dat <- dplyr::select(dat_1, gene_name, start, end)

# Calculate gene length as end - start - 1
all_data_from_dat$geneLength <- (all_data_from_dat$end - all_data_from_dat$start) - 1

gene_name_length <- dplyr::select(all_data_from_dat, gene_name, geneLength)

# Step 2: Read and process raw gene counts data
data0_csv <- read.csv("WGCNA_gene_count_7_June_2023.csv", header = TRUE)

# Rename first column for merging
colnames(data0_csv)[1] <- "gene_name"

# Merge gene length data with counts (geneLength retained for reference)
data0 <- merge(data0_csv, gene_name_length, by = "gene_name", all.x = TRUE)

# Set rownames as gene names and remove gene_name column
colnames(data0)[1] <- ""
rownames(data0) <- data0[, 1]
data0 <- data0[, -1]

# Convert numeric columns to integer counts (raw counts)
data0[] <- lapply(data0, function(x) if(is.numeric(x)) as.integer(x) else x)

# Step 3: Load sample metadata
sample_metadata <- read.csv(file = "WGCNA_Heat_Call_Experimental_Design.csv")
colnames(sample_metadata) <- c('sample_ID', 'Sex', 'Condition')
rownames(sample_metadata) <- sample_metadata$sample_ID

# Normalization and Quality Control
# Step 4: Normalize raw counts using DESeq2 VST (Variance Stabilizing Transformation)

# Remove geneLength column before normalization
counts_for_norm <- data0[, !colnames(data0) %in% "geneLength"]

# Ensure sample names match exactly between counts and metadata
stopifnot(setequal(colnames(counts_for_norm), rownames(sample_metadata)))

# Create DESeq2 dataset (design ~1 for normalization only)
dds <- DESeqDataSetFromMatrix(countData = counts_for_norm,
                              colData = sample_metadata,
                              design = ~ 1)

# Apply variance stabilizing transformation (blind = TRUE)
vsd <- vst(dds, blind = TRUE)

# Extract normalized expression matrix and transpose for WGCNA
datExpr <- assay(vsd)
datExpr <- t(datExpr)

# Step 5: Sample Quality Control with Outlier Detection
sampleTree <- hclust(dist(datExpr), method = "average")

pdf("sample_clustering_outliers.pdf", width = 12, height = 6)
plot(sampleTree, main = "Sample Clustering for Outlier Detection")
dev.off()

# Optional: Remove outliers identified visually
# outlier_samples <- c("Sample1", "Sample2")
# datExpr <- datExpr[!rownames(datExpr) %in% outlier_samples, ]

# Network Construction and Module Detection
# Step 6: Soft Threshold Diagnostics for Network Construction
powers <- c(1:20, seq(22, 30, 2))

# Disable parallel processing to avoid forking issues
disableWGCNAThreads()

sft <- pickSoftThreshold(
  datExpr,
  powerVector = powers,
  verbose = 5,
  networkType = "signed",
  blockSize = ncol(datExpr)
)

pdf("soft_threshold_diagnostics.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Fit (R²)",
     main = "Scale Independence", type = "n")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, col = "red")
abline(h = 0.85, col = "red", lty = 2)

plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     main = "Mean Connectivity", type = "n")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, col = "red")
dev.off()

# Select optimal power
if (is.na(sft$powerEstimate)) {
  optimal_power <- sft$fitIndices$Power[which.max(sft$fitIndices$SFT.R.sq)]
  message("No power reached R² > 0.85. Using max R² at power: ", optimal_power)
} else {
  optimal_power <- sft$powerEstimate
}

save.image("06_July_2025_Pre_Parameter.Rdata")

# Step 7: Parameter Sensitivity Analysis for Module Detection
parameter_grid <- expand.grid(
  minModuleSize = c(20, 30, 40),
  mergeCutHeight = c(0.15, 0.25, 0.35)
)

module_results <- list()

for (i in 1:nrow(parameter_grid)) {
  params <- parameter_grid[i, ]
  
  assign("cor", WGCNA::cor, envir = .GlobalEnv)
  
  net <- blockwiseModules(
    datExpr,
    power = optimal_power,
    TOMType = "signed",
    minModuleSize = params$minModuleSize,
    mergeCutHeight = params$mergeCutHeight,
    numericLabels = TRUE,
    saveTOMs = FALSE,
    verbose = 3
  )
  
  assign("cor", stats::cor, envir = .GlobalEnv)
  
  module_results[[i]] <- list(
    params = params,
    n_modules = length(unique(net$colors)),
    module_sizes = table(net$colors)
  )
  
  mergedColors <- labels2colors(net$colors)
  pdf(paste0("module_dendrogram_minSize", params$minModuleSize, 
             "_mergeCut", params$mergeCutHeight, ".pdf"), 
      width = 10, height = 6)
  plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                      "Module colors", dendroLabels = FALSE)
  dev.off()
}

parameter_summary <- do.call(rbind, lapply(module_results, function(x) {
  data.frame(
    minModuleSize = x$params$minModuleSize,
    mergeCutHeight = x$params$mergeCutHeight,
    n_modules = x$n_modules,
    min_module_size = min(x$module_sizes),
    max_module_size = max(x$module_sizes)
  )
}))

print(parameter_summary)

final_params <- list(
  minModuleSize = 30,
  mergeCutHeight = 0.25
)

save.image("06_July_2025_pre_network_WGCNA.RData")

# Step 8: Final Network Construction
enableWGCNAThreads(nThreads = 20)

assign("cor", WGCNA::cor, envir = .GlobalEnv)

net <- blockwiseModules(
  datExpr,
  power = optimal_power,
  TOMType = "signed",
  minModuleSize = final_params$minModuleSize,
  mergeCutHeight = final_params$mergeCutHeight,
  numericLabels = TRUE,
  saveTOMs = TRUE,
  saveTOMFileBase = "network_TOM",
  verbose = 3,
  maxBlockSize = 10000
)

assign("cor", stats::cor, envir = .GlobalEnv)

moduleColors <- labels2colors(net$colors)
write.csv(data.frame(Gene = colnames(datExpr), Module = moduleColors),
          "final_module_assignments.csv", row.names = FALSE)

pdf("final_module_dendrogram.pdf", width = 15, height = 10)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors", dendroLabels = FALSE)
dev.off()

save.image("06_July_2025_post_network_WGCNA.RData")

#Module-Trait Relationships
# Step 9: Module-Trait Relationships
traits <- sample_metadata[rownames(datExpr), c("Sex", "Condition")]
traits$Sex <- as.numeric(factor(traits$Sex))
traits$Condition <- as.numeric(factor(traits$Condition))

MEs <- moduleEigengenes(datExpr, moduleColors)$eigengenes

moduleTraitCor <- cor(MEs, traits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

pdf("module_trait_relationships.pdf", width = 10, height = 8)
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(traits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1, 1),
               main = "Module-Trait Relationships")
dev.off()

# Gene Ontology Enrichment
library(org.Gg.eg.db, lib.loc = "/project/cugbf/software/R/4.4.1/library")

module_genes <- split(colnames(datExpr), moduleColors)

perform_GO_analysis <- function(module_name, genes) {
  if (length(genes) < 10) {
    message("Skipping GO for ", module_name, " (only ", length(genes), " genes)")
    return(NULL)
  }
  
  genes_entrez <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", 
                       OrgDb = org.Gg.eg.db, drop = TRUE)
  
  if (nrow(genes_entrez) < 5) {
    message("Skipping GO for ", module_name, " (only ", nrow(genes_entrez), " ENTREZ IDs)")
    return(NULL)
  }
  
  enrichGO(
    gene = genes_entrez$ENTREZID,
    OrgDb = org.Gg.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )
}

go_results <- lapply(names(module_genes), function(module) {
  result <- tryCatch({
    perform_GO_analysis(module, module_genes[[module]])
  }, error = function(e) {
    message("GO failed for ", module, ": ", conditionMessage(e))
    return(NULL)
  })
  
  if (!is.null(result) && nrow(result) > 0) {
    pdf(paste0("GO_enrichment_", module, ".pdf"), width = 10, height = 8)
    print(barplot(result, showCategory = 20, title = paste("GO Enrichment -", module)))
    dev.off()
  }
  
  return(result)
})
names(go_results) <- names(module_genes)

saveRDS(go_results, "06_July_2025_go_enrichment_results.rds")

# Hub Gene Identification and Validation
ADJ <- adjacency(datExpr, power = optimal_power, type = "signed")
kIN <- intramodularConnectivity(ADJ, moduleColors, scaleByMax = TRUE)

hub_genes <- lapply(unique(moduleColors), function(mod) {
  modGenes <- (moduleColors == mod)
  kIN_mod <- kIN[modGenes, ]
  head(kIN_mod[order(kIN_mod$kWithin, decreasing = TRUE), ], 10)
})
names(hub_genes) <- unique(moduleColors)

saveRDS(hub_genes, "hub_genes_per_module.rds")

pdf("hub_gene_validation.pdf", width = 12, height = 8)
for (mod in names(hub_genes)) {
  top_gene <- rownames(hub_genes[[mod]])[1]
  plot(datExpr[, top_gene], 
       MEs[, paste0("ME", mod)], 
       main = paste("Hub Gene Validation -", mod),
       xlab = top_gene,
       ylab = "Module Eigengene")
}
dev.off()

save.image("07_july_2025_WGCNA_analysis_final.RData")

# Network Visualization 
library(igraph)
library(ggraph)
library(tidygraph)
library(dplyr)
library(ggplot2)

module_of_interest <- "green"

module_genes <- colnames(datExpr)[moduleColors == module_of_interest]

ADJ <- adjacency(datExpr, power = optimal_power, type = "signed")

module_adj <- ADJ[module_genes, module_genes]

edge_threshold <- 0.1
module_adj[module_adj < edge_threshold] <- 0

graph <- graph_from_adjacency_matrix(module_adj, mode = "undirected", weighted = TRUE, diag = FALSE)

kIN_module <- kIN[moduleColors == module_of_interest, ]

top_hubs <- rownames(kIN_module)[order(kIN_module$kWithin, decreasing = TRUE)[1:10]]

V(graph)$hub <- V(graph)$name %in% top_hubs

V(graph)$color <- ifelse(V(graph)$hub, "red", "green")

V(graph)$size <- ifelse(V(graph)$hub, 8, 4)

E(graph)$width <- E(graph)$weight * 5

p <- ggraph(graph, layout = "fr") +
  geom_edge_link(aes(width = weight), alpha = 0.3, color = "grey50") +
  geom_node_point(aes(color = color, size = size)) +
  geom_node_text(aes(label = ifelse(hub, name, "")), repel = TRUE, size = 4, color = "black", fontface = "bold") +
  scale_color_identity() +
  scale_size_identity() +
  theme_void() +
  ggtitle("WGCNA Green Module Network with Hub Genes Highlighted") +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"))

print(p)

library(WGCNA)
exportNetworkToCytoscape(
  adjMat = module_adj,
  edgeFile = "CytoscapeInput-edges-green.txt",
  nodeFile = "CytoscapeInput-nodes-green.txt",
  weighted = TRUE,
  threshold = 0.2,
  nodeNames = module_genes,
  nodeAttr = moduleColors[moduleColors == module_of_interest],
  includeColNames = TRUE
)

# Less busy plot with top 5 hubs and their neighbors
top_hubs <- rownames(kIN_module)[order(kIN_module$kWithin, decreasing = TRUE)[1:5]]

edge_threshold <- 0.4

neighbors <- unique(unlist(lapply(top_hubs, function(hub) {
  hub_edges <- module_adj[hub, ]
  top_neighbors <- names(sort(hub_edges, decreasing = TRUE)[1:5])
  top_neighbors[hub_edges[top_neighbors] > edge_threshold]
})))
genes_to_plot <- unique(c(top_hubs, neighbors))

sub_adj <- module_adj[genes_to_plot, genes_to_plot]
sub_graph <- graph_from_adjacency_matrix(sub_adj, mode = "undirected", weighted = TRUE, diag = FALSE)

V(sub_graph)$hub <- V(sub_graph)$name %in% top_hubs
V(sub_graph)$color <- ifelse(V(sub_graph)$hub, "red", "green")
V(sub_graph)$size <- ifelse(V(sub_graph)$hub, 8, 4)

p_sub <- ggraph(sub_graph, layout = "fr", niter = 4000, area = 40000, repulserad = 20000) +
  geom_edge_link(aes(width = weight), alpha = 0.3, color = "grey50") +
  geom_node_point(aes(color = color, size = size)) +
  geom_node_text(aes(label = ifelse(hub, name, "")), repel = TRUE, size = 4, color = "black", fontface = "bold", max.overlaps = Inf) +
  scale_color_identity() +
  scale_size_identity() +
  theme_void() +
  ggtitle("Green Module: Hub Genes and Their Top Neighbors") +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"))

print(p_sub)

hub_color <- "#0072B2"
neighbor_color <- "#D55E00"

V(sub_graph)$color <- ifelse(V(sub_graph)$hub, hub_color, neighbor_color)
V(sub_graph)$size <- ifelse(V(sub_graph)$hub, 10, 5)
V(sub_graph)$frame.color <- "black"

p_kk <- ggraph(sub_graph, layout = "kk") +
  geom_edge_link(aes(width = weight), alpha = 0.25, color = "grey70") +
  geom_node_point(aes(color = color, size = size), shape = 21, stroke = 1.2) +
  geom_node_text(aes(label = ifelse(hub, name, "")), repel = TRUE, size = 4, color = "black", fontface = "bold", bg.color = "white", bg.r = 0.15, max.overlaps = Inf) +
  scale_color_identity() +
  scale_size_identity() +
  theme_void() +
  ggtitle("Green Module: Hub Genes and Their Top Neighbors (Kamada-Kawai)") +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"), legend.position = "none")

ggsave("green_module_hub_subnetwork_nature.pdf", p_kk, width = 12, height = 12, dpi = 600)
ggsave("green_module_hub_subnetwork_nature.png", p_kk, width = 12, height = 12, dpi = 600)

print(p_kk)

# Cell-Type Enrichment Analysis
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pheatmap)
library(RColorBrewer)

# Load cell-type marker genes (converted to gene symbols)
top_markers_converted <- read.csv("hypomap_cell_markers_gene_id.csv")

# Create a named list of marker genes per cell type
cell_markers_list <- split(top_markers_converted$gene, top_markers_converted$cluster)

# Load module assignments
module_gene_df <- read.csv("final_module_assignments.csv")  # columns: Gene, Module

# Create a named list of genes per module
module_genes <- split(module_gene_df$Gene, module_gene_df$Module)

# Define the background gene universe (all genes assigned to modules)
all_genes <- unique(module_gene_df$Gene)

# Perform hypergeometric enrichment test with FDR correction
enrichment_pvals <- matrix(NA, nrow = length(module_genes), ncol = length(cell_markers_list))
rownames(enrichment_pvals) <- names(module_genes)
colnames(enrichment_pvals) <- names(cell_markers_list)

for (mod in names(module_genes)) {
  for (ct in names(cell_markers_list)) {
    overlap <- intersect(module_genes[[mod]], cell_markers_list[[ct]])
    enrichment_pvals[mod, ct] <- phyper(
      length(overlap) - 1,
      length(cell_markers_list[[ct]]),
      length(all_genes) - length(cell_markers_list[[ct]]),
      length(module_genes[[mod]]),
      lower.tail = FALSE
    )
  }
}

# Multiple testing correction (FDR)
enrichment_fdr <- matrix(
  p.adjust(as.vector(enrichment_pvals), method = "BH"),
  nrow = nrow(enrichment_pvals),
  dimnames = dimnames(enrichment_pvals)
)

# Annotate FDR values with asterisks for significance
annotations <- matrix("", nrow = nrow(enrichment_fdr), ncol = ncol(enrichment_fdr))
annotations[enrichment_fdr < 0.001] <- "***"
annotations[enrichment_fdr < 0.01 & enrichment_fdr >= 0.001] <- "**"
annotations[enrichment_fdr < 0.05 & enrichment_fdr >= 0.01] <- "*"

display_numbers <- matrix("", nrow = nrow(enrichment_fdr), ncol = ncol(enrichment_fdr))
for (i in seq_len(nrow(enrichment_fdr))) {
  for (j in seq_len(ncol(enrichment_fdr))) {
    if (enrichment_fdr[i, j] < 0.001) {
      display_numbers[i, j] <- paste0("<0.001", annotations[i, j])
    } else {
      display_numbers[i, j] <- paste0(sprintf("%.3f", enrichment_fdr[i, j]), annotations[i, j])
    }
  }
}

# Define color palette: low FDR (strong enrichment) = dark red, high FDR = white
color_palette <- colorRampPalette(c("darkred", "white"))(100)

# Plot publication-quality heatmap
pheatmap(
  enrichment_fdr,
  color = color_palette,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  fontsize_row = 10,
  fontsize_col = 10,
  angle_col = 45,
  treeheight_row = 0,
  treeheight_col = 0,
  legend = TRUE,
  display_numbers = display_numbers,
  number_color = "black",
  main = "Module Enrichment in Cell Types (FDR)",
  legend_title = "FDR values",
  width = 12,
  height = 10,
  filename = "module_cell_type_enrichment_FDR_reversed.pdf"
)
