# Figure 1c ---------------------------------------------------------------
# Adhesome DEG heatmap
library(tidyverse)
library(ComplexHeatmap)

# Import adhesome gene data
adhesome <- read.csv("~/chen-origer-xist-adhesome/data/adhesome-components.csv")

# Get DESeq2 results from PMID:39546568 Dataset S01
# Import as res_X7 and res_X9
res_X7 <- read.csv("~/chen-origer-xist-adhesome/data/deseq-res-X7.csv")
res_X9 <- read.csv("~/chen-origer-xist-adhesome/data/deseq-res-X9.csv")

res_X7_degs <- res_X7 %>%
  # Label DEGs
  dplyr::mutate(X7_deg = ifelse(abs(log2FoldChange) > 0.58 & padj < 0.05,
                yes = T, no = F)) %>%
  # Remove unnecessary columns
  dplyr::select(c("gene","log2FoldChange","X7_deg")) %>%
  # Rename in preparation for merge
  dplyr::rename(X7 = log2FoldChange)

res_X9_degs <- res_X9 %>%
  # Label DEGs
  dplyr::mutate(X9_deg = ifelse(abs(log2FoldChange) > 0.58 & padj < 0.05,
                                yes = T, no = F)) %>%
  # Remove unnecessary columns
  dplyr::select(c("gene","log2FoldChange","X9_deg")) %>%
  # Rename in preparation for merge
  dplyr::rename(X9 = log2FoldChange)

# Merge deg dataframes
xist_degs <- merge(res_X7_degs, res_X9_degs, by = "gene", all = T)

xist_degs <- xist_degs %>%
  # Filter for adhesome genes
  dplyr::filter(gene %in% adhesome$gene) %>%
  # Filter for genes that are DEGs in at least one clone
  dplyr::filter(X7_deg | X9_deg) %>%
  # Add clone DEG label for sorting
  dplyr::mutate(clone = ifelse(X7_deg & X9_deg,
                               yes = "Overlap",
                               no = ifelse(X7_deg,
                                           yes = "X7-only",
                                           no = "X9-only"))) %>%
  # Remove unnecessary columns
  dplyr::select(c("gene", "X7", "X9", "clone")) %>%
  # Sort for plotting
  dplyr::arrange(clone, desc(X7))

# Create matrix of log2FoldChange values with gene names as rownames for plotting
xist_degs_mat <- xist_degs %>% 
  dplyr::select(c("gene","X7","X9")) %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames(var = "gene") %>% 
  as.matrix()

# Create matrix of adhesome category with gene names as rownames for plotting
xist_cat_mat <- adhesome[match(xist_degs$gene, adhesome$gene),] %>% 
  dplyr::select(c("gene","cad","int")) %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames(var = "gene") %>% 
  as.matrix()

# Adjust columns and order and make numerical
colnames(xist_cat_mat) <- c("Cad.", "Int.")
xist_cat_mat <- xist_cat_mat[,c(2,1)]
xist_cat_mat <- ifelse(xist_cat_mat, yes = 1, no = 0)

ht1 <- ComplexHeatmap::Heatmap(xist_degs_mat, 
                               na_col = "grey", 
                               name = "Log\nFold\nChange", 
                               cluster_rows = F,
                               show_row_dend = F, 
                               cluster_columns = F, 
                               row_names_side = "left", 
                               row_names_rot = 0, 
                               column_names_rot = 0,
                               column_names_centered = T,
                               width = unit(0.75, "npc"))

ht2 <- ComplexHeatmap::Heatmap(xist_cat_mat, 
                               col = c("white", "black"), 
                               cluster_rows = F,
                               cluster_columns = F, 
                               row_names_side = "left", 
                               row_names_rot = 0, 
                               column_names_rot = 0,
                               column_names_centered = T,
                               show_heatmap_legend = F,
                               width = unit(0.4, "npc"))

draw(ht1+ht2)

# Figure 1d ---------------------------------------------------------------
# GSEA plots
library(ggplot2)
library(tidyverse)
library(tidytext)
library(gprofiler2)

# Get DESeq2 results from PMID:39546568 Dataset S01
# Import as res_X7 and res_X9
res_X7 <- read.csv("~/chen-origer-xist-adhesome/data/deseq-res-X7.csv")
res_X9 <- read.csv("~/chen-origer-xist-adhesome/data/deseq-res-X9.csv")

# Filter for significant adjust p value and >1.5 fold change
upDEG_X7 <- dplyr::filter(res_X7, log2FoldChange > 0.58 & padj < 0.05)
upDEG_X9 <- dplyr::filter(res_X9, log2FoldChange > 0.58 & padj < 0.05)

# Get all genes in either X7 or X9 object ("background")
#background_genes <- unique(res_X7$gene)
background_genes <- unique(res_X9$gene)

# Select the DEG list to test
genes <- upDEG_X7$gene

# Get gene list functional enrichment
gost_res <- gprofiler2::gost(
  query = genes,
  organism = "hsapiens",
  sources = c("GO:BP"),
  custom_bg = background_genes,
  correction_method = "fdr"
)

# Subset significant results
sig_terms <- gost_res$result %>%
  dplyr::filter(p_value < 0.05)

# Get the top 20 results sorted by p value
top20 <- sig_terms %>%
  dplyr::arrange(p_value) %>%
  dplyr::slice(1:20)

# Clean results for export/downstream plots
top20_clean <- top20 %>%
  dplyr::mutate(parents = sapply(parents, paste, collapse = ";")) %>%
  dplyr::select(term_id,
                term_name,
                source,
                p_value,
                intersection_size,
                term_size,
                precision,
                recall)

write.csv(top20_clean,
          "~/chen-origer-xist-adhesome/data/fig1d-GOBP-upDEG-X9.csv",
          row.names = FALSE)

# Import saved GSEA results
up_X7_GOBP <- read.csv("~/chen-origer-xist-adhesome/data/fig1d-GOBP-upDEG-X7.csv") %>%
  dplyr::mutate(category = "up_X7")
up_X9_GOBP <- read.csv("~/chen-origer-xist-adhesome/data/fig1d-GOBP-upDEG-X9.csv") %>%
  dplyr::mutate(category = "up_X9")

# Merge into single object
res_GOBP <- rbind(up_X7_GOBP, up_X9_GOBP)
rm(up_X7_GOBP, up_X9_GOBP)

# Capitalize terms
res_GOBP$term_name <- stringr::str_to_title(res_GOBP$term_name)

# Sort within categories and add -log(p_value) as a column
res_GOBP <- res_GOBP %>%
  dplyr::group_by(category) %>%
  dplyr::arrange(p_value, .by_group = T) %>%
  dplyr::mutate(neg_log_p = -log10(p_value))

# Create dot plot
ggplot(res_GOBP, 
       aes(x = recall, 
           y = tidytext::reorder_within(term_name, recall, category), 
           color = neg_log_p, 
           size = intersection_size)) +
  geom_point(alpha = 0.8) +
  scale_size_continuous(name = "Gene Count") +
  scale_color_gradient(low = "#f0907b", high = "#cb181d") +
  tidytext::scale_y_reordered() +
  theme_minimal(base_size = 14) +
  labs(title = NULL,
       x = NULL,
       y = NULL,
       color = "-Log P",
       size = "Count") +
  facet_grid(category ~ .,
             scales = "free_y")

# Supp Figure 1a/b ---------------------------------------------------------
# DEG volcano plots
library(ggplot2)
library(ggrepel)

# Import DEG data derived from DESeq2 results PMID:39546568 Dataset S01
xist_volcano <- read.csv("~/chen-origer-xist-adhesome/data/supp-fig1ab-volcano.csv")

ggplot(xist_volcano, aes(x = log2FoldChange, y = -log10(padj), col = deg, label = delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel(size = 3) +
  scale_color_manual(values = c("downregulated" = "blue", 
                                "none" = "gray", 
                                "upregulated" = "red")) +
  geom_vline(xintercept = c(-0.58, 0.58), color = "black", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), color = "black", linetype = "dashed") +
  theme(panel.spacing = unit(0.1, "lines"),
        plot.title = element_text(size = 7, hjust = 0.5, face = "plain"),
        plot.subtitle = element_text(size = 7, hjust = 0.5, face = "plain"),
        legend.position = "none",
        axis.title = element_text(size = 7),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.text = element_text(size = 7, color = "black"),
        axis.line = element_line(color = 'black', linewidth = 0.75),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        strip.text.x = element_text(size = 7)) +
  labs(title = NULL, x = "Log Fold Change", y = "-Log P") +
  facet_grid(. ~ cell_line)


# Supp Figure 1c/d --------------------------------------------------------
# Adhesome gene category proportions
library(ggplot2)
library(tidyverse)

# Import adhesome gene data
adhesome <- read.csv("~/chen-origer-xist-adhesome/data/adhesome-components.csv")

# Import DEG data derived from DESeq2 results PMID:39546568 Dataset S01
xist_degs <- read.csv("~/chen-origer-xist-adhesome/data/supp-fig1ab-volcano.csv")

xist_degs <- merge(xist_degs, 
                   dplyr::select(adhesome, "gene", "cat"),
                   by = "gene")

xist_degs_both <- xist_degs %>% dplyr::filter(cat == "Both")
xist_degs_cad <- xist_degs_both %>% dplyr::mutate(cat = "Cadherin")
xist_degs_int <- xist_degs_both %>% dplyr::mutate(cat = "Integrin")

xist_degs <- xist_degs %>% 
  dplyr::filter(cat != "Both") %>% 
  rbind(xist_degs_cad, xist_degs_int)

xist_degs$deg <- stringr::str_to_title(xist_degs$deg)
xist_degs$cat <- factor(xist_degs$cat, levels = c("Integrin", "Cadherin"))

ggplot(data = filter(xist_degs, deg != "None" & cat != "Both"), 
       aes(x = deg, fill = cat)) +
  geom_bar(position = "fill", stat = "count") + 
  theme_classic() +
  scale_fill_manual(values = c("Integrin" = "red", 
                                "Cadherin" = "pink")) +
  theme(plot.title = element_text(size = 7, hjust = 0.5, face = "plain"),
        plot.subtitle = element_text(size = 7, hjust = 0.5, face = "plain"),
        legend.position = "none",
        axis.title = element_text(size = 7),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.text = element_text(size = 7, color = "black"),
        axis.line = element_line(color = 'black', linewidth = 0.75),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        strip.text.x = element_text(size = 7)) +
  labs(title = NULL, x = NULL, y = "Proportion") +
  facet_grid(. ~ cell_line)

# Supp Figure 1e ----------------------------------------------------------
# XIST knockdown lines PCA
library(DESeq2)
library(ggplot2)

count_data <- read.csv("~/chen-origer-xist-adhesome/data/deseq-countData.csv",
                       row.names = 1)
col_data <- read.csv("~/chen-origer-xist-adhesome/data/deseq-colData.csv",
                     row.names = 1)

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = count_data,
  colData = col_data,
  design = ~ group
)

# Prefilter low counts (recommended)
dds <- dds[rowSums(counts(dds)) > 10, ]

# Perform DESeq
dds <- DESeq2::DESeq(dds)

# Variance stabilizing transformation (for PCA)
vsd <- DESeq2::vst(dds, blind = TRUE)

pca_data <- DESeq2::plotPCA(vsd, intgroup = "group", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

ggplot(pca_data, aes(PC1, PC2, color = group)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percent_var[1], "% Variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% Variance")) +
  theme_classic(base_size = 14) +
  scale_color_manual(name = "Group",
                     values = c("Ctrl" = "#85b22c", 
                                "X7" = "#ffbc1f", 
                                "X9" = "#990f0f")) +
  coord_fixed(ratio = 3/2)

# Figure 2a ---------------------------------------------------------------
# lncRNA-adhesome correlations
library(tidyverse)
library(WGCNA)
library(ggridges)
library(ggplot2)

adhesome_expression_df <- readRDS("~/rds/gtex-adhesome-exp.rds") %>% 
  dplyr::select(!RNAss)
lncrna_expression_df <- readRDS("~/rds/gtex-lncrna-exp.rds") %>% 
  dplyr::select(!RNAss)

# Initialize matrices for bicor values and p-values
correlation_results <- matrix(NA, nrow = ncol(adhesome_expression_df), ncol = ncol(lncrna_expression_df))
p_value_results <- matrix(NA, nrow = ncol(adhesome_expression_df), ncol = ncol(lncrna_expression_df))

rownames(correlation_results) <- colnames(adhesome_expression_df)
colnames(correlation_results) <- colnames(lncrna_expression_df)
rownames(p_value_results) <- colnames(adhesome_expression_df)
colnames(p_value_results) <- colnames(lncrna_expression_df)

# Loop through each lncRNA and adhesome gene pair
for (lncRNA in colnames(lncrna_expression_df)) {
  for (gene in colnames(adhesome_expression_df)) {
    result <- WGCNA::bicorAndPvalue(lncrna_expression_df[[lncRNA]], adhesome_expression_df[[gene]])
    correlation_results[gene, lncRNA] <- result$bicor
    p_value_results[gene, lncRNA] <- result$p
  }
}

correlation_results_df <- as.data.frame(correlation_results)
p_value_df <- as.data.frame(p_value_results)

# then test if XIST is significant left shifted when comparing to other lncRNA

# Extract the correlations for XIST
xist_correlations <- correlation_results_df[ , "XIST"]

# Extract the correlations for all other lncRNAs
lncrna_names <- colnames(correlation_results_df)
lncrna_names <- lncrna_names[lncrna_names != "XIST"]

# Initialize vectors to store results
p_values <- numeric(length(lncrna_names))
names(p_values) <- lncrna_names

# Loop through each lncRNA name
for (lncrna in lncrna_names) {
  # Extract the correlation vector for the current lncRNA
  other_vector <- correlation_results_df[ , lncrna]
  
  # Perform the Wilcoxon test to check if XIST correlations are greater than the other lncRNAs
  test_result <- wilcox.test(xist_correlations, other_vector, alternative = "greater")
  
  # Store the p-value
  p_values[lncrna] <- test_result$p.value
}

# Print raw p-values
print(p_values)

lnc_corrs <- correlation_results_df %>%
  tibble::rownames_to_column(var = "adhesome_gene") %>%
  tidyr::pivot_longer(cols = -adhesome_gene, 
                      names_to = "lncRNA", 
                      values_to = "corr") %>%
  dplyr::mutate(condition = "GTEx")

# Load saved results
lnc_corrs <- read.csv("~/chen-origer-xist-adhesome/data/fig2a-lncRNA-adhesome-corrs.csv")

ggplot(lnc_corrs, aes(x = corr, y = lncRNA, fill = condition)) +
  geom_density_ridges(scale = 2, rel_min_height = 0.01) +
  scale_fill_manual(values = c("GTEx" = "#66a182",
                               "TCGA" = "#dd6e42")) +
  labs(title = "lncRNA-adhesome correlation in Xena", x = "Biweight Midcorrelation", y = NULL) +
  theme_classic() +
  theme(panel.spacing = unit(0.1, "lines"),
        plot.title = element_text(size = 7, hjust = 0.5, face = "plain"),
        plot.subtitle = element_text(size = 7, hjust = 0.5, face = "plain"),
        legend.position = "none",
        axis.title = element_text(size = 7),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.text = element_text(size = 7, color = "black"),
        axis.line = element_line(color = 'black', linewidth = 0.75),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        strip.text.x = element_text(size = 7))

# Figure 2b ---------------------------------------------------------------
# Z-score difference of adhesome genes in GTEx vs TCGA
library(ggplot2)
library(ggrepel)

# Load saved results
rtz <- read.csv("~/chen-origer-xist-adhesome/data/r-to-z-corrs.csv")

# Create an artificial maximum for plotting
rtz$neg_log_p <- ifelse(rtz$padj_diff <= 1e-282,
                        yes = 300,
                        no = -log10(rtz$padj_diff))

ggplot(rtz, aes(x = z_diff, y = neg_log_p)) +
  geom_point(aes(color = z_diff), size = 2) + # Add points colored by correlation difference
  scale_color_gradient2(low = "#DD6E42", mid = "#FFFFBF", high = "#66A182",
                        name = "Z-score Difference\n(GTEx - TCGA)",
                        midpoint = 0) +
  labs(title = NULL,
       x = "Z-score Difference",
       y = "-Log P") +
  theme_minimal() +
  theme(legend.position = "right",
        axis.title = element_text(size = 7),
        axis.text = element_text(size = 7, color = "black"),
        strip.text.y = element_text(size = 7),
        aspect.ratio = 1) +
  ggrepel::geom_text_repel(data = rbind(head(rtz[order(rtz$z_diff, decreasing = F),]),
                                        tail(rtz[order(rtz$z_diff, decreasing = F),],10)), 
                           aes(label = gene_name),
                           size = 3, 
                           box.padding = 0.5,
                           max.overlaps = 100)

# Figure 2c ---------------------------------------------------------------
# XIST-adhesome correlation by tissue type
library(tidyverse)
library(WGCNA)
library(ggplot2)

# Example data prep with gtex
# See "GTEx Data Preparation" for generating expression matrices
gtex <- readRDS("~/big-data/gtex-gene-tpm-female.rds")

# Initialize list to store correlation results for each tissue type
tissue_correlation_results <- list()

# Get the unique tissue types from the SMTS column
tissue_types <- unique(sample_df_female_RNAseq$SMTS)

# Initialize a list to store correlation and p-value results for each tissue type
tissue_correlation_results <- list()
tissue_pvalue_results <- list()

# Loop through each tissue type
for (tissue in tissue_types) {
  
  # Extract SAMPID for the current tissue type
  tissue_samples <- sample_df_female_RNAseq %>%
    dplyr::filter(SMTS == tissue) %>%
    dplyr::pull(SAMPID)
  
  # Extract the gene expression data for the current tissue from GTEx
  tissue_expression <- GTEx_gene_tpm_female_RNAseq %>%
    dplyr::select(Description, all_of(tissue_samples))  # Ensure to select relevant samples
  
  # Extract expression data for lncRNAs
  lncrna_expression <- tissue_expression %>%
    dplyr::filter(Description %in% lncRNAs)
  lncrna_expression_df <- as.data.frame(t(lncrna_expression[-1]))
  colnames(lncrna_expression_df) <- lncrna_expression$Description
  
  # Extract expression data for adhesome genes
  adhesome_expression <- tissue_expression %>%
    dplyr::filter(Description %in% Adhesome_genes_updated$gene_name)
  adhesome_expression_df <- as.data.frame(t(adhesome_expression[-1]))
  colnames(adhesome_expression_df) <- adhesome_expression$Description
  
  # Replace zeros with NA in lncRNA_expression_df
  lncrna_expression_df[lncrna_expression_df == 0] <- NA
  
  # Replace zeros with NA in adhesome_expression_df
  adhesome_expression_df[adhesome_expression_df == 0] <- NA
  
  # Initialize matrices to store correlation and p-value results
  correlation_results <- matrix(NA, nrow = ncol(adhesome_expression_df), ncol = ncol(lncrna_expression_df))
  pvalue_results <- matrix(NA, nrow = ncol(adhesome_expression_df), ncol = ncol(lncrna_expression_df))
  rownames(correlation_results) <- colnames(adhesome_expression_df)
  colnames(correlation_results) <- colnames(lncrna_expression_df)
  rownames(pvalue_results) <- colnames(adhesome_expression_df)
  colnames(pvalue_results) <- colnames(lncrna_expression_df)
  
  # Calculate biweight midcorrelation and p-value for each lncRNA and adhesome gene pair
  for (lncRNA in colnames(lncrna_expression_df)) {
    for (gene in colnames(adhesome_expression_df)) {
      result <- WGCNA::bicorAndPvalue(lncrna_expression_df[[lncRNA]], adhesome_expression_df[[gene]], use = "pairwise.complete.obs")
      correlation_results[gene, lncRNA] <- result$bicor
      pvalue_results[gene, lncRNA] <- result$p
    }
  }
  
  # Convert results to data frames for storage
  correlation_results_df <- as.data.frame(correlation_results)
  pvalue_results_df <- as.data.frame(pvalue_results)
  
  # Store the correlation and p-value results for the current tissue type
  tissue_correlation_results[[tissue]] <- correlation_results_df
  tissue_pvalue_results[[tissue]] <- pvalue_results_df
}

# Load saved results
tissue_corrs <- read.csv("~/chen-origer-xist-adhesome/data/fig2c-xist-adhesome-tissue-corrs.csv")

tissue_corrs$tissue <- factor(tissue_corrs$tissue, levels = rev(levels(tissue_corrs$tissue)))
ggplot(tissue_corrs, aes(x = corr, y = tissue, fill = study)) +
  geom_density_ridges(scale = 2, rel_min_height = 0.01) +
  scale_fill_manual(values = c("GTEx" = "#66A182",
                               "TCGA" = "#DD6E42")) +
  labs(title = NULL, x = NULL, y = "Biweight Midcorrelation") +
  theme_classic() +
  theme(panel.spacing = unit(0.1, "lines"),
        plot.title = element_text(size = 7, hjust = 0.5, face = "plain"),
        plot.subtitle = element_text(size = 7, hjust = 0.5, face = "plain"),
        legend.position = "none",
        axis.title = element_text(size = 7),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.text = element_text(size = 7, color = "black"),
        axis.line = element_line(color = 'black', linewidth = 0.75),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        strip.text.x = element_text(size = 7))

# Supp Figure 2a/b --------------------------------------------------------
# XIST expression by tissue type
library(tidyverse)
library(ggplot2)
library(ggridges)

# Import XIST expression data and reformat
gtex <- read.csv("~/chen-origer-xist-adhesome/data/supp-fig2a-gtex-xist.csv") %>%
  dplyr::select(any_of(c("SAMPID","TPM_log","SMTS"))) %>%
  dplyr::rename(sample_id = SAMPID,
                tissue = SMTS) %>%
  dplyr::mutate(study = "GTEx")

# Import XIST expression data and reformat
tcga <- read.csv("~/chen-origer-xist-adhesome/data/supp-fig2b-tcga-xist.csv") %>%
  dplyr::select(any_of(c("sample","TPM_log","X_primary_site"))) %>%
  dplyr::rename(sample_id = sample,
                tissue = X_primary_site) %>%
  dplyr::mutate(study = "TCGA")

# Make names consistent
tcga$tissue <- gsub("Thyroid Gland", "Thyroid", tcga$tissue)
gtex$tissue <- gsub("Cervix Uteri", "Cervix", gtex$tissue)
tcga$tissue <- gsub("Adrenal gland", "Adrenal Gland", tcga$tissue)

gtex$tissue <- factor(gtex$tissue, levels = rownames(gtex %>% 
                                                       dplyr::group_by(tissue) %>% 
                                                       dplyr::summarise(mean = mean(TPM_log)) %>% 
                                                       tibble::column_to_rownames(var = "tissue") %>%
                                                       dplyr::arrange(mean)))

ggplot(gtex, aes(x = TPM_log, y = tissue, fill = study)) +
  geom_density_ridges(scale = 2, rel_min_height = 0.01) +
  #geom_boxplot() +
  scale_fill_manual(values = c("GTEx" = "#66A182",
                               "TCGA" = "#dd6e42")) +
  labs(title = "XIST expression in GTEx", x = "log(TPM)", y = NULL) +
  theme_classic() +
  theme(panel.spacing = unit(0.1, "lines"),
        plot.title = element_text(size = 7, hjust = 0.5, face = "plain"),
        plot.subtitle = element_text(size = 7, hjust = 0.5, face = "plain"),
        legend.position = "none",
        axis.title = element_text(size = 7),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.text = element_text(size = 7, color = "black"),
        axis.line = element_line(color = 'black', linewidth = 0.75),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        strip.text.x = element_text(size = 7)) +
  xlim(-3, 10)

tcga$tissue <- factor(tcga$tissue, levels = rownames(tcga %>% 
                                                       dplyr::group_by(tissue) %>% 
                                                       dplyr::summarise(mean = mean(TPM_log)) %>% 
                                                       tibble::column_to_rownames(var = "tissue") %>%
                                                       dplyr::arrange(mean)))

ggplot(tcga, aes(x = TPM_log, y = tissue, fill = study)) +
  geom_density_ridges(scale = 2, rel_min_height = 0.01) +
  #geom_boxplot() +
  scale_fill_manual(values = c("GTEx" = "#66A182",
                               "TCGA" = "#dd6e42")) +
  labs(title = "XIST expression in TCGA", x = "log(TPM)", y = NULL) +
  theme_classic() +
  theme(panel.spacing = unit(0.1, "lines"),
        plot.title = element_text(size = 7, hjust = 0.5, face = "plain"),
        plot.subtitle = element_text(size = 7, hjust = 0.5, face = "plain"),
        legend.position = "none",
        axis.title = element_text(size = 7),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.text = element_text(size = 7, color = "black"),
        axis.line = element_line(color = 'black', linewidth = 0.75),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        strip.text.x = element_text(size = 7)) +
  xlim(-4, 15)

ggplot(filter(rbind(gtex,tcga), tissue %in% (gtex$tissue[gtex$tissue %in% tcga$tissue])), 
       aes(x = TPM_log, y = tissue, fill = study)) +
  geom_density_ridges(scale = 2, rel_min_height = 0.01) +
  #geom_boxplot() +
  scale_fill_manual(values = c("GTEx" = "#66A182",
                               "TCGA" = "#dd6e42")) +
  labs(title = "XIST expression in TCGA", x = "log(TPM)", y = NULL) +
  theme_classic() +
  theme(panel.spacing = unit(0.1, "lines"),
        plot.title = element_text(size = 7, hjust = 0.5, face = "plain"),
        plot.subtitle = element_text(size = 7, hjust = 0.5, face = "plain"),
        legend.position = "none",
        axis.title = element_text(size = 7),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.text = element_text(size = 7, color = "black"),
        axis.line = element_line(color = 'black', linewidth = 0.75),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        strip.text.x = element_text(size = 7)) +
  xlim(-2, 12)

# Get statistics
ggplot(filter(rbind(gtex,tcga), tissue %in% (gtex$tissue[gtex$tissue %in% tcga$tissue])), 
       aes(x = TPM_log, y = study, fill = study)) +
  #geom_density_ridges(scale = 2, rel_min_height = 0.01) +
  geom_boxplot() +
  scale_fill_manual(values = c("GTEx" = "#66A182",
                               "TCGA" = "#dd6e42")) +
  labs(title = "XIST expression in TCGA", x = "log(TPM)", y = NULL) +
  theme_classic() +
  theme(panel.spacing = unit(0.1, "lines"),
        plot.title = element_text(size = 7, hjust = 0.5, face = "plain"),
        plot.subtitle = element_text(size = 7, hjust = 0.5, face = "plain"),
        legend.position = "none",
        axis.title = element_text(size = 7),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.text = element_text(size = 7, color = "black"),
        axis.line = element_line(color = 'black', linewidth = 0.75),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        strip.text.x = element_text(size = 7)) +
  ggpubr::stat_compare_means(comparisons = list(c("GTEx", "TCGA")),
                             label = "p.signif") +
  xlim(-2, 12) +
  facet_wrap(. ~ tissue)


# Supp Figure 2c/d --------------------------------------------------------
# XIST correlation with other gene sets
library(ggplot2)
library(tidyverse)
library(WGCNA)

set.seed(123)

# See "GTEx Data Preparation" for generating expression matrices
gtex <- readRDS("~/chen-origer-xist-adhesome/data/gtex-gene-tpm-female.rds")

# Import adhesome gene data
adhesome <- read.csv("~/chen-origer-xist-adhesome/data/adhesome-components.csv")

# Convert dataset to numeric matrix
expr_matrix <- as.matrix(gtex[ , -1])
rownames(expr_matrix) <- gtex$Description

# Replace zeros with NA
expr_matrix[expr_matrix == 0] <- NA

# Extract XIST expression
xist_vector <- expr_matrix["XIST", ]

# Define adhesome gene set
adhesome_genes <- intersect(adhesome$gene,
                            rownames(expr_matrix))

# Get adhesome gene expression matrix
adhesome_matrix <- expr_matrix[adhesome_genes, ]

# Remove zero variance genes
adhesome_matrix <- adhesome_matrix[apply(adhesome_matrix, 1, var, na.rm=TRUE) > 0, ]

# Compute correlations in one shot
adhesome_cor <- WGCNA::bicor(t(adhesome_matrix), xist_vector,
                             use = "pairwise.complete.obs")

# Store median adhesome correlation for later plotting
observed_median <- median(adhesome_cor, na.rm=TRUE)

# Define background genes (not XIST)
# Any gene location
background_genes <- setdiff(rownames(expr_matrix),
                            c("XIST"))

# Only X-linked genes
library(biomaRt)
ensembl <- biomaRt::useMart("ensembl", 
                            dataset = "hsapiens_gene_ensembl")
gene_locations <- biomaRt::getBM(attributes = c("hgnc_symbol", 
                                                "chromosome_name"),
                                 filters = "hgnc_symbol",
                                 values = background_genes,
                                 mart = ensembl)
background_x_genes <- gene_locations %>% 
  dplyr::filter(chromosome_name == "X") %>%
  dplyr::pull(hgnc_symbol)

# Define number of permutations and empty vector
n_perm <- 2000
null_medians <- numeric(n_perm)

# Run permutations
for (i in 1:n_perm) {
  
  random_genes <- sample(background_genes, length(adhesome_genes))
  
  random_matrix <- expr_matrix[random_genes, ]
  
  random_matrix <- random_matrix[
    apply(random_matrix, 1, var, na.rm=TRUE) > 0, ]
  
  random_cor <- bicor(t(random_matrix), xist_vector,
                      use = "pairwise.complete.obs")
  
  null_medians[i] <- median(random_cor, na.rm=TRUE)
}

# Empirical p-value
empirical_p <- mean(null_medians >= observed_median) #fraction of sets where random median ≥ observed

# Only empirical p % of random sets were as extreme as your observed
message("Adhesome median: ", observed_median)
message("Empirical p-value: ", empirical_p)

# Optional: adjust empirical p-value if zero
n_perm <- length(null_medians)
empirical_p <- (1 + sum(null_medians >= observed_median)) / (1 + n_perm)

# Plot
ggplot(data.frame(null_medians = null_medians,
                  study = rep("GTEx", length(null_medians))), 
       aes(x = null_medians, fill = study)) + 
  geom_density(color = "black") +
  scale_fill_manual(values = c("GTEx" = "#66A182")) +
  theme_classic() +
  theme(panel.spacing = unit(0.1, "lines"),
        plot.title = element_text(size = 7, hjust = 0.5, face = "plain"),
        plot.subtitle = element_text(size = 7, hjust = 0.5, face = "plain"),
        legend.position = "none",
        axis.title = element_text(size = 7),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.text = element_text(size = 7, color = "black"),
        axis.line = element_line(color = 'black', linewidth = 0.75),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        strip.text.x = element_text(size = 7)) +
  labs(title = NULL, x = "XIST-random gene set correlations in GTEx", y = "Frequency") +
  geom_vline(xintercept = observed_median,
             color = "black", linetype = "dashed", linewidth = 0.75)

# Supp Figure 3a/b --------------------------------------------------------
# lncRNA-adhesome correlation ridgeplots
library(tidyverse)
library(ggplot2)
library(ggridges)

# See Figure 2a for data prep
# Load saved results
lnc_corrs <- read.csv("~/chen-origer-xist-adhesome/data/fig2a-lncRNA-adhesome-corrs.csv")

gtex_long <- dplyr::filter(lnc_corrs, condition == "GTEx")

gtex_long$lncRNA <- factor(gtex_long$lncRNA, levels = rownames(gtex_long %>% 
                                                                 dplyr::group_by(lncRNA) %>% 
                                                                 dplyr::summarise(mean = mean(corr)) %>% 
                                                                 tibble::column_to_rownames(var = "lncRNA") %>%
                                                                 dplyr::arrange(mean)))

ggplot(gtex_long, aes(x = corr, y = lncRNA, fill = condition)) + #fill = after_stat(x)
  #geom_density_ridges_gradient(scale = 2, rel_min_height = 0.01) +
  geom_density_ridges(scale = 2, rel_min_height = 0.01) +
  scale_fill_manual(values = c("GTEx" = "#66A182",
                               "TCGA" = "#dd6e42")) +
  #scale_fill_viridis(name = "Bicor", option = "C") +
  labs(title = "lncRNA-adhesome correlation in GTEx", x = "Biweight Midcorrelation", y = NULL) +
  theme_classic() +
  theme(panel.spacing = unit(0.1, "lines"),
        plot.title = element_text(size = 7, hjust = 0.5, face = "plain"),
        plot.subtitle = element_text(size = 7, hjust = 0.5, face = "plain"),
        legend.position = "none",
        axis.title = element_text(size = 7),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.text = element_text(size = 7, color = "black"),
        axis.line = element_line(color = 'black', linewidth = 0.75),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        strip.text.x = element_text(size = 7))

tcga_long <- dplyr::filter(lnc_corrs, condition == "TCGA")

tcga_long$lncRNA <- factor(tcga_long$lncRNA, levels = rownames(tcga_long %>% 
                                                                 dplyr::group_by(lncRNA) %>% 
                                                                 dplyr::summarise(mean = mean(corr)) %>% 
                                                                 tibble::column_to_rownames(var = "lncRNA") %>%
                                                                 dplyr::arrange(mean)))

ggplot(tcga_long, aes(x = corr, y = lncRNA, fill = condition)) + #fill = after_stat(x)
  #geom_density_ridges_gradient(scale = 2, rel_min_height = 0.01) +
  geom_density_ridges(scale = 2, rel_min_height = 0.01) +
  scale_fill_manual(values = c("GTEx" = "#66A182",
                               "TCGA" = "#dd6e42")) +
  #scale_fill_viridis(name = "Bicor", option = "C") +
  labs(title = "lncRNA-adhesome correlation in TCGA", x = "Biweight Midcorrelation", y = NULL) +
  theme_classic() +
  theme(panel.spacing = unit(0.1, "lines"),
        plot.title = element_text(size = 7, hjust = 0.5, face = "plain"),
        plot.subtitle = element_text(size = 7, hjust = 0.5, face = "plain"),
        legend.position = "none",
        axis.title = element_text(size = 7),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.text = element_text(size = 7, color = "black"),
        axis.line = element_line(color = 'black', linewidth = 0.75),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        strip.text.x = element_text(size = 7))

# Supp Figure 3c/d ---------------------------------------------------------
# XIST-adhesome correlations by tissue type
library(tidyverse)
library(ggplot2)

# See "Figure 2c" for data prep example
# Load saved results
gtex <- read.csv("~/chen-origer-xist-adhesome/data/supp-fig3c-xist-adhesome-tissue-corrs-gtex.csv")
tcga <- read.csv("~/chen-origer-xist-adhesome/data/supp-fig3d-xist-adhesome-tissue-corrs-tcga.csv")

# Make names more consistent
tcga$tissue <- gsub("Thyroid Gland", "Thyroid", tcga$tissue)
gtex$tissue <- gsub("Cervix Uteri", "Cervix", gtex$tissue)
tcga$tissue <- gsub("Adrenal gland", "Adrenal Gland", tcga$tissue)

gtex$tissue <- factor(gtex$tissue, levels = rownames(gtex %>% 
                                                       dplyr::group_by(tissue) %>% 
                                                       dplyr::summarise(mean = mean(corr)) %>% 
                                                       tibble::column_to_rownames(var = "tissue") %>%
                                                       dplyr::arrange(mean)))

ggplot(gtex, aes(x = corr, y = tissue, fill = study)) +
  geom_density_ridges(scale = 2, rel_min_height = 0.01) +
  scale_fill_manual(values = c("GTEx" = "#66A182",
                               "TCGA" = "#dd6e42")) +
  labs(title = "XIST-adhesome correlation in GTEx", x = "Biweight Midcorrelation", y = NULL) +
  theme_classic() +
  theme(panel.spacing = unit(0.1, "lines"),
        plot.title = element_text(size = 7, hjust = 0.5, face = "plain"),
        plot.subtitle = element_text(size = 7, hjust = 0.5, face = "plain"),
        legend.position = "none",
        axis.title = element_text(size = 7),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.text = element_text(size = 7, color = "black"),
        axis.line = element_line(color = 'black', linewidth = 0.75),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        strip.text.x = element_text(size = 7))

tcga$tissue <- factor(tcga$tissue, levels = rownames(tcga %>% 
                                                       dplyr::group_by(tissue) %>% 
                                                       dplyr::summarise(mean = mean(corr)) %>% 
                                                       tibble::column_to_rownames(var = "tissue") %>%
                                                       dplyr::arrange(mean)))
tcga <- tcga[complete.cases(tcga),]
ggplot(tcga, aes(x = corr, y = tissue, fill = study)) +
  geom_density_ridges(scale = 2, rel_min_height = 0.01) +
  scale_fill_manual(values = c("GTEx" = "#66A182",
                               "TCGA" = "#dd6e42")) +
  labs(title = "XIST-adhesome correlation in TCGA", x = "Biweight Midcorrelation", y = NULL) +
  theme_classic() +
  theme(panel.spacing = unit(0.1, "lines"),
        plot.title = element_text(size = 7, hjust = 0.5, face = "plain"),
        plot.subtitle = element_text(size = 7, hjust = 0.5, face = "plain"),
        legend.position = "none",
        axis.title = element_text(size = 7),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.text = element_text(size = 7, color = "black"),
        axis.line = element_line(color = 'black', linewidth = 0.75),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        strip.text.x = element_text(size = 7))


# Supp Figure 3e ----------------------------------------------------------
# lncRNA-adhesome correlations in mouse data
library(tidyverse)
library(ggplot2)

# Load saved data
mgi <- read.csv("~/chen-origer-xist-adhesome/data/supp-fig3e-lncRNA-adhesome-corrs-mouse.csv") %>% 
  dplyr::mutate(study = "MGI")

# Sort
mgi$lncRNA <- factor(mgi$lncRNA, levels = unique(mgi$lncRNA[order(mgi$mean_bicor)]))

ggplot(mgi, aes(x = corr, y = lncRNA, fill = study)) +
  geom_density_ridges(scale = 2, rel_min_height = 0.01) +
  scale_fill_manual(values = c("MGI" = "#66A182")) +
  labs(title = "lncRNA-adhesome correlation in MGI", x = "Biweight Midcorrelation", y = NULL) +
  theme_classic() +
  theme(panel.spacing = unit(0.1, "lines"),
        plot.title = element_text(size = 7, hjust = 0.5, face = "plain"),
        plot.subtitle = element_text(size = 7, hjust = 0.5, face = "plain"),
        legend.position = "none",
        axis.title = element_text(size = 7),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.text = element_text(size = 7, color = "black"),
        axis.line = element_line(color = 'black', linewidth = 0.75),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        strip.text.x = element_text(size = 7))

# Check p values
ggplot(mgi, aes(x = lncRNA, y = corr, fill = study)) +
  geom_violin() +
  scale_fill_manual(values = c("MGI" = "#66A182")) +
  theme_classic() +
  theme(panel.spacing = unit(0.1, "lines"),
        plot.title = element_text(size = 7, hjust = 0.5, face = "plain"),
        plot.subtitle = element_text(size = 7, hjust = 0.5, face = "plain"),
        legend.position = "none",
        axis.title = element_text(size = 7),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.text = element_text(size = 7, color = "black"),
        axis.line = element_line(color = 'black', linewidth = 0.75),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        strip.text.x = element_text(size = 7)) +
  ggpubr::stat_compare_means(comparisons = list(c("Xist","Neat1"),
                                                c("Xist","Hottip"),
                                                c("Xist","Malat1"),
                                                c("Xist","Hotair"),
                                                c("Xist","H19"),
                                                c("Xist","Pvt1"),
                                                c("Xist","Meg3"),
                                                c("Xist","Airn"),
                                                c("Xist","Gas5"),
                                                c("Xist","Kcnq1ot1")),
                             method = "wilcox.test")

# Figure 3a ---------------------------------------------------------------
# Cellular compartment stratified XIST-adhesome correlation
library(ggplot2)

# Load data
nucl_corrs <- read.csv("~/chen-origer-xist-adhesome/data/fig3a-xist-adhesome-corrs-nuclear.csv")

# Factor for plots
nucl_corrs$sample_type <- factor(nucl_corrs$sample_type, levels = c("Normal Tissue", "Primary Tumor"))
nucl_corrs$nucl_bin <- factor(nucl_corrs$nucl_bin, levels = c("Nuclear-localized", "Non-nuclear-localized"))
ggplot(nucl_corrs, aes(x = nucl_bin, y = bicor, fill = sample_type)) + 
  geom_violin(color = "black") + 
  scale_fill_manual(values = c("Normal Tissue" = "#66A182",
                               "Primary Tumor" = "#DD6E42")) +
  theme_classic() +
  theme(panel.spacing = unit(0.1, "lines"),
        plot.title = element_text(size = 7, hjust = 0.5, face = "plain"),
        plot.subtitle = element_text(size = 7, hjust = 0.5, face = "plain"),
        legend.position = "none",
        axis.title = element_text(size = 7),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.text = element_text(size = 7, color = "black"),
        axis.line = element_line(color = 'black', linewidth = 0.75),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        strip.text.x = element_text(size = 7)) +
  stat_summary(fun.data = median_hilow, fun.args = list(conf.int = 0.5),
               geom = "pointrange", color = "black") +
  labs(title = NULL, x = NULL, y = "Biweight midcorrelation") +
  facet_grid(. ~ sample_type) +
  ggpubr::stat_compare_means(comparisons = list(c("Nuclear-localized", "Non-nuclear-localized")),
                             label = "p.signif")

# Figure 3c ---------------------------------------------------------------
# Cellular compartment stratified hazard ratios
library(tidyverse)
library(ggplot2)

# Load data
nucl_hr <- read.csv("~/chen-origer-xist-adhesome/data/fig3c-adhesome-coxph-nuclear.csv")
nucl_hr <- nucl_hr %>% dplyr::mutate(sample_type = "Primary Tumor")

# Factor for plots
nucl_hr$nucl_bin <- factor(nucl_hr$nucl_bin, levels = c("Nuclear-localized", "Non-nuclear-localized"))

ggplot(nucl_hr, aes(x = nucl_bin, y = abs_log2_ratio, fill = sample_type)) + 
  geom_violin(color = "black") + 
  scale_fill_manual(values = c("Normal Tissue" = "#66A182",
                               "Primary Tumor" = "#DD6E42")) +
  theme_classic() +
  theme(panel.spacing = unit(0.1, "lines"),
        plot.title = element_text(size = 7, hjust = 0.5, face = "plain"),
        plot.subtitle = element_text(size = 7, hjust = 0.5, face = "plain"),
        legend.position = "none",
        axis.title = element_text(size = 7),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.text = element_text(size = 7, color = "black"),
        axis.line = element_line(color = 'black', linewidth = 0.75),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        strip.text.x = element_text(size = 7)) +
  stat_summary(fun.data = median_hilow, fun.args = list(conf.int = 0.5),
               geom = "pointrange", color = "black") +
  labs(title = NULL, x = NULL, y = "Abs Log2 Hazard Ratio") +
  facet_grid(. ~ sample_type) +
  ggpubr::stat_compare_means(comparisons = list(c("Nuclear-localized", "Non-nuclear-localized")),
                             label = "p.signif")

# Figure 3d ---------------------------------------------------------------
# XIST-binned hazard ratio forest plots
library(tidyverse)
library(survival)
library(survminer)
library(ggplot2)

# Import adhesome gene data
adhesome <- read.csv("~/chen-origer-xist-adhesome/data/adhesome-components.csv")

# Load survival data from Xena
xena_survival <- read.delim("~/rds/xena-tcga-gtex-survival.txt")

# Load expression data from Xena
xena <- readRDS("~/rds/xena-tcga-gtex_some.RDS")

get_dat_wide <- function(ct, genes) {
  # Subset expression/meta for one cancer with filters
  dat <- xena %>%
    dplyr::filter(
      cancer_type %in% ct,
      gene %in% genes,
      if (!is.null(restrict_gender)) gender == restrict_gender else TRUE,
      if (!is.null(restrict_sample_type)) sample_type == restrict_sample_type else TRUE
    ) %>%
    dplyr::select(sample, cancer_type, gene, expression_value)
  
  if (nrow(dat) == 0) return(NULL)
  
  # If some genes are entirely missing in this cancer, we still proceed for the present ones.
  present_genes <- intersect(genes, unique(dat$gene))
  if (length(present_genes) == 0) return(NULL)
  
  # Per-gene quartiles (Q1 and Q3) within this cancer
  th <- dat %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(
      q1 = quantile(expression_value, low_quantile,  na.rm = TRUE),
      q3 = quantile(expression_value, high_quantile, na.rm = TRUE),
      .groups = "drop"
    )
  q1 <- setNames(th$q1, th$gene)
  q3 <- setNames(th$q3, th$gene)
  
  # Pivot to wide to align samples across genes
  dat_wide <- dat %>%
    dplyr::distinct(sample, cancer_type, gene, expression_value) %>%
    tidyr::pivot_wider(names_from = gene, values_from = expression_value)
  
  # Construct per-gene High/Low bins (middle = NA)
  for (g in present_genes) {
    dat_wide[[paste0("bin_", g)]] <- ifelse(
      dat_wide[[g]] >= q3[g], paste0("High ", g),
      ifelse(dat_wide[[g]] <= q1[g], paste0("Low ", g), NA)
    )
  }
  
  # Attach survival
  dat_wide <- dat_wide %>%
    dplyr::left_join(select(xena_survival, sample, OS, OS.time), by = "sample") %>%
    dplyr::filter(!is.na(OS), !is.na(OS.time))
  
  if (nrow(dat_wide) == 0) return(NULL)
  dat_wide
}

# Network genes
genes_of_interest <- c("XIST","CTNNB1","CSNK1E","VEZT","SORBS3","ABL1","DOCK1","PLCG1","INPPL1","MAPK8","FER","FLNB","SH2B1","NISCH","ITGB3BP")

# Filters for get_dat_wide
{
  restrict_gender      <- "Female"        # set to NULL to disable
  restrict_sample_type <- "Primary Tumor" # set to NULL to disable
  
  # Quartile thresholds for binning
  low_quantile  <- 0.25
  high_quantile <- 0.75
}

cancers <- sort(unique(xena$cancer_type))

xena.subset <- get_dat_wide(cancers, c(adhesome$gene,"XIST")) #genes_of_interest
xena.subset <- xena.subset %>% dplyr::select(bin_XIST, OS, OS.time, tidyselect::any_of(c(adhesome$gene,"XIST")))

# Create cox models for high and low XIST bins
{
  cox.res.high <- survival::coxph(formula = Surv(OS.time, OS) ~ XIST+ACTN1+CFL1+CORO1B+CTTN+FLNA+KEAP1+LASP1+ENAH+NEXN+SVIL+VASP+CORO2A+ACTB+SORBS2+BCAR1+CAV1+SMPX+SH3KBP1+CRK+CRKL+EZR+FHL2+GAB1+GRB2+GRB7+HAX1+NEDD9+CASS4+TGFB1I1+ITGB1BP1+FERMT1+FERMT2+FERMT3+LPXN+PPFIA1+LPP+NF2+FBLIM1+MSN+NCK2+PALLD+PARVA+PARVB+PXN+LIMS1+LIMS2+RDX+OSTF1+NUDT16L1+SYNM+SDCBP+TLN1+TNS1+TES+TRIP6+VCL+SORBS3+LDB3+ZYX+NDEL1+SH2B1+ZFYVE21+SLC3A2+KTN1+LRP1+PVR+SDC4+ITGA1+ITGA2+ITGA3+ITGA4+ITGA5+ITGA6+ITGA7+ITGA8+ITGA9+ITGA10+ITGA11+ITGAD+ITGAE+ITGAL+ITGAM+ITGAV+ITGAX+ITGB1+ITGB2+ITGB3+ITGB4+ITGB5+ITGB6+ITGB7+ITGB8+NRP1+NRP2+CD151+PDE4D+PKD1+TRPM7+CALR+ASAP3+GIT1+GIT2+ARHGAP26+ASAP2+ARHGAP24+DLC1+AGAP2+STARD13+RASA1+DEF6+DOCK1+ELMO1+ARHGEF6+ARHGEF7+DNM2+RHOU+PLCG1+INPP5D+INPPL1+ITGB3BP+RAVER1+STAT3+SPTLC1+ILK+PAK1+PDPK1+PRKCA+PPM1M+PPM1F+PPP2CA+ABL1+CSK+PTK2+PTK2B+SRC+PEAK1+PTPRF+PTPN12+PTPRA+PTPN6+PTPN11+PRNP+MYH9+MACF1+ARPC2+MARCKS+PFN1+ABI1+ABI2 +ABI3+ANKRD28+CSRP1+IRS1+MAPK8IP3+PLEC+SORBS1+SHC1+MYOM1+TSPAN1+TUBA1B+VIM+SHARPIN+FABP3+NISCH+CIB1+CIB2+SRCIN1+ADAM12+CEACAM1+ENG+CD47+LAYN+SIRPA+THY1+PLAUR+KCNH2+SLC16A3+SLC9A1+HSPB1+HSPA2+CBL+RNF5+RNF185+ARHGAP5+ARHGAP32+BCAR3+RAPGEF1+SOS1+TIAM1+TRIO+VAV1+VAV2+VAV3+ARHGEF12+ARHGEF2+CYTH2+ARF1+HRAS+RAC1+RHOA+PLD1+CAPN1+CAPN2+CASP8+MMP14+PIK3CA+PIP5K1C+PTEN+PABPC1+SSH1+AKT1+MAPK1+MAPK8+LIMK1+PRKACA+ROCK1+ILKAP+INSR+FYN+LYN+SYK+TESK1+PTPN1+PTPRO+PTPRH+PTPN2+ACTN4+ADD1+AKAP5+CTNNA1+CTNNA2+CTNNA3+DBN1+EPB41+FLNB+FSCN1+LIMA1+SHROOM3+TJP1+TRIOBP+ABI2+ANK3+ARVCF+CD2AP+CGNL1+CLIP1+CTNNB1+CTNND1+CTNND2+DLG1+DLG5+IQGAP1+JUP+KRIT1+LIN7A+LMO7+MAGI1+MPP7+NUMB+PARD3+PDZD2+PKP4+PLEKHA7+SCRIB+SPTBN1+SSX2IP+TJP2+WIPF2+ACTR3+CAPZA1+DIAPH1+WASF2+WASL+CDH1+CDH10+CDH11+CDH12+CDH15+CDH18+CDH19+CDH2+CDH20+CDH22+CDH24+CDH3+CDH4+CDH5+CDH6+CDH7+CDH8+CDH9+FAT1+SDC1+AJAP1+HEG1+VEZT+TRPC4+TRPV4+DYNC1H1+MYH10+MYO6+MYO7A+ARF6+CDC42+RAP1A+SIPA1L3+RAPGEF2+ARHGAP12+ARHGAP21+ARHGAP35+RACGAP1+SH3BP1+ECT2+CSNK1E+CSNK2A1+PRKCD+PRKD1+EGFR+MET+FER+ACP1+PTPN14+PTPRJ+PTPRK+PTPRM+PTPRT+PTPRU+ADAM10+ADAM9+CASP3+PSEN1+CBLL1+GNA12+EXOC3+NME1, 
                        data = as.data.frame(filter(xena.subset, bin_XIST == "High XIST"))) # Evaluate among cancer types
  cox.high <- cox.res.high %>%
    broom::tidy() %>% 
    dplyr::mutate(hr = exp(estimate)) %>% 
    dplyr::select(term,p.value,hr)
  colnames(cox.high) <- paste0(colnames(cox.high), "_high")
  colnames(cox.high)[1] <- "gene"
  
  cox.res.low <- survival::coxph(formula = Surv(OS.time, OS) ~ XIST+ACTN1+CFL1+CORO1B+CTTN+FLNA+KEAP1+LASP1+ENAH+NEXN+SVIL+VASP+CORO2A+ACTB+SORBS2+BCAR1+CAV1+SMPX+SH3KBP1+CRK+CRKL+EZR+FHL2+GAB1+GRB2+GRB7+HAX1+NEDD9+CASS4+TGFB1I1+ITGB1BP1+FERMT1+FERMT2+FERMT3+LPXN+PPFIA1+LPP+NF2+FBLIM1+MSN+NCK2+PALLD+PARVA+PARVB+PXN+LIMS1+LIMS2+RDX+OSTF1+NUDT16L1+SYNM+SDCBP+TLN1+TNS1+TES+TRIP6+VCL+SORBS3+LDB3+ZYX+NDEL1+SH2B1+ZFYVE21+SLC3A2+KTN1+LRP1+PVR+SDC4+ITGA1+ITGA2+ITGA3+ITGA4+ITGA5+ITGA6+ITGA7+ITGA8+ITGA9+ITGA10+ITGA11+ITGAD+ITGAE+ITGAL+ITGAM+ITGAV+ITGAX+ITGB1+ITGB2+ITGB3+ITGB4+ITGB5+ITGB6+ITGB7+ITGB8+NRP1+NRP2+CD151+PDE4D+PKD1+TRPM7+CALR+ASAP3+GIT1+GIT2+ARHGAP26+ASAP2+ARHGAP24+DLC1+AGAP2+STARD13+RASA1+DEF6+DOCK1+ELMO1+ARHGEF6+ARHGEF7+DNM2+RHOU+PLCG1+INPP5D+INPPL1+ITGB3BP+RAVER1+STAT3+SPTLC1+ILK+PAK1+PDPK1+PRKCA+PPM1M+PPM1F+PPP2CA+ABL1+CSK+PTK2+PTK2B+SRC+PEAK1+PTPRF+PTPN12+PTPRA+PTPN6+PTPN11+PRNP+MYH9+MACF1+ARPC2+MARCKS+PFN1+ABI1+ABI2 +ABI3+ANKRD28+CSRP1+IRS1+MAPK8IP3+PLEC+SORBS1+SHC1+MYOM1+TSPAN1+TUBA1B+VIM+SHARPIN+FABP3+NISCH+CIB1+CIB2+SRCIN1+ADAM12+CEACAM1+ENG+CD47+LAYN+SIRPA+THY1+PLAUR+KCNH2+SLC16A3+SLC9A1+HSPB1+HSPA2+CBL+RNF5+RNF185+ARHGAP5+ARHGAP32+BCAR3+RAPGEF1+SOS1+TIAM1+TRIO+VAV1+VAV2+VAV3+ARHGEF12+ARHGEF2+CYTH2+ARF1+HRAS+RAC1+RHOA+PLD1+CAPN1+CAPN2+CASP8+MMP14+PIK3CA+PIP5K1C+PTEN+PABPC1+SSH1+AKT1+MAPK1+MAPK8+LIMK1+PRKACA+ROCK1+ILKAP+INSR+FYN+LYN+SYK+TESK1+PTPN1+PTPRO+PTPRH+PTPN2+ACTN4+ADD1+AKAP5+CTNNA1+CTNNA2+CTNNA3+DBN1+EPB41+FLNB+FSCN1+LIMA1+SHROOM3+TJP1+TRIOBP+ABI2+ANK3+ARVCF+CD2AP+CGNL1+CLIP1+CTNNB1+CTNND1+CTNND2+DLG1+DLG5+IQGAP1+JUP+KRIT1+LIN7A+LMO7+MAGI1+MPP7+NUMB+PARD3+PDZD2+PKP4+PLEKHA7+SCRIB+SPTBN1+SSX2IP+TJP2+WIPF2+ACTR3+CAPZA1+DIAPH1+WASF2+WASL+CDH1+CDH10+CDH11+CDH12+CDH15+CDH18+CDH19+CDH2+CDH20+CDH22+CDH24+CDH3+CDH4+CDH5+CDH6+CDH7+CDH8+CDH9+FAT1+SDC1+AJAP1+HEG1+VEZT+TRPC4+TRPV4+DYNC1H1+MYH10+MYO6+MYO7A+ARF6+CDC42+RAP1A+SIPA1L3+RAPGEF2+ARHGAP12+ARHGAP21+ARHGAP35+RACGAP1+SH3BP1+ECT2+CSNK1E+CSNK2A1+PRKCD+PRKD1+EGFR+MET+FER+ACP1+PTPN14+PTPRJ+PTPRK+PTPRM+PTPRT+PTPRU+ADAM10+ADAM9+CASP3+PSEN1+CBLL1+GNA12+EXOC3+NME1, 
                       data = as.data.frame(filter(xena.subset, bin_XIST == "Low XIST"))) # Evaluate among cancer types
  cox.low <- cox.res.low %>%
    broom::tidy() %>% 
    dplyr::mutate(hr = exp(estimate)) %>% 
    dplyr::select(term,p.value,hr)
  colnames(cox.low) <- paste0(colnames(cox.low), "_low")
  colnames(cox.low)[1] <- "gene"
}

# Figure 4a ---------------------------------------------------------------
# XIST/RNAss scatter plot by sample type
library(tidyverse)
library(ggplot2)

gtex <- readRDS("~/rds/gtex-xist-adhesome-exp.rds")
tcga <- readRDS("~/rds/tcga-xist-adhesome-exp.rds")
mets <- readRDS("~/rds/mets-xist-adhesome-exp.rds")

# Merge data for plotting
xist_rnass <- rbind(gtex %>% 
                      dplyr::select("XIST", "RNAss") %>% 
                      dplyr::mutate(sample_type = "Normal Tissue"),
                    tcga %>% 
                      dplyr::select("XIST", "RNAss") %>% 
                      dplyr::mutate(sample_type = "Primary Tumor"),
                    mets %>% 
                      dplyr::select("XIST", "RNAss") %>% 
                      dplyr::mutate(sample_type = "Metastatic"))

# Load data
xist_rnass <- read.csv("~/chen-origer-xist-adhesome/data/fig4ab-xist-rnass.csv")

# Factor for plots
xist_rnass$sample_type <- factor(xist_rnass$sample_type, levels = c("Normal Tissue", "Primary Tumor", "Metastatic"))

# Plot
ggplot(xist_rnass, aes(x = XIST, y = RNAss, color = sample_type)) + 
  geom_point() + 
  geom_smooth(method = "gam", color = "black", formula = y ~ x) +
  scale_color_manual(values = c("Normal Tissue" = "#66A182",
                               "Primary Tumor" = "#DD6E42",
                               "Metastatic" = "#F2E863")) +
  theme_classic() +
  theme(panel.spacing = unit(0.1, "lines"),
        plot.title = element_text(size = 7, hjust = 0.5, face = "plain"),
        plot.subtitle = element_text(size = 7, hjust = 0.5, face = "plain"),
        legend.position = "none",
        axis.title = element_text(size = 7),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.text = element_text(size = 7, color = "black"),
        axis.line = element_line(color = 'black', linewidth = 0.75),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        strip.text.x = element_text(size = 7)) +
  labs(title = NULL, x = "XIST Expression (Log TPM)", y = "RNAss", color = "Sample Type") +
  facet_grid(. ~ sample_type, scales = "free")

# Figure 4b ---------------------------------------------------------------
# XIST binned RNAss violin plot by sample type
library(tidyverse)
library(ggplot2)

gtex <- readRDS("~/rds/gtex-xist-adhesome-exp.rds")
tcga <- readRDS("~/rds/tcga-xist-adhesome-exp.rds")
mets <- readRDS("~/rds/mets-xist-adhesome-exp.rds")

# Merge data for plotting
xist_rnass <- rbind(gtex %>% 
                      dplyr::select("XIST", "RNAss") %>% 
                      dplyr::mutate(sample_type = "Normal Tissue"),
                    tcga %>% 
                      dplyr::select("XIST", "RNAss") %>% 
                      dplyr::mutate(sample_type = "Primary Tumor"),
                    mets %>% 
                      dplyr::select("XIST", "RNAss") %>% 
                      dplyr::mutate(sample_type = "Metastatic"))

# Define 25% and 75% quantiles
q_low  <- quantile(xist_rnass$XIST, 0.25, na.rm = TRUE)
q_high <- quantile(xist_rnass$XIST, 0.75, na.rm = TRUE)

xist_rnass$bin_XIST <- ifelse(xist_rnass$XIST > q_high, "High XIST",
                              ifelse(xist_rnass$XIST < q_low, "Low XIST", NA))

# Load data
xist_rnass <- read.csv("~/chen-origer-xist-adhesome/data/fig4ab-xist-rnass.csv")

# Factor for plots
xist_rnass$sample_type <- factor(xist_rnass$sample_type, levels = c("Normal Tissue", "Primary Tumor", "Metastatic"))
xist_rnass$XIST_group <- factor(xist_rnass$bin_XIST, levels = c("Low XIST", "High XIST"))

# Plot
ggplot(filter(xist_rnass, !is.na(bin_XIST)), aes(x = bin_XIST, y = RNAss, fill = sample_type)) + 
  geom_violin(color = "black") + 
  scale_fill_manual(values = c("Normal Tissue" = "#66A182",
                               "Primary Tumor" = "#DD6E42",
                               "Metastatic" = "#F2E863")) +
  theme_classic() +
  theme(panel.spacing = unit(0.1, "lines"),
        plot.title = element_text(size = 7, hjust = 0.5, face = "plain"),
        plot.subtitle = element_text(size = 7, hjust = 0.5, face = "plain"),
        legend.position = "none",
        axis.title = element_text(size = 7),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.text = element_text(size = 7, color = "black"),
        axis.line = element_line(color = 'black', linewidth = 0.75),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        strip.text.x = element_text(size = 7)) +
  stat_summary(fun.data = median_hilow, fun.args = list(conf.int = 0.5),
               geom = "pointrange", color = "black") +
  labs(title = NULL, x = "XIST", y = "RNAss") +
  facet_grid(. ~ sample_type) +
  ggpubr::stat_compare_means(comparisons = list(c("Low XIST", "High XIST")),
                             label = "p.signif")

# Figure 4c ---------------------------------------------------------------
# Adhesome/RNAss correlation density plot by sample type
library(tidyverse)
library(ggplot2)

gtex <- readRDS("~/rds/gtex-xist-adhesome-exp.rds")
tcga <- readRDS("~/rds/tcga-xist-adhesome-exp.rds")
mets <- readRDS("~/rds/mets-xist-adhesome-exp.rds")

# Calculate bicor between each gene and RNAss
genes_to_test <- c(adhesome$gene, "XIST")

# Switch to gtex or mets
expr_mat <- tcga

results <- lapply(genes_to_test, function(g) {
  if (! g %in% colnames(expr_mat)) return(NULL)
  r <- bicor(expr_mat[, g], expr_mat[, "RNAss"], use = "pairwise.complete.obs")
  p <- corPvalueStudent(r, nSamples = nrow(expr_mat))
  data.frame(gene = g, bicor = r, pval = p)
})

results_df <- do.call(rbind, results)

# Density plot of adhesome and stemness score correlation
bicor_adhesome_df <- results_df %>%
  filter(gene %in% adhesome$gene)

# Load saved results
adhesome_rnass <- read.csv("~/chen-origer-xist-adhesome/data/fig4c-adhesome-rnass.csv")

adhesome_rnass_medians <- adhesome_rnass %>%
  group_by(sample_type) %>%
  summarize(median = median(bicor))

# Factor for plots
adhesome_rnass$sample_type <- factor(adhesome_rnass$sample_type, levels = c("Normal Tissue", "Primary Tumor", "Metastatic"))

# Plot
ggplot(adhesome_rnass, aes(x = bicor, fill = sample_type)) + 
  geom_density(color = "black") +
  scale_fill_manual(values = c("Normal Tissue" = "#66A182",
                               "Primary Tumor" = "#DD6E42",
                               "Metastatic" = "#F2E863")) +
  theme_classic() +
  theme(panel.spacing = unit(0.1, "lines"),
        plot.title = element_text(size = 7, hjust = 0.5, face = "plain"),
        plot.subtitle = element_text(size = 7, hjust = 0.5, face = "plain"),
        legend.position = "none",
        axis.title = element_text(size = 7),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.text = element_text(size = 7, color = "black"),
        axis.line = element_line(color = 'black', linewidth = 0.75),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        strip.text.x = element_text(size = 7)) +
  labs(title = NULL, x = "Biweight midcorrelation", y = "Density") +
  facet_grid(. ~ sample_type) + 
  geom_vline(data = adhesome_rnass_medians, aes(xintercept = median),
             color = "black", linetype = "dashed", linewidth = 0.75)


# Figure 4d ---------------------------------------------------------------
# XIST-binned adhesome/RNAss correlation violin plot
library(ggplot2)

# Load saved results
xist_adhesome_rnass <- read.csv("~/chen-origer-xist-adhesome/data/fig4d-xist-adhesome-rnass.csv")

# Factor for plots
xist_adhesome_rnass$sample_type <- factor(xist_adhesome_rnass$sample_type, levels = c("Normal Tissue", "Primary Tumor", "Metastatic"))
xist_adhesome_rnass$XIST_group <- factor(xist_adhesome_rnass$XIST_group, levels = c("Low XIST", "High XIST"))
ggplot(xist_adhesome_rnass, aes(x = XIST_group, y = bicor, fill = sample_type)) + 
  geom_violin(color = "black") + 
  scale_fill_manual(values = c("Normal Tissue" = "#66A182",
                               "Primary Tumor" = "#DD6E42",
                               "Metastatic" = "#F2E863")) +
  theme_classic() +
  theme(panel.spacing = unit(0.1, "lines"),
        plot.title = element_text(size = 7, hjust = 0.5, face = "plain"),
        plot.subtitle = element_text(size = 7, hjust = 0.5, face = "plain"),
        legend.position = "none",
        axis.title = element_text(size = 7),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.text = element_text(size = 7, color = "black"),
        axis.line = element_line(color = 'black', linewidth = 0.75),
        axis.ticks = element_line(colour = "black", linewidth = 0.75),
        strip.text.x = element_text(size = 7)) +
  stat_summary(fun.data = median_hilow, fun.args = list(conf.int = 0.5),
               geom = "pointrange", color = "black") +
  labs(title = NULL, x = NULL, y = "Biweight midcorrelation", color = "Sample Type") +
  facet_grid(. ~ sample_type) +
  ggpubr::stat_compare_means(comparisons = list(c("Low XIST", "High XIST")),
                             label = "p.signif")
