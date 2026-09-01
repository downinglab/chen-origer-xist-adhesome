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
# Adhesome dysregulation during XIST knockdown

library(patchwork)
library(tidyverse)

# Get DESeq2 results from PMID:39546568 Dataset S01
res_X7 <- read.csv("~/chen-origer-xist-adhesome/data/deseq-res-X7.csv")
res_X9 <- read.csv("~/chen-origer-xist-adhesome/data/deseq-res-X9.csv")

# Import adhesome gene data
adhesome <- read.csv("~/chen-origer-xist-adhesome/data/adhesome-components.csv")

# Function for expression-matched permutation
run_xist_permutation <- function(res_df, label_name, seed = 42) {
  set.seed(seed)
  
  # Clean the input data and calculate baseline log2 expression
  clean_df <- res_df %>%
    filter(!is.na(gene), !is.na(baseMean), !is.na(log2FoldChange)) %>%
    mutate(log2_baseline = log2(baseMean + 1)) %>%
    mutate(group = ifelse(gene %in% adhesome$gene, "Adhesome", "Global Background"))
  
  # Define the exact adhesome genes present in this dataset
  adhesome_present <- clean_df %>% filter(group == "Adhesome")
  background_pool <- clean_df %>% filter(group == "Global Background")
  
  message("\n--- Running Permutations for ", label_name, " ---")
  message("Adhesome genes analyzed: ", nrow(adhesome_present))
  message("Background genes available: ", nrow(background_pool))
  
  # Divide the genome into 5 expression bins based on baseline expression
  quant_breaks <- quantile(clean_df$log2_baseline, probs = seq(0, 1, 0.2), na.rm = TRUE)
  unique_breaks <- unique(quant_breaks)
  
  clean_df$quintile <- cut(
    clean_df$log2_baseline, 
    breaks = unique_breaks, 
    include.lowest = TRUE, 
    labels = FALSE
  )
  
  # Re-split with the newly added quintile column
  adhesome_mapped <- clean_df %>% filter(group == "Adhesome")
  background_mapped <- clean_df %>% filter(group == "Global Background")
  
  
  # Map how many adhesome genes land in each quintile
  adhesome_bin_counts <- table(adhesome_mapped$quintile)
  background_by_bin <- split(background_mapped$gene, background_mapped$quintile)
  
  # Run the 2,000-iteration permutation loop
  n_perm <- 2000
  perm_results <- replicate(n_perm, {
    # Sample background genes matching the exact baseline profile of the adhesome
    sampled_genes <- lapply(names(adhesome_bin_counts), function(bin_num) {
      bin_pool <- background_by_bin[[bin_num]]
      sample(bin_pool, size = adhesome_bin_counts[bin_num], replace = FALSE)
    })
    random_gene_set <- unlist(sampled_genes)
    
    # Subset results for these sampled background genes
    random_subset <- clean_df %>% filter(gene %in% random_gene_set)
    
    c(
      mean_lfc = mean(random_subset$log2FoldChange, na.rm = TRUE),
      mean_abs_lfc = mean(abs(random_subset$log2FoldChange), na.rm = TRUE)
    )
  })
  
  perm_df <- as.data.frame(t(perm_results))
  
  # Calculate observed statistics for the real adhesome
  obs_mean_abs_lfc <- mean(abs(adhesome_mapped$log2FoldChange), na.rm = TRUE)
  
  # Calculate Empirical P-Values
  # Check our hypothesis that direct targets go UP (positive LFC)
  p_val_volatility <- (1 + sum(perm_df$mean_abs_lfc >= obs_mean_abs_lfc)) / (1 + n_perm)
  
  # Print stats to console
  message("Observed Adhesome Mean |LFC|: ", round(obs_mean_abs_lfc, 4))
  message("Observed Adhesome Median |LFC|: ", round(median(abs(adhesome_mapped$log2FoldChange)), 4))
  message("Adhesome genes with |LFC| > 2: ", sum(abs(adhesome_mapped$log2FoldChange) > 2.0))
  message("Matched Background Mean |LFC|: ", round(mean(perm_df$mean_abs_lfc), 4))
  message("Volatility Empirical p-value: ", round(p_val_volatility, 4))
  
  # Generate Plots
  p1 <- ggplot(perm_df, aes(x = mean_abs_lfc)) +
    geom_density(fill = "gray95", color = "slategray", linewidth = 0.8) +
    geom_vline(xintercept = obs_mean_abs_lfc, color = "darkred", linetype = "dashed", linewidth = 1) +
    annotate("text", x = obs_mean_abs_lfc, y = Inf, 
             label = paste0("Adhesome (|LFC| = ", round(obs_mean_abs_lfc, 3), ")"), 
             color = "darkred", fontface = "bold", size = 2.5, vjust = 1.5, hjust = 1.1) +
    theme_classic() +
    coord_cartesian(clip = "off") +
    theme(plot.title = element_text(size = 9, face = "bold"),
          axis.title = element_text(size = 8),
          axis.text = element_text(size = 7)) +
    labs(title = paste(label_name),
         x = "Mean Absolute Log2 Fold Change", y = "Density")
  
  return(list(plot = p1, p_val = p_val_volatility))
}

# Run permutation test on X7
results_X7 <- run_xist_permutation(res_X7, "X7")
results_X7$plot

# Run permutation test on X9
results_X9 <- run_xist_permutation(res_X9, "X9")
results_X9$plot

# Combine the plots (X7 vs X9)
combined_volatility_plot <- results_X7$plot + results_X9$plot +
  plot_annotation(
    title = "XIST Knockout DEG Results",
    subtitle = "Comparing Adhesome Absolute |Log2 Fold Change| to Expression-Matched Background Genes",
    theme = theme(plot.title = element_text(face="bold", size=12),
                  plot.subtitle = element_text(face="italic", size=10))
  )

# Render the plots
print(combined_volatility_plot)

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
							 
							 