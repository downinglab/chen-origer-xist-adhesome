# lncRNA Adhesome Correlation ---------------------------------------------

library(dplyr)
library(tidyverse)
library(WGCNA)

# Define list of representative lncRNAs
lncRNAs <- c("XIST", "MALAT1", "AIRN", "HOTAIR", "H19", "NEAT1", "HOTTIP", "GAS5", "MEG3", "PVT1", "KCNQ1OT1", "CCAT1", "PCAT1", "UCA1")

# Specify the path to GCT file downloaded from GTEx
file_path2 <- "~/data/GTEx_gene_tpm.gct"

# Read the GCT file
GTEx_gene_tpm <- read.table(file_path2, sep = "\t", header = TRUE, skip = 2)  # Skip the first two lines which contain metadata

# each subject is one patient with multiple tissue data in each sample id; 
# sample id is started with subject id (10 characters for 'GTEX-1; 9 character for 'GTEX-'; 5 character for 'K-')
sample.df <- GTEx_Analysis_v8_Annotations_SampleAttributesDS
subject.df <- GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS
# Add SUBJID to sample.df
sample.df <- sample.df %>%
  mutate(SUBJID = case_when(
    grepl("^GTEX-1", SAMPID) ~ substr(SAMPID, 1, 10),
    grepl("^GTEX-", SAMPID) ~ substr(SAMPID, 1, 9),
    grepl("^K-", SAMPID) ~ substr(SAMPID, 1, 5),
    TRUE ~ NA_character_  # In case there are other patterns, handle them appropriately
  ))
merged_sample_subject_df <- sample.df %>%
  inner_join(subject.df, by = "SUBJID")

# filter for female sample only
merged_sample_subject_df_female <- filter(merged_sample_subject_df, SEX == '2')

# SELECT SAMPID with RNASEQ data
merged_sample_subject_df_female_RNAseq <- filter(merged_sample_subject_df_female, SMAFRZE == 'RNASEQ')

# Modify SAMPID in the sample info dataframe to match the format in the gene expression dataframe
merged_sample_subject_df_female_RNAseq <- merged_sample_subject_df_female_RNAseq %>%
  mutate(SAMPID = gsub("-", ".", SAMPID))

# Check the common sample id from female RNAseq sample and gene_tpm
common_samples <- merged_sample_subject_df_female_RNAseq$SAMPID[merged_sample_subject_df_female_RNAseq$SAMPID %in% colnames(GTEx_gene_tpm)]
print(common_samples)

# Extract the female RNAseq dataframe from gene_tpm and prepare female sample information dataframe
GTEx_gene_tpm_filtered <- GTEx_gene_tpm %>%
  select(Description, all_of(common_samples))
sample_df_female_RNAseq <- merged_sample_subject_df_female_RNAseq

# Extract expression data for lncRNAs and adhesome genes (GTEx_gene_tpm_filtered is female data generated from GTEx_gene_tpm.gct)
lncrna_expression <- GTEx_gene_tpm_filtered %>%
  filter(Description %in% lncRNAs)
lncrna_expression_df <- as.data.frame(t(lncrna_expression[-1]))
colnames(lncrna_expression_df) <- lncrna_expression$Description

adhesome_expression <- GTEx_gene_tpm_filtered %>%
  filter(Description %in% Adhesome_genes_updated$gene_name)
adhesome_expression_df <- as.data.frame(t(adhesome_expression[-1]))
colnames(adhesome_expression_df) <- adhesome_expression$Description

# Replace zeros with NA in lncRNA_expression_df
lncrna_expression_df[lncrna_expression_df == 0] <- NA
# Filter columns with non-zero variance
valid_columns <- apply(lncrna_expression_df, 2, function(x) var(x, na.rm = TRUE) > 0)
lncrna_expression_df <- lncrna_expression_df[, valid_columns]

# Replace zeros with NA in adhesome_expression_df
adhesome_expression_df[adhesome_expression_df == 0] <- NA
# Filter columns with non-zero variance
valid_columns <- apply(adhesome_expression_df, 2, function(x) var(x, na.rm = TRUE) > 0)
adhesome_expression_df <- adhesome_expression_df[, valid_columns]

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
    result <- bicorAndPvalue(lncrna_expression_df[[lncRNA]], adhesome_expression_df[[gene]])
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


# Analyze the XIST-adhesome correlation across different tissue types
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
    filter(SMTS == tissue) %>%
    pull(SAMPID)
  
  # Extract the gene expression data for the current tissue from GTEx
  tissue_expression <- GTEx_gene_tpm_female_RNAseq %>%
    select(Description, all_of(tissue_samples))  # Ensure to select relevant samples
  
  # Extract expression data for lncRNAs
  lncrna_expression <- tissue_expression %>%
    filter(Description %in% lncRNAs)
  lncrna_expression_df <- as.data.frame(t(lncrna_expression[-1]))
  colnames(lncrna_expression_df) <- lncrna_expression$Description
  
  # Extract expression data for adhesome genes
  adhesome_expression <- tissue_expression %>%
    filter(Description %in% Adhesome_genes_updated$gene_name)
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
      result <- bicorAndPvalue(lncrna_expression_df[[lncRNA]], adhesome_expression_df[[gene]], use = "pairwise.complete.obs")
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


# GTEx RNAss Calculation --------------------------------------------------

library(TCGAbiolinks)
library(ggplot2)
library(dplyr)

# Follow the section 'lncRNA Adhesome Correlation', this script is used to calculate stemness score in normal tissues (using log2(TPM+0.001) instead of raw TPM).
gtex_data <- GTEx_gene_tpm_filtered # use female RNAseq dataframe
# Remove duplicated genes based on 'Description' and set rownames
gtex_data <- gtex_data[!duplicated(gtex_data$Description), ]
rownames(gtex_data) <- gtex_data$Description
gtex_data <- gtex_data[, -1]  # remove Description column

# Convert TPM to log2(TPM + 0.001)
gtex_data <- log2(gtex_data + 0.001)

# Quantile filtering (optional, keeps top 75% expressed genes)
dataFilt <- TCGAanalyze_Filtering(
  tabDF = gtex_data,
  method = "quantile",
  qnt.cut = 0.25
)

# Load stemness signature
data(SC_PCBC_stemSig)
print(SC_PCBC_stemSig)

# Calculate stemness scores using log2 TPM
Stemness_score <- TCGAanalyze_Stemness(
  stemSig = SC_PCBC_stemSig,
  dataGE = dataFilt
)

# GTEx XIST Adhesome RNAss ------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(WGCNA)
library(purrr)
library(car)

# Follow the section 'lncRNA Adhesome Correlation' to generate female RNAseq dataframe
gtex_data <- GTEx_gene_tpm_filtered
# Load Stemness_score calculated using script 'Calculate RNAss using GTEx'
Stemness_score <- GTEX_Stemness_Score_log2TPM

# Filter gtex_data for XIST and adhesome genes
genes_of_interest <- c("XIST", adhesome_genes)
gtex_filtered <- gtex_data %>%
  filter(Description %in% genes_of_interest)

# Convert data to long format
gtex_long <- gtex_filtered %>%
  pivot_longer(
    cols = -Description,
    names_to = "Sample",
    values_to = "expression_value"
  )

# Join stemness score
gtex_long <- gtex_long %>%
  inner_join(Stemness_score, by = c("Sample" = "Sample"))

# Log2 transform TPM
gtex_long <- gtex_long %>%
  mutate(log2_expr = log2(expression_value + 0.001))

#### Calculate the XIST-stemness score and adhesome-stemness score correlation
bicor_results <- gtex_long %>%
  group_by(Description) %>%
  summarize(
    tmp = list(bicorAndPvalue(log2_expr, stemness_score, use = "pairwise.complete.obs")),
    .groups = "drop"
  ) %>%
  mutate(
    bicor_value = map_dbl(tmp, ~ .x$bicor),
    p_value     = map_dbl(tmp, ~ .x$p)
  ) %>%
  select(Description, bicor_value, p_value)

# Scatter plot to show the XIST-adhesome correlation
xist_data <- gtex_long %>%
  filter(Description == "XIST")

ggplot(xist_data, aes(x = log2_expr, y = stemness_score)) +
  geom_point(alpha = 0.6, color = "#33a02c", size = 2) +  # points
  geom_smooth(method = "lm", color = "black", se = TRUE) + # regression line
  theme_minimal() +
  labs(
    x = "XIST Expression (log2(TPM+0.001))",
    y = "Stemness Score",
    title = NULL
  ) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold"),
    panel.grid = element_blank()   # removes all grid lines
  )


# Density plot to show the adhesome-stemness correlation
bicor_filtered <- bicor_results %>%
  filter(Description %in% adhesome_genes)

ggplot(bicor_filtered, aes(x = bicor_value)) +
  geom_density(fill = "#33a02c", alpha = 0.6) +
  theme_minimal() +
  labs(
    x = "Bicor",
    y = "Density",
    title = "Adhesome–Stemness Score Bicor"
  ) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold")
  )


# Check normality first
shapiro_test <- shapiro.test(bicor_filtered$bicor_value)
shapiro_test$p.value
# One-sample test: null hypothesis bicor = 0
if(shapiro_test$p.value > 0.05){
  # t-test
  test_result <- t.test(bicor_filtered$bicor_value, mu = 0, alternative = "less") # test if mean < 0
} else {
  # Wilcoxon signed-rank test
  test_result <- wilcox.test(bicor_filtered$bicor_value, mu = 0, alternative = "less")
}

test_result


#### Plot stemness score in low vs high XIST samples

# Calculate 25th and 75th quantiles for XIST expression
xist_quantiles <- quantile(xist_data$log2_expr, probs = c(0.25, 0.75), na.rm = TRUE)

# Assign each sample to Low, High, or Intermediate XIST group
xist_data <- xist_data %>%
  mutate(XIST_group = case_when(
    log2_expr <= xist_quantiles[1] ~ "Low XIST",
    log2_expr >= xist_quantiles[2] ~ "High XIST",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(XIST_group))  # keep only Low and High groups

# Make XIST_group a factor with desired order
xist_data$XIST_group <- factor(xist_data$XIST_group, levels = c("Low XIST", "High XIST"))

# Step 3: Plot stemness score by XIST group
ggplot(xist_data, aes(x = XIST_group, y = stemness_score, fill = XIST_group)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, color = "black", fill = "yellow") +
  scale_fill_manual(values = c("Low XIST" = "#33a02c", "High XIST" = "#33a02c")) +
  theme_minimal() +
  labs(
    x = "XIST Expression Group",
    y = "Stemness Score",
    title = NULL
  ) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, face = "bold")
  )

# Step 4: Optional: Statistical comparison
# Check normality
shapiro_low <- shapiro.test(xist_data$stemness_score[xist_data$XIST_group == "Low XIST"])
shapiro_high <- shapiro.test(xist_data$stemness_score[xist_data$XIST_group == "High XIST"])

# Choose test based on normality
if(shapiro_low$p.value > 0.05 & shapiro_high$p.value > 0.05){
  test_result <- t.test(stemness_score ~ XIST_group, data = xist_data, var.equal = TRUE)
} else {
  test_result <- wilcox.test(stemness_score ~ XIST_group, data = xist_data)
}

print(test_result)


#### Compare adhesome-stemness score cor in low vs high XIST samples
# Extract XIST expression
xist_expr <- gtex_long %>%
  filter(Description == "XIST") %>%
  select(Sample, log2_expr)

# Merge with stemness score
xist_expr <- xist_expr %>%
  inner_join(GTEX_Stemnes_Score, by = "Sample")

# Define Low/High XIST based on 25% and 75% quantiles
quantiles <- quantile(xist_expr$log2_expr, probs = c(0.25, 0.75), na.rm = TRUE)
xist_expr <- xist_expr %>%
  mutate(XIST_group = case_when(
    log2_expr <= quantiles[1] ~ "Low XIST",
    log2_expr >= quantiles[2] ~ "High XIST",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(XIST_group))

# Prepare adhesome expression
adhesome_expr <- gtex_long %>%
  filter(Description %in% adhesome_genes) %>%
  select(Sample, Description, log2_expr)

# Merge adhesome expression with XIST group and stemness score
adhesome_long <- adhesome_expr %>%
  inner_join(xist_expr %>% select(Sample, XIST_group, stemness_score), by = "Sample")

# Calculate bicor for each adhesome gene and stemness within each XIST group
bicor_results <- adhesome_long %>%
  group_by(Description, XIST_group) %>%
  summarize(
    tmp = list(bicorAndPvalue(log2_expr, stemness_score, use = "pairwise.complete.obs")),
    .groups = "drop"
  ) %>%
  mutate(
    bicor_value = map_dbl(tmp, ~ .x$bicor),
    p_value     = map_dbl(tmp, ~ .x$p)
  ) %>%
  select(Description, XIST_group, bicor_value, p_value)

# Plot distribution of bicor values by XIST group
bicor_results$XIST_group <- factor(bicor_results$XIST_group, levels = c("Low XIST", "High XIST"))

ggplot(bicor_results, aes(x = XIST_group, y = bicor_value, fill = XIST_group)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, color = "black", fill = "yellow") +
  scale_fill_manual(values = c("Low XIST" = "#33a02c", "High XIST" = "#33a02c")) +
  theme_minimal() +
  labs(y = "Adhesome Bicor with Stemness") +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  )

# Normality test and appropriate statistical comparison
shapiro_low <- shapiro.test(bicor_results$bicor_value[bicor_results$XIST_group == "Low XIST"])
shapiro_high <- shapiro.test(bicor_results$bicor_value[bicor_results$XIST_group == "High XIST"])

# If both pass normality, check for equal variance and use t-test; otherwise, Wilcoxon test
if(shapiro_low$p.value > 0.05 & shapiro_high$p.value > 0.05){
  # Check variance equality
  library(car)
  levene_result <- leveneTest(bicor_value ~ XIST_group, data = bicor_results)
  equal_var <- ifelse(levene_result$`Pr(>F)`[1] > 0.05, TRUE, FALSE)
  
  test_result <- t.test(bicor_value ~ XIST_group, data = bicor_results, var.equal = equal_var)
} else {
  test_result <- wilcox.test(bicor_value ~ XIST_group, data = bicor_results)
}

print(test_result)


# TCGA XIST Adhesome RNAss ------------------------------------------------


library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)
library(WGCNA)
library(tidyr)
library(ggplot2)

# Download 'TcgaTargetGtex_rsem_gene_tpm.RDS' from Xena and load it into R
my_data <- readRDS("~/data/TcgaTargetGtex_rsem_gene_tpm.RDS")
# Clean Ensembl IDs (remove version numbers)
ensembl_ids <- gsub("\\..*", "", rownames(my_data))
# Map Ensembl IDs to HGNC gene symbols
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = ensembl_ids,
                       column = "SYMBOL",   
                       keytype = "ENSEMBL", 
                       multiVals = "first") 
my_data_with_names <- cbind(gene_name = gene_symbols, my_data)

# Filter samples by Primary Tumor or Metastatic before next step
filtered_samples <- TCGA_phenotype_combined %>%
  filter(X_sample_type == "Primary Tumor",
         X_gender == "Female",
         X_study == "TCGA") %>%
  pull(sample)  # extract the sample column as a vector

# Keep 'gene_name' column + filtered sample columns
my_data_filtered <- my_data_with_names %>%
  select(gene_name, all_of(filtered_samples))

# Combine adhesome and lncRNA gene names into one vector
genes_of_interest <- c(Adhesome_genes_updated$gene_name, "XIST")

# Filter my_data_filtered to keep only those genes
my_data_subset <- my_data_filtered %>%
  filter(gene_name %in% genes_of_interest)

# prepare RNAss for matching samples; TCGA phenotype downloaded from Xena
RNAss_df <- TCGA_phenotype_combined %>%
  select(sample, RNAss)

common_samples <- intersect(colnames(my_data_subset)[-1], RNAss_df$sample)

expr_mat <- my_data_subset %>%
  select(gene_name, all_of(common_samples))

RNAss_vec <- RNAss_df %>%
  filter(sample %in% common_samples) %>%
  arrange(match(sample, common_samples)) %>%
  pull(RNAss)

# append RNAss as a new row
expr_with_RNAss <- bind_rows(
  expr_mat,
  data.frame(gene_name = "RNAss", t(RNAss_vec)) %>% 
    setNames(colnames(expr_mat))
)

# Transpose to genes as columns
# Remove existing rownames
rownames(expr_with_RNAss) <- NULL

# Convert column to rownames and transpose
expr_long <- expr_with_RNAss %>%
  column_to_rownames(var = "gene_name") %>%
  t()

# Calculate bicor between each gene and RNAss
genes_to_test <- c(Adhesome_genes_updated$gene_name, "XIST")

results <- lapply(genes_to_test, function(g) {
  if (! g %in% colnames(expr_long)) return(NULL)
  r <- bicor(expr_long[, g], expr_long[, "RNAss"], use = "pairwise.complete.obs")
  p <- corPvalueStudent(r, nSamples = nrow(expr_long))
  data.frame(gene = g, bicor = r, pval = p)
})

results_df <- do.call(rbind, results)



## scatter plot of XIST and stemness score correlation
plot_df <- data.frame(
  XIST = expr_long[, "XIST"],
  RNAss = expr_long[, "RNAss"]
)

ggplot(plot_df, aes(x = XIST, y = RNAss)) +
  geom_point(color = "#DA7049", alpha = 0.6) +
  geom_smooth(method = "lm", color = "black", se = FALSE) +
  labs(title = NULL,
       x = "XIST expression", y = "RNAss") +
  theme_classic()




## Density plot of adhesome and stemness score correlation
bicor_adhesome_df <- results_df %>%
  filter(gene %in% Adhesome_genes_updated$gene_name)

ggplot(bicor_adhesome_df, aes(x = bicor)) +
  geom_density(fill = "#DA7049", alpha = 0.5) +
  labs(title = "Density of adhesome/RNAss bicor",
       x = "Bicor", y = "Density") +
  theme_classic()

# Check normality first
shapiro_test <- shapiro.test(bicor_adhesome_df$bicor)
shapiro_test$p.value
# One-sample test: null hypothesis bicor = 0
if(shapiro_test$p.value > 0.05){
  # t-test
  test_result <- t.test(bicor_adhesome_df$bicor, mu = 0, alternative = "less") # test if mean < 0
} else {
  # Wilcoxon signed-rank test
  test_result <- wilcox.test(bicor_adhesome_df$bicor, mu = 0, alternative = "less")
}

test_result



## Compare stemness score in low vs high XIST samples
XIST_expr <- expr_long[, "XIST"]
# Define 25% and 75% quantiles
q_low  <- quantile(XIST_expr, 0.25, na.rm = TRUE)
q_high <- quantile(XIST_expr, 0.75, na.rm = TRUE)

XIST_group <- ifelse(XIST_expr > q_high, "High_XIST",
                     ifelse(XIST_expr < q_low, "Low_XIST", NA))

expr_long_df <- as.data.frame(expr_long)
expr_long_df$XIST_group <- XIST_group

expr_long_df_filtered <- expr_long_df %>% filter(!is.na(XIST_group))

# Calculate mean and SD for each group
stats_df <- expr_long_df_filtered %>%
  group_by(XIST_group) %>%
  summarise(
    mean_RNAss = mean(RNAss, na.rm = TRUE),
    sd_RNAss = sd(RNAss, na.rm = TRUE)
  )

# Convert XIST_group to factor with desired order
expr_long_df_filtered$XIST_group <- factor(expr_long_df_filtered$XIST_group,
                                           levels = c("Low_XIST", "High_XIST"))

# Violin plot with boxplot + mean ± SD
ggplot(expr_long_df_filtered, aes(x = XIST_group, y = RNAss, fill = XIST_group)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, color = "black", fill = "yellow") +
  scale_fill_manual(values = c("Low_XIST" = "#DA7049", "High_XIST" = "#DA7049")) +
  theme_minimal() +
  labs(
    x = NULL,
    y = "RNAss",
    title = NULL
  ) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16)
  )

# Normality test
shapiro_low  <- shapiro.test(expr_long_df_filtered$RNAss[expr_long_df_filtered$XIST_group == "Low_XIST"])
shapiro_high <- shapiro.test(expr_long_df_filtered$RNAss[expr_long_df_filtered$XIST_group == "High_XIST"])

shapiro_low$p.value
shapiro_high$p.value

# Choose test based on normality
if(shapiro_low$p.value > 0.05 & shapiro_high$p.value > 0.05){
  test_result <- t.test(RNAss ~ XIST_group, data = expr_long_df_filtered)
} else {
  test_result <- wilcox.test(RNAss ~ XIST_group, data = expr_long_df_filtered)
}

test_result



## compare adhesome and stemness score correlation in low vs high XIST samples
low_samples  <- expr_long_df_filtered %>% filter(XIST_group == "Low_XIST") %>% rownames()
high_samples <- expr_long_df_filtered %>% filter(XIST_group == "High_XIST") %>% rownames()

bicor_by_group <- lapply(c("Low_XIST", "High_XIST"), function(group_name) {
  samples <- if(group_name == "Low_XIST") low_samples else high_samples
  bicor_values <- sapply(Adhesome_genes_updated$gene_name, function(gene) {
    bicor(expr_long[samples, gene], expr_long[samples, "RNAss"], use = "pairwise.complete.obs")
  })
  data.frame(
    XIST_group = group_name,
    gene = names(bicor_values),
    bicor = bicor_values
  )
})

bicor_df <- do.call(rbind, bicor_by_group)

# Violin plot
ggplot(bicor_df, aes(x = factor(XIST_group, levels = c("Low_XIST", "High_XIST")), y = bicor, fill = XIST_group)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, color = "black", fill = "yellow") +
  scale_fill_manual(values = c("Low_XIST" = "#DA7049", "High_XIST" = "#DA7049")) +
  theme_minimal() +
  labs(
    x = NULL,
    y = "Adhesome Bicor with RNAss",
    title = NULL
  ) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16)
  )

# Normality and statistical test
shapiro_low  <- shapiro.test(bicor_df$bicor[bicor_df$XIST_group == "Low_XIST"])
shapiro_high <- shapiro.test(bicor_df$bicor[bicor_df$XIST_group == "High_XIST"])

if(shapiro_low$p.value > 0.05 & shapiro_high$p.value > 0.05){
  test_result <- t.test(bicor ~ XIST_group, data = bicor_df)
} else {
  test_result <- wilcox.test(bicor ~ XIST_group, data = bicor_df)
}

test_result
