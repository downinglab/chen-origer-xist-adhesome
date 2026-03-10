# GTEx Data Preparation ---------------------------------------------------

library(dplyr)
library(tidyverse)
library(WGCNA)

# Import adhesome gene data
adhesome <- read.csv("~/chen-origer-xist-adhesome/data/adhesome-components.csv")

# Define list of representative lncRNAs
lncRNAs <- c("XIST", "MALAT1", "AIRN", "HOTAIR", "H19", "NEAT1", "HOTTIP", "GAS5", "MEG3", "PVT1", "KCNQ1OT1", "CCAT1", "PCAT1", "UCA1")

# Read GTEx GCT file
# Genes as rows and samples as columns
GTEx_gene_tpm <- read.table("~/big-data/GTEx_gene_tpm.gct", sep = "\t", header = TRUE, skip = 2)  # Skip the first two lines which contain metadata

# Each subject is one patient with multiple tissue data in each sample id; 
# Sample id is started with subject id (10 characters for 'GTEX-1; 9 character for 'GTEX-'; 5 character for 'K-')
sample_df <- read.csv("~/big-data/GTEx_Analysis_v11_Annotations_SampleAttributesDS.csv", header = TRUE)
subject_df <- read.table("~/big-data/GTEx_Analysis_v11_Annotations_SubjectPhenotypesDS.txt", sep = "\t", header = TRUE)

# Add SUBJID to sample_df
sample_df <- sample_df %>%
  dplyr::mutate(SUBJID = dplyr::case_when(
    grepl("^GTEX-1", SAMPID) ~ substr(SAMPID, 1, 10),
    grepl("^GTEX-", SAMPID) ~ substr(SAMPID, 1, 9),
    grepl("^K-", SAMPID) ~ substr(SAMPID, 1, 5),
    TRUE ~ NA_character_  # In case there are other patterns, handle them appropriately
  ))
merged_sample_subject_df <- sample_df %>%
  dplyr::inner_join(subject_df, by = "SUBJID")

# Filter for female samples only
merged_sample_subject_df_female <- dplyr::filter(merged_sample_subject_df, SEX == '2')

# Filter for samples with RNAseq data
merged_sample_subject_df_female_RNAseq <- dplyr::filter(merged_sample_subject_df_female, SMAFRZE == 'RNASEQ')

# Modify SAMPID in the sample info dataframe to match the format in the gene expression dataframe
merged_sample_subject_df_female_RNAseq <- merged_sample_subject_df_female_RNAseq %>%
  dplyr::mutate(SAMPID = gsub("-", ".", SAMPID))

# Check the common sample id from female RNAseq sample and gene_tpm
common_samples <- merged_sample_subject_df_female_RNAseq$SAMPID[merged_sample_subject_df_female_RNAseq$SAMPID %in% colnames(GTEx_gene_tpm)]
print(common_samples)

# Extract the female RNAseq dataframe from gene_tpm and prepare female sample information dataframe
GTEx_gene_tpm_filtered <- GTEx_gene_tpm %>%
  dplyr::select(Description, all_of(common_samples))
sample_df_female_RNAseq <- merged_sample_subject_df_female_RNAseq

gtex_data <- GTEx_gene_tpm_filtered # use female RNAseq dataframe
# Remove duplicated genes based on 'Description' and set rownames
gtex_data <- gtex_data[!duplicated(gtex_data$Description), ]
rownames(gtex_data) <- gtex_data$Description
gtex_data <- gtex_data[, -1]  # remove Description column

# Convert TPM to log2(TPM + 0.001)
gtex_data <- log2(gtex_data + 0.001)

library(TCGAbiolinks)
# Quantile filtering (optional, keeps top 75% expressed genes)
dataFilt <- TCGAbiolinks::TCGAanalyze_Filtering(
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
Stemness_score <- Stemness_score %>% 
  dplyr::select(stemness_score) %>% 
  dplyr::rename(RNAss = stemness_score)

saveRDS(Stemness_score, "~/rds/gtex-stemness-score.rds")

# Extract expression data for lncRNAs and adhesome genes (GTEx_gene_tpm_filtered is female data generated from GTEx_gene_tpm.gct)
lncrna_expression <- GTEx_gene_tpm_filtered %>%
  dplyr::filter(Description %in% lncRNAs)
lncrna_expression_df <- as.data.frame(t(lncrna_expression[-1]))
colnames(lncrna_expression_df) <- lncrna_expression$Description

adhesome_expression <- GTEx_gene_tpm_filtered %>%
  dplyr::filter(Description %in% adhesome$gene)
adhesome_expression_df <- as.data.frame(t(adhesome_expression[-1]))
colnames(adhesome_expression_df) <- adhesome_expression$Description

# Extract expression data for XIST and adhesome genes
xist_adhesome_expression <- GTEx_gene_tpm_filtered %>%
  dplyr::filter(Description %in% c(adhesome$gene, "XIST"))
xist_adhesome_expression_df <- as.data.frame(t(xist_adhesome_expression[-1]))
colnames(xist_adhesome_expression_df) <- xist_adhesome_expression$Description

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

# Replace zeros with NA in xist_adhesome_expression_df
xist_adhesome_expression_df[xist_adhesome_expression_df == 0] <- NA
# Filter columns with non-zero variance
valid_columns <- apply(xist_adhesome_expression_df, 2, function(x) var(x, na.rm = TRUE) > 0)
xist_adhesome_expression_df <- xist_adhesome_expression_df[, valid_columns]

xist_adhesome_expression_df <- merge(xist_adhesome_expression_df, 
                                     Stemness_score, 
                                     by = "row.names") %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames(var = "Row.names")

adhesome_expression_df <- merge(adhesome_expression_df, 
                                Stemness_score, 
                                by = "row.names") %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames(var = "Row.names")

lncrna_expression_df <- merge(lncrna_expression_df, 
                              Stemness_score, 
                              by = "row.names") %>% 
  tibble::remove_rownames() %>% 
  tibble::column_to_rownames(var = "Row.names")

saveRDS(xist_adhesome_expression_df, "~/rds/gtex-xist-adhesome-exp.rds")
saveRDS(adhesome_expression_df, "~/rds/gtex-adhesome-exp.rds")
saveRDS(lncrna_expression_df, "~/rds/gtex-lncrna-exp.rds")
saveRDS(sample_df_female_RNAseq, "~/rds/gtex-female-metadata.rds")

# TCGA Data Preparation ---------------------------------------------------

library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)
library(WGCNA)
library(tidyr)
library(ggplot2)

# Import adhesome gene data
adhesome <- read.csv("~/chen-origer-xist-adhesome/data/adhesome-components.csv")

# Download 'TcgaTargetGtex_rsem_gene_tpm.txt' from Xena and load it into R
tcga <- read.delim("~/big-data/TcgaTargetGtex_rsem_gene_tpm.txt")

# Make Ensembl IDs rownames
tcga <- tcga %>% tibble::remove_rownames() %>% tibble::column_to_rownames(var = "sample")

# Clean Ensembl IDs (remove version numbers)
ensembl_ids <- gsub("\\..*", "", rownames(tcga))

# Map Ensembl IDs to HGNC gene symbols
gene_symbols <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                      keys = ensembl_ids,
                                      column = "SYMBOL",   
                                      keytype = "ENSEMBL", 
                                      multiVals = "first") 
tcga_anno <- cbind(gene = gene_symbols, tcga)

# Filter samples by Primary Tumor or Metastatic before next step
tcga_meta <- readRDS("~/rds/xena-tcga-gtex-meta.rds")
filtered_samples <- tcga_meta %>%
  dplyr::filter(sample_type == "Primary Tumor", # or Metastatic
         gender == "Female",
         study == "TCGA") %>%
  dplyr::pull(sample)  # extract the sample column as a vector


# tcga_anno colnames: "TCGA.S9.A7J2.01"
# tcga_anno sample names: "TCGA-S9-A7J2-01"
# Replace tcga_anno colnames "." with "-"
colnames(tcga_anno) <- stringr::str_replace_all(colnames(tcga_anno),
                                                pattern = "\\.",
                                                replacement = "-")

# Keep 'gene_name' column + filtered sample columns
tcga_filtered <- tcga_anno %>%
  dplyr::select(gene, any_of(filtered_samples))

# Adhesome genes and XIST
genes_of_interest <- c(adhesome$gene, "XIST")

# Define list of representative lncRNAs
lncRNAs <- c("XIST", "MALAT1", "AIRN", "HOTAIR", "H19", "NEAT1", "HOTTIP", "GAS5", "MEG3", "PVT1", "KCNQ1OT1", "CCAT1", "PCAT1", "UCA1")

# Filter my_data_filtered to keep only those genes
tcga_subset <- tcga_filtered %>%
  dplyr::filter(gene %in% genes_of_interest)

# Prepare RNAss for matching samples; TCGA phenotype downloaded from Xena
RNAss_df <- tcga_meta %>%
  dplyr::select(sample, RNAss)

common_samples <- intersect(colnames(tcga_subset)[-1], RNAss_df$sample)

expr_mat <- tcga_subset %>%
  dplyr::select(gene, all_of(common_samples))

RNAss_vec <- RNAss_df %>%
  dplyr::filter(sample %in% common_samples) %>%
  dplyr::arrange(match(sample, common_samples)) %>%
  dplyr::pull(RNAss)

# append RNAss as a new row
expr_with_RNAss <- dplyr::bind_rows(
  expr_mat,
  data.frame(gene_name = "RNAss", t(RNAss_vec)) %>% 
    setNames(colnames(expr_mat))
)

# Convert column to rownames and transpose
expr_long <- expr_with_RNAss %>%
  tibble::remove_rownames() %>%
  tibble::column_to_rownames(var = "gene") %>%
  t() %>% 
  as.data.frame()


saveRDS(expr_long, "~/rds/tcga-xist-adhesome-exp.rds")
