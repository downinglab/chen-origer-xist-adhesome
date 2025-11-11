# Libraries ---------------------------------------------------------------
library(ggplot2)
library(tidyverse)

# Figure 1c ---------------------------------------------------------------

xist.kd.degs <- read.csv("~/data/xist-kd-degs.csv", header = T)

xist.kd.degs <- merge(xist.kd.degs, adhesome, by = "gene", sort = F)
xist.deg.mat <- xist.kd.degs %>% select(c("gene","X7","X9")) %>% remove_rownames() %>% column_to_rownames(var = "gene") %>% as.matrix()
xist.cat.mat <- xist.kd.degs %>% select(c("gene","cad","int")) %>% remove_rownames() %>% column_to_rownames(var = "gene") %>% as.matrix()
colnames(xist.cat.mat) <- c("Cad.", "Int.")
xist.cat.mat <- xist.cat.mat[,c(2,1)]
xist.cat.mat <- ifelse(xist.cat.mat, yes = 1, no = 0)

library(ComplexHeatmap)
ht1 <- ComplexHeatmap::Heatmap(xist.deg.mat, 
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

ht2 <- ComplexHeatmap::Heatmap(xist.cat.mat, 
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

up_X7_GOBP <- read.delim("~/data/up_X7_GOBP.txt")
up_X9_GOBP <- read.delim("~/data/up_X9_GOBP.txt")
down_X9_GOBP <- read.delim("~/data/down_X9_GOBP.txt")

# Select top 20 rows from each dataset
up_X7_GOBP_top20 <- head(up_X7_GOBP, 20)
up_X9_GOBP_top20 <- head(up_X9_GOBP, 20)
down_X9_GOBP_top20 <- head(down_X9_GOBP, 20)

# Use gsub() to remove the pattern from 'Term'
up_X7_GOBP_top20$Term <- gsub("^GO[^~]*~", "", up_X7_GOBP_top20$Term)

# Capitalize terms
up_X7_GOBP_top20$Term <- stringr::str_to_title(up_X7_GOBP_top20$Term)

# Sort
up_X7_GOBP_top20 <- up_X7_GOBP_top20[order(up_X7_GOBP_top20$Fold.Enrichment, decreasing = T),]

# Ensure 'Term' is a factor with levels in the original order
up_X7_GOBP_top20$Term <- factor(up_X7_GOBP_top20$Term, levels = rev(up_X7_GOBP_top20$Term))

# Apply log transformation for better color scaling
up_X7_GOBP_top20$logPValue <- -log10(up_X7_GOBP_top20$PValue)

# Remove outlier
up_X7_GOBP_top20 <- up_X7_GOBP_top20[-1,]

# Create dot plot
ggplot(up_X7_GOBP_top20, aes(x = Term, y = Fold.Enrichment, color = logPValue, size = Count)) +
  geom_point() +
  scale_color_gradient(low = "white", high = "red", breaks = c(2:10), limits = c(2,10)) +  # Improved color contrast
  scale_size_continuous(breaks = c(50, 100, 150), limits = c(0, 150)) +
  scale_y_continuous(limit = c(1,9),
                     breaks = c(1,3,5,7,9)) + 
  labs(title = "Up DEG in X7",
       x = NULL,
       y = NULL,
       color = "-log10(P-value)",
       size = "Count") +
  coord_flip() +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        plot.title = element_text(size = 7))

# Use gsub() to remove the pattern from 'Term'
up_X9_GOBP_top20$Term <- gsub("^GO[^~]*~", "", up_X9_GOBP_top20$Term)

# Capitalize terms
up_X9_GOBP_top20$Term <- stringr::str_to_title(up_X9_GOBP_top20$Term)

# Sort
up_X9_GOBP_top20 <- up_X9_GOBP_top20[order(up_X9_GOBP_top20$Fold.Enrichment, decreasing = T),]

# Ensure 'Term' is a factor with levels in the original order
up_X9_GOBP_top20$Term <- factor(up_X9_GOBP_top20$Term, levels = rev(up_X9_GOBP_top20$Term))

# Apply log transformation for better color scaling
up_X9_GOBP_top20$logPValue <- -log10(up_X9_GOBP_top20$PValue)

# Create dot plot
ggplot(up_X9_GOBP_top20, aes(x = Term, y = Fold.Enrichment, color = logPValue, size = Count)) +
  geom_point() +
  scale_color_gradient(low = "white", high = "red", breaks = c(2:10), limits = c(2,10)) +  # Improved color contrast
  scale_size_continuous(breaks = c(50, 100, 150), limits = c(0, 150)) +
  scale_y_continuous(limit = c(1,9),
                     breaks = c(1,3,5,7,9)) + 
  labs(title = "Up DEG in X9",
       x = NULL,
       y = NULL,
       color = "-log10(P Value)",
       size = "Count") +
  coord_flip() +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        plot.title = element_text(size = 7))

# Use gsub() to remove the pattern from 'Term'
down_X9_GOBP_top20$Term <- gsub("^GO[^~]*~", "", down_X9_GOBP_top20$Term)

# Capitalize terms
down_X9_GOBP_top20$Term <- stringr::str_to_title(down_X9_GOBP_top20$Term)

# Sort
down_X9_GOBP_top20 <- down_X9_GOBP_top20[order(down_X9_GOBP_top20$Fold.Enrichment, decreasing = T),]

# Ensure 'Term' is a factor with levels in the original order
down_X9_GOBP_top20$Term <- factor(down_X9_GOBP_top20$Term, levels = rev(down_X9_GOBP_top20$Term))

# Apply log transformation for better color scaling
down_X9_GOBP_top20$logPValue <- -log10(down_X9_GOBP_top20$PValue)

# Create dot plot
ggplot(down_X9_GOBP_top20, aes(x = Term, y = Fold.Enrichment, color = logPValue, size = Count)) +
  geom_point() +
  scale_color_gradient(low = "white", high = "red", breaks = c(2:10), limits = c(2,10)) +  # Improved color contrast
  scale_size_continuous(breaks = c(50, 100, 150), limits = c(0, 150)) +
  scale_y_continuous(limit = c(1,9),
                     breaks = c(1,3,5,7,9)) + 
  labs(title = "down DEG in X9",
       x = NULL,
       y = NULL,
       color = "-log10(P Value)",
       size = "Count") +
  coord_flip() +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        plot.title = element_text(size = 7))

up_X7_GOBP_top20$Condition <- rep("up_X7", length(up_X7_GOBP_top20$Category))
up_X9_GOBP_top20$Condition <- rep("up_X9", length(up_X9_GOBP_top20$Category))
down_X9_GOBP_top20$Condition <- rep("down_X9", length(down_X9_GOBP_top20$Category))

up_X7_GOBP_top20 <- up_X7_GOBP_top20[order(up_X7_GOBP_top20$Fold.Enrichment, decreasing = T),]
up_X7_GOBP_top20$Term <- factor(up_X7_GOBP_top20$Term, levels = up_X7_GOBP_top20$Term)

up_X9_GOBP_top20 <- up_X9_GOBP_top20[order(up_X9_GOBP_top20$Fold.Enrichment, decreasing = T),]
up_X9_GOBP_top20$Term <- factor(up_X9_GOBP_top20$Term, levels = up_X9_GOBP_top20$Term)

down_X9_GOBP_top20 <- down_X9_GOBP_top20[order(down_X9_GOBP_top20$Fold.Enrichment, decreasing = T),]
down_X9_GOBP_top20$Term <- factor(down_X9_GOBP_top20$Term, levels = down_X9_GOBP_top20$Term)

combined_GOBP <- rbind(up_X7_GOBP_top20, up_X9_GOBP_top20, down_X9_GOBP_top20)

combined_GOBP$Condition <- factor(combined_GOBP$Condition, levels = c("up_X7", "up_X9", "down_X9"))

combined_GOBP$Term <- varhandle::unfactor(combined_GOBP$Term)

# Create dot plot
ggplot(combined_GOBP, aes(x = Fold.Enrichment, 
                          y = tidytext::reorder_within(Term, Fold.Enrichment, Condition), 
                          color = logPValue, size = Count)) +
  geom_point() +
  scale_color_gradient(low = "#fee5d9", high = "#FF0000", breaks = c(2:10), limits = c(1.5,10.5)) +  # Improved color contrast
  scale_size_continuous(breaks = c(50, 100, 150), limits = c(0, 150), range = c(1,4)) +
  scale_x_continuous(limit = c(1,9),
                     breaks = c(1,3,5,7,9)) + 
  labs(title = "Combined",
       x = NULL,
       y = NULL,
       color = "-log10(P Value)",
       size = "Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        plot.title = element_text(size = 7)) +
  facet_grid(Condition ~ ., 
             scales = "free_y")


# Supp Figure 1a/b ---------------------------------------------------------

xist.volcano <- read.csv("~/data/xist-kd-volcano.csv", row.names = F)

library(ggrepel)
ggplot(xist.volcano, aes(x = log2FoldChange, y = -log10(padj), col = deg, label = delabel)) +
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

# Load adhesome gene data
adhesome <- read.csv("~/data/adhesome-components.csv", header = T)

xist.kd.degs <- read.csv("~/data/xist-kd-volcano.csv", row.names = F)
colnames(xist.kd.degs)[1] <- "gene"
xist.kd.degs <- merge(xist.kd.degs, 
                      select(adhesome, "gene", "cat"),
                      by = "gene")

xist.kd.degs.both <- filter(xist.kd.degs, cat == "Both")
xist.kd.degs.cad <- xist.kd.degs.both %>% mutate(cat = "Cadherin")
xist.kd.degs.int <- xist.kd.degs.both %>% mutate(cat = "Integrin")

xist.kd.degs <- xist.kd.degs %>% filter(cat != "Both") %>% rbind(xist.kd.degs.cad, xist.kd.degs.int)

xist.kd.degs$deg <- stringr::str_to_title(xist.kd.degs$deg)
xist.kd.degs$cat <- factor(xist.kd.degs$cat, levels = c("Integrin", "Cadherin"))
ggplot(data = filter(xist.kd.degs, deg != "None" & cat != "Both"), 
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

# Figure 2a ---------------------------------------------------------------

library(ggridges)

lnc.corrs <- read.csv("~/data/lncRNA_adhesome_corrs.csv")

ggplot(lnc.corrs, aes(x = corr, y = lncRNA, fill = condition)) +
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

rtz <- read.csv("~/data/r-to-z-corrs.csv")

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

tissue.corrs <- read.csv("~/data/xist-adhesome-tissue-corrs.csv")

tissue.corrs$tissue <- factor(tissue.corrs$tissue, levels = rev(levels(tissue.corrs$tissue)))
ggplot(tissue.corrs, aes(x = corr, y = tissue, fill = study)) +
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

lnc.corrs <- read.csv("~/data/lncRNA-adhesome-corrs.csv")

gtex.long <- filter(lnc.corrs, condition == "GTEx")

gtex.long$lncRNA <- factor(gtex.long$lncRNA, levels = rownames(gtex.long %>% 
                                                                 dplyr::group_by(lncRNA) %>% 
                                                                 dplyr::summarise(mean = mean(corr)) %>% 
                                                                 tibble::column_to_rownames(var = "lncRNA") %>%
                                                                 dplyr::arrange(mean)))

ggplot(gtex.long, aes(x = corr, y = lncRNA, fill = condition)) + #fill = after_stat(x)
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

tcga.long <- filter(lnc.corrs, condition == "TCGA")

tcga.long$lncRNA <- factor(tcga.long$lncRNA, levels = rownames(tcga.long %>% 
                                                                 dplyr::group_by(lncRNA) %>% 
                                                                 dplyr::summarise(mean = mean(corr)) %>% 
                                                                 tibble::column_to_rownames(var = "lncRNA") %>%
                                                                 dplyr::arrange(mean)))

ggplot(tcga.long, aes(x = corr, y = lncRNA, fill = condition)) + #fill = after_stat(x)
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


# Supp Figure c/d ---------------------------------------------------------

gtex <- read.csv("~/data/xist-adhesome-tissue-corrs-gtex.csv")
tcga <- read.csv("~/data/xist-adhesome-tissue-corrs-tcga.csv")

# Make more consistent
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


# Supp Figure 2e ----------------------------------------------------------

# Data import and cleanup
mgi <- read.csv("~/data/lncRNA-adhesome-corrs-mouse.csv")
mgi <- mgi %>% mutate(study = "MGI")

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

nucl.corrs <- read.csv("~/data/xist-adhesome-corrs-nuclear.csv")

# Factor for plots
combined$sample_type <- factor(combined$sample_type, levels = c("Normal Tissue", "Primary Tumor"))
combined$nucl_bin <- factor(combined$nucl_bin, levels = c("Nuclear-localized", "Non-nuclear-localized"))
ggplot(combined, aes(x = nucl_bin, y = bicor, fill = sample_type)) + 
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

# Load adhesome-XIST correlations from Figure 2a
wald <- read.csv("~/data/adhesome-coxph-nuclear.csv")
wald <- wald %>% mutate(sample_type = "Primary Tumor")

# Factor for plots
wald$nucl_bin <- factor(wald$nucl_bin, levels = c("Nuclear-localized", "Non-nuclear-localized"))

ggplot(wald, aes(x = nucl_bin, y = abs_log2_ratio, fill = sample_type)) + 
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

library(survival)
library(survminer)

# Load survival data
xena.survival <- read.delim("~/data/xena-survival.txt")
xena <- readRDS("~/data/xena-subset.rds")

get_dat_wide <- function(ct, genes) {
  # Subset expression/meta for one cancer with filters
  dat <- xena %>%
    filter(
      cancer_type %in% ct,
      gene %in% genes,
      if (!is.null(restrict_gender)) gender == restrict_gender else TRUE,
      if (!is.null(restrict_sample_type)) sample_type == restrict_sample_type else TRUE
    ) %>%
    select(sample, cancer_type, gene, expression_value)
  
  if (nrow(dat) == 0) return(NULL)
  
  # If some genes are entirely missing in this cancer, we still proceed for the present ones.
  present_genes <- intersect(genes, unique(dat$gene))
  if (length(present_genes) == 0) return(NULL)
  
  # Per-gene quartiles (Q1 and Q3) within this cancer
  th <- dat %>%
    group_by(gene) %>%
    summarise(
      q1 = quantile(expression_value, low_quantile,  na.rm = TRUE),
      q3 = quantile(expression_value, high_quantile, na.rm = TRUE),
      .groups = "drop"
    )
  q1 <- setNames(th$q1, th$gene)
  q3 <- setNames(th$q3, th$gene)
  
  # Pivot to wide to align samples across genes
  dat_wide <- dat %>%
    distinct(sample, cancer_type, gene, expression_value) %>%
    pivot_wider(names_from = gene, values_from = expression_value)
  
  # Construct per-gene High/Low bins (middle = NA)
  for (g in present_genes) {
    dat_wide[[paste0("bin_", g)]] <- ifelse(
      dat_wide[[g]] >= q3[g], paste0("High ", g),
      ifelse(dat_wide[[g]] <= q1[g], paste0("Low ", g), NA)
    )
  }
  
  # Attach survival
  dat_wide <- dat_wide %>%
    left_join(select(xena.survival, sample, OS, OS.time), by = "sample") %>%
    filter(!is.na(OS), !is.na(OS.time))
  
  if (nrow(dat_wide) == 0) return(NULL)
  dat_wide
}

# Network genes
genes_of_interest <- c("XIST","CTNNB1","CSNK1E","VEZT","SORBS3","ABL1","DOCK1","PLCG1","INPPL1","MAPK8","FER","FLNB","SH2B1","NISCH","ITGB3BP")

# Load adhesome gene data
adhesome <- read.csv("~/data/adhesome-components.csv", header = T)

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
xena.subset <- xena.subset %>% select(bin_XIST, OS, OS.time, any_of(c(adhesome$gene,"XIST")))

# Create cox models for high and low XIST bins
{
  cox.res.high <- coxph(formula = Surv(OS.time, OS) ~ XIST+ACTN1+CFL1+CORO1B+CTTN+FLNA+KEAP1+LASP1+ENAH+NEXN+SVIL+VASP+CORO2A+ACTB+SORBS2+BCAR1+CAV1+SMPX+SH3KBP1+CRK+CRKL+EZR+FHL2+GAB1+GRB2+GRB7+HAX1+NEDD9+CASS4+TGFB1I1+ITGB1BP1+FERMT1+FERMT2+FERMT3+LPXN+PPFIA1+LPP+NF2+FBLIM1+MSN+NCK2+PALLD+PARVA+PARVB+PXN+LIMS1+LIMS2+RDX+OSTF1+NUDT16L1+SYNM+SDCBP+TLN1+TNS1+TES+TRIP6+VCL+SORBS3+LDB3+ZYX+NDEL1+SH2B1+ZFYVE21+SLC3A2+KTN1+LRP1+PVR+SDC4+ITGA1+ITGA2+ITGA3+ITGA4+ITGA5+ITGA6+ITGA7+ITGA8+ITGA9+ITGA10+ITGA11+ITGAD+ITGAE+ITGAL+ITGAM+ITGAV+ITGAX+ITGB1+ITGB2+ITGB3+ITGB4+ITGB5+ITGB6+ITGB7+ITGB8+NRP1+NRP2+CD151+PDE4D+PKD1+TRPM7+CALR+ASAP3+GIT1+GIT2+ARHGAP26+ASAP2+ARHGAP24+DLC1+AGAP2+STARD13+RASA1+DEF6+DOCK1+ELMO1+ARHGEF6+ARHGEF7+DNM2+RHOU+PLCG1+INPP5D+INPPL1+ITGB3BP+RAVER1+STAT3+SPTLC1+ILK+PAK1+PDPK1+PRKCA+PPM1M+PPM1F+PPP2CA+ABL1+CSK+PTK2+PTK2B+SRC+PEAK1+PTPRF+PTPN12+PTPRA+PTPN6+PTPN11+PRNP+MYH9+MACF1+ARPC2+MARCKS+PFN1+ABI1+ABI2 +ABI3+ANKRD28+CSRP1+IRS1+MAPK8IP3+PLEC+SORBS1+SHC1+MYOM1+TSPAN1+TUBA1B+VIM+SHARPIN+FABP3+NISCH+CIB1+CIB2+SRCIN1+ADAM12+CEACAM1+ENG+CD47+LAYN+SIRPA+THY1+PLAUR+KCNH2+SLC16A3+SLC9A1+HSPB1+HSPA2+CBL+RNF5+RNF185+ARHGAP5+ARHGAP32+BCAR3+RAPGEF1+SOS1+TIAM1+TRIO+VAV1+VAV2+VAV3+ARHGEF12+ARHGEF2+CYTH2+ARF1+HRAS+RAC1+RHOA+PLD1+CAPN1+CAPN2+CASP8+MMP14+PIK3CA+PIP5K1C+PTEN+PABPC1+SSH1+AKT1+MAPK1+MAPK8+LIMK1+PRKACA+ROCK1+ILKAP+INSR+FYN+LYN+SYK+TESK1+PTPN1+PTPRO+PTPRH+PTPN2+ACTN4+ADD1+AKAP5+CTNNA1+CTNNA2+CTNNA3+DBN1+EPB41+FLNB+FSCN1+LIMA1+SHROOM3+TJP1+TRIOBP+ABI2+ANK3+ARVCF+CD2AP+CGNL1+CLIP1+CTNNB1+CTNND1+CTNND2+DLG1+DLG5+IQGAP1+JUP+KRIT1+LIN7A+LMO7+MAGI1+MPP7+NUMB+PARD3+PDZD2+PKP4+PLEKHA7+SCRIB+SPTBN1+SSX2IP+TJP2+WIPF2+ACTR3+CAPZA1+DIAPH1+WASF2+WASL+CDH1+CDH10+CDH11+CDH12+CDH15+CDH18+CDH19+CDH2+CDH20+CDH22+CDH24+CDH3+CDH4+CDH5+CDH6+CDH7+CDH8+CDH9+FAT1+SDC1+AJAP1+HEG1+VEZT+TRPC4+TRPV4+DYNC1H1+MYH10+MYO6+MYO7A+ARF6+CDC42+RAP1A+SIPA1L3+RAPGEF2+ARHGAP12+ARHGAP21+ARHGAP35+RACGAP1+SH3BP1+ECT2+CSNK1E+CSNK2A1+PRKCD+PRKD1+EGFR+MET+FER+ACP1+PTPN14+PTPRJ+PTPRK+PTPRM+PTPRT+PTPRU+ADAM10+ADAM9+CASP3+PSEN1+CBLL1+GNA12+EXOC3+NME1, 
                        data = as.data.frame(filter(xena.subset, bin_XIST == "High XIST"))) # Evaluate among cancer types
  cox.high <- cox.res.high %>%
    broom::tidy() %>% 
    mutate(hr = exp(estimate)) %>% 
    select(term,p.value,hr)
  colnames(cox.high) <- paste0(colnames(cox.high), "_high")
  colnames(cox.high)[1] <- "gene"
  
  cox.res.low <- coxph(formula = Surv(OS.time, OS) ~ XIST+ACTN1+CFL1+CORO1B+CTTN+FLNA+KEAP1+LASP1+ENAH+NEXN+SVIL+VASP+CORO2A+ACTB+SORBS2+BCAR1+CAV1+SMPX+SH3KBP1+CRK+CRKL+EZR+FHL2+GAB1+GRB2+GRB7+HAX1+NEDD9+CASS4+TGFB1I1+ITGB1BP1+FERMT1+FERMT2+FERMT3+LPXN+PPFIA1+LPP+NF2+FBLIM1+MSN+NCK2+PALLD+PARVA+PARVB+PXN+LIMS1+LIMS2+RDX+OSTF1+NUDT16L1+SYNM+SDCBP+TLN1+TNS1+TES+TRIP6+VCL+SORBS3+LDB3+ZYX+NDEL1+SH2B1+ZFYVE21+SLC3A2+KTN1+LRP1+PVR+SDC4+ITGA1+ITGA2+ITGA3+ITGA4+ITGA5+ITGA6+ITGA7+ITGA8+ITGA9+ITGA10+ITGA11+ITGAD+ITGAE+ITGAL+ITGAM+ITGAV+ITGAX+ITGB1+ITGB2+ITGB3+ITGB4+ITGB5+ITGB6+ITGB7+ITGB8+NRP1+NRP2+CD151+PDE4D+PKD1+TRPM7+CALR+ASAP3+GIT1+GIT2+ARHGAP26+ASAP2+ARHGAP24+DLC1+AGAP2+STARD13+RASA1+DEF6+DOCK1+ELMO1+ARHGEF6+ARHGEF7+DNM2+RHOU+PLCG1+INPP5D+INPPL1+ITGB3BP+RAVER1+STAT3+SPTLC1+ILK+PAK1+PDPK1+PRKCA+PPM1M+PPM1F+PPP2CA+ABL1+CSK+PTK2+PTK2B+SRC+PEAK1+PTPRF+PTPN12+PTPRA+PTPN6+PTPN11+PRNP+MYH9+MACF1+ARPC2+MARCKS+PFN1+ABI1+ABI2 +ABI3+ANKRD28+CSRP1+IRS1+MAPK8IP3+PLEC+SORBS1+SHC1+MYOM1+TSPAN1+TUBA1B+VIM+SHARPIN+FABP3+NISCH+CIB1+CIB2+SRCIN1+ADAM12+CEACAM1+ENG+CD47+LAYN+SIRPA+THY1+PLAUR+KCNH2+SLC16A3+SLC9A1+HSPB1+HSPA2+CBL+RNF5+RNF185+ARHGAP5+ARHGAP32+BCAR3+RAPGEF1+SOS1+TIAM1+TRIO+VAV1+VAV2+VAV3+ARHGEF12+ARHGEF2+CYTH2+ARF1+HRAS+RAC1+RHOA+PLD1+CAPN1+CAPN2+CASP8+MMP14+PIK3CA+PIP5K1C+PTEN+PABPC1+SSH1+AKT1+MAPK1+MAPK8+LIMK1+PRKACA+ROCK1+ILKAP+INSR+FYN+LYN+SYK+TESK1+PTPN1+PTPRO+PTPRH+PTPN2+ACTN4+ADD1+AKAP5+CTNNA1+CTNNA2+CTNNA3+DBN1+EPB41+FLNB+FSCN1+LIMA1+SHROOM3+TJP1+TRIOBP+ABI2+ANK3+ARVCF+CD2AP+CGNL1+CLIP1+CTNNB1+CTNND1+CTNND2+DLG1+DLG5+IQGAP1+JUP+KRIT1+LIN7A+LMO7+MAGI1+MPP7+NUMB+PARD3+PDZD2+PKP4+PLEKHA7+SCRIB+SPTBN1+SSX2IP+TJP2+WIPF2+ACTR3+CAPZA1+DIAPH1+WASF2+WASL+CDH1+CDH10+CDH11+CDH12+CDH15+CDH18+CDH19+CDH2+CDH20+CDH22+CDH24+CDH3+CDH4+CDH5+CDH6+CDH7+CDH8+CDH9+FAT1+SDC1+AJAP1+HEG1+VEZT+TRPC4+TRPV4+DYNC1H1+MYH10+MYO6+MYO7A+ARF6+CDC42+RAP1A+SIPA1L3+RAPGEF2+ARHGAP12+ARHGAP21+ARHGAP35+RACGAP1+SH3BP1+ECT2+CSNK1E+CSNK2A1+PRKCD+PRKD1+EGFR+MET+FER+ACP1+PTPN14+PTPRJ+PTPRK+PTPRM+PTPRT+PTPRU+ADAM10+ADAM9+CASP3+PSEN1+CBLL1+GNA12+EXOC3+NME1, 
                       data = as.data.frame(filter(xena.subset, bin_XIST == "Low XIST"))) # Evaluate among cancer types
  cox.low <- cox.res.low %>%
    broom::tidy() %>% 
    mutate(hr = exp(estimate)) %>% 
    select(term,p.value,hr)
  colnames(cox.low) <- paste0(colnames(cox.low), "_low")
  colnames(cox.low)[1] <- "gene"
}

# Figure 4a ---------------------------------------------------------------

# Load data
xist.rnass <- read.csv("~/data/xist-rnass.csv")

# Factor for plots
xist.rnass$sample_type <- factor(xist.rnass$sample_type, levels = c("Normal Tissue", "Primary Tumor", "Metastatic"))
ggplot(xist.rnass, aes(x = XIST, y = RNAss, color = sample_type)) + 
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

# build model then plot
mod.normal <- lm(RNAss ~ XIST, data = filter(xist.rnass, sample_type == "Normal Tissue"))
summary(mod.normal)$r.squared
mod.cancer <- lm(RNAss ~ XIST, data = filter(xist.rnass, sample_type == "Primary Tumor"))
summary(mod.cancer)$r.squared
mod.met <- lm(RNAss ~ XIST, data = filter(xist.rnass, sample_type == "Metastatic"))
summary(mod.met)$r.squared

# Figure 4b ---------------------------------------------------------------

xist.bin.rnass <- read.csv("~/data/xist-bin-rnass.csv")

# Factor for plots
xist.bin.rnass$sample_type <- factor(xist.bin.rnass$sample_type, levels = c("Normal Tissue", "Primary Tumor", "Metastatic"))
xist.bin.rnass$XIST_group <- factor(xist.bin.rnass$XIST_group, levels = c("Low XIST", "High XIST"))
ggplot(xist.bin.rnass, aes(x = XIST_group, y = RNAss, fill = sample_type)) + 
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
  labs(title = NULL, x = "XIST", y = "RNAss", color = "Sample Type") +
  facet_grid(. ~ sample_type) +
  ggpubr::stat_compare_means(comparisons = list(c("Low XIST", "High XIST")),
                             label = "p.signif")

# Figure 4c ---------------------------------------------------------------

adhesome.rnass <- read.csv("~/data/adhesome-rnass.csv")

adhesome.rnass.medians <- adhesome.rnass %>%
  group_by(sample_type) %>%
  summarize(median = median(bicor))

# Factor for plots
adhesome.rnass$sample_type <- factor(adhesome.rnass$sample_type, levels = c("Normal Tissue", "Primary Tumor", "Metastatic"))
ggplot(adhesome.rnass, aes(x = bicor, fill = sample_type)) + 
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
  geom_vline(data = adhesome.rnass.medians, aes(xintercept = median),
             color = "black", linetype = "dashed", linewidth = 0.75)


# Figure 4d ---------------------------------------------------------------

xist.bin.adhesome.rnass <- read.csv("~/data/xist-bin-adhesome-rnass.csv")

# Factor for plots
xist.bin.adhesome.rnass$sample_type <- factor(xist.bin.adhesome.rnass$sample_type, levels = c("Normal Tissue", "Primary Tumor", "Metastatic"))
xist.bin.adhesome.rnass$XIST_group <- factor(xist.bin.adhesome.rnass$XIST_group, levels = c("Low XIST", "High XIST"))
ggplot(xist.bin.adhesome.rnass, aes(x = XIST_group, y = bicor, fill = sample_type)) + 
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
