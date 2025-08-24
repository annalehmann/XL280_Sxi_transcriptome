library(tidyverse)
library(dplyr)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)
library(ggtext)
library(VennDiagram)
library(matrixStats)
library(RColorBrewer)

setwd("")

# bring in counts files 
counts <- list.files(pattern = "featurecounts.txt")

# simplify sample names
col_name <- gsub("(JH-V8_72h_|_featurecounts.txt)", "", counts)

# empty list to store df names
df_list <- list()

# read in counts files and keep gene ids and count values only
for (i in seq_along(counts)) {
  df <- read.delim(counts[i], 
                   skip = 1,
                   sep = "\t",
                   header = TRUE) %>%
    select(1,7) 
  colnames(df) <- c("gene_ID", col_name[i])
  df$gene_ID <- sub("gene-", "", df$gene_ID)
  df_list[[i]] <- df
}
names(df_list) <- col_name

# get gff info for later tpm calcs
gff <- read.delim("PATH/XL280_liftmerged.gff", 
                  header = FALSE, 
                  col.names = c("contig", "source", "cat", "start", "end", ".", "strand", ".", "ID"),
                  sep = "\t", 
                  skip = 3)
gff_genes <- gff[gff$cat == "gene", ]
gff_genes$gene <- gsub("(ID=gene-|;.*)", "", gff_genes$ID)
gene_lengths <- gff_genes[, c("gene", "start", "end")]
gene_lengths$length <- gene_lengths$end - gene_lengths$start
row.names(gene_lengths) <- gene_lengths$gene

# MATalpha unisex colnames
alpha_unisex <- c("XL280alpha_rep1", "XL280alpha_rep2", "XL280alpha_rep3", 
                  "JHG4_rep1", "JHG4_rep2", "JHG4_rep3")
## data wrangling
### merge counts dfs
alpha_unisex_counts = Reduce(function(x, y) 
  merge(x, y, by = "gene_ID", all = TRUE), 
  df_list[names(df_list) %in% alpha_unisex])

### make gene ids into row names
rownames(alpha_unisex_counts) <- alpha_unisex_counts$gene_ID
alpha_unisex_counts <- alpha_unisex_counts[,-1]

## TPM filter
#gene_counts_by_ID$gene <- row.names(gene_counts_by_ID)
alpha_unisex_counts$gene <- row.names(alpha_unisex_counts)
tpm_counts_alpha_unisex <- merge(alpha_unisex_counts, gene_lengths[ , c("gene", "length")], by = "gene")
## calculate tpms
### convert length from bp to kbp
tpm_counts_alpha_unisex$length <- tpm_counts_alpha_unisex$length/1000
### df/character vector to hold tpm values
tpm_calc_alpha_unisex <- tpm_counts_alpha_unisex[, c("gene"), drop = FALSE]
### get cols except gene
tpm_cols_alpha_unisex <- setdiff(colnames(tpm_counts_alpha_unisex), c("gene", "length"))
### get tpm
for (gene in tpm_cols_alpha_unisex) {
  rpk <- tpm_counts_alpha_unisex[[gene]] / tpm_counts_alpha_unisex$length
  tpm <- rpk / sum(rpk) * 1e6
  tpm_calc_alpha_unisex[[gene]] <- tpm
}
## average tpm for technical reps
tpm_calc_alpha_unisex <- tpm_calc_alpha_unisex %>%
  pivot_longer(
    cols = -gene,
    names_to = "replicate",
    values_to = "tpm"
  )

tpm_calc_alpha_unisex_long <- tpm_calc_alpha_unisex %>%
  mutate(sample = str_extract(replicate, "^[^_]+_[^_]+"))

tpm_avg_alpha_unisex <- tpm_calc_alpha_unisex_long %>%
  group_by(gene, sample) %>%
  summarise(avg = mean(tpm), .groups = "drop") %>%
  pivot_wider(names_from = "sample", values_from = "avg")

tpm_filtered_genes_alpha_unisex <- tpm_avg_alpha_unisex %>%
  filter(if_any(-gene, ~ . >= 1)) %>%
  pull(gene)

alpha_unisex_counts_tpmfiltered <- alpha_unisex_counts[rownames(alpha_unisex_counts) %in% tpm_filtered_genes_alpha_unisex, ]
# remove last column of gene ids
alpha_unisex_counts_tpmfiltered <- alpha_unisex_counts_tpmfiltered[ ,1:6]

## DEG analysis
### make df with sample info
sample_unisex_alpha <- c(colnames(alpha_unisex_counts_tpmfiltered))
condition_solo <- c('mut', 'mut', 'mut', 'wt', 'wt', 'wt')
colData_unisex_alpha <- data.frame(sample_unisex_alpha, condition_solo)

### run deseq
dds_unisex_alpha <- DESeqDataSetFromMatrix(countData = alpha_unisex_counts_tpmfiltered,
                                           colData = colData_unisex_alpha,
                                           design = ~ condition_solo)

dds_unisex_alpha <- DESeq(dds_unisex_alpha)

#### EDA - PCA
trans_unisex_alpha <- vst(dds_unisex_alpha, blind = FALSE)
##### extract matrix and get variances
genevar_unisex_alpha <- rowVars(assay(trans_unisex_alpha))

genevar_unisex_alpha_plot <- plot(sort(genevar_unisex_alpha, decreasing = TRUE), 
                                  type = "l", 
                                  ylab = "Variance", 
                                  xlab = "Gene")

# plot on screen just to see variances
plotPCA(trans_unisex_alpha, intgroup = c("condition_solo")) # PC1: 70%; PC2: 21%

# plot for fig
pca_unisex_alpha_data <- plotPCA(trans_unisex_alpha, intgroup = c("condition_solo"), returnData=TRUE)
pca_unisex_alpha <- ggplot(pca_unisex_alpha_data, 
                           aes(PC1, PC2, color=condition_solo, shape=condition_solo)) + 
  geom_richtext(data = pca_unisex_alpha_data, 
                aes(label = condition_solo),
                color = "black") +
  geom_point(size = 1.5, shape = c(22, 22, 22, 21, 21, 21),  fill = "white", stroke = 1) +
  scale_color_manual(values=c("#edaaa8","#d8d8d8")) +
  theme(panel.grid.major = element_line(color = "#eeeeee", linewidth = 0.3), 
        panel.grid.minor = element_line(color = "#eeeeee", linewidth = 0.1), 
        panel.background = element_blank(), 
        axis.line = element_line(color = "black"),
        panel.border = element_rect(color = "black", 
                                    fill=NA, 
                                    linewidth = 0.75),
        legend.position = "none") +
  xlab("PC1: 70%") +
  ylab("PC2: 21%") 

ggsave(file="pca_unisex_alpha.svg", plot=pca_unisex_alpha, width=2, height=2, units= "in")

# get vst genes for heatmap
mat_unisex_alpha <- assay(trans_unisex_alpha)
matrix_unisex_alpha <- SummarizedExperiment::assay(trans_unisex_alpha)
topVarGenes_unisex_alpha <- head(order(rowVars(matrix_unisex_alpha), decreasing = TRUE),500)
topGenesMatrix_unisex_alpha <- matrix_unisex_alpha[topVarGenes_unisex_alpha, ]

heatmap_unisex_alpha <- pheatmap(topGenesMatrix_unisex_alpha,
                                 scale = "row",
                                 cluster_rows = TRUE,
                                 cluster_cols = TRUE,
                                 show_rownames = FALSE) # JHG4_rep3 clusters away from other 2 reps

ggsave(file="heatmap_unisex_alpha.svg", plot=heatmap_unisex_alpha, width=5, height=5)

#### dispersion plot
dispersion_unisex_alpha <- plotDispEsts(dds_unisex_alpha)

### remove JHG4_rep3
dds_unisex_alpha_cleaned <- dds_unisex_alpha[, colnames(dds_unisex_alpha) !="JHG4_rep3"]
dds_unisex_alpha_cleaned$condition_solo <- droplevels(dds_unisex_alpha_cleaned$condition_solo)
#### re-run deseq
dds_unisex_alpha_cleaned <- DESeq(dds_unisex_alpha_cleaned)

#### PCA
trans_unisex_alpha_cleaned <- vst(dds_unisex_alpha_cleaned, blind = FALSE)
##### extract matrix and get variances
genevar_unisex_alpha_cleaned <- rowVars(assay(trans_unisex_alpha_cleaned))

plot(sort(genevar_unisex_alpha_cleaned, decreasing = TRUE), 
     type = "l", 
     ylab = "Variance", 
     xlab = "Gene")

pca_unisex_alpha_cleaned_data <- plotPCA(trans_unisex_alpha_cleaned, intgroup = c("condition_solo"), returnData=TRUE) # PC1: 76%; PC2: 14%
pca_unisex_alpha_cleaned <- ggplot(pca_unisex_alpha_cleaned_data, 
                                   aes(PC1, PC2, color=condition_solo, shape=condition_solo)) + 
  geom_richtext(data = pca_unisex_alpha_cleaned_data, 
                aes(label = condition_solo),
                color = "black") +
  geom_point(size = 1.5, shape = c(22, 22, 21, 21, 21),  fill = "white", stroke = 1) +
  scale_color_manual(values=c("#edaaa8","#d8d8d8")) +
  theme(panel.grid.major = element_line(color = "#eeeeee", linewidth = 0.3), 
        panel.grid.minor = element_line(color = "#eeeeee", linewidth = 0.1), 
        panel.background = element_blank(), 
        axis.line = element_line(color = "black"),
        panel.border = element_rect(color = "black", 
                                    fill=NA, 
                                    linewidth = 0.75),
        legend.position = "none") +
  xlab("PC1: 76%") +
  ylab("PC2: 14%") 
ggsave(file="pca_unisex_alpha_cleaned.svg", plot=pca_unisex_alpha_cleaned, width=2, height=2, units= "in")

mat_unisex_alpha_cleaned <- assay(trans_unisex_alpha_cleaned)
matrix_unisex_alpha_cleaned <- SummarizedExperiment::assay(trans_unisex_alpha_cleaned)
topVarGenes_unisex_alpha_cleaned <- head(order(rowVars(matrix_unisex_alpha_cleaned), decreasing = TRUE),500)
topGenesMatrix_unisex_alpha_cleaned <- matrix_unisex_alpha_cleaned[topVarGenes_unisex_alpha_cleaned, ]

heatmap_unisex_alpha_cleaned <- pheatmap(topGenesMatrix_unisex_alpha_cleaned,
                                         scale = "row",
                                         cluster_rows = TRUE,
                                         cluster_cols = TRUE,
                                         show_rownames = FALSE,
) 
ggsave(file="heatmap_unisex_alpha_cleaned.png", plot=heatmap_unisex_alpha_cleaned, width=5, height=5)

#### dispersion plot
dispersion_unisex_alpha_cleaned <- plotDispEsts(dds_unisex_alpha_cleaned)

### compare mut to wt
#dds_unisex_alpha_cleaned$condition_solo <- relevel(dds_unisex_alpha_cleaned$condition_solo, ref = "wt")
#dds_unisex_alpha_cleaned <- DESeq(dds_unisex_alpha)

res_unisex_alpha <- results(dds_unisex_alpha_cleaned,
                            contrast = c("condition_solo", "mut",
                                         "wt"),
                            alpha = 0.05)

res_unisex_alpha <- as.data.frame(res_unisex_alpha)
# remove degs that are artifacts (SXI1, FAO1, SXI2)
artifacts <- c("CND05950", "CND05960", "AF542530.2:10466..12990")
res_unisex_alpha <- res_unisex_alpha[!rownames(res_unisex_alpha) %in% artifacts, ]

# save data
write.csv(res_unisex_alpha, "unisex_alpha.csv")

### make volcano plot (fig. 5A)
unisex_alpha_volcano <- EnhancedVolcano(res_unisex_alpha,
                                        lab = rep(NA, nrow(res_unisex_alpha)),
                                        x = 'log2FoldChange',
                                        y = 'padj',
                                        pCutoff = 0.05,
                                        FCcutoff = 1,
                                        cutoffLineType = 'dashed',
                                        title = NULL,
                                        subtitle = NULL,
                                        caption = NULL,
                                        legendLabels = NULL,
                                        legendPosition = 'none',
                                        axisLabSize = 6,
                                        pointSize = 0.5,
                                        gridlines.major = 0.03,
                                        gridlines.minor = 0.01,
                                        col=c('#dad9d9', '#dad9d9', 'grey', '#edaaa8'),
                                        colAlpha = 1,
                                        xlim = c(-5, 5),
                                        ylim = c(-0.5, 250))

ggsave(file="unisex_alpha_volcano.svg", plot=unisex_alpha_volcano, width=2, height=2, units = "in")


# MATa unisex
a_unisex <- c("XL280a_rep1", "XL280a_rep2", "XL280a_rep3", 
              "JHG6_rep1", "JHG6_rep2", "JHG6_rep3")
## data wrangling
### merge counts dfs
a_unisex_counts = Reduce(function(x, y) 
  merge(x, y, by = "gene_ID", all = TRUE), 
  df_list[names(df_list) %in% a_unisex])

### make gene ids into row names
rownames(a_unisex_counts) <- a_unisex_counts$gene_ID
a_unisex_counts <- a_unisex_counts[,-1]

## TPM filter
#gene_counts_by_ID$gene <- row.names(gene_counts_by_ID)
a_unisex_counts$gene <- row.names(a_unisex_counts)
# get rid of "rna-" from row names
a_unisex_counts$gene <- gsub("rna-", "", a_unisex_counts$gene)
rownames(a_unisex_counts) <- a_unisex_counts$gene
tpm_counts_a_unisex <- merge(a_unisex_counts, gene_lengths[ , c("gene", "length")], by = "gene")
## calculate tpms
### convert length from bp to kbp
tpm_counts_a_unisex$length <- tpm_counts_a_unisex$length/1000
### df/character vector to hold tpm values
tpm_calc_a_unisex <- tpm_counts_a_unisex[, c("gene"), drop = FALSE]
### get cols except gene
tpm_cols_a_unisex <- setdiff(colnames(tpm_counts_a_unisex), c("gene", "length"))
### get tpm
for (gene in tpm_cols_a_unisex) {
  rpk <- tpm_counts_a_unisex[[gene]] / tpm_counts_a_unisex$length
  tpm <- rpk / sum(rpk) * 1e6
  tpm_calc_a_unisex[[gene]] <- tpm
}
## average tpm for technical reps
tpm_calc_a_unisex <- tpm_calc_a_unisex %>%
  pivot_longer(
    cols = -gene,
    names_to = "replicate",
    values_to = "tpm"
  )

tpm_calc_a_unisex_long <- tpm_calc_a_unisex %>%
  mutate(sample = str_extract(replicate, "^[^_]+_[^_]+"))

tpm_avg_a_unisex <- tpm_calc_a_unisex_long %>%
  group_by(gene, sample) %>%
  summarise(avg = mean(tpm), .groups = "drop") %>%
  pivot_wider(names_from = "sample", values_from = "avg")

tpm_filtered_genes_a_unisex <- tpm_avg_a_unisex %>%
  filter(if_any(-gene, ~ . >= 1)) %>%
  pull(gene)

a_unisex_counts_tpmfiltered <- a_unisex_counts[rownames(a_unisex_counts) %in% tpm_filtered_genes_a_unisex, ]
# remove last column of gene ids
a_unisex_counts_tpmfiltered <- a_unisex_counts_tpmfiltered[ ,1:6]

## DEG analysis
### make df with appropriate sample info
sample_unisex_a <- c(colnames(a_unisex_counts_tpmfiltered))

colData_unisex_a <- data.frame(sample_unisex_a, condition_solo)

### run deseq
dds_unisex_a <- DESeqDataSetFromMatrix(countData = a_unisex_counts_tpmfiltered,
                                       colData = colData_unisex_a,
                                       design = ~ condition_solo)

#### check quality
trans_unisex_a <- vst(dds_unisex_a, blind = FALSE)
plotPCA(trans_unisex_a, intgroup = c("condition_solo")) # PC1: 53%; PC2: 28%
pca_unisex_a_data <- plotPCA(trans_unisex_a, intgroup = c("condition_solo"), returnData=TRUE) 
pca_unisex_a <- ggplot(pca_unisex_a_data, 
                       aes(PC1, PC2, color=condition_solo, shape=condition_solo)) + 
  geom_richtext(data = pca_unisex_a_data, 
                aes(label = condition_solo),
                color = "black") +
  geom_point(size = 1.5, shape = c(24, 24, 24, 21, 21, 21),  fill = "white", stroke = 1) +
  scale_color_manual(values=c("#b9eefb","#b2b4b6")) +
  theme(panel.grid.major = element_line(color = "#eeeeee", linewidth = 0.3), 
        panel.grid.minor = element_line(color = "#eeeeee", linewidth = 0.1), 
        panel.background = element_blank(), 
        axis.line = element_line(color = "black"),
        panel.border = element_rect(color = "black", 
                                    fill=NA, 
                                    linewidth = 0.75),
        legend.position = "none") +
  xlab("PC1: 53%") +
  ylab("PC2: 28%") 
ggsave(file="pca_unisex_a.svg", plot=pca_unisex_a, width=2, height=2, units = "in")

mat_unisex_a <- assay(trans_unisex_a)
matrix_unisex_a <- SummarizedExperiment::assay(trans_unisex_a)
topVarGenes_unisex_a <- head(order(rowVars(matrix_unisex_a), decreasing = TRUE),500)
topGenesMatrix_unisex_a <- matrix_unisex_a[topVarGenes_unisex_a, ]

heatmap_unisex_a <- pheatmap(topGenesMatrix_unisex_a,
                             scale = "row",
                             cluster_rows = TRUE,
                             cluster_cols = TRUE,
                             show_rownames = FALSE,
) # clustering checks out
ggsave(file="heatmap_unisex_a.png", plot=heatmap_unisex_a, width=5, height=5)

### compare mut to wt
dds_unisex_a$condition_solo <- relevel(dds_unisex_a$condition_solo, ref = "wt")
dds_unisex_a <- DESeq(dds_unisex_a)

res_unisex_a <- results(dds_unisex_a,
                        contrast = c("condition_solo", "mut",
                                     "wt"),
                        alpha = 0.05)

res_unisex_a <- as.data.frame(res_unisex_a)
res_unisex_a <- res_unisex_a[!rownames(res_unisex_a) %in% artifacts, ]
# save data
write.csv(res_unisex_a, "unisex_a.csv")

### make volcano plot (fig 5B)
unisex_a_volcano <- EnhancedVolcano(res_unisex_a,
                                    lab = rep(NA, nrow(res_unisex_a)),
                                    x = 'log2FoldChange',
                                    y = 'padj',
                                    pCutoff = 0.05,
                                    FCcutoff = 1,
                                    cutoffLineType = 'dashed',
                                    title = NULL,
                                    subtitle = NULL,
                                    caption = NULL,
                                    legendLabels = NULL,
                                    legendPosition = 'none',
                                    axisLabSize = 6,
                                    pointSize = 0.5,
                                    gridlines.major = 0.03,
                                    gridlines.minor = 0.01,
                                    col=c('#dad9d9', '#dad9d9', 'grey', '#b9eefb'),
                                    colAlpha = 1,
                                    xlim = c(-5, 5),
                                    ylim = c(-0.5, 250))

ggsave(file="unisex_a_volcano.svg", plot=unisex_a_volcano, width=2, height=2, units = "in")



# analysis of unisex degs
unisex_alpha_degs_full <- subset(res_unisex_alpha, res_unisex_alpha$padj <= 0.05 & abs(res_unisex_alpha$log2FoldChange) >= 1)
write.csv(unisex_alpha_degs_full, "unisex_alpha_upreg.csv")
unisex_alpha_degs <- row.names(unisex_alpha_degs_full)

unisex_a_degs_full <- subset(res_unisex_a, res_unisex_a$padj <= 0.05 & abs(res_unisex_a$log2FoldChange) >= 1)
unisex_a_upreg <- subset(res_unisex_a, res_unisex_a$padj <= 0.05 & res_unisex_a$log2FoldChange >= 1)
write.csv(unisex_a_upreg, "unisex_a_upreg.csv")
unisex_a_downreg <- subset(res_unisex_a, res_unisex_a$padj <= 0.05 & res_unisex_a$log2FoldChange <= -1)
write.csv(unisex_a_downreg, "unisex_a_downreg.csv")
unisex_a_degs <- row.names(unisex_a_degs_full)

# make venn diagram (fig 5C)
venn.diagram(
  x = list(unisex_alpha_degs, unisex_a_degs),
  disable.logging = TRUE,
  category.names = c("" , ""),
  resolution = 500, 
  imagetype = "tiff",
  filename = 'unisex_venn.tif',
  output=FALSE,
  col = c("black", "black"),
  fill = c("#edaaa8", "#b9eefb"),
  cat.fontfamily = "Helvetica",  
  cat.cex = 1.8,                   
  fontfamily = "Helvetica",      
  cex = 0,
  scaled = TRUE,  # Proportional scaling of circles
  margin = 0.15,   # Margin around the diagram
  circle.padding = 0.05
)

# alphaxa crosses
bisex <- c("XL280alpha-a_rep1", "XL280alpha-a_rep2", "XL280alpha-a_rep3", 
           "JHG4-XL280a_rep1", "JHG4-XL280a_rep2", "JHG4-XL280a_rep3",
           "XL280alpha-JHG6_rep1", "XL280alpha-JHG6_rep2", "XL280alpha-JHG6_rep3", 
           "JHG4-JHG6_rep1", "JHG4-JHG6_rep2", "JHG4-JHG6_rep3")
## data wrangling
### merge counts dfs
bisex_counts = Reduce(function(x, y) 
  merge(x, y, by = "gene_ID", all = TRUE), 
  df_list[names(df_list) %in% bisex])

### make gene ids into row names
rownames(bisex_counts) <- bisex_counts$gene_ID
bisex_counts <- bisex_counts[,-1]

## TPM filter
#gene_counts_by_ID$gene <- row.names(gene_counts_by_ID)
bisex_counts$gene <- row.names(bisex_counts)
# get rid of "rna-" from row names
bisex_counts$gene <- gsub("rna-", "", bisex_counts$gene)
rownames(bisex_counts) <- bisex_counts$gene
tpm_counts_bisex <- merge(bisex_counts, gene_lengths[ , c("gene", "length")], by = "gene")
## calculate tpms
### convert length from bp to kbp
tpm_counts_bisex$length <- tpm_counts_bisex$length/1000
### df/character vector to hold tpm values
tpm_calc_bisex <- tpm_counts_bisex[, c("gene"), drop = FALSE]
### get cols except gene
tpm_cols_bisex <- setdiff(colnames(tpm_counts_bisex), c("gene", "length"))
### get tpm
for (gene in tpm_cols_bisex) {
  rpk <- tpm_counts_bisex[[gene]] / tpm_counts_bisex$length
  tpm <- rpk / sum(rpk) * 1e6
  tpm_calc_bisex[[gene]] <- tpm
}
## average tpm for technical reps
tpm_calc_bisex <- tpm_calc_bisex %>%
  pivot_longer(
    cols = -gene,
    names_to = "replicate",
    values_to = "tpm"
  )

tpm_calc_bisex_long <- tpm_calc_bisex %>%
  mutate(sample = str_extract(replicate, "^[^_]+_[^_]+"))

tpm_avg_bisex <- tpm_calc_bisex_long %>%
  group_by(gene, sample) %>%
  summarise(avg = mean(tpm), .groups = "drop") %>%
  pivot_wider(names_from = "sample", values_from = "avg")

tpm_filtered_genes_bisex <- tpm_avg_bisex %>%
  filter(if_any(-gene, ~ . >= 1)) %>%
  pull(gene)

bisex_counts_tpmfiltered <- bisex_counts[rownames(bisex_counts) %in% tpm_filtered_genes_bisex, ]
# remove last column of gene ids
bisex_counts_tpmfiltered <- bisex_counts_tpmfiltered[ ,1:12]

## DEG analysis
### make df with appropriate sample info
sample_bisex <- c(colnames(bisex_counts_tpmfiltered))
condition_bisex <- c('mutalpha_muta', 'mutalpha_muta', 'mutalpha_muta', 
                     'mutalpha_wta', 'mutalpha_wta', 'mutalpha_wta',
                     'wtalpha_wta', 'wtalpha_wta', 'wtalpha_wta',
                     'wtalpha_muta', 'wtalpha_muta', 'wtalpha_muta'
) 
colData_bisex <- data.frame(sample_bisex, condition_bisex)

### run deseq
dds_bisex <- DESeqDataSetFromMatrix(countData = bisex_counts_tpmfiltered,
                                    colData = colData_bisex,
                                    design = ~ condition_bisex)

#### check quality
trans_bisex <- vst(dds_bisex, blind = FALSE)
plotPCA(trans_bisex, intgroup = c("condition_bisex")) # PC1: 62%; PC2: 25%
pca_bisex_data <- plotPCA(trans_bisex, intgroup = c("condition_bisex"), returnData=TRUE) 
pca_bisex <- ggplot(pca_bisex_data, 
                    aes(PC1, PC2, color=condition_bisex, shape=condition_bisex)) + 
  geom_richtext(data = pca_bisex_data, 
                aes(label = condition_bisex),
                color = "black") +
  geom_point(size = 1.5, shape = c(23, 23, 23, 22, 22, 22, 21, 21, 21, 24, 24, 24),  fill = "white", stroke = 1) +
  scale_color_manual(values=c("#d2ceb3","#edaaa8", "#b9eefb","#d8d8d8")) + 
  theme(panel.grid.major = element_line(color = "#eeeeee", linewidth = 0.3), 
        panel.grid.minor = element_line(color = "#eeeeee", linewidth = 0.1),  
        panel.background = element_blank(), 
        axis.line = element_line(color = "black"),
        panel.border = element_rect(color = "black", 
                                    fill=NA, 
                                    linewidth = 0.75),
        legend.position = "none") +
  xlab("PC1: 62%") +
  ylab("PC2: 25%")
ggsave(file="pca_bisex.svg", plot=pca_bisex, width=2, height=2, units = "in")

mat_bisex <- assay(trans_bisex)
matrix_bisex <- SummarizedExperiment::assay(trans_bisex)
topVarGenes_bisex <- head(order(rowVars(matrix_bisex), decreasing = TRUE),500)
topGenesMatrix_bisex <- matrix_bisex[topVarGenes_bisex, ]

heatmap_bisex <- pheatmap(topGenesMatrix_bisex,
                          scale = "row",
                          cluster_rows = TRUE,
                          cluster_cols = TRUE,
                          show_rownames = FALSE,
) # clustering checks out mostly

ggsave(file="heatmap_bisex.png", plot=heatmap_bisex, width=5, height=5)

### compare muts to wt
dds_bisex$condition_bisex <- relevel(dds_bisex$condition_bisex, ref = "wtalpha_wta")
dds_bisex <- DESeq(dds_bisex)


res_bisex_JHG4_wta <- results(dds_bisex,
                             contrast = c("condition_bisex", "mutalpha_wta",
                                          "wtalpha_wta"),
                             alpha = 0.05)

res_bisex_JHG4_wta <- as.data.frame(res_bisex_JHG4_wta)
res_bisex_JHG4_wta <- res_bisex_JHG4_wta[!rownames(res_bisex_JHG4_wta) %in% artifacts, ]
# save data
write.csv(res_bisex_JHG4_wta, "bisex_JHG4_wta.csv")


### make volcano plot (Fig 5D)
bisex_JHG4_wta_volcano <- EnhancedVolcano(res_bisex_JHG4_wta,
                                          lab = rep(NA, nrow(res_bisex_JHG4_wta)),
                                          x = 'log2FoldChange',
                                          y = 'padj',
                                          pCutoff = 0.05,
                                          FCcutoff = 1,
                                          cutoffLineType = 'dashed',
                                          title = NULL,
                                          subtitle = NULL,
                                          caption = NULL,
                                          legendLabels = NULL,
                                          legendPosition = 'none',
                                          axisLabSize = 6,
                                          pointSize = 0.5,
                                          gridlines.major = 0.03,
                                          gridlines.minor = 0.01,
                                          col=c('#dad9d9', '#dad9d9', 'grey', '#edaaa8'),
                                          colAlpha = 1,
                                          xlim = c(-5, 5),
                                          ylim = c(-0.5, 250))

ggsave(file="bisex_JHG4_wta_volcano.svg", plot=bisex_JHG4_wta_volcano, width=2, height=2, units = "in")


# common genes with JHG4 solo
res_bisex_JHG4_wta_degs_full <- subset(res_bisex_JHG4_wta, res_bisex_JHG4_wta$padj <= 0.05 & abs(res_bisex_JHG4_wta$log2FoldChange) >= 1)
bisex_JHG4_wta_upreg <- subset(res_bisex_JHG4_wta, res_bisex_JHG4_wta$padj <= 0.05 & res_bisex_JHG4_wta$log2FoldChange >= 1)
write.csv(bisex_JHG4_wta_upreg, "bisex_JHG4_wta_upreg.csv")
bisex_JHG4_wta_downreg <- subset(res_bisex_JHG4_wta, res_bisex_JHG4_wta$padj <= 0.05 & res_bisex_JHG4_wta$log2FoldChange <= -1)
write.csv(bisex_JHG4_wta_downreg, "bisex_JHG4_wta_downreg.csv")
bisex_JHG4_wta_degs <- row.names(res_bisex_JHG4_wta_degs_full)

## wtalpha v. muta
res_bisex_wtalpha_JHG6 <- results(dds_bisex,
                                  contrast = c("condition_bisex", "wtalpha_muta",
                                               "wtalpha_wta"),
                                  alpha = 0.05)

res_bisex_wtalpha_JHG6 <- as.data.frame(res_bisex_wtalpha_JHG6)
res_bisex_wtalpha_JHG6 <- res_bisex_wtalpha_JHG6[!rownames(res_bisex_wtalpha_JHG6) %in% artifacts, ]
# save data
write.csv(res_bisex_wtalpha_JHG6, "bisex_wtalpha_JGH6.csv")

# fig 5E
bisex_wtalpha_JHG6_volcano <- EnhancedVolcano(res_bisex_wtalpha_JHG6,
                                              lab = rep(NA, nrow(res_bisex_wtalpha_JHG6)),
                                              x = 'log2FoldChange',
                                              y = 'padj',
                                              pCutoff = 0.05,
                                              FCcutoff = 1,
                                              cutoffLineType = 'dashed',
                                              title = NULL,
                                              subtitle = NULL,
                                              caption = NULL,
                                              legendLabels = NULL,
                                              legendPosition = 'none',
                                              axisLabSize = 6,
                                              pointSize = 0.5,
                                              gridlines.major = 0.03,
                                              gridlines.minor = 0.01,
                                              col=c('#dad9d9', '#dad9d9', 'grey', '#b9eefb'),
                                              colAlpha = 1,
                                              xlim = c(-5, 5),
                                              ylim = c(-0.5, 250))

ggsave(file="bisex_wtalpha_JHG6_volcano.svg", plot=bisex_wtalpha_JHG6_volcano, width=2, height=2, units = "in")

bisex_wtalpha_JHG6_degs_full <- subset(res_bisex_wtalpha_JHG6, res_bisex_wtalpha_JHG6$padj <= 0.05 & abs(res_bisex_wtalpha_JHG6$log2FoldChange) >= 1)
bisex_wtalpha_JHG6_upreg <- subset(res_bisex_wtalpha_JHG6, res_bisex_wtalpha_JHG6$padj <= 0.05 & res_bisex_wtalpha_JHG6$log2FoldChange >= 1)
write.csv(bisex_wtalpha_JHG6_upreg, "bisex_wtalpha_JHG6_upreg.csv")

bisex_wtalpha_JHG6_downreg <- subset(res_bisex_wtalpha_JHG6, res_bisex_wtalpha_JHG6$padj <= 0.05 & res_bisex_wtalpha_JHG6$log2FoldChange <= -1)
write.csv(bisex_wtalpha_JHG6_downreg, "bisex_wtalpha_JHG6_downreg.csv")
bisex_wtalpha_JHG6_degs <- row.names(bisex_wtalpha_JHG6_degs_full)

## mutalpha v. muta
res_bisex_JHG4_JHG6 <- results(dds_bisex,
                               contrast = c("condition_bisex", "mutalpha_muta",
                                            "wtalpha_wta"),
                               alpha = 0.05)

res_bisex_JHG4_JHG6 <- as.data.frame(res_bisex_JHG4_JHG6)
res_bisex_JHG4_JHG6 <- res_bisex_JHG4_JHG6[!rownames(res_bisex_JHG4_JHG6) %in% artifacts, ]
# save data
write.csv(res_bisex_JHG4_JHG6, "bisex_JHG4_JGH6.csv")
bisex_JHG4_JHG6_degs_full <- subset(res_bisex_JHG4_JHG6, res_bisex_JHG4_JHG6$padj <= 0.05 & abs(res_bisex_JHG4_JHG6$log2FoldChange) >= 1)
bisex_JHG4_JHG6_upreg <- subset(res_bisex_JHG4_JHG6, res_bisex_JHG4_JHG6$padj <= 0.05 & res_bisex_JHG4_JHG6$log2FoldChange >= 1)
write.csv(bisex_JHG4_JHG6_upreg, "bisex_JHG4_JHG6_upreg.csv")
bisex_JHG4_JHG6_downreg <- subset(res_bisex_JHG4_JHG6, res_bisex_JHG4_JHG6$padj <= 0.05 & res_bisex_JHG4_JHG6$log2FoldChange <= -1)
write.csv(bisex_JHG4_JHG6_downreg, "bisex_JHG4_JHG6_downreg.csv")

### make volcano plot (Fig 5F)
bisex_JHG4_JHG6_volcano <- EnhancedVolcano(res_bisex_JHG4_JHG6,
                                           lab = rep(NA, nrow(res_bisex_JHG4_JHG6)),
                                           x = 'log2FoldChange',
                                           y = 'padj',
                                           pCutoff = 0.05,
                                           FCcutoff = 1,
                                           cutoffLineType = 'dashed',
                                           title = NULL,
                                           subtitle = NULL,
                                           caption = NULL,
                                           legendLabels = NULL,
                                           legendPosition = 'none',
                                           axisLabSize = 6,
                                           pointSize = 0.5,
                                           gridlines.major = 0.03,
                                           gridlines.minor = 0.01,
                                           col=c('#dad9d9', '#dad9d9', 'grey', '#d2ceb3'),
                                           colAlpha = 1,
                                           xlim = c(-5, 5),
                                           ylim = c(-0.5, 250))

ggsave(file="bisex_JHG4_JHG6_volcano.svg", plot=bisex_JHG4_JHG6_volcano, width=2, height=2, units = "in")

# compare sxi1∆ DEGs between unisex and bisex (Fig. 5G)
venn.diagram(
  x = list(unisex_alpha_degs, bisex_JHG4_wta_degs),
  disable.logging = TRUE,
  category.names = c("" , ""),
  resolution = 500, 
  imagetype = "tiff",
  filename = 'alpha_unisex_bisex_venn.tif',
  output=FALSE,
  col = c("black", "black"),
  fill = c("#edaaa8", "#bf4040"),
  cat.fontfamily = "Helvetica",  
  cat.cex = 1.8,                   
  fontfamily = "Helvetica",      
  cex = 0,
  scaled = TRUE,  # Proportional scaling of circles
  margin = 0.15,   # Margin around the diagram
  circle.padding = 0.05
)
# get list of shared DEGs
Reduce(intersect, list(unisex_alpha_degs, bisex_JHG4_wta_degs))

# compare sxi2∆ DEGs between unisex and bisex (Fig. 5H)
venn.diagram(
  x = list(unisex_a_degs, bisex_wtalpha_JHG6_degs),
  disable.logging = TRUE,
  category.names = c("" , ""),
  resolution = 500, 
  imagetype = "tiff",
  filename = 'a_unisex_bisex_venn.tif',
  output=FALSE,
  col = c("black", "black"),
  fill = c("#b9eefb", "#4979ad"),
  cat.fontfamily = "Helvetica",  
  cat.cex = 1.8,                   
  fontfamily = "Helvetica",      
  cex = 0,
  scaled = TRUE,  # Proportional scaling of circles
  margin = 0.15,   # Margin around the diagram
  circle.padding = 0.05
)
# get list of shared DEGs
Reduce(intersect, list(unisex_a_degs, bisex_wtalpha_JHG6_degs))

## MAT loci heatmaps (Fig. 5J and K)
# make merged df with 3 bisexual crosses 
bisex_dfs <- list(
  JH4_wta = res_bisex_JH4_wta,
  wtalpha_JH6 = res_bisex_wtalpha_JHG6,
  JH4_JH6 = res_bisex_JHG4_JHG6
)

bisex_dfs_short <- list()

for (df in names(bisex_dfs)) {
  d <- bisex_dfs[[df]]
  d$gene <- row.names(d)
  df_simple <- d[, c("gene", "log2FoldChange", "padj")]
  colnames(df_simple)[2:3] <- c(paste0("log2FC_", df), paste0("padj_", df))
  bisex_dfs_short[[df]] <- df_simple
}


bisex_merged_df <- Reduce(function(x,y) merge(x, y, by = "gene"), bisex_dfs_short)

## pull out mat alpha genes
### coords: Chr4 1541390-1646179
matalpha_genes <- c("CND05680", "CND05690", "CND05700", "CND05710", "CND05720",
                    "CND05740", "CND05750", "CND05760", "CND05770", "CND05780",
                    "CND05790", "CND05800", "CND05810", "CND05820", "CND05830",
                    "CND05840", "CND05850", "CND05860", "CND05870", "CND05880",
                    "CND05890", "CND05910", "CND05920", "CND05940")

bisex_genes_matalpha <- subset(bisex_merged_df, bisex_merged_df$gene %in% matalpha_genes)
### set log2FC genes with padj > 0.05 = 0
bisex_samples <- c("JH4_wta", "wtalpha_JH6", "JH4_JH6")
for (s in bisex_samples) {
  log2FC_col <- paste0("log2FC_", s)
  padj_col <- paste0("padj_", s)
  bisex_genes_matalpha[[log2FC_col]][bisex_genes_matalpha[[padj_col]] > 0.05] <- 0
}

bisex_genes_matalpha_forheatmap <- subset(bisex_genes_matalpha, select = colnames(bisex_genes_matalpha) %in% c("gene", "log2FC_JH4_wta", "log2FC_wtalpha_JH6", "log2FC_JH4_JH6"))
rownames(bisex_genes_matalpha_forheatmap) <- bisex_genes_matalpha_forheatmap$gene
bisex_genes_matalpha_forheatmap <- bisex_genes_matalpha_forheatmap[,2:4]

## manually define colors and scale for heatmaps
breaks <- seq(-6, 6, length.out = 101)
colors_heatmap <- colorRampPalette(c("#000c91", "#6b74d0", "#e0e2ff", "white", "#e79536", "#cb432a", "#ba202d"))(length(breaks) - 1)

bisex_matalpha_heatmap <- pheatmap(bisex_genes_matalpha_forheatmap,
                                   color = colors_heatmap,
                                   breaks = breaks)

ggsave(file="bisex_matalpha_heatmap.svg", plot=bisex_matalpha_heatmap, width=5, height=8)

## mata
mata_genes <- c("AF542530.2:14556..20135", "AF542530.2:22185..26193", "AF542530.2:26786..29505", "AF542530.2:31113..31587",
                "AF542530.2:34695..36940", "AF542530.2:45822..47844", "AF542530.2:48200..51121", "AF542530.2:51957..55098", "AF542530.2:56006..62001",
                "AF542530.2:65419..68176", "AF542530.2:71749..72350", "AF542530.2:74911..77572", "AF542530.2:77703..81239", "AF542530.2:82207..86105",
                "AF542530.2:87889..88017", "AF542530.2:89248..91842", "AF542530.2:93853..95179", "AF542530.2:98365..98493", "AF542530.2:106876..107004",
                "AF542530.2:112926..115296", "AF542530.2:117554..123345", "AF542530.2:124039..126330", "AF542530.2:127675..129090")

bisex_genes_mata <- subset(bisex_merged_df, bisex_merged_df$gene %in% mata_genes)

for (s in bisex_samples) {
  log2FC_col <- paste0("log2FC_", s)
  padj_col <- paste0("padj_", s)
  bisex_genes_mata[[log2FC_col]][bisex_genes_mata[[padj_col]] > 0.05] <- 0
}

bisex_genes_mata_forheatmap <- subset(bisex_genes_mata, select = colnames(bisex_genes_mata) %in% c("gene", "log2FC_JH4_wta", "log2FC_wtalpha_JH6", "log2FC_JH4_JH6"))
rownames(bisex_genes_mata_forheatmap) <- bisex_genes_mata_forheatmap$gene
bisex_genes_mata_forheatmap <- bisex_genes_mata_forheatmap[,2:4]

bisex_mata_heatmap <- pheatmap(bisex_genes_mata_forheatmap,
                               color = colors_heatmap,
                               breaks = breaks)
ggsave(file="bisex_mata_heatmap.svg", plot=bisex_mata_heatmap, width=6.5, height=8)

# all bisex DEGs
bisex_DEGs_all <- Reduce(union, list(bisex_JHG4_wta_degs, bisex_wtalpha_JHG6_degs, bisex_JHG4_JHG6_degs))
bisex_DEGs_df <- subset(bisex_merged_df, bisex_merged_df$gene %in% bisex_DEGs_all)
for (s in bisex_samples) {
  log2FC_col <- paste0("log2FC_", s)
  padj_col <- paste0("padj_", s)
  bisex_DEGs_df[[log2FC_col]][bisex_DEGs_df[[padj_col]] > 0.05] <- 0
}

bisex_DEGs_forheatmap <- subset(bisex_DEGs_df, select = colnames(bisex_DEGs_df) %in% c("gene", "log2FC_JH4_wta", "log2FC_wtalpha_JH6", "log2FC_JH4_JH6"))
rownames(bisex_DEGs_forheatmap) <- bisex_DEGs_forheatmap$gene
bisex_DEGs_forheatmap <- bisex_DEGs_forheatmap[,2:4]

breaks2 <- seq(-7.5, 7.5, length.out = 101)
colors_heatmap2 <- colorRampPalette(c("#000c91", "#6b74d0", "#e0e2ff", "white", "#e79536", "#cb432a", "#ba202d"))(length(breaks2) - 1)
bisex_DEGs_heatmap <- pheatmap(bisex_DEGs_forheatmap,
                               color = colors_heatmap2,
                               breaks = breaks2)

# get DEGs for different regulatory scenarios (Fig. 5I)
## divide DEGs into upreg/downreg genes
bisex_JHG4_wta_up <- rownames(res_bisex_JHG4_wta_degs_full[res_bisex_JHG4_wta_degs_full$log2FoldChange > 1 ,]) # 147
bisex_JHG4_wta_down <- rownames(res_bisex_JHG4_wta_degs_full[res_bisex_JHG4_wta_degs_full$log2FoldChange < -1 ,]) # 189

bisex_wtalpha_JHG6_up <- rownames(bisex_wtalpha_JHG6_degs_full[bisex_wtalpha_JHG6_degs_full$log2FoldChange > 1 ,]) # 83
bisex_wtalpha_JHG6_down <- rownames(bisex_wtalpha_JHG6_degs_full[bisex_wtalpha_JHG6_degs_full$log2FoldChange < -1 ,]) # 107

bisex_JHG4_JHG6_up <- rownames(bisex_JHG4_JHG6_degs_full[bisex_JHG4_JHG6_degs_full$log2FoldChange > 1 ,]) # 99
bisex_JHG4_JHG6_down <- rownames(bisex_JHG4_JHG6_degs_full[bisex_JHG4_JHG6_degs_full$log2FoldChange < -1 ,]) # 103

## Case 1: Sxi1/Sxi2 heterodimer regulates
### 1A: Heterodimer repressed; increased expression in all comparisons
Case1A_DEGs <- Reduce(intersect, list(bisex_JHG4_wta_up, bisex_wtalpha_JHG6_up, bisex_JHG4_JHG6_up)) # 57

### 1B: Heterodimer activated; decreased expression in all comparisons
Case1B_DEGs <- Reduce(intersect, list(bisex_JHG4_wta_down, bisex_wtalpha_JHG6_down, bisex_JHG4_JHG6_down)) # 65

## Case 2: Independently regulated by Sxi1
### 2A: Sxi1 represses; up in sxi1∆/wt BUT NOT wt/sxi2∆
Case2A_DEGs <- setdiff(bisex_JHG4_wta_up, bisex_wtalpha_JHG6_up) # 78

### 2B: Sxi1 activates; down in sxi1∆/wt BUT NOT wt/sxi2∆
Case2B_DEGs <- setdiff(bisex_JHG4_wta_down, bisex_wtalpha_JHG6_down) # 109

## Case 3: Independently regulated by Sxi2
### 3A: Sxi2 represses; up in wt/sxi2∆ BUT NOT sxi1∆/wt
Case3A_DEGs <- setdiff(bisex_wtalpha_JHG6_up, bisex_JHG4_wta_up) # 14

### 3B: Sxi2 activates; down in wt/sxi2∆ BUT NOT sxi1∆/wt
Case3B_DEGs <- setdiff(bisex_wtalpha_JHG6_down, bisex_JHG4_wta_down) # 27
