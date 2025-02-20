# GSE75192 datasets 
# Merve KUTAY
#28.11.24

################################################################################
# Load necessary library
library(readxl)
library(dplyr)
library(DESeq2)
library(limma)
library(org.Mm.eg.db)
library(clusterProfiler)
library(ggplot2)
library(data.table)
library(AUCell)
library(writexl)
library(edgeR)



data <- read_excel("GSE75192_counts_liver.xls")

data<- data[,-1]

groups <- factor(c(rep("2months",4), rep("9months",5), rep("15months",5), rep("24months",5), rep("30months",5)))

# Step 1: Remove duplicate genes, keeping the first occurrence
data_unique <- data %>%
  group_by(external_gene_id) %>%
  slice_max(rowSums(across(where(is.numeric))), n = 1, with_ties = FALSE) %>%
  ungroup()

# Step 2: Remove genes starting with "Gm" followed by numbers (predicted genes)
data_filtered <- data_unique %>%
  filter(!grepl("^Gm[0-9]+$", external_gene_id))

# Step 3: Remove genes with all zero values or non-zero values in at most 2 samples
data_cleaned <- data_filtered %>%
  filter(rowSums(select(., -c(external_gene_id)) != 0) > 2)


# Save the cleaned data
data_gse75192 <- data_cleaned

write_xlsx(data_cleaned, "data_gse75192.xlsx")

################################################################################
################################################################################

# HDS scoring
################################################################################

# Delete the first column, gene_id
data_gse75192 <- data_gse75192[, -1]

# Assign the second column, gene_symbol, as rownames and remove it from the column

rownames(data_gse75192) <- data_cleaned$external_gene_id

#######################################################################################
#######################################################################################
# Step 1: Load libraries
  # Assuming AUCell is required for the calculation, install if not already

# Step 2: Load the necessary data
hd_genes <- read.csv("HDAG.csv", stringsAsFactors = FALSE)  # Load the HDAG gene list

# Step 3: Source the SharedFunctions.R to load the specific function
source("SharedFunctions.R")

HDAG <- hd_genes

# Step 5: Calculate the damage score using the function
# Assuming the function is like `calc_damage_score()`, you need to adapt based on the function's exact signature

# Here is a pseudo-function usage, replace with the actual function details
HDS.scores <- DS_calc.func(exprMatrices = data_gse75192,
                           DSignature = HDAG,
                           ntop = 42, 
                           ceilThrsh = 0.2,
                           wghtd = TRUE, 
                           progStat = FALSE)



# Step 6: Save the results
write.table(HDS.scores, file = "damage_score_results.csv", sep = ",", row.names = TRUE, quote = FALSE)


scores <- as.data.frame(HDS.scores)
sample_group2 <- rownames(scores)
hds<- scores$HDS.scores
sample_group <- factor(c(rep("2months",4), rep("9months",5), rep("15months",5), rep("24months",5), rep("30months",5)))

dataplot <- data.frame(hds, sample_group)
library(ggplot2)

ggplot(dataplot, aes(x = sample_group2, y = hds, color = sample_group)) +
  geom_point(size = 3) + 
  labs(title = "HDS Scoring for All Samples", x = "Samples", y = "Scores") +
  theme_minimal() +  
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 


ggplot(dataplot, aes(x = sample_group, y = hds, fill = sample_group)) +
  geom_violin(trim = FALSE, na.rm = TRUE) +  
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) + 
  labs(title = "HDS Scoring for Young and Old Groups", x = "Groups", y = "Scores") + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) 

################################################################################
################################################################################


# DESeq2 dataset 
dds <- DESeqDataSetFromMatrix(countData = data_gse75192,
                              colData = data.frame(condition = factor(c(rep("2months",4), rep("9months",5), rep("15months",5), rep("24months",5), rep("30months",5)))),
                              design = ~ condition)

# normalize data
vsd <- vst(dds, blind = TRUE)

# PCA analyse
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))


# PCA graph
ggplot(pcaData, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw()
################################################################################
data_gse75192_oy<-cbind.data.frame(data_gse75192[,5:9], data_gse75192[,20:24]) # it means old vs young
rownames(data_gse75192_oy)<- rownames(data_gse75192)

dds <- DESeqDataSetFromMatrix(countData = data_gse75192_oy,
                              colData = data.frame(condition = factor(c( rep("9months",5), rep("30months",5)))),
                              design = ~ condition)

# normalize
vsd <- vst(dds, blind = TRUE)

# PCA 
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))


# PCA graph
ggplot(pcaData, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw()

# DESeq2  for old vs young
dds <- DESeq(dds)
res_deseq <- results(dds)
res_deseq <- as.data.frame(res_deseq[order(res_deseq$padj), ])


res_deseq$gene_id <- rownames(res_deseq)
res_deseq <- merge(res_deseq, data_cleaned, by.x = "gene_id", by.y = "external_gene_id", all.x = TRUE)

# DESeq2 results order
res_deseq <- res_deseq[order(res_deseq$padj), ]

########################################################################################################################
# DGEA - limma for old vs young
design <- model.matrix(~ 0 + factor(c( rep("Young",5), rep("Old",5))))
colnames(design) <- c("Young", "Old")
v <- voom(data_gse75192_oy, design, plot = TRUE)
fit <- lmFit(v, design)
contrast <- makeContrasts(Old-Young, levels = design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)
res_limma <- topTable(fit2, adjust = "fdr", number = Inf)


res_limma$gene_id <- rownames(res_limma)
res_limma <- merge(res_limma, data_cleaned, by.x = "gene_id", by.y = "external_gene_id", all.x = TRUE)

# limma results order
res_limma <- res_limma[order(res_limma$adj.P.Val), ]

################################################################################
# data preperation
# 
group <- factor(c( rep("Young",5), rep("Old",5)))
y <- DGEList(counts = data_gse75192_oy, group = group)

# Filter genes
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes = FALSE]

# Normalize data
y <- calcNormFactors(y)

# produce matrix
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

# prediction disperition
y <- estimateDisp(y, design)

#fit to model
fit <- glmFit(y, design)

# produce contrast matrix
contrast <- makeContrasts(Old - Young, levels = design)
fit2 <- glmLRT(fit, contrast = contrast)

# Results
res_edgeR <- topTags(fit2, n = Inf)$table
res_edgeR$gene_id <- rownames(res_edgeR)

# order
res_edgeR <- res_edgeR[order(res_edgeR$FDR), ]

###########################################################################################
###########################################################################################

# Significant genes for 3 algorithms for 9m vs 30m comparison
sig_genes_deseq <- res_deseq$gene_id[which(res_deseq$padj < 0.05)]
sig_genes_limma <- res_limma$gene_id[which(res_limma$adj.P.Val < 0.05)]
sig_genes_edgeR <- res_edgeR$gene_id[which(res_edgeR$FDR < 0.05)]
# NA adj pvalue elimination
sig_genes_deseq <- na.omit(sig_genes_deseq)
sig_genes_limma <- na.omit(sig_genes_limma)
sig_genes_edgeR <- na.omit(sig_genes_edgeR)

# save result
write.csv(res_edgeR, file = "edgeR_results.csv", row.names = FALSE)
write.csv(sig_genes_edgeR, file = "edgeR_significant_genes.csv", row.names = FALSE)
write.csv(res_deseq, "DESeq2_results.csv", row.names = FALSE)
write.csv(res_limma, "limma_results.csv", row.names = FALSE)

####################################################################################################
####################################################################################################

# Comparison of significantly changed genes for 3 algorithms

library(VennDiagram)


gen_list1 <- sig_genes_deseq
gen_list2 <- sig_genes_limma
gene_list3 <- sig_genes_edgeR
# find common genes
common_genes <- intersect(gen_list1, gen_list2)

venn.plot <- venn.diagram(
  x = list("deseq2" = gen_list1, "limma" = gen_list2, "edgeR" = gene_list3),
  filename = NULL,
  output = TRUE,
  fill = c("lightblue", "lightgreen", "pink"),
  alpha = 0.5,
  cex = 2,
  cat.cex = 2,
  cat.pos = 1,
  cat.dist = 0,1,
  main = "Comparison of Limma, Deseq2, and edgeR results"
)

# Venn 
grid.newpage()

grid.draw(venn.plot)

# Saving this gene lists for ranking section
# It is related with another study

max_length <- max(length(gen_list1), length(gen_list2), length(common_genes))
gen_list1 <- c(gen_list1, rep(NA, max_length - length(gen_list1)))
gen_list2 <- c(gen_list2, rep(NA, max_length - length(gen_list2)))
common_genes <- c(common_genes, rep(NA, max_length - length(common_genes)))

output <- data.frame(Gen_List1 = gen_list1, Gen_List2 = gen_list2, Common_Genes = common_genes)

writexl::write_xlsx(output, "GSE75192_siggenes.xlsx")

#######################################################################################
#######################################################################################

#GO Analysis for only deseq2 and limma
# GO - DESeq2
ego_deseq <- enrichGO(gene = sig_genes_deseq,
                      OrgDb = org.Mm.eg.db,
                      keyType = "SYMBOL",
                      ont = "BP",
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.05,
                      readable = TRUE)
go_results_deseq <- as.data.frame(ego_deseq)

# GO - limma
ego_limma <- enrichGO(gene = sig_genes_limma,
                      OrgDb = org.Mm.eg.db,
                      keyType = "SYMBOL",
                      ont = "BP",
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.05,
                      readable = TRUE)
go_results_limma <- as.data.frame(ego_limma)

#writing results
write.csv(go_results_deseq, "GO_results_DESeq2.csv")
write.csv(go_results_limma, "GO_results_limma.csv")
