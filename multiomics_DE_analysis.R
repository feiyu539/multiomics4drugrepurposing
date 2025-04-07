# Multi-Omics Differential Expression Analysis
# Includes: RNA-seq (mRNA), Metabolome, and Methylation (RRBS)
# Author: ldy | Date: 2025-02-07

# Set working directory
setwd("your/project/path")  # Update to project root

# ==== 1. RNA-seq Differential Expression using edgeR ====
if (!require(edgeR)) install.packages("edgeR")
if (!require(limma)) install.packages("limma")
if (!require(ggplot2)) install.packages("ggplot2")

library(edgeR)
library(limma)
library(ggplot2)

count <- read.csv("data/count.csv")
group_info <- read.csv("data/infor.csv")
rownames(count) <- count$X
count <- count[, -1]

high.control <- group_info[1:62, ]
low.treat <- group_info[63:121, ]
order <- c(high.control$Subject_ID, low.treat$Subject_ID)
count2 <- count[, order]
group <- c(rep("high", 62), rep("low", 59))

dgelist <- DGEList(counts = count2, group = group)
keep <- rowSums(cpm(dgelist) > 1) >= 2
dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]
dgelist_norm <- calcNormFactors(dgelist, method = "TMM")

design <- model.matrix(~group)
dge <- estimateDisp(dgelist_norm, design, robust = TRUE)
fit <- glmFit(dge, design, robust = TRUE)
lrt <- topTags(glmLRT(fit), n = nrow(dgelist$counts))
write.table(lrt, file = "output/control_treat.glmLRT.txt", sep = "\t", col.names = NA, quote = FALSE)

fit <- glmQLFit(dge, design, robust = TRUE)
lrt <- topTags(glmQLFTest(fit), n = nrow(dgelist$counts))
write.table(lrt, file = "output/control_treat.glmQLFit.txt", sep = "\t", col.names = NA, quote = FALSE)

gene_diff <- read.delim("output/control_treat.glmLRT.txt", row.names = 1, sep = "\t", check.names = FALSE)
gene_diff <- gene_diff[order(gene_diff$FDR, gene_diff$logFC, decreasing = c(FALSE, TRUE)), ]

PValue_cutoff <- 3.422313e-06
logFC_cutoff <- 1.5
gene_diff$sig <- "none"
gene_diff[gene_diff$logFC >= logFC_cutoff & gene_diff$PValue < PValue_cutoff, "sig"] <- "up"
gene_diff[gene_diff$logFC <= -logFC_cutoff & gene_diff$PValue < PValue_cutoff, "sig"] <- "down"

write.csv(subset(gene_diff, sig %in% c("up", "down")), file = "output/control_treat.glmQLFit.select.csv")
write.csv(subset(gene_diff, sig == "up"), file = "output/control_treat.glmQLFit.up.csv")
write.csv(subset(gene_diff, sig == "down"), file = "output/control_treat.glmQLFit.down.csv")

ggplot(gene_diff, aes(x = logFC, y = -log10(PValue), color = sig)) +
  geom_point(size = 1) +
  scale_color_manual(values = c("red", "gray", "green")) +
  labs(x = "log2 Fold Change", y = "-log10 P-value", title = "control vs treat") +
  theme_minimal() +
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), lty = 3) +
  geom_hline(yintercept = -log10(PValue_cutoff), lty = 3) +
  xlim(-10, 10) + ylim(0, 10)
ggsave("output/volcano_plot.png")


# ==== 2. Metabolome Analysis ====
if (!require(plyr)) install.packages("plyr")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require(mixOmics)) BiocManager::install("mixOmics")

library(plyr)
library(mixOmics)

pos <- read.csv("data/pos.csv"); rownames(pos) <- pos$row.ID; pos <- pos[, -1]
neg <- read.csv("data/neg.csv"); rownames(neg) <- neg$row.ID; neg <- neg[, -1]
group <- read.csv("data/group.csv", header = FALSE)[, 1]

design <- model.matrix(~ 0 + factor(group))
colnames(design) <- c("control", "disease")
rownames(design) <- colnames(pos)
contrast.matrix <- makeContrasts(control - disease, levels = design)

# POS
fit1 <- lmFit(pos, design)
fit1 <- eBayes(contrasts.fit(fit1, contrast.matrix))
nrDEG1 <- na.omit(topTable(fit1, coef = 1, n = Inf, adjust.method = "BH"))
write.csv(nrDEG1[nrDEG1$P.Value < 0.05, ], file = "output/DM_pos.csv")

# NEG
fit2 <- lmFit(neg, design)
fit2 <- eBayes(contrasts.fit(fit2, contrast.matrix))
nrDEG2 <- na.omit(topTable(fit2, coef = 1, n = Inf, adjust.method = "BH"))
write.csv(nrDEG2[nrDEG2$P.Value < 0.05, ], file = "output/DM_neg.csv")

# T-tests
group_df <- read.csv("data/group.csv", header = FALSE)
rownames(group_df) <- colnames(pos)
control <- rownames(group_df)[group_df$V1 == 0]
disease <- rownames(group_df)[group_df$V1 == 1]

# t-test for pos
poscontrol <- pos[, control]
posdisease <- pos[, disease]
a <- cbind(poscontrol, posdisease)
meta1 <- ldply(lapply(1:nrow(a), function(i) {
  r <- t.test(a[i, 1:62], a[i, 63:121])
  data.frame(Pvalue = r$p.value, FDR = p.adjust(r$p.value, "BH"))
}), data.frame)
rownames(meta1) <- rownames(pos)
write.csv(meta1[meta1$Pvalue < 0.05, ], file = "output/t_test_pos.csv")

# t-test for neg
negcontrol <- neg[, control]
negdisease <- neg[, disease]
b <- cbind(negcontrol, negdisease)
meta2 <- ldply(lapply(1:nrow(b), function(i) {
  r <- t.test(b[i, 1:62], b[i, 63:121])
  data.frame(Pvalue = r$p.value, FDR = p.adjust(r$p.value, "BH"))
}), data.frame)
rownames(meta2) <- rownames(neg)
write.csv(meta2[meta2$Pvalue < 0.05, ], file = "output/t_test_neg.csv")

# PLS-DA
pos_t <- t(pos)
group_f <- factor(group_df$V1)
pls_model_pos <- plsda(pos_t, group_f, ncomp = 2)
vip_pos <- vip(pls_model_pos)
write.csv(intersect(rownames(vip_pos)[vip_pos$comp1 > 1.5], rownames(vip_pos)[vip_pos$comp2 > 1.5]),
          file = "output/plsda_pos.csv")

neg_t <- t(neg)
pls_model_neg <- plsda(neg_t, group_f, ncomp = 2)
vip_neg <- vip(pls_model_neg)
write.csv(intersect(rownames(vip_neg)[vip_neg$comp1 > 1.5], rownames(vip_neg)[vip_neg$comp2 > 1.5]),
          file = "output/plsda_neg.csv")


# ==== 3. RRBS Methylation Analysis ====
BiocManager::install(c("methylKit", "methylSig", "DSS", "annotatr", "genomation"), ask = FALSE)

library(methylKit)
library(methylSig)
library(annotatr)
library(genomation)

folders <- list.files("data/RRBS_NEW")
fold.list <- as.list(file.path("data/RRBS_NEW", folders))
sampleid <- as.list(gsub("_CpG.txt", "", folders))
group <- read.csv("data/group.csv", header = FALSE)[, 1]

myobj <- methRead(fold.list,
                  sample.id = sampleid,
                  assembly = "hg19",
                  treatment = group,
                  context = "CpG")

meth <- unite(myobj, destrand = FALSE)
myDiff <- calculateDiffMeth(meth)
all <- getMethylDiff(myDiff, difference = 5, qvalue = 0.05, type = "all", chunk.size = 121)

myDiff25p.hyper <- getMethylDiff(myDiff, difference = 5, qvalue = 0.01, type = "hyper")
myDiff25p.hypo  <- getMethylDiff(myDiff, difference = 5, qvalue = 0.01, type = "hypo")

regions <- tileMethylCounts(myobj, win.size = 1000, step.size = 1000)
meth1 <- unite(regions, destrand = FALSE)
myDiff1 <- calculateDiffMeth(meth1)
all1 <- getMethylDiff(myDiff1, difference = 5, qvalue = 0.05, type = "all", chunk.size = 121)

gene.obj <- readTranscriptFeatures(system.file("extdata", "hg19Tables.txt", package = "methylKit"))

diffAnn  <- annotateWithGeneParts(as(meth, "GRanges"), gene.obj)
diffAnn1 <- annotateWithGeneParts(as(meth1, "GRanges"), gene.obj)

gene  <- diffAnn@dist.to.TSS$feature.name
gene1 <- diffAnn1@dist.to.TSS$feature.name
write.csv(unique(c(gene, gene1)), file = "output/gene_combind.csv")

diffgene  <- unique(annotateWithGeneParts(as(all, "GRanges"), gene.obj)@dist.to.TSS$feature.name)
diffgene1 <- unique(annotateWithGeneParts(as(all1, "GRanges"), gene.obj)@dist.to.TSS$feature.name)

meth_data <- data.frame(chr = meth$chr, base = meth$start, end = meth$end)
write.csv(meth_data, file = "output/meth_data.csv")

anno_total <- diffAnn@dist.to.TSS
chr <- data.frame(chr = meth$chr, base = meth$start)
anno <- chr[anno_total$target.row, ]
anno_methylation <- cbind(anno, anno_total)

freq.C <- read.csv("data/freq_C.csv")[, -1]
names(meth_data) <- c("chr", "base", "end")
methy.expr <- merge(anno_methylation, freq.C, by = c("chr", "base"))
expr <- methy.expr[, -c(1:4, 6)]

NM <- expr$feature.name
nm <- gsub(".$", "", NM)
nm1 <- gsub("[[:punct:]]", "", nm)
expr$id <- nm1
expr <- expr[, -1]

refgene <- read.csv("data/refGene1.csv", header = FALSE, stringsAsFactors = FALSE)[, c(2, 13)]
colnames(refgene) <- c("id", "symbol")
gene_symbol <- merge(refgene, expr, by = "id")
gene_symbol <- gene_symbol[, -1]
symbol <- aggregate(. ~ symbol, data = gene_symbol, mean)
rownames(symbol) <- symbol$symbol
symbol <- symbol[, -1]
symbol <- round(symbol, digits = 3)

write.csv(symbol, file = "output/expression.csv")
write.csv(diffgene, file = "output/Diff_sitegene.csv")
write.csv(unique(gene), file = "output/total_sitegene.csv")
