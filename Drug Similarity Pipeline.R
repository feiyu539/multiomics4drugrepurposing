# Drug Repositioning and Similarity Integration Pipeline
# GitHub version: all paths cleaned, 'cell-line' naming removed for generalization

# --- Setup ---
install.packages("tidyverse")
install.packages("limma")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("GeneExpressionSignature")

library(tidyverse)
library(limma)
library(GeneExpressionSignature)
library(tidyr)
library(dplyr)

# --- Load LINCS Data ---
AlllogFC <- read.csv("./data/level4.csv", stringsAsFactors = FALSE)
info <- read.csv("./data/inst_info.csv", header = TRUE, stringsAsFactors = FALSE)

info <- info %>% mutate(inst_id = gsub(":", ".", inst_id))
AlllogFC$x <- 1:nrow(AlllogFC)

Case_IDs <- info %>% filter(pert_id != "DMSO") %>% pull(inst_id)
Control_IDs <- info %>% filter(pert_id == "DMSO") %>% pull(inst_id)

Allrank <- do.call(cbind, lapply(Case_IDs, function(i) {
  sub <- AlllogFC[, c(i, "x")] %>% arrange(desc(.[[1]]))
  sub$x
}))
Allrank <- as.data.frame(Allrank)
colnames(Allrank) <- 1:ncol(Allrank)

rownames(info) <- info$inst_id
pdata <- data.frame(pert_id = info[Case_IDs, ]$pert_id)
rownames(pdata) <- colnames(Allrank)
pdata$pert_id <- as.factor(pdata$pert_id)

metadata <- data.frame(labelDescription = c("drug name"), row.names = c("pert_id"))
phenodata <- new("AnnotatedDataFrame", data = pdata, varMetadata = metadata)

Allrank <- as.matrix(Allrank)
data <- ExpressionSet(assayData = Allrank, phenoData = phenodata)

MergingSet <- RankMerging(data, "Spearman", weighted = TRUE)
ds <- ScoreGSEA(MergingSet, 250, "avg") %>% as.data.frame()
write.csv(ds, file = "./output/drug_expression_similarity.csv")

# --- Drug Similarity Matrix Construction ---
inter <- read.csv("./data/chemical_chemical.csv", stringsAsFactors = FALSE)
name <- read.csv("./data/infoFinal.csv", stringsAsFactors = FALSE, header = FALSE)[-1, -1]

textming <- inter[, c(1,2,5)]
structure <- na.omit(inter[, 1:3])
NCI <- inter[, c(1,2,4)]
overlap <- intersect(inter$chemical1, inter$chemical2)

expression <- read.csv("./data/drug_expr_similarity.csv")
rownames(expression) <- expression$X
expression <- expression[, -1]
colnames(expression) <- rownames(expression)

expre2 <- expression[name$V2, name$V2]

cidname <- name$V10
newcid <- sprintf("CIDs%08d", as.integer(cidname))
rownames(expre2) <- newcid
colnames(expre2) <- newcid
cidname <- newcid
lap <- intersect(overlap, cidname)

str.over <- expand.grid(cidname, cidname) %>% mutate(across(everything(), as.character))
names(structure) <- c("Var1", "Var2", "Value")
str.new <- left_join(str.over, structure, by = c("Var1", "Var2")) %>% mutate(Value = ifelse(Var1 == Var2, 1000, Value))
str.new1 <- matrix(str.new$Value, ncol = length(cidname), byrow = TRUE) %>% as.data.frame()
colnames(str.new1) <- rownames(str.new1) <- cidname
str.new1[is.na(str.new1)] <- 0

names(NCI) <- c("Var1", "Var2", "Value")
expe.new <- left_join(str.over, NCI, by = c("Var1", "Var2")) %>% mutate(Value = ifelse(Var1 == Var2, 1000, Value))
expe.new1 <- matrix(expe.new$Value, ncol = length(cidname), byrow = TRUE) %>% as.data.frame()
colnames(expe.new1) <- rownames(expe.new1) <- cidname
expe.new1[is.na(expe.new1)] <- 0

names(textming) <- c("Var1", "Var2", "Value")
text <- left_join(str.over, textming, by = c("Var1", "Var2")) %>% mutate(Value = ifelse(Var1 == Var2, 1000, Value))
text1 <- matrix(text$Value, ncol = length(cidname), byrow = TRUE) %>% as.data.frame()
colnames(text1) <- rownames(text1) <- cidname
text1[is.na(text1)] <- 0

expression.1000 <- expre2 * 1000
expr.list <- expression.1000 %>%
  as_tibble(rownames = "Var1") %>%
  pivot_longer(-Var1, names_to = "Var2", values_to = "Value")

expr.1 <- expand.grid(cidname, cidname) %>% mutate(across(everything(), as.character))
expr1 <- left_join(expr.1, expr.list, by = c("Var1", "Var2")) %>% mutate(Value = ifelse(Var1 == Var2, 1000, Value))
expr2 <- matrix(expr1$Value, ncol = length(cidname), byrow = TRUE) %>% as.data.frame()
colnames(expr2) <- rownames(expr2) <- cidname
expr2[is.na(expr2)] <- 0

# --- SNF integration ---
library(SNFtool)
W2 <- SNF(list(as.matrix(str.new1), as.matrix(expe.new1), as.matrix(text1), as.matrix(expr2)), 20, 20)
W2 <- as.data.frame(W2)
colnames(W2) <- rownames(W2) <- name$V4
write.csv(W2, file = "./output/final_drug_similarity.csv")
