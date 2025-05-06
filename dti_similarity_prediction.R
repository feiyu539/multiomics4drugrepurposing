# Drug-Target Interaction (DTI) Similarity & DT-Hybrid Prediction Pipeline
# Author: ldy | Cleaned for GitHub release

# === Setup ===
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GOSemSim", "org.Hs.eg.db", "clusterProfiler"), ask = FALSE)
install.packages("SNFtool")
install.packages("tidyverse")

library(data.table)
library(GOSemSim)
library(org.Hs.eg.db)
library(clusterProfiler)
library(SNFtool)
library(tidyverse)

# === Load DTI Data ===
dti_file <- "data/DTI_symbol.csv"  # Update with your path
DTI <- fread(dti_file)[, -1]
target_symbols <- unique(DTI$V2)

# === Convert SYMBOL to ENTREZ ===
gene.df <- bitr(target_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# === GO Semantic Similarity ===
d_mf <- godata("org.Hs.eg.db", keytype = "ENTREZID", ont = "MF", computeIC = FALSE)
d_bp <- godata("org.Hs.eg.db", keytype = "ENTREZID", ont = "BP", computeIC = FALSE)
d_cc <- godata("org.Hs.eg.db", keytype = "ENTREZID", ont = "CC", computeIC = FALSE)

sim_mf <- mgeneSim(gene.df$ENTREZID, semData = d_mf, measure = "Wang")
sim_bp <- mgeneSim(gene.df$ENTREZID, semData = d_bp, measure = "Wang")
sim_cc <- mgeneSim(gene.df$ENTREZID, semData = d_cc, measure = "Wang")

# === Convert Similarity Matrix to Full Symmetric Matrix ===
make_full_matrix <- function(sim, genes) {
  sim <- as.data.frame(sim)
  sim_list <- sim %>% as_tibble() %>%
    add_column(Var1 = rownames(sim)) %>%
    pivot_longer(-Var1, names_to = "Var2", values_to = "Value")
  pairs <- expand.grid(genes, genes) %>%
    transform(Var1 = as.character(Var1), Var2 = as.character(Var2))
  sim_full <- left_join(pairs, sim_list, by = c("Var1", "Var2")) %>%
    mutate(Value = ifelse(Var1 == Var2, 1, Value))
  mat <- matrix(sim_full$Value, ncol = length(genes), byrow = TRUE)
  colnames(mat) <- rownames(mat) <- genes
  mat[is.na(mat)] <- 0
  return(as.matrix(mat))
}

mat_bp <- make_full_matrix(sim_bp, gene.df$ENTREZID)
mat_cc <- make_full_matrix(sim_cc, gene.df$ENTREZID)
mat_mf <- make_full_matrix(sim_mf, gene.df$ENTREZID)

# === Fuse Similarity Using SNF ===
W1 <- SNF(list(mat_bp, mat_cc, mat_mf), 20, 20)
W1 <- as.data.frame(W1)
colnames(W1) <- rownames(W1) <- gene.df$ENTREZID

# === Map to Standard Target ID ===
num_map <- fread("data/all_num.csv")  # expects columns: SYMBOL, num
symbol_map <- gene.df %>% left_join(num_map, by = c("SYMBOL" = "all"))

# Remove duplicates manually if needed (e.g., TEC)
W1 <- W1[rownames(W1) != "100124696", colnames(W1) != "100124696"]
final_ids <- symbol_map$num[symbol_map$ENTREZID %in% rownames(W1)]
colnames(W1) <- rownames(W1) <- final_ids
write.csv(W1, file = "output/target_similarity.csv")

# === Prepare Drug-Target Matrix and Drug Similarity ===
library(DTHybrid)
DTI_mat <- fread("data/DTIMatrix.csv")
drug_ids <- DTI_mat[1, -1] %>% as.character()
target_ids <- DTI_mat[-1, 1] %>% as.character()

# Clean drug format (e.g., CID padding)
targetdrug <- sapply(drug_ids, function(i) paste0("CIDs", str_pad(i, 8, pad = "0")))
DTI_matrix <- DTI_mat[-1, -1]
colnames(DTI_matrix) <- targetdrug
rownames(DTI_matrix) <- target_ids

# Load drug similarity matrix
drug_sim <- fread("data/drug_similarity.csv")
drug_sim <- as.data.frame(drug_sim[,-1])
rownames(drug_sim) <- colnames(drug_sim) <- fread("data/drug_similarity.csv")[[1]][-1]

# Align drug/target names
target_sim <- fread("output/target_similarity.csv")
target_order <- target_sim$V1[-1]
target_sim <- as.data.frame(target_sim[-1,-1])
rownames(target_sim) <- colnames(target_sim) <- target_order

DTI_matrix <- DTI_matrix[target_order, intersect(colnames(DTI_matrix), colnames(drug_sim))]
missing <- setdiff(colnames(drug_sim), colnames(DTI_matrix))
DTI_matrix <- cbind(DTI_matrix, matrix(0, nrow = nrow(DTI_matrix), ncol = length(missing)))
colnames(DTI_matrix)[(ncol(DTI_matrix)-length(missing)+1):ncol(DTI_matrix)] <- missing
DTI_matrix <- DTI_matrix[, colnames(drug_sim)]

# === Run DT-Hybrid ===
A <- t(as.matrix(DTI_matrix))
res <- computeRecommendation(DTI_matrix, lambda = 0.5, alpha = 0.5, target_sim, drug_sim, NA)
write.csv(res, file = "output/DT-Hybrid_results.csv")

# === Filter Off-target Interactions ===
res_df <- as.data.frame(res)
drug_ids <- colnames(res_df)
target_ids <- rownames(res_df)

predicted <- bind_rows(lapply(drug_ids, function(drug) {
  targets <- target_ids[res_df[, drug] > 0.1]
  data.frame(Drug = drug, Target = targets)
}))

original_DTI <- fread("data/DTI_Num.csv")
original_DTI$Drug <- paste0("CIDs", str_pad(original_DTI$Drug, 8, pad = "0"))
overlap <- semi_join(predicted, original_DTI, by = c("Drug", "Target"))
off_target <- anti_join(predicted, original_DTI, by = c("Drug", "Target"))

write.csv(predicted, file = "output/all_interactions.csv", row.names = FALSE)
write.csv(off_target, file = "output/off_target.csv", row.names = FALSE)
