# GSEA + TE Scoring Pipeline
# Cleaned and combined version for GitHub

# === Load Required Packages ===
library(data.table)

# === GSEA Enrichment Score Function ===
GSEA.EnrichmentScore <- function(gene.list, gene.set, weighted.score.type = 1, correl.vector = NULL) {
  tag.indicator <- sign(match(gene.list, gene.set, nomatch = 0))
  no.tag.indicator <- 1 - tag.indicator
  N <- length(gene.list)
  Nh <- length(gene.set)
  Nm <- N - Nh
  if (weighted.score.type == 0) correl.vector <- rep(1, N)
  alpha <- weighted.score.type
  correl.vector <- abs(correl.vector^alpha)
  sum.correl.tag <- sum(correl.vector[tag.indicator == 1])
  norm.tag <- 1 / sum.correl.tag
  norm.no.tag <- 1 / Nm
  RES <- cumsum(tag.indicator * correl.vector * norm.tag - no.tag.indicator * norm.no.tag)
  max.ES <- max(RES)
  min.ES <- min(RES)
  if (max.ES > -min.ES) {
    ES <- signif(max.ES, digits = 5)
    arg.ES <- which.max(RES)
  } else {
    ES <- signif(min.ES, digits = 5)
    arg.ES <- which.min(RES)
  }
  return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator))
}

# === Load Data ===
mcf_expr <- fread("data/drug_expression.csv")  # Gene expression matrix, columns: drugs, rows: genes
mcf_expr <- na.omit(mcf_expr)
druggene <- as.character(mcf_expr[[1]])
mcf_expr <- mcf_expr[, -1]
colnames(mcf_expr) <- make.unique(colnames(mcf_expr))
druglist <- colnames(mcf_expr)

# Create ranked gene lists per drug
list_muscMCF <- lapply(1:ncol(mcf_expr), function(i) {
  gene_order <- order(mcf_expr[[i]], decreasing = TRUE)
  drug_ordered <- druggene[gene_order]
  return(drug_ordered)
})
names(list_muscMCF) <- druglist

# === Load Driver Network ===
driver_data <- fread("data/network_detail.csv")
driver_names <- driver_data[[1]]
driver_data <- as.data.frame(driver_data[, -1])
rownames(driver_data) <- driver_names

# === Compute ES matrix ===
gene.list <- as.data.frame(list_muscMCF)
gene.set <- t(driver_data)

drug_work_ES <- matrix(NA, nrow = length(druglist), ncol = ncol(gene.set))
rownames(drug_work_ES) <- druglist
colnames(drug_work_ES) <- colnames(gene.set)

for (d in 1:length(druglist)) {
  genelist <- gene.list[[d]]
  for (g in 1:ncol(gene.set)) {
    geneset <- gene.set[, g]
    geneset <- geneset[geneset != 0]
    es <- GSEA.EnrichmentScore(genelist, geneset, weighted.score.type = 0, correl.vector = NULL)$ES
    drug_work_ES[d, g] <- es
  }
}
write.csv(drug_work_ES, file = "output/ES_results.csv")

# === Compute TE Scores ===
es_matrix <- fread("output/ES_results.csv")
es_drugs <- es_matrix[[1]]
es_matrix <- as.data.frame(es_matrix[, -1])
rownames(es_matrix) <- es_drugs

eita_matrix <- fread("data/eita.csv")
eita_drugs <- eita_matrix[[1]]
eita_matrix <- as.data.frame(eita_matrix[, -1])
rownames(eita_matrix) <- eita_drugs

common_drugs <- intersect(rownames(es_matrix), rownames(eita_matrix))
network_list <- colnames(eita_matrix)

es_eita <- matrix(NA, nrow = length(common_drugs), ncol = length(network_list))
rownames(es_eita) <- common_drugs
colnames(es_eita) <- network_list

for (j in common_drugs) {
  for (i in network_list) {
    es_eita[j, i] <- eita_matrix[j, i] * es_matrix[j, i]
  }
}

# Sum across networks to get TE score
TE <- rowSums(es_eita, na.rm = TRUE)
TE <- data.frame(Drug = common_drugs, TE = TE)
write.csv(TE, file = "output/TE_scores.csv", row.names = FALSE)

# Merge with drug info (if available)
if (file.exists("data/infoFinal.csv")) {
  lincsinfo <- fread("data/infoFinal.csv")
  final_results <- merge(lincsinfo, TE, by.x = "pert_iname", by.y = "Drug", all.y = TRUE)
  write.csv(final_results, file = "output/FINAL_RESULTS.csv", row.names = FALSE)
}
