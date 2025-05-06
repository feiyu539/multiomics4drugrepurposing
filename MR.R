# MR Pipeline: Drug Target Mapping, GWAS Integration, and MR Analysis
# Author: ldy-GitHub-ready version (personal paths and IDs generalized)
# Requirements: data.table, TwoSampleMR, devtools, ieugwasr

# === Setup ===
library(data.table)
library(tidyverse)
library(TwoSampleMR)
library(ieugwasr)

setwd("./6-MR")  # Use relative path for portability

# === Step 1: Prepare Drug-Target Pairs ===
results_df <- fread("FINAL_RESULTS.CSV")
thresh <- nrow(results_df) * 0.05
selected_drugs <- results_df[1:507, ]$eitadrugname %>% as.character()

DTI <- fread("Drugbank_DTI.csv")
DTI <- DTI[, -1]
colnames(DTI) <- c("drug", "target")

intersected_drugs <- intersect(selected_drugs, DTI$drug)
targets_df <- DTI[drug %in% intersected_drugs]
write.csv(targets_df, "needMR.csv", row.names = FALSE)

# === Step 2: Map Targets to Gene Symbols ===
mapping <- fread("idmapping_2023_12_12.tsv")
targets_df <- fread("needMR.csv")
targets_filtered <- targets_df[target %in% mapping$From]
mapped_targets <- left_join(targets_filtered, mapping, by = c("target" = "From"))

# === Step 3: Load GWAS and MAF Data ===
load("beta_se_data.RData")  # Contains FDR, MAF2
FDR <- as.data.frame(FDR)
MAF2 <- as.data.frame(MAF2)

genes <- unique(mapped_targets$To)
FDR_selected <- FDR[FDR$GeneSymbol %in% genes, ] %>% na.omit()
MAF_selected <- MAF2[unique(FDR_selected$SNP), ] %>% na.omit()

# === Step 4: Compute Beta and SE ===
beta_se_list <- vector("list", nrow(FDR_selected))
for (i in 1:nrow(FDR_selected)) {
  snp <- as.character(FDR_selected[i, "SNP"])
  z <- as.numeric(FDR_selected[i, "Zscore"])
  n <- as.numeric(FDR_selected[i, "N"])
  af <- as.numeric(MAF_selected[snp, 9])
  beta <- z / sqrt(2 * af * (1 - af) * (n + z^2))
  se <- 1 / sqrt(2 * af * (1 - af) * (n + z^2))
  beta_se_list[[i]] <- data.frame(snp = snp, beta = beta, se = se)
}
beta_se_df <- do.call(rbind, beta_se_list)
GWAS <- cbind(FDR_selected, beta_se_df)
write.csv(GWAS, "GWAS_processed.csv", row.names = FALSE)

# === Step 5: Format for TwoSampleMR ===
eQTL <- fread("GWAS_processed.csv")
IDs <- fread("BMDGWASID.csv")$GWASIDs

DBID <- fread("MRtargetdata.csv")
DBID <- data.frame(drug = DBID$drug, target = DBID$To)
unique_genes <- intersect(unique(DBID$target), unique(eQTL$GeneSymbol))

eQTL_filtered <- eQTL %>% filter(GeneSymbol %in% unique_genes)
DBID_filtered <- DBID %>% filter(target %in% unique_genes)

# === Step 6: Run MR ===
all_drugs <- unique(DBID_filtered$drug)
results_all <- list()

for (drug in all_drugs) {
  targets <- DBID_filtered %>% filter(drug == !!drug) %>% pull(target)
  gene_data <- eQTL_filtered %>% filter(GeneSymbol %in% targets)

  if (nrow(gene_data) == 0) next

  exposure_dat <- format_data(
    gene_data,
    type = "exposure",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "AssessedAllele",
    other_allele_col = "OtherAllele",
    eaf_col = "AlleleB_all",
    pval_col = "PValue"
  )
  exposure_dat <- clump_data(exposure_dat, clump_r2 = 0.8)

  for (id in IDs) {
    outcome_dat <- extract_outcome_data(snps = exposure_dat$SNP, outcomes = id)
    dat <- harmonise_data(exposure_dat = exposure_dat, outcome_dat = outcome_dat)
    res <- mr(dat, method_list = c("mr_ivw"))
    hetero <- mr_heterogeneity(dat, method_list = c("mr_ivw"))
    pleio <- mr_pleiotropy_test(dat)

    res <- cbind(
      res,
      heterogeneity.Q = hetero$Q,
      heterogeneity.Q_pval = hetero$Q_pval,
      pleiotropy.pval = pleio$pval,
      exposure = drug,
      id.exposure = drug
    )
    results_all[[length(results_all) + 1]] <- res
  }
  cat("Completed:", drug, "\n")
}

results_df <- do.call(rbind, results_all)
write.csv(results_df, "results.csv", row.names = FALSE)

# === Step 7: Filter and Annotate Results ===
threshold <- 0.05 / (length(IDs) * length(all_drugs))
filtered <- results_df %>%
  filter(pval < threshold, heterogeneity.Q_pval > 0.05, pleiotropy.pval > 0.05)

# Annotate with drug names
drug_info <- fread("infoFinal.csv")
drug_info <- drug_info %>% mutate(pubchem_cid = as.character(pubchem_cid))
filtered$DrugName <- drug_info[match(filtered$id.exposure, drug_info$pubchem_cid), ]$pert_iname

write.csv(filtered, "MRpositiveresults.csv", row.names = FALSE)
write.csv(unique(filtered$DrugName), "drugresults.csv", row.names = FALSE)
