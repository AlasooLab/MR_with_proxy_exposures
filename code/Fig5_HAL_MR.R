library("MungeSumstats")
library("locuszoomr")
library("dplyr")
library("EnsDb.Hsapiens.v86")
library("tidyr")
library("mrlocus")

#Read all HAL missense betas
dat = readr::read_tsv("data/glucose_pyruvate_histidine_missense.tsv")
HAL_variants = dplyr::filter(dat, gene_name == "HAL") %>% 
  dplyr::select(variant, chromosome, ref, alt, His_beta, His_se, VitD_beta, VitD_se)


#Estimating the causal effect of HAL activity on VitD levels (including the skin eQTL variant + 3rd missense)
res = list(
  beta_hat_a = c(0.827, 0.492, 0.0137711, 0.664),
  beta_hat_b = c(-0.12841, -0.10174,-0.0414, 0.0222),
  sd_a = c(0.0177, 0.0203, 0.00288072, 0.0350),
  sd_b = c(0.015054, 0.017242,  0.00245, 0.0303),
  alleles = data.frame(id = c("chr12_95977953_C_T","chr12_95986106_C_T", "chr12_95984993_C_T"),
                       ref = c("C","C","C"),
                       eff = c("T","T","T"))
)
hal_coding_fit = mrlocus::fitSlope(res, iter = 100000)
pdf("figures/HAL/HAL_missense_3+eQTL_MRLocus.pdf", width = 5, height = 5)
plotMrlocus(hal_coding_fit, a = "His", b = "VitD", label = "Variant effect on", legend = T)
dev.off()

#Estimating the causal effect of HAL activity on VitD levels (3 missense)
res = list(
  beta_hat_a = c(0.827, 0.492, 0.664),
  beta_hat_b = c(-0.12841, -0.10174, 0.0222),
  sd_a = c(0.0177, 0.0203, 0.0350),
  sd_b = c(0.015054, 0.017242, 0.0303),
  alleles = data.frame(id = c("chr12_95977953_C_T","chr12_95986106_C_T", "chr12_95984993_C_T"),
                       ref = c("C","C","C"),
                       eff = c("T","T","T"))
)
hal_coding_fit = mrlocus::fitSlope(res, iter = 100000)
pdf("figures/HAL/HAL_missense_3.pdf", width = 5, height = 5)
plotMrlocus(hal_coding_fit, a = "His", b = "VitD", label = "Variant effect on", legend = T)
dev.off()


#Estimating the causal effect of HAL activity on VitD levels (including the skin eQTL variant)
res = list(
  beta_hat_a = c(0.827, 0.492),
  beta_hat_b = c(-0.12841, -0.10174),
  sd_a = c(0.0177, 0.0203),
  sd_b = c(0.015054, 0.017242),
  alleles = data.frame(id = c("chr12_95977953_C_T","chr12_95986106_C_T", "chr12_95984993_C_T"),
                       ref = c("C","C","C"),
                       eff = c("T","T","T"))
)
hal_coding_fit = mrlocus::fitSlope(res, iter = 100000)
pdf("figures/HAL/HAL_missense_2_MRLocus.pdf", width = 5, height = 5)
plotMrlocus(hal_coding_fit, a = "His", b = "VitD", label = "Variant effect on", legend = T)
dev.off()


#Replicate MR using betas from EstBB
his_estbb = arrow::read_parquet("big_data/His_EstBB.parquet")
dplyr::filter(his_estbb, CHROM == 12, GENPOS %in% c(95977953, 95986106, 95994812))

#Estimating the causal effect of HAL activity on VitD levels (including the skin eQTL variant)
res = list(
  beta_hat_a = c(0.816, 0.474),
  beta_hat_b = c(-0.12841, -0.10174),
  sd_a = c(0.0496, 0.0329),
  sd_b = c(0.015054, 0.017242),
  alleles = data.frame(id = c("chr12_95977953_C_T","chr12_95986106_C_T", "chr12_95984993_C_T"),
                       ref = c("C","C","C"),
                       eff = c("T","T","T"))
)
hal_coding_fit = mrlocus::fitSlope(res, iter = 100000)
pdf("figures/HAL/HAL_EstBB_missense_MRLocus.pdf", width = 5, height = 5)
plotMrlocus(hal_coding_fit, a = "His", b = "VitD", label = "Variant effect on", legend = T)
dev.off()


#Extract betas for His from our larger UKB_EUR GWAS:
ukb_eur = readr::read_delim("big_data/EUR_assoc_chr12_His.regenie.gz", delim =" ")
dplyr::filter(ukb_eur, CHROM == 12, GENPOS %in% c(95977953, 95986106, 95994812))


#Estimate the causal effect of HAL disruption on skin cancer (UKB + MVP + FinnGen meta-analysis)
res = list(
  beta_hat_a = c(0.816, 0.474),
  beta_hat_b = c(-0.155, -0.118),
  sd_a = c(0.0496, 0.0329),
  sd_b = c(0.0316, 0.0338),
  alleles = data.frame(id = c("chr12_95977953_C_T","chr12_95986106_C_T", "chr12_95984993_C_T"),
                       ref = c("C","C","C"),
                       eff = c("T","T","T"))
)
hal_coding_fit = mrlocus::fitSlope(res, iter = 100000)
pdf("figures/HAL/HAL_His_vs_skin_cancer_missense_MRLocus.pdf", width = 5, height = 5)
plotMrlocus(hal_coding_fit, a = "HAL function \n(poxied by reduction in plasma histidine)", b = "skin cancer", label = "Variant effect on", legend = T)
dev.off()

#Wald ratio based on the eQTL variant
-0.0339/0.49


#OLD!Estimating the causal effect of HAL activity on VitD levels (including the skin eQTL variant)
res = list(
  beta_hat_a = c(0.866349, 0.55, 0.0137711),
  beta_hat_b = c(-0.12841, -0.10174,-0.0414),
  sd_a = c(0.0283073, 0.03, 0.00288072),
  sd_b = c(0.015054, 0.017242,  0.00245),
  alleles = data.frame(id = c("chr12_95977953_C_T","chr12_95986106_C_T", "chr12_95984993_C_T"),
                       ref = c("C","C","C"),
                       eff = c("T","T","T"))
)
hal_coding_fit = mrlocus::fitSlope(res, iter = 100000)
pdf("figures/HAL/HAL_eQTL+missense_MRLocus.pdf", width = 5, height = 5)
plotMrlocus(hal_coding_fit, a = "His", b = "VitD", label = "Variant effect on", legend = T)
dev.off()

#OLD!Estimating the causal effect of HAL activity on VitD levels (including the skin eQTL variant)
res = list(
  beta_hat_a = c(0.866349, 0.55),
  beta_hat_b = c(-0.12841, -0.10174),
  sd_a = c(0.0283073, 0.03),
  sd_b = c(0.015054, 0.017242),
  alleles = data.frame(id = c("chr12_95977953_C_T","chr12_95986106_C_T"),
                       ref = c("C","C"),
                       eff = c("T","T"))
)
hal_coding_fit = mrlocus::fitSlope(res, iter = 100000)
pdf("figures/HAL/HAL_missense_MRLocus.pdf", width = 5, height = 5)
plotMrlocus(hal_coding_fit, a = "His", b = "VitD", label = "Variant effect on", legend = T)
dev.off()