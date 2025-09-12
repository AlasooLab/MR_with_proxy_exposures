library("MungeSumstats")
library("locuszoomr")
library("dplyr")
library("EnsDb.Hsapiens.v86")
library("tidyr")
library("mrlocus")

#Import skin cancer sumstats
sumstats = readr::read_tsv("big_data/skin cancer/GCST90137411_buildGRCh37.tsv.gz")
munged = dplyr::transmute(sumstats, SNP = variant_id, CHR = chromosome, BP = base_pair_location, 
                          A1 = other_allele, A2 = effect_allele, BETA = beta, SE = standard_error, p_value)
skin_cancer = MungeSumstats::liftover(munged, ref_genome = "GRCh37", convert_ref_genome = "GRCh38") 

#Make LocusZoom
HAL_skin_cancer = dplyr::transmute(skin_cancer, chrom = CHR, pos = BP, rsid = SNP, other_allele = A1, effect_allele = A2, p = p_value, beta = BETA, se = SE, variant = SNP) %>%
  as.data.frame()
loc <- locus(data = HAL_skin_cancer, gene = 'HAL', flank = 1e5, ens_db = "EnsDb.Hsapiens.v86")
summary(loc)
pdf("figures/HAL/HAL_skin_cancer.pdf", width = 5.5, height = 4.5)
locus_plot(loc)
dev.off()

#Import trans-urocanate sumstats
sumstats2 = readr::read_tsv("big_data/HAL_coloc/trans_urocanate_607_upload_filtered.tsv.gz")
munged2 = dplyr::transmute(sumstats2, SNP, CHR = chr, BP = position, A1 = coded_all, A2 = noncoded_all, BETA = beta, SE, p_value = pval)
trans_urocanate = MungeSumstats::liftover(munged2, ref_genome = "GRCh37", convert_ref_genome = "GRCh38") 

#Make LocusZoom
HAL_trans_urocanate = dplyr::transmute(trans_urocanate, chrom = CHR, pos = BP, rsid = SNP, other_allele = A1, effect_allele = A2, p = p_value, beta = BETA, se = SE, variant = SNP) %>%
  as.data.frame()
loc <- locus(data = HAL_trans_urocanate, gene = 'HAL', flank = 1e5, ens_db = "EnsDb.Hsapiens.v86", index_snp = "rs3819817")
summary(loc)
pdf("figures/HAL/HAL_transurocanate.pdf", width = 5.5, height = 4.5)
locus_plot(loc)
dev.off()

#Import sunburns
sumstats = readr::read_tsv("big_data/HAL_coloc/childhood_sunburns_1737.gwas.imputed_v3.both_sexes.tsv.bgz")
split = tidyr::separate(sumstats, variant, c("CHR","BP", "A1", "A2"), sep = ":", remove = FALSE)
munged3 = dplyr::transmute(split, SNP = variant, CHR, BP, A1, A2, BETA = beta, SE = se, p_value = pval)
sunburns = MungeSumstats::liftover(munged3, ref_genome = "GRCh37", convert_ref_genome = "GRCh38") %>%
  dplyr::mutate(SNP = paste0("chr", CHR, "_", BP,"_",A1,"_", A2))

#Make LocusZoom
HAL_sunburns = dplyr::transmute(sunburns, chrom = CHR, pos = BP, rsid = SNP, other_allele = A1, effect_allele = A2, p = p_value, beta = BETA, se = SE, variant = SNP) %>%
  as.data.frame()
loc <- locus(data = HAL_sunburns, gene = 'HAL', flank = 1e5, ens_db = "EnsDb.Hsapiens.v86", index_snp = "chr12_95984993_C_T")
summary(loc)
pdf("figures/HAL/HAL_sunburns.pdf", width = 5.5, height = 4.5)
locus_plot(loc)
dev.off()

#Vitamin D
vitD = readr::read_tsv("big_data/masa_sumstats/UKBB.VitD.all.tsv.gz")

#Make LocusZoom
HAL_vitD = dplyr::transmute(vitD, chrom = chromosome, pos = position, rsid = variant, other_allele = ref, effect_allele = alt, p = pvalue, beta = beta, se = se, variant = variant) %>%
  as.data.frame()
loc <- locus(data = HAL_vitD, gene = 'HAL', flank = 1e5, ens_db = "EnsDb.Hsapiens.v86", index_snp = "chr12_95984993_C_T")
summary(loc)
pdf("figures/HAL/HAL_vitD.pdf", width = 5.5, height = 4.5)
locus_plot(loc)
dev.off()

#Histidine
sumstats4 = readr::read_tsv("big_data/UKBB_sumstats/His_full_regenie.tsv")
munged4 = dplyr::transmute(sumstats4, SNP = ID, CHR = CHROM, BP = GENPOS, A1 = ALLELE0, A2 = ALLELE1, BETA = BETA, SE = SE, p_value = 10^-LOG10P)
His = MungeSumstats::liftover(munged4, ref_genome = "GRCh37", convert_ref_genome = "GRCh38") 

HAL_His = dplyr::transmute(His, chrom = CHR, pos = BP, rsid = SNP, other_allele = A1, effect_allele = A2, p = p_value, beta = BETA, se = SE, variant = SNP) %>%
  as.data.frame()
loc <- locus(data = HAL_His, gene = 'HAL', flank = 1e5, ens_db = "EnsDb.Hsapiens.v86", index_snp = "rs3819817")
summary(loc)
pdf("figures/HAL/HAL_His.pdf", width = 5.5, height = 4.5)
locus_plot(loc)
dev.off()

#HAL eQTL in TwinsUK
HAL = readr::read_tsv("big_data/HAL_coloc/QTD000544.cc.tsv") %>%
  dplyr::filter(molecular_trait_id == "ENSG00000084110")

HAL_eQTL = dplyr::transmute(HAL, chrom = chromosome, pos = position, rsid = rsid, other_allele = ref, effect_allele = alt, p = pvalue, beta = beta, se = se, variant = variant) %>%
  as.data.frame()
loc <- locus(data = HAL_eQTL, gene = 'HAL', flank = 1e5, ens_db = "EnsDb.Hsapiens.v86", index_snp = "rs3819817")
summary(loc)
pdf("figures/HAL/HAL_TwinsUK_eQTL.pdf", width = 5.5, height = 4.5)
locus_plot(loc)
dev.off()

#Make a forest plot of of the eQTL variant effect sizes
library("ggforestplot")
data = readr::read_tsv("data/HAL_eQTL_betas.tsv")
data = dplyr::mutate(data, trait = factor(trait, levels = c("skin cancer", "childhood sunburns", "vitamin D", "trans-urocanate", "histidine", "HAL expression (skin)"))) %>%
  dplyr::arrange(desc(trait))
plot = ggforestplot::forestplot(df = data, name = trait, estimate = beta, se = se)
ggsave("figures/HAL/HAL_eQTL_forest.pdf", plot = plot, width = 4, height = 3.5)

