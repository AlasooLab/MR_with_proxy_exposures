library("dplyr")
library("ggplot2")
library("ggrepel")
library("mrlocus")

old_beta_df = readr::read_tsv("data/glucose_pyruvate_missense_betas.tsv") %>%
  dplyr::mutate(missense_trait = c("Pyruvate", "Glucose","Glucose","Glucose","Glucose","Glucose","Pyruvate","Pyruvate","Pyruvate","Pyruvate","Pyruvate","Pyruvate","Pyruvate")) %>%
  dplyr::mutate(glucolysis_pathway = c(0,0,0,0,0,0,1,1,0,0,1,1,0)) %>%
  dplyr::select(variant, gene_name, missense_trait, glucolysis_pathway)

beta_df = readr::read_tsv("data/glucose_pyruvate_missense.tsv") %>%
  dplyr::left_join(old_beta_df, by = "variant")
  

#Pyruvate missense variants
pyruvate_df = dplyr::filter(beta_df, missense_trait == "Pyruvate") %>%
  dplyr::mutate(HbA1c_beta = HbA1c_beta*sign(pyruvate_beta)) %>%
  dplyr::mutate(Hb_beta = Hb_beta*sign(pyruvate_beta)) %>%
  dplyr::mutate(Ht_beta = Ht_beta*sign(pyruvate_beta)) %>%
  dplyr::mutate(T2D_beta = T2D_beta*sign(pyruvate_beta)) %>%
  dplyr::mutate(RBC_beta = RBC_beta*sign(pyruvate_beta)) %>%
  dplyr::mutate(glucose_beta = RBC_beta*sign(glucose_beta)) %>%
  dplyr::mutate(pyruvate_beta = pyruvate_beta*sign(pyruvate_beta))


###### Pyruvate vs RBC #####
pyr_plot = ggplot(pyruvate_df, aes(x = pyruvate_beta, y = RBC_beta, label = gene_name)) +
  geom_point() +
  geom_errorbar(aes(xmin=pyruvate_beta-pyruvate_se, xmax=pyruvate_beta+pyruvate_se)) +
  geom_errorbar(aes(ymin = RBC_beta-RBC_se, ymax = RBC_beta+RBC_se)) +
  xlim(0,0.5) +
  ylim(-0.1,0.25) + 
  geom_label_repel()
ggsave("figures/pyruvate_vs_RBC_all.pdf", plot = pyr_plot, width = 4, height = 4)

#Pyruvate vs RBC (glucolysis pathway)
pyr_glycolysis = dplyr::filter(pyruvate_df, glucolysis_pathway == 1)
pyr_plot2 = ggplot(pyr_glycolysis, aes(x = pyruvate_beta, y = RBC_beta, label = gene_name)) +
  geom_point() +
  geom_errorbar(aes(xmin=pyruvate_beta-pyruvate_se, xmax=pyruvate_beta+pyruvate_se)) +
  geom_errorbar(aes(ymin = RBC_beta-RBC_se, ymax = RBC_beta+RBC_se)) +
  xlim(0,0.5) +
  ylim(-0.1,0.25) + 
  geom_label_repel()
ggsave("figures/pyruvate_vs_RBC_glyc.pdf", plot = pyr_plot2, width = 4, height = 4)

#Run MRLocus on all pyruvate missense variants
mrlocus_object = list(
  beta_hat_a = pyruvate_df$pyruvate_beta,
  beta_hat_b = pyruvate_df$RBC_beta,
  sd_a = pyruvate_df$pyruvate_se,
  sd_b = pyruvate_df$RBC_se,
  alleles = data.frame(id = pyruvate_df$variant,
                       ref = pyruvate_df$ref,
                       eff = pyruvate_df$alt)
)
mrlocus_fit = mrlocus::fitSlope(mrlocus_object, iter = 100000)
pdf("figures/pyruvate_RBC_all_variants_MRLocus.pdf", width = 4.5, height = 4.5)
plotMrlocus(mrlocus_fit, a = "pyruvate", b = "RBC", label = "Variant effect on", legend = T)
dev.off()

#Run MRLocus on glyc pathway
mrlocus_object = list(
  beta_hat_a = pyr_glycolysis$pyruvate_beta,
  beta_hat_b = pyr_glycolysis$RBC_beta,
  sd_a = pyr_glycolysis$pyruvate_se,
  sd_b = pyr_glycolysis$RBC_se,
  alleles = data.frame(id = pyr_glycolysis$variant,
                       ref = pyr_glycolysis$ref,
                       eff = pyr_glycolysis$alt)
)
mrlocus_fit = mrlocus::fitSlope(mrlocus_object, iter = 100000)
pdf("figures/pyruvate_glycolysis_MRLocus.pdf", width = 4.5, height = 4.5)
plotMrlocus(mrlocus_fit, a = "glyocolysis", b = "RBC", label = "Variant effect on", legend = T)
dev.off()

#Now repeat this for non-pathway genes
pyr_not_glycolysis = dplyr::filter(pyruvate_df, glucolysis_pathway == 0)

mrlocus_object = list(
  beta_hat_a = pyr_not_glycolysis$pyruvate_beta,
  beta_hat_b = pyr_not_glycolysis$RBC_beta,
  sd_a = pyr_not_glycolysis$pyruvate_se,
  sd_b = pyr_not_glycolysis$RBC_se,
  alleles = data.frame(id = pyr_not_glycolysis$variant,
                       ref = pyr_not_glycolysis$ref,
                       eff = pyr_not_glycolysis$alt)
)
mrlocus_fit = mrlocus::fitSlope(mrlocus_object, iter = 100000)
pdf("figures/pyruvate_not_glycolysis_MRLocus.pdf", width = 5.5, height = 5.5)
plotMrlocus(mrlocus_fit, a = "pyruvate", b = "RBC", label = "Variant effect on", legend = T)
dev.off()


#Replicate the glyc pathway analysis using metabolite data from the EstBB
pyr_glycolysis
estbb_pyruvate = arrow::read_parquet("big_data/Pyruvate_EstBB.parquet")
extracted_variants = dplyr::filter(estbb_pyruvate, GENPOS %in% c(155291845, 155291918, 48118502, 44326728))

mrlocus_object

mrlocus_object = list(
  beta_hat_a = c(0.247, 0.112),
  beta_hat_b = c(0.102870, 0.027256),
  sd_a = c(0.0550, 0.00409),
  sd_b = c(0.02811820, 0.00363329),
  alleles = data.frame(id = pyr_glycolysis$variant,
                       ref = pyr_glycolysis$ref,
                       eff = pyr_glycolysis$alt)
)
mrlocus_fit = mrlocus::fitSlope(mrlocus_object, iter = 100000)
pdf("figures/pyruvate_RBC_EstBB_MRLocus.pdf", width = 5, height = 5)
plotMrlocus(mrlocus_fit, a = "pyruvate", b = "RBC", label = "Variant effect on", legend = T)
dev.off()


#Pyruvate vs Hb (glucolysis pathway)
pyr_plot2 = ggplot(dplyr::filter(pyruvate_df, glucolysis_pathway == 1), aes(x = pyruvate_beta, y = Hb_beta, label = gene_name)) +
  geom_point() +
  geom_errorbar(aes(xmin=pyruvate_beta-pyruvate_se, xmax=pyruvate_beta+pyruvate_se)) +
  geom_errorbar(aes(ymin = Hb_beta-Hb_se, ymax = Hb_beta+Hb_se)) +
  xlim(0,0.5) +
  ylim(-0.1,0.25) + 
  geom_label_repel()
ggsave("figures/pyruvate_vs_Hb_glyc.pdf", plot = pyr_plot2, width = 4, height = 4)


#Pyruvate vs HbA1c_beta (glucolysis pathway)
pyr_plot2 = ggplot(dplyr::filter(pyruvate_df, glucolysis_pathway == 1), aes(x = pyruvate_beta, y = HbA1c_beta, label = gene_name)) +
  geom_point() +
  geom_errorbar(aes(xmin=pyruvate_beta-pyruvate_se, xmax=pyruvate_beta+pyruvate_se)) +
  geom_errorbar(aes(ymin = HbA1c_beta-HbA1c_se, ymax = HbA1c_beta+HbA1c_se)) +
  xlim(0,0.5) +
  ylim(-0.25,0.1) + 
  geom_label_repel()
ggsave("figures/pyruvate_vs_HbA1c_glyc.pdf", plot = pyr_plot2, width = 4, height = 4)


#Pyruvate vs HbA1c_beta (glucolysis pathway)
pyr_plot2 = ggplot(dplyr::filter(pyruvate_df, glucolysis_pathway == 1), aes(x = pyruvate_beta, y = T2D_beta, label = gene_name)) +
  geom_point() +
  geom_errorbar(aes(xmin=pyruvate_beta-pyruvate_se, xmax=pyruvate_beta+pyruvate_se)) +
  geom_errorbar(aes(ymin = T2D_beta-T2D_se, ymax = T2D_beta+T2D_se)) +
  xlim(0,0.5) +
  ylim(-0.5,0.5) + 
  geom_label_repel()
ggsave("figures/pyruvate_vs_T2D_glyc.pdf", plot = pyr_plot2, width = 4, height = 4)


#Pyruvate vs Glucose (glucolysis pathway)
pyr_plot2 = ggplot(dplyr::filter(pyruvate_df, glucolysis_pathway == 1), aes(x = pyruvate_beta, y = glucose_beta, label = gene_name)) +
  geom_point() +
  geom_errorbar(aes(xmin=pyruvate_beta-pyruvate_se, xmax=pyruvate_beta+pyruvate_se)) +
  geom_errorbar(aes(ymin = glucose_beta-glucose_se, ymax = glucose_beta+glucose_se)) +
  xlim(0,0.5) +
  ylim(-0.5,0.5) + 
  geom_label_repel()


#Replication in EstBB
pyruvate_estbb = arrow::read_parquet("big_data/Pyruvate_EstBB.parquet")
selected_vars = dplyr::filter(pyruvate_estbb, GENPOS %in% c(27508073, 155291845,155291918,137205865,10487356,48118502,44326728,24503382))

