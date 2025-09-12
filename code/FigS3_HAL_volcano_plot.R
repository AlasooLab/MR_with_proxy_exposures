library("dplyr")
library("ggplot2")
library("ggrepel")

#Import dataset metadata
dataset_metadata = readr::read_tsv("https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/data_tables/dataset_metadata.tsv")

#HAL skin eQTL volcano plot
betas = readr::read_tsv("data/HAL_skin_eQTL.tsv") %>%
  dplyr::left_join(dplyr::select(dataset_metadata, dataset_id, study_label, tissue_label), by = "dataset_id")

label_df = dplyr::filter(betas, pvalue < 1e-5) %>% 
  dplyr::mutate(plot_label = paste(study_label, tissue_label))

hal_plot = ggplot(betas, aes(x = beta, y = -log(pvalue, 10))) + 
  geom_point() + 
  geom_label_repel(aes(label = plot_label), data = label_df) + 
  theme_light() +
#  ggtitle("Effect of chr12_95984993_C_T on HAL expression") +
  ylab("-log10 pvalue")
ggsave("figures/HAL_skin_eQTL.pdf", hal_plot, width = 2.3, height = 2.3)
