library("dplyr")
library("eQTLUtils")

#Create a mapping table between UT and NG UKBB sample ids
ng_samples = readr::read_delim("private_data/ng_to_ut_sample_maps/ng_sample_file.sample", delim = " ")
ut_samples = readr::read_delim("private_data/ng_to_ut_sample_maps/ut_sample_file.sample", delim = " ")

match = dplyr::tibble(eid_30418 = ng_samples$ID_1, eid = ut_samples$ID_1) %>%
  dplyr::mutate(eid_30418 = as.character(eid_30418), eid = as.character(eid))

#Import sample lists on the X chromosome (smallest subset of them all)
gt_samples = readr::read_delim("private_data/ukb22828_cX_b0_v3_s486535.sample", delim = " ") %>% 
  dplyr::mutate(ID_1 = as.character(ID_1))

#Import metabolite data
phase1 = readr::read_tsv("private_data/ukb-phase1__results.tsv") %>% 
  dplyr::mutate(eid_30418 = as.character(eid_30418))
phase2 = readr::read_tsv("private_data/ukb-phase2__results.tsv") %>% 
  dplyr::mutate(eid_30418 = as.character(eid_30418))

#Extract Main Phase
metabolites = dplyr::bind_rows(phase1,phase2) %>% 
  dplyr::filter(visit == "Main Phase") %>% 
  dplyr::left_join(match, by = c("eid_30418")) %>%
  dplyr::select(-eid_30418) %>%
  dplyr::select(eid, everything()) %>%
  dplyr::arrange(eid)

#Remove blind duplicates
unique_samples = dplyr::select(metabolites, eid, sample_id) %>% 
  dplyr::group_by(eid) %>%
  dplyr::filter(row_number() == 1) %>%
  dplyr::ungroup()
metabolites_uniq = dplyr::semi_join(metabolites, unique_samples, by = c("eid","sample_id")) %>%
  dplyr::group_by(eid) %>%
  dplyr::filter(row_number() == 1) %>%
  dplyr::ungroup()
  
#Remove samples without GT data
metabolites_gt = dplyr::filter(metabolites_uniq, eid %in% gt_samples$ID_1)

#Convert to matrix
metabolite_matrix = as.matrix(dplyr::select(metabolites_gt, -eid, -sample_id, -plate_id, -plate_position, -visit, -spectrometer))
rownames(metabolite_matrix) = metabolites_gt$eid

#Remove metabolites with more than 8000 missing values
metabolite_matrix = metabolite_matrix[,!colSums(is.na(metabolite_matrix)) > 8000]
#Remove individuals with more than 5 missing values
metabolite_matrix = metabolite_matrix[!(rowSums(is.na(metabolite_matrix)) > 5), ]

#Replace missing values with median values for the same metabolite
log_matrix = log(metabolite_matrix + 1)
imputed_matrix = zoo::na.aggregate(log_matrix, FUN = median)

#Perform INT on the columns of the matrix
metabolites_norm = eQTLUtils::quantileNormaliseCols(imputed_matrix)

#Export data frame
metabolites_df = dplyr::as_tibble(round(metabolites_norm, 4)) %>%
  dplyr::mutate(eid = rownames(metabolites_norm)) %>%
  dplyr::select(eid, everything())
write.table(metabolites_df, "private_data/ukbb_nightingale_NMR_phase1+2_INT_270423.tsv", row.names = F, quote = F, sep = "\t")
