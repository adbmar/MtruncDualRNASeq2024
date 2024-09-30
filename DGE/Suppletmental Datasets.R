# Must run GO Enrichment Analysis.R before running this script

dir_out_datasets <- file.path(dir_out, "datasets")
dir.create(dir_out_datasets, showWarnings = FALSE)
setwd(dir_out_datasets)
### Dataset_0S1 Gene expression quantification data
counts_out <- bind_rows(medicago_counts %>% mutate(organism = "Medicago truncatula"),
                        rhizob21_counts %>% mutate(organism = "Ensifer meliloti strain USDA1021"),
                        rhizob22_counts %>% mutate(organism = "Ensifer meliloti strain WSM1022"),
                    nematode_counts %>% mutate(organism = "Meloidogyne hapla"))
write.csv(counts_out, file = "Dataset_S01.csv", row.names = FALSE)

### Dataset_S02 Expression status of genes from all organisms in all sample types
exp_out <- full_join(
  counts_out %>% select("gene", "organism"),
  t %>% select("gene", starts_with("exp_") & !starts_with("exp_s")),
  by = "gene")
write.csv(exp_out, file = "Dataset_S02.csv", row.names = FALSE)

### Dataset_S03 Differential expression statistics for Medicago genes
med_DE <- t %>% filter(organism == "Medicago") %>% select("gene", starts_with("med_"), "NCR") %>%
  select(!ends_with(c(".LFC", ".SE", ".sig", ".significance", ".sigu", ".sigd", "sigu_ns", "sigd_ns")))
write.csv(med_DE, file = "Dataset_S03.csv", row.names = FALSE)

### Dataset_S04 Differential expression statistics of Ensifer meliloti genes in all sample types
rhi_DE <- t %>% filter(organism == "Rhizob21" | organism == "Rhizob22") %>% select("gene", starts_with("r2")) %>%
  select(!ends_with(c(".LFC", ".SE", ".sig", ".significance", ".sigu", ".sigd", "sigu_ns", "sigd_ns")))
write.csv(rhi_DE, file = "Dataset_S04.csv", row.names = FALSE)

### Dataset_S05. Differential expression statistics of Meloidogyne hapla genes in all sample types
nem_DE <-t %>% filter(organism == "Nematode") %>% select("gene", starts_with("nem_")) %>%
  select(!ends_with(c(".LFC", ".SE", ".sig", ".significance", ".sigu", ".sigd", "sigu_ns", "sigd_ns")))
write.csv(nem_DE, file = "Dataset_S05.csv", row.names = FALSE)

### Dataset_S06. All gene ontology over representation analysis results
enrich_go_summation <- go_gall_r %>% mutate(analysis = "Response to rhizobial difference in galls")
enrich_go_summation <- rbind(enrich_go_summation, go_gall_r_up   %>% mutate(analysis = "Response to rhizobial differences in galls -- upregulated in USDA1021"))
enrich_go_summation <- rbind(enrich_go_summation, go_gall_r_down %>% mutate(analysis = "Response to rhizobial differences in galls -- upregulated in WSM1022"))
enrich_go_summation <- rbind(enrich_go_summation, go_nodp_r      %>% mutate(analysis = "Response to rhizobial differences in nodules from hosts with parasites"))
enrich_go_summation <- rbind(enrich_go_summation, go_nodp_r_up   %>% mutate(analysis = "Response to rhizobial differences in nodules from hosts with parasites -- upregulated in USDA1021"))
enrich_go_summation <- rbind(enrich_go_summation, go_nodp_r_down %>% mutate(analysis = "Response to rhizobial differences in nodules from hosts with parasites -- upregulated in WSM1022"))
enrich_go_summation <- rbind(enrich_go_summation, go_nods_r21_n      %>% mutate(analysis = "Response to nematode infection in nodules with USDA1021"))
enrich_go_summation <- rbind(enrich_go_summation, go_nods_r21_n_up   %>% mutate(analysis = "Response to nematode infection in nodules with USDA1021 -- upregulation in infected hosts"))
enrich_go_summation <- rbind(enrich_go_summation, go_nods_r21_n_down %>% mutate(analysis = "Response to nematode infection in nodules with USDA1021 -- upregulated in uninfected hosts"))
enrich_go_summation <- rbind(enrich_go_summation, go_nods_r22_n      %>% mutate(analysis = "Response to nematode infection in nodules with WSM1022"))
enrich_go_summation <- rbind(enrich_go_summation, go_nods_r22_n_up   %>% mutate(analysis = "Response to nematode infection in nodules with WSM1022 -- upregulated in infected hosts"))
enrich_go_summation <- rbind(enrich_go_summation, go_nods_r22_n_down %>% mutate(analysis = "Response to nematode infection in nodules with WSM1022 -- upregulated in uninfected hosts"))
enrich_go_summation <- rbind(enrich_go_summation, go_nods_r      %>% mutate(analysis = "Similar response to rhizobial differences in nodules from infected and uninfected hosts"))
enrich_go_summation <- rbind(enrich_go_summation, go_nods_r_up   %>% mutate(analysis = "Similar response to rhizobial differences in nodules from infected and uninfected hosts -- upregulated in USDA1021"))
enrich_go_summation <- rbind(enrich_go_summation, go_nods_r_down %>% mutate(analysis = "Similar response to rhizobial differences in nodules from infected and uninfected hosts -- upregulated in WSM1022"))
enrich_go_summation <- rbind(enrich_go_summation, go_nods_n      %>% mutate(analysis = "Similar response to nematode infection in nodules with USDA1021 and nodules with WSM1022"))
enrich_go_summation <- rbind(enrich_go_summation, go_nods_n_up   %>% mutate(analysis = "Similar response to nematode infection in nodules with USDA1021 and nodules with WSM1022 -- upregaulted in USDA1021"))
enrich_go_summation <- rbind(enrich_go_summation, go_nods_n_down %>% mutate(analysis = "Similar response to nematode infection in nodules with USDA1021 and nodules with WSM1022 -- upregulated in WSM1022"))
enrich_go_summation <- rbind(enrich_go_summation, go_nods_i      %>% mutate(analysis = "Interation effect of rhizobial strain and nematode infection status in nodules"))
enrich_go_summation <- rbind(enrich_go_summation, go_gano_r      %>% mutate(analysis = "Similar response to rhizobial differences in galls and nodules"))
enrich_go_summation <- rbind(enrich_go_summation, go_gano_r_up   %>% mutate(analysis = "Similar response to rhizobial differences in galls and nodules -- upregulated in USDA1021"))
enrich_go_summation <- rbind(enrich_go_summation, go_gano_r_down %>% mutate(analysis = "Similar response to rhizobial differences in galls and nodules -- upregulated in WSM1022"))
enrich_go_summation <- rbind(enrich_go_summation, go_gano_i      %>% mutate(analysis = "Interavtion effect of rhizobial strain and organ type in galls and nodules"))
enrich_go_summation <- rbind(enrich_go_summation, go_garo_r      %>% mutate(analysis = "Similar response to rhizobial differences in galls and roots"))
enrich_go_summation <- rbind(enrich_go_summation, go_garo_r_up   %>% mutate(analysis = "Similar response to rhizobial differences in galls and roots -- upregulated in USDA1021"))
enrich_go_summation <- rbind(enrich_go_summation, go_garo_r_down %>% mutate(analysis = "Similar response to rhizobial differences in galls and roots -- upregulated in WSM1022"))
enrich_go_summation <- rbind(enrich_go_summation, go_garo_i %>% mutate(analysis = "Interaction effect of rhizobial strain and organ type in galls and roots"))
GOEnr <- enrich_go_summation
write.csv(GOEnr, file = "Dataset_S06.csv", row.names = FALSE)


med_genome <- t %>% filter(organism == "Medicago") %>% pull(gene) %>% length()
t %>% filter(organism == "Medicago") %>% select("gene", ends_with("_norm")) %>%
  rowwise() %>%
  mutate(cv_gall_1 = 100*mean(s13_norm, s17_norm, s21_norm, na.rm = TRUE)/sd(c(s13_norm, s17_norm, s21_norm), na.rm = TRUE),
         cv_gall_2 = 100*mean(s14_norm, s18_norm, s22_norm, na.rm = TRUE)/sd(c(s14_norm, s18_norm, s22_norm), na.rm = TRUE),
         cv_nodp_1 = 100*mean(s25_norm, s29_norm, s33_norm, na.rm = TRUE)/sd(c(s25_norm, s29_norm, s33_norm), na.rm = TRUE),
         cv_nodp_2 = 100*mean(s26_norm, s30_norm, s34_norm, na.rm = TRUE)/sd(c(s26_norm, s30_norm, s34_norm), na.rm = TRUE),
         cv_nodm_1 = 100*mean(s37_norm, s41_norm, s45_norm, na.rm = TRUE)/sd(c(s37_norm, s41_norm, s45_norm), na.rm = TRUE),
         cv_nodm_2 = 100*mean(s38_norm, s42_norm, s46_norm, na.rm = TRUE)/sd(c(s38_norm, s42_norm, s46_norm), na.rm = TRUE),
         cv_root_1 = 100*mean(s1_norm, s5_norm, s9_norm, na.rm = TRUE)/sd(c(s1_norm, s5_norm, s9_norm)),
         cv_root_2 = 100*mean(s2_norm, s6_norm, s10_norm, na.rm = TRUE)/sd(c(s2_norm, s6_norm, s10_norm)),
         cv_all = 100*mean(s1_norm, s2_norm, s5_norm, s6_norm, s9_norm, s10_norm,
                           s13_norm, s14_norm, s17_norm, s18_norm, s21_norm, s22_norm,
                           s25_norm, s26_norm, s29_norm, s30_norm, s33_norm, s34_norm,
                           s37_norm, s38_norm, s41_norm, s42_norm, s45_norm, s46_norm, na.rm = TRUE)/
           sd(c(s1_norm, s2_norm, s5_norm, s6_norm, s9_norm, s10_norm,
              s13_norm, s14_norm, s17_norm, s18_norm, s21_norm, s22_norm,
              s25_norm, s26_norm, s29_norm, s30_norm, s33_norm, s34_norm,
              s37_norm, s38_norm, s41_norm, s42_norm, s45_norm, s46_norm))) %>%
  mutate(cv_greater_gall_1 = if_else(cv_all > cv_gall_1, 1, 0),
         cv_greater_gall_2 = if_else(cv_all > cv_gall_2, 1, 0),
         cv_greater_nodp_1 = if_else(cv_all > cv_nodp_1, 1, 0),
         cv_greater_nodp_2 = if_else(cv_all > cv_nodp_2, 1, 0),
         cv_greater_nodm_1 = if_else(cv_all > cv_nodm_1, 1, 0),
         cv_greater_nodm_2 = if_else(cv_all > cv_nodm_2, 1, 0),
         cv_greater_root_1 = if_else(cv_all > cv_root_1, 1, 0),
         cv_greater_root_2 = if_else(cv_all > cv_root_2, 1, 0)) %>% ungroup() %>%
  summarise(cv_sum_gall_1 = sum(cv_greater_gall_1, na.rm = TRUE)/med_genome,
            cv_sum_gall_2 = sum(cv_greater_gall_2, na.rm = TRUE)/med_genome,
            cv_sum_nodp_1 = sum(cv_greater_nodp_1, na.rm = TRUE)/med_genome,
            cv_sum_nodp_2 = sum(cv_greater_nodp_2, na.rm = TRUE)/med_genome,
            cv_sum_nodm_1 = sum(cv_greater_nodm_1, na.rm = TRUE)/med_genome,
            cv_sum_nodm_2 = sum(cv_greater_nodm_2, na.rm = TRUE)/med_genome,
            cv_sum_root_1 = sum(cv_greater_root_1, na.rm = TRUE)/med_genome,
            cv_sum_root_2 = sum(cv_greater_root_2, na.rm = TRUE)/med_genome)


nem_genome <- t %>% filter(organism == "Nematode") %>% pull(gene) %>% length()
t %>% filter(organism == "Nematode") %>% select("gene", ends_with("_norm")) %>%
  rowwise() %>%
  mutate(cv_gall_1 = 100*mean(s13_norm, s17_norm, s21_norm, na.rm = TRUE)/sd(c(s13_norm, s17_norm, s21_norm), na.rm = TRUE),
         cv_gall_2 = 100*mean(s14_norm, s18_norm, s22_norm, na.rm = TRUE)/sd(c(s14_norm, s18_norm, s22_norm), na.rm = TRUE),
         cv_all = 100*mean(s13_norm, s14_norm, s17_norm, s18_norm, s21_norm, s22_norm, na.rm = TRUE)/
           sd(c(s13_norm, s14_norm, s17_norm, s18_norm, s21_norm, s22_norm, na.rm = TRUE))) %>%
  mutate(cv_greater_gall_1 = if_else(cv_all > cv_gall_1, 1, 0),
         cv_greater_gall_2 = if_else(cv_all > cv_gall_2, 1, 0)) %>% ungroup() %>%
  summarise(cv_sum_gall_1 = sum(cv_greater_gall_1, na.rm = TRUE)/nem_genome,
            cv_sum_gall_2 = sum(cv_greater_gall_2, na.rm = TRUE)/nem_genome)
         


r21_genome <- t %>% filter(organism == "Rhizob21") %>% pull(gene) %>% length()
t %>% filter(organism == "Rhizob21") %>% select("gene", ends_with("_norm")) %>%
  rowwise() %>%
  mutate(cv_nodp_1 = 100*mean(s25_norm, s29_norm, s33_norm, na.rm = TRUE)/sd(c(s25_norm, s29_norm, s33_norm), na.rm = TRUE),
         cv_nodm_1 = 100*mean(s37_norm, s41_norm, s45_norm, na.rm = TRUE)/sd(c(s37_norm, s41_norm, s45_norm), na.rm = TRUE),
         cv_all = 100*mean(s25_norm, s29_norm, s33_norm, s37_norm, s41_norm, s45_norm, na.rm = TRUE)/
           sd(c(s25_norm, s29_norm, s33_norm, s37_norm, s41_norm, s45_norm, na.rm = TRUE))) %>%
  mutate(cv_greater_nodp_1 = if_else(cv_all > cv_nodp_1, 1, 0),
         cv_greater_nodm_1 = if_else(cv_all > cv_nodm_1, 1, 0)) %>% ungroup() %>%
  summarise(cv_sum_nodp_1 = sum(cv_greater_nodp_1, na.rm = TRUE)/r21_genome,
            cv_sum_nodm_1 = sum(cv_greater_nodm_1, na.rm = TRUE)/r21_genome)



r22_genome <- t %>% filter(organism == "Rhizob22") %>% pull(gene) %>% length()
t %>% filter(organism == "Rhizob22") %>% select("gene", ends_with("_norm")) %>%
  rowwise() %>%
  mutate(cv_nodp_2 = 100*mean(s26_norm, s30_norm, s34_norm, na.rm = TRUE)/sd(c(s26_norm, s30_norm, s34_norm), na.rm = TRUE),
         cv_nodm_2 = 100*mean(s38_norm, s42_norm, s46_norm, na.rm = TRUE)/sd(c(s38_norm, s42_norm, s46_norm), na.rm = TRUE),
         cv_all = 100*mean(s26_norm, s30_norm, s34_norm, s38_norm, s42_norm, s46_norm, na.rm = TRUE)/
           sd(c(s26_norm, s30_norm, s34_norm, s38_norm, s42_norm, s46_norm, na.rm = TRUE))) %>%
  mutate(cv_greater_nodp_2 = if_else(cv_all > cv_nodp_2, 1, 0),
         cv_greater_nodm_2 = if_else(cv_all > cv_nodm_2, 1, 0)) %>% ungroup() %>%
  summarise(cv_sum_nodp_2 = sum(cv_greater_nodp_2, na.rm = TRUE)/r22_genome,
            cv_sum_nodm_2 = sum(cv_greater_nodm_2, na.rm = TRUE)/r22_genome)
