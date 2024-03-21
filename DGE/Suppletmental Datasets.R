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
med_DE <- t %>% filter(organism == "Medicago") %>% select("gene", starts_with("med_"), "NCR") %>% drop_na() %>%
  select(!ends_with(c(".LFC", ".SE", ".sig", ".significance", ".sigu", ".sigd", "sigu_ns", "sigd_ns")))
write.csv(med_DE, file = "Dataset_S03.csv", row.names = FALSE)

### Dataset_S04 Differential expression statistics of Ensifer meliloti genes in all sample types
rhi_DE <-t %>% filter(organism == "Rhizob21" | organism == "Rhizob22") %>% select("gene", starts_with("r21_")) %>% drop_na() %>%
  select(!ends_with(c(".LFC", ".SE", ".sig", ".significance", ".sigu", ".sigd", "sigu_ns", "sigd_ns")))
write.csv(rhi_DE, file = "Dataset_S04.csv", row.names = FALSE)

### Dataset_S05. Differential expression statistics of Meloidogyne hapla genes in all sample types
nem_DE <-t %>% filter(organism == "Nematode") %>% select("gene", starts_with("nem_")) %>% drop_na() %>%
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
