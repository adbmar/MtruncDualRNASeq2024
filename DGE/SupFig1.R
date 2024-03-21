library(ggplot2)
library(gridExtra)
library(ggplotify)
library(ggpubr)
library(ggpp)


### Generating useful dataframe for upset plotting

t_upset <- t %>% filter(organism == "Medicago") %>%
  mutate(exp_stat_G = ifelse(exp_gall, ifelse(exp_gall21, ifelse(exp_gall22, 3, 1), 2), 0)) %>%
  mutate(exp_stat_N = ifelse(exp_nodp, ifelse(exp_nodp21, ifelse(exp_nodp22, 3, 1), 2), 0)) %>%
  mutate(exp_stat_n = ifelse(exp_nodm, ifelse(exp_nodm21, ifelse(exp_nodm22, 3, 1), 2), 0)) %>%
  mutate(exp_stat_R = ifelse(exp_root, ifelse(exp_root21, ifelse(exp_root22, 3, 1), 2), 0)) %>%
  mutate(exp_stat_1_n = ifelse(exp_nods21, ifelse(exp_nodp21, ifelse(exp_nodm21, 3, 1), 2), 0)) %>%
  mutate(exp_stat_2_n = ifelse(exp_nods22, ifelse(exp_nodp22, ifelse(exp_nodm22, 3, 1), 2), 0))

t_upset <- t_upset %>% rowwise() %>% 
  mutate(exp_stat = ifelse(exp_gall & !exp_nodm & !exp_nodp & !exp_root,
                           ifelse(exp_stat_G > 2, "Expressed in all sample types", ifelse(exp_stat_G > 1, "Expression in Em1022 Only", ifelse(exp_stat_G > 0, "Expression in Em1021 Only", "Error"))), 
                           ifelse(!exp_gall & exp_nodm & !exp_nodp & !exp_root,
                                  ifelse(exp_stat_n > 2, "Expressed in all sample types", ifelse(exp_stat_n > 1, "Expression in Em1022 Only", ifelse(exp_stat_n > 0, "Expression in Em1021 Only", "Error"))),
                                  ifelse(!exp_gall & !exp_nodm & exp_nodp & !exp_root,
                                         ifelse(exp_stat_N > 2, "Expressed in all sample types", ifelse(exp_stat_N > 1, "Expression in Em1022 Only", ifelse(exp_stat_N > 0, "Expression in Em1021 Only", "Error"))),
                                         ifelse(!exp_gall & !exp_nodm & !exp_nodp & exp_root, 
                                                ifelse(exp_stat_R > 2, "Expressed in all sample types", ifelse(exp_stat_R > 1, "Expression in Em1022 Only", ifelse(exp_stat_R > 0, "Expression in Em1021 Only", "Error"))),
                                                ifelse(exp_gall & exp_nodm & !exp_nodp & !exp_root, 
                                                       ifelse(exp_stat_G == exp_stat_n, ifelse(exp_stat_G == 3, "Expressed in all sample types", ifelse(exp_stat_G == 2, "Expression in Em1022 Only", ifelse(exp_stat_G == 1, "Expression in Em1021 Only", "Error"))), "Expression differs across rhizobia and tissue combinations"),
                                                       ifelse(exp_gall & !exp_nodm & exp_nodp & !exp_root, 
                                                              ifelse(exp_stat_G == exp_stat_N, ifelse(exp_stat_G == 3, "Expressed in all sample types", ifelse(exp_stat_G == 2, "Expression in Em1022 Only", ifelse(exp_stat_G == 1, "Expression in Em1021 Only", "Error"))), "Expression differs across rhizobia and tissue combinations"),
                                                              ifelse(exp_gall & !exp_nodm & !exp_nodp & exp_root,
                                                                     ifelse(exp_stat_G == exp_stat_R, ifelse(exp_stat_G == 3, "Expressed in all sample types", ifelse(exp_stat_G == 2, "Expression in Em1022 Only", ifelse(exp_stat_G == 1, "Expression in Em1021 Only", "Error"))), "Expression differs across rhizobia and tissue combinations"),
                                                                     ifelse(!exp_gall & exp_nodm & exp_nodp & !exp_root,
                                                                            ifelse(exp_stat_N == exp_stat_n, ifelse(exp_stat_N == 3, "Expressed in all sample types", ifelse(exp_stat_N == 2, "Expression in Em1022 Only", ifelse(exp_stat_N == 1, "Expression in Em1021 Only", "Error"))), "Expression differs across rhizobia and tissue combinations"),
                                                                            ifelse(!exp_gall & exp_nodm & !exp_nodp & exp_root,
                                                                                   ifelse(exp_stat_R == exp_stat_n, ifelse(exp_stat_R == 3, "Expressed in all sample types", ifelse(exp_stat_R == 2, "Expression in Em1022 Only", ifelse(exp_stat_R == 1, "Expression in Em1021 Only", "Error"))), "Expression differs across rhizobia and tissue combinations"),
                                                                                   ifelse(!exp_gall & !exp_nodm & exp_nodp & exp_root, 
                                                                                          ifelse(exp_stat_R == exp_stat_N, ifelse(exp_stat_R == 3, "Expressed in all sample types", ifelse(exp_stat_R == 2, "Expression in Em1022 Only", ifelse(exp_stat_R == 1, "Expression in Em1021 Only", "Error"))), "Expression differs across rhizobia and tissue combinations"),
                                                                                          ifelse(exp_gall & exp_nodm & exp_nodp & !exp_root,
                                                                                                 ifelse(exp_stat_G == exp_stat_N & exp_stat_G == exp_stat_n, ifelse(exp_stat_G == 3, "Expressed in all sample types", ifelse(exp_stat_G == 2, "Expression in Em1022 Only", ifelse(exp_stat_G == 1, "Expression in Em1021 Only", "Error"))), "Expression differs across rhizobia and tissue combinations"),
                                                                                                 ifelse(exp_gall & exp_nodm & !exp_nodp & exp_root, 
                                                                                                        ifelse(exp_stat_G == exp_stat_n & exp_stat_G== exp_stat_R, ifelse(exp_stat_G == 3, "Expressed in all sample types", ifelse(exp_stat_G == 2, "Expression in Em1022 Only", ifelse(exp_stat_G == 1, "Expression in Em1021 Only", "Error"))), "Expression differs across rhizobia and tissue combinations"),
                                                                                                        ifelse(exp_gall & !exp_nodm & exp_nodp & exp_root,
                                                                                                               ifelse(exp_stat_G == exp_stat_N & exp_stat_G == exp_stat_R, ifelse(exp_stat_G == 3, "Expressed in all sample types", ifelse(exp_stat_G == 2, "Expression in Em1022 Only", ifelse(exp_stat_G == 1, "Expression in Em1021 Only", "Error"))), "Expression differs across rhizobia and tissue combinations"),
                                                                                                               ifelse(!exp_gall & exp_nodm & exp_nodp & exp_root, 
                                                                                                                      ifelse(exp_stat_N == exp_stat_n & exp_stat_N == exp_stat_R, ifelse(exp_stat_N == 3, "Expressed in all sample types", ifelse(exp_stat_N == 2, "Expression in Em1022 Only", ifelse(exp_stat_N == 1, "Expression in Em1021 Only", "Error"))), "Expression differs across rhizobia and tissue combinations"),
                                                                                                                      ifelse(exp_gall & exp_nodm & exp_nodp & exp_root, 
                                                                                                                             ifelse(exp_stat_G == exp_stat_n & exp_stat_G == exp_stat_N & exp_stat_G == exp_stat_R, ifelse(exp_stat_G == 3, "Expressed in all sample types", ifelse(exp_stat_G == 2, "Expression in Em1022 Only", ifelse(exp_stat_G == 1, "Expression in Em1021 Only", "Error"))), "Expression differs across rhizobia and tissue combinations"),
                                                                                                                             "No Expression")))))))))))))))) %>%
  
  mutate(exp_statn = case_when(exp_nodp21 & exp_nodp22 & exp_nodm21 & exp_nodm22 ~ "Expressed in all sample types",
                               !exp_nodp21 & !exp_nodp22 & !exp_nodm21 & !exp_nodm22 ~ "Not expressed",
                               exp_nodp21 & !exp_nodp22 & !exp_nodm21 & !exp_nodm22 ~ "Expressed in nodules with parasites only",
                               !exp_nodp21 & exp_nodp22 & !exp_nodm21 & !exp_nodm22 ~ "Expressed in nodules with parasites only",
                               !exp_nodp21 & !exp_nodp22 & exp_nodm21 & !exp_nodm22 ~ "Expressed in nodules without parasites only",
                               !exp_nodp21 & !exp_nodp22 & !exp_nodm21 & exp_nodm22 ~ "Expressed in nodules without parasites only",
                               exp_nodp21 & exp_nodp22 & !exp_nodm21 & !exp_nodm22 ~ "Expressed in nodules with parasites only",
                               exp_nodp21 & !exp_nodp22 & exp_nodm21 & !exp_nodm22 ~ "Expressed in all sample types",
                               exp_nodp21 & !exp_nodp22 & !exp_nodm21 & exp_nodm22 ~ "Opposite expression pattern across parasite status",
                               !exp_nodp21 & exp_nodp22 & exp_nodm21 & !exp_nodm22 ~ "Opposite expression pattern across parasite status",
                               !exp_nodp21 & exp_nodp22 & !exp_nodm21 & exp_nodm22 ~ "Expressed in all sample types",
                               !exp_nodp21 & !exp_nodp22 & exp_nodm21 & exp_nodm22 ~ "Expressed in nodules without parasites only",
                               exp_nodp21 & exp_nodp22 & exp_nodm21 & !exp_nodm22 ~ "Expression differs across rhizobia and parasite combinations",
                               exp_nodp21 & exp_nodp22 & !exp_nodm21 & exp_nodm22 ~ "Expression differs across rhizobia and parasite combinations",
                               exp_nodp21 & !exp_nodp22 & exp_nodm21 & exp_nodm22 ~ "Expression differs across rhizobia and parasite combinations",
                               !exp_nodp21 & exp_nodp22 & exp_nodm21 & exp_nodm22 ~ "Expression differs across rhizobia and parasite combinations",
                               .default = "error"))
    
  #   exp_stat_1_n == 3 & exp_stat_2_n == 3, "Expressed in all sample types",
  #                            ifelse((exp_stat_1_n == 1 & exp_stat_2_n == 2) | (exp_stat_1_n == 2 & exp_stat_2_n == 1), "Opposite expression pattern with parasites",
  #                                   ifelse(exp_stat_1_n == 1 & exp_stat_2_n == 1, "Expression in hosts without parasites only",
  #                                          ifelse(exp_stat_1_n == 2 & exp_stat_2_n == 2, "Expression in hosts with parasites only", 
  #                                                 ifelse((exp_stat_1_n == 0 & exp_stat_2_n > 0) | (exp_stat_1_n > 0 & exp_stat_2_n == 0), "Expression differs across rhizobia and parasite combinations", "other")
  #                                          )
  #                                   )
  #                            )
  # )
  # )


t_upset_long <- t_upset %>% select("gene", starts_with("exp_stat_")) %>% pivot_longer(cols = starts_with("exp_stat_"))
t_upset_long$value <- as.character(t_upset_long$value)

t_upset_long_n <-
  rbind(
    t_upset %>% select("gene", "exp_nods21", "exp_nods22", "exp_statn") %>% pivot_longer(cols = c("exp_statn")) %>% filter(exp_nods21) %>% mutate(strain = ifelse(exp_nods21, "Em1021", "error")) %>% select("gene", "strain", "name", "value"),
    t_upset %>% select("gene", "exp_nods21", "exp_nods22", "exp_statn") %>% pivot_longer(cols = c("exp_statn")) %>% filter(exp_nods22) %>% mutate(strain = ifelse(exp_nods22, "Em1022", "error")) %>% select("gene", "strain", "name", "value"))

t_upset_long_n_2 <- rbind(
    t %>% filter(organism == "Medicago") %>% filter(exp_nods21) %>% mutate(exp_stat = ifelse(exp_nodm21 & exp_nodp21, "Expressed in all sample types", ifelse(exp_nodm21 & !exp_nodp21, "Expressed in nodules without parasites", ifelse(!exp_nodm21 & exp_nodp21, "Expressed in nodules with parasites", "error")))) %>% mutate(strain = "Em1021") %>% select("gene", "strain", "exp_stat"),
    t %>% filter(organism == "Medicago") %>% filter(exp_nods22) %>% mutate(exp_stat = ifelse(exp_nodm22 & exp_nodp22, "Expressed in all sample types", ifelse(exp_nodm22 & !exp_nodp22, "Expressed in nodules without parasites", ifelse(!exp_nodm22 & exp_nodp22, "Expressed in nodules with parasites", "error")))) %>% mutate(strain = "Em1022") %>% select("gene", "strain", "exp_stat")) 

c(
  "exp_gall"="Expressed in galls",
  "exp_root"="Expressed in roots",
  "exp_nodm" = "Expressed in nodules from\nhosts without parasites",
  "exp_nodp" = "Expressed in nodules from\nhosts with parasites",
  "exp_nodm"="Expressed in nodules\nwithout parasites",
  "exp_nods21"="Expressed in nodules\nwith Em1021",
  "exp_nods22" = "Expressed in nodules\nwith Em1022",
  "exp_nodp" = "Expressed in nodules\nwith parasites",
  "med_gano_r.sig" = "Similar differential expression\nacross rhizobia strain\nin galls and nodules",
  "med_garo_r.sig" = "Similar differential expression\nacross rhizobia strain\nin galls and roots",
  "med_noro_r.sig" = "Similar differential expression\nacross rhizobia strain\nin nodules and roots",
  "med_gall_r.sig" = "Differential expression\nin galls",
  "med_nodp_r.sig" = "Differential expression\nin nodules",
  "med_root_r.sig" = "Differential expression\nin roots"
  )


med_total_genes <- t %>% filter(organism == "Medicago") %>% pull(gene) %>% length()
nem_total_genes <- t %>% filter(organism == "Nematode") %>% pull(gene) %>% length()
r21_total_genes <- t %>% filter(organism == "Rhizob21") %>% pull(gene) %>% length()
r22_total_genes <- t %>% filter(organism == "Rhizob22") %>% pull(gene) %>% length()

















### Supplemental figure 1: Expression in different tissue types as it relates to rhizobial identity

sup_p1_1 <- ComplexUpset::upset(t_upset, c("exp_gall", "exp_nodp", "exp_nodm", "exp_root"),
                              sort_intersections_by=c('degree', 'cardinality'),
                              name = "Expression intersections",
                              height_ratio = 1:3,
                              set_sizes = FALSE,
                              base_annotation = list("Gene count" =
                                                       ComplexUpset::intersection_size() +
                                                       geom_hline(yintercept = med_total_genes, color = "dark gray") +
                                                       annotate("text", y = med_total_genes, x = 0.4, hjust = 1, vjust = 0, color = "dark gray", cex = 3,
                                                                label = paste0("Total genes in genome (", format(med_total_genes, big.mark = ","), " genes)")) +
                                                       coord_cartesian(clip = "off")),
                              annotations = list(
                                "Expression status" = (
                                  ggplot(mapping = aes(fill = exp_stat)) +
                                    geom_bar(stat="count", position = "fill") +
                                    scale_y_continuous(labels = scales::percent_format()) +
                                    ylab("Proportion of intersection\nshowing different\nexpression patterns") +
                                    scale_fill_manual(name = "Expression pattern",
                                                      values = c("Expressed in all sample types" = "blue",
                                                                 "Expression differs across rhizobia and tissue combinations" = "orange",
                                                                 "Expression in Em1021 Only" = "cyan",
                                                                 "Expression in Em1022 Only" = "dark blue",
                                                                 "No Expression" = "gray"),
                                                      guide = FALSE))),
                              labeller =ggplot2::as_labeller(c(
                                "exp_gall"="Expressed in galls",
                                "exp_root"="Expressed in roots",
                                "exp_nodm" = "Expressed in nodules from\nhosts without parasites",
                                "exp_nodp" = "Expressed in nodules from\nhosts with parasites")))


sup_p1_legend <- get_legend(ggplot(t_upset %>% select(c(gene, exp_stat)), aes(x = exp_stat, fill = exp_stat)) +
  geom_bar(stat="count") +
  scale_fill_manual(name = "Expression pattern",
                    values = c("Expressed in all sample types" = "blue",
                               "Expression differs across rhizobia and tissue combinations" = "orange",
                               "Expression in Em1021 Only" = "cyan",
                               "Expression in Em1022 Only" = "dark blue",
                               "No Expression" = "gray")))


sup_p1_2 <- ggplot(t_upset_long %>% filter(value != "0" & name != "exp_stat_1_n" & name !=  "exp_stat_2_n")) +
  aes(x = name, fill = value) + 
  geom_bar(position = "stack") +
  ylab("Number of genes") +
  xlab("Sample type") +
  scale_x_discrete(labels = c("exp_stat_G" = "Galls",
                              "exp_stat_n" = "Nodules from hosts\nwithout parasites",
                              "exp_stat_N" = "Nodules from hosts\nwith parasites",
                              "exp_stat_R" = "Roots from hosts\nwith parasites")) +
  scale_fill_manual(name = "Expression status across rhizobia strain samples",
                    values = c("0" = "gray",
                               "1" = "cyan",
                               "2" = "dark blue", 
                               "3" = "blue"),
                    labels = c("0" = "Not expressed",
                               "1" = "Expressed when hosts are with Em1021 only", 
                               "2" = "Expressed when hosts are with Em1022 only", 
                               "3" = "Expressed when hosts are with either strain")) +
  geom_text(aes(label = format((..count..), format = "f", big.mark = ",", digits = 1)),
            stat = "count", color = "white", angle = 270,
            position = position_stack(vjust = 0.5)) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank()) +
  coord_flip()




sup_p1 <- as.ggplot(grid.arrange(grobs = list(as.ggplot(sup_p1_1), as.ggplot(sup_p1_2), as.ggplot(sup_p1_legend)),
                                  row = 2, ncol = 2, heights = c(200,105), widths = c(5,2),
                                  layout_matrix = rbind(c(1,3), c(1,2)))) +
  ggtitle("Medicgo gene expression across different organs") + 
  theme(panel.background = element_rect(fill = 'white', color = 'white'))




### Supplemental figure 2: Expression in different nodule samples


sup_p2_1 <- ComplexUpset::upset(t_upset, c("exp_nods21", "exp_nods22"),
                                sort_intersections_by=c('degree', 'cardinality'),
                                name = "Expression intersections",
                                height_ratio = 1:3,
                                set_sizes = FALSE,
                                labeller = ggplot2::as_labeller(c("exp_nods21" = "Expressed in nodules with Em1021", "exp_nods22" = "Expressed in nodules with Em1022")), 
                                base_annotation = list("Gene count" =
                                                         ComplexUpset::intersection_size() +
                                                         geom_hline(yintercept = med_total_genes, color = "dark gray") +
                                                         annotate("text", y = med_total_genes, x = 0.4, hjust = 1, vjust = 0, color = "dark gray", cex = 3,
                                                                  label = paste0("Total genes in genome (", format(med_total_genes, big.mark = ","), " genes)")) +
                                                         coord_cartesian(clip = "off")),
                                annotations = list(
                                  "Expression status" = (
                                    ggplot(mapping = aes(fill = exp_statn)) +
                                      geom_bar(stat="count", position = "fill") +
                                      scale_fill_manual(name = "Expression pattern",
                                                        values = c("Not expressed" = "light gray",
                                                                   "Expressed in all sample types" = "black",
                                                                   "Expressed in nodules with parasites only" = "red",
                                                                   "Expressed in nodules without parasites only" = "blue",
                                                                   "Expression differs across rhizobia and parasite combinations" = "purple",
                                                                   "Opposite expression pattern across parasite status" = "orange"),
                                                        guide = "none") +
                                      scale_y_continuous(labels = scales::percent_format()))))


sup_p2_legend <- get_legend(ggplot(t_upset %>% select(c(gene, exp_statn)), aes(x = exp_statn, fill = exp_statn)) +
                              geom_bar(stat="count") +
                              scale_fill_manual(name = "Expression pattern",
                                                values = c("Not expressed" = "light gray",
                                                           "Expressed in all sample types" = "black",
                                                           "Expressed in nodules with parasites only" = "red",
                                                           "Expressed in nodules without parasites only" = "blue",
                                                           "Expression differs across rhizobia and parasite combinations" = "purple",
                                                           "Opposite expression pattern across parasite status" = "orange")))

sup_p2_2 <- ggplot(t_upset_long_n_2 %>% filter(exp_stat != "Not Expressed")) +
  aes(x = strain, fill = exp_stat) + geom_bar(position = "stack") +
  ylab("Number of genes") +
  xlab("Sample type") +
  geom_text(aes(label = format((..count..), format = "f", big.mark = ",", digits = 1)),
            stat = "count", color = "white", angle = 270,
            position = position_stack(vjust = 0.5)) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank()) +
  coord_flip() +
  scale_fill_manual(name = "Expression pattern",
                    values = c("Expressed in all sample types" = "black",
                               "Expressed in nodules with parasites" = "red",
                               "Expressed in nodules without parasites" = "blue"),
                    guide = "none") 




sup_p2 <- as.ggplot(grid.arrange(grobs = list(as.ggplot(sup_p2_1), as.ggplot(sup_p2_2), as.ggplot(sup_p2_legend)),
                                 row = 2, ncol = 2, heights = c(200,105), widths = c(5,2),
                                 layout_matrix = rbind(c(1,3), c(1,2)))) +
  ggtitle("Medicgo gene expression in nodules across different rhizobia strains") +
  theme(panel.background = element_rect(fill = 'white', color = 'white'))
  


### Supplemental figure 3: Differential expression in galls and nodules

sup_p3_a <- ComplexUpset::upset(t %>% filter(organism == "Medicago") %>% filter(med_gall_r.sig | med_nodp_r.sig | med_gano_r.sig),
                    c("med_gall_r.sig", "med_nodp_r.sig", "med_gano_r.sig"),
                    name = "Expression intersections",
                    set_sizes = ComplexUpset::upset_set_size(position = "right") + ylab("Number of genes"),
                    labeller = ggplot2::as_labeller(c("med_gano_r.sig" = "Similar differential expression\naross rhizobia strain\nin galls and nodules",
                                                      "med_gall_r.sig" = "Differential expression\nacross rhizobia\nin galls",
                                                      "med_nodp_r.sig" = "Differential expression\nacross rhizobia\nin nodules")),
                    annotation = list(
                      "Log fold change across rhizobia status in galls" = (
                        ggplot() +
                          geom_hline(yintercept = 0) +
                          geom_linerange(alpha = 0.2,
                                         aes(y = med_gall_r.ns_LFC,
                                             ymin = med_gall_r.ns_LFC - med_gall_r.ns_SE,
                                             ymax = med_gall_r.ns_LFC + med_gall_r.ns_SE),
                                         position = position_jitternudge(seed = 1,
                                                                         width = 0.4,
                                                                         x = 0,
                                                                         nudge.from = "original.y")) +
                          geom_point(aes(y = med_gall_r.ns_LFC),
                                     shape = 21,  alpha = 0.7, fill = "red",
                                     position = position_jitternudge(seed = 1,
                                                                     width = 0.4,
                                                                     x = 0,
                                                                     nudge.from = "original.y")) +
                          ylab("Log fold change\nin galls")),
                      "Log fold change across rhizobia status in nodules" = (
                        ggplot() +
                          geom_hline(yintercept = 0) +
                          geom_linerange(alpha = 0.2,
                                         aes(y = med_nodp_r.ns_LFC,
                                             ymin = med_nodp_r.ns_LFC - med_nodp_r.ns_SE,
                                             ymax = med_nodp_r.ns_LFC + med_nodp_r.ns_SE),
                                         position = position_jitternudge(seed = 1,
                                                                         width = 0.4,
                                                                         x = 0,
                                                                         nudge.from = "original.y")) +
                          geom_point(aes(y = med_nodp_r.ns_LFC),
                                     shape = 21,  alpha = 0.7, fill = "blue",
                                     position = position_jitternudge(seed = 1,
                                                                     width = 0.4,
                                                                     x = 0,
                                                                     nudge.from = "original.y")) +
                          ylab("Log fold change\nin nodules")),
                      "Comparison of log fold changes in nodules and galls" = (
                        ggplot() + geom_hline(yintercept = 0) +
                          geom_violin(aes(y = med_nodp_r.ns_LFC), 
                                      fill = "blue", width = 0.5,
                                      position = position_nudge(x = -0.15),
                                      alpha = 0.5,
                                      draw_quantiles = c(0.25, 0.5, 0.75),) +
                          geom_violin(aes(y = med_gall_r.ns_LFC), 
                                      fill = "red", width = 0.5, 
                                      position = position_nudge(x = 0.15),
                                      alpha = 0.5,
                                      draw_quantiles = c(0.25, 0.5, 0.75),) +
                          ylab("Log fold change distribution\nin nodules and galls"))
                      )) %>% as.ggplot() +
  ggtitle("Differential expression in galls and nodules across rhizobial status")



### Supplemental figure 4: Differential expression in galls and roots

sup_p3_b <- ComplexUpset::upset(t %>% filter(organism == "Medicago") %>% filter(med_gall_r.sig | med_root_r.sig | med_garo_r.sig),
                    c("med_gall_r.sig", "med_root_r.sig", "med_garo_r.sig"),
                    name = "Expression intersections",
                    set_sizes = ComplexUpset::upset_set_size(position = "right") + ylab("Number of genes"),
                    labeller = ggplot2::as_labeller(c("med_garo_r.sig" = "Similar differential expression\naross rhizobia strain\nin galls and roots",
                                                      "med_gall_r.sig" = "Differential expression\nacross rhizobia\nin galls",
                                                      "med_root_r.sig" = "Differential expression\nacross rhizobia\nin roots")),
                    annotation = list(
                      "Log fold change across rhizobia status in galls" = (
                        ggplot() +
                          geom_hline(yintercept = 0) +
                          geom_linerange(alpha = 0.2,
                                         aes(y = med_gall_r.ns_LFC,
                                             ymin = med_gall_r.ns_LFC - med_gall_r.ns_SE,
                                             ymax = med_gall_r.ns_LFC + med_gall_r.ns_SE),
                                         position = position_jitternudge(seed = 1,
                                                                         width = 0.4,
                                                                         x = 0,
                                                                         nudge.from = "original.y")) +
                          geom_point(aes(y = med_gall_r.ns_LFC),
                                     shape = 21,  alpha = 0.7, fill = "red",
                                     position = position_jitternudge(seed = 1,
                                                                     width = 0.4,
                                                                     x = 0,
                                                                     nudge.from = "original.y")) +
                          ylab("Log fold change\nin galls")),
                      "Log fold change across rhizobia status in roots" = (
                        ggplot() +
                          geom_hline(yintercept = 0) +
                          geom_linerange(alpha = 0.2,
                                         aes(y = med_root_r.ns_LFC,
                                             ymin = med_root_r.ns_LFC - med_root_r.ns_SE,
                                             ymax = med_root_r.ns_LFC + med_root_r.ns_SE),
                                         position = position_jitternudge(seed = 1,
                                                                         width = 0.4,
                                                                         x = 0,
                                                                         nudge.from = "original.y")) +
                          geom_point(aes(y = med_root_r.ns_LFC),
                                     shape = 21,  alpha = 0.7, fill = "tan",
                                     position = position_jitternudge(seed = 1,
                                                                     width = 0.4,
                                                                     x = 0,
                                                                    nudge.from = "original.y")) +
                          ylab("Log fold change\nin roots")),
                      "Comparison of log fold changes in galls and roots" = (
                        ggplot() + geom_hline(yintercept = 0) +
                          geom_violin(aes(y = med_root_r.ns_LFC), 
                                      fill = "tan", width = 0.5,
                                      position = position_nudge(x = -0.15),
                                      alpha = 0.5,
                                      draw_quantiles = c(0.25, 0.5, 0.75),) +
                          geom_violin(aes(y = med_gall_r.ns_LFC), 
                                      fill = "red", width = 0.5, 
                                      position = position_nudge(x = 0.15),
                                      alpha = 0.5,
                                      draw_quantiles = c(0.25, 0.5, 0.75),) +
                          ylab("Log fold change distribution\nin nodules and galls"))
                    )) %>% as.ggplot() +
  ggtitle("Differential expression in galls and roots across rhizobial status")



### Supplemental figure 5: Differential expression in nodules and roots

sup_p3_c <- ComplexUpset::upset(t %>% filter(organism == "Medicago") %>% filter(med_root_r.sig | med_nodp_r.sig | med_noro_r.sig),
                    c("med_root_r.sig", "med_nodp_r.sig", "med_noro_r.sig"),
                    name = "Expression intersections",
                    set_sizes = ComplexUpset::upset_set_size(position = "right") + ylab("Number of genes"),
                    
                    labeller = ggplot2::as_labeller(c("med_noro_r.sig" = "Similar differential expression\naross rhizobia strain\nin nodules and roots",
                                                      "med_root_r.sig" = "Differential expression\nacross rhizobia\nin roots",
                                                      "med_nodp_r.sig" = "Differential expression\nacross rhizobia\nin nodules")),
                    annotation = list(
                      "Log fold change across rhizobia status in roots" = (
                        ggplot() +
                          geom_hline(yintercept = 0) +
                          geom_linerange(alpha = 0.2,
                                         aes(y = med_root_r.ns_LFC,
                                             ymin = med_root_r.ns_LFC - med_root_r.ns_SE,
                                             ymax = med_root_r.ns_LFC + med_root_r.ns_SE),
                                         position = position_jitternudge(seed = 1,
                                                                         width = 0.4,
                                                                         x = 0,
                                                                         nudge.from = "original.y")) +
                          geom_point(aes(y = med_root_r.ns_LFC),
                                     shape = 21,  alpha = 0.7, fill = "tan",
                                     position = position_jitternudge(seed = 1,
                                                                     width = 0.4,
                                                                     x = 0,
                                                                     nudge.from = "original.y")) +
                          ylab("Log fold change in roots")),
                      "Log fold change across rhizobia status in nodules" = (
                        ggplot() +
                          geom_hline(yintercept = 0) +
                          geom_linerange(alpha = 0.2,
                                         aes(y = med_nodp_r.ns_LFC,
                                             ymin = med_nodp_r.ns_LFC - med_nodp_r.ns_SE,
                                             ymax = med_nodp_r.ns_LFC + med_nodp_r.ns_SE),
                                         position = position_jitternudge(seed = 1,
                                                                         width = 0.4,
                                                                         x = 0,
                                                                         nudge.from = "original.y")) +
                          geom_point(aes(y = med_nodp_r.ns_LFC),
                                     shape = 21,  alpha = 0.7, fill = "blue",
                                     position = position_jitternudge(seed = 1,
                                                                     width = 0.4,
                                                                     x = 0,
                                                                     nudge.from = "original.y")) +
                          ylab("Log fold change in nodules")),
                      "Comparison of log fold changes in nodules and roots" = (
                        ggplot() + geom_hline(yintercept = 0) +
                          geom_violin(aes(y = med_nodp_r.ns_LFC), 
                                      fill = "blue", width = 0.5,
                                      position = position_nudge(x = -0.15),
                                      alpha = 0.5,
                                      draw_quantiles = c(0.25, 0.5, 0.75),) +
                          geom_violin(aes(y = med_root_r.ns_LFC), 
                                      fill = "tan", width = 0.5, 
                                      position = position_nudge(x = 0.15),
                                      alpha = 0.5,
                                      draw_quantiles = c(0.25, 0.5, 0.75),) +
                          ylab("Log fold change distribution\nin nodules and roots"))
                    )) %>% as.ggplot() +
  ggtitle("Differential expression in roots and nodules across rhizobial status")


sup_p3 <- as.ggplot(grid.arrange(grobs = list(as.ggplot(sup_p3_a), as.ggplot(sup_p3_b), as.ggplot(sup_p3_c)),
                                 row = 1, ncol = 3)) +
  theme(panel.background = element_rect(fill = 'white', color = 'white'))

### Supplemental figure 6: Showing resolution to seeming statistical paradox
sup_p4 <- ggplot(
  rbind(
    t %>% filter(med_gall_r.sig) %>% select(med_gall_r.ns_LFC) %>% pivot_longer(cols = med_gall_r.ns_LFC) %>% mutate(model = "DEGs from model with just gall samples") %>% mutate(LFC = "Log fold change\nfrom model with\njust gall samples") %>% mutate(Sample_size = "One organ type"),
    t %>% filter(med_gano_r.sig) %>% select(med_gall_r.ns_LFC) %>% pivot_longer(cols = med_gall_r.ns_LFC) %>% mutate(model = "DEGs from model with gall and nodule samples") %>% mutate(LFC = "Log fold change\nfrom model with\njust gall samples") %>% mutate(Sample_size = "Two organ types"),
    t %>% filter(med_garo_r.sig) %>% select(med_gall_r.ns_LFC) %>% pivot_longer(cols = med_gall_r.ns_LFC) %>% mutate(model = "DEGs from model with gall and root samples") %>% mutate(LFC = "Log fold change\nfrom model with\njust gall samples") %>% mutate(Sample_size = "Two organ types"),
    
    t %>% filter(med_nodp_r.sig) %>% select(med_nodp_r.ns_LFC) %>% pivot_longer(cols = med_nodp_r.ns_LFC) %>% mutate(model = "DEGs from model with just nodule samples") %>% mutate(LFC = "Log fold change\nfrom model with\njust nodule samples") %>% mutate(Sample_size = "One organ type"),
    t %>% filter(med_gano_r.sig) %>% select(med_nodp_r.ns_LFC) %>% pivot_longer(cols = med_nodp_r.ns_LFC) %>% mutate(model = "DEGs from model with gall and nodule samples") %>% mutate(LFC = "Log fold change\nfrom model with\njust nodule samples") %>% mutate(Sample_size = "Two organ types"),
    t %>% filter(med_noro_r.sig) %>% select(med_nodp_r.ns_LFC) %>% pivot_longer(cols = med_nodp_r.ns_LFC) %>% mutate(model = "DEGs from model with nodule and root samples") %>% mutate(LFC = "Log fold change\nfrom model with\njust nodule samples") %>% mutate(Sample_size = "Two organ types"),
    
    t %>% filter(med_root_r.sig) %>% select(med_root_r.ns_LFC) %>% pivot_longer(cols = med_root_r.ns_LFC) %>% mutate(model = "DEGs from model with just root samples") %>% mutate(LFC = "Log fold change\nfrom model with\njust root samples") %>% mutate(Sample_size = "One organ type"),
    t %>% filter(med_garo_r.sig) %>% select(med_root_r.ns_LFC) %>% pivot_longer(cols = med_root_r.ns_LFC) %>% mutate(model = "DEGs from model with gall and root samples") %>% mutate(LFC = "Log fold change\nfrom model with\njust root samples") %>% mutate(Sample_size = "Two organ types"),
    t %>% filter(med_noro_r.sig) %>% select(med_root_r.ns_LFC) %>% pivot_longer(cols = med_root_r.ns_LFC) %>% mutate(model = "DEGs from model with nodule and root samples") %>% mutate(LFC = "Log fold change\nfrom model with\njust root samples") %>% mutate(Sample_size = "Two organ types")
  )) + 
  aes(x = abs(value), y =LFC , fill = model, color = Sample_size) + 
  geom_boxplot(outlier.shape = 21, linewidth = 1.2) +
  scale_fill_manual(name = "Organ types provided to\nmodel to identify differentially\nexpressed genes",
    values = c("DEGs from model with just gall samples" = "red",
               "DEGs from model with just nodule samples" = "blue",
               "DEGs from model with just root samples" = "yellow",
               "DEGs from model with gall and nodule samples" = "purple",
               "DEGs from model with gall and root samples" = "orange",
               "DEGs from model with nodule and root samples" = "green"),
    labels = c(
      "DEGs from model with just gall samples" = "Just galls",
      "DEGs from model with just nodule samples" = "Just nodules with parasites",
      "DEGs from model with just root samples" = "Just roots",
      "DEGs from model with gall and nodule samples" = "Galls and nodules with parasites",
      "DEGs from model with gall and root samples" = "Galls and roots",
      "DEGs from model with nodule and root samples" = "Nodules with parasites and roots")) +
  scale_color_manual(name = "Number of organ types", 
    values = c("Two organ types" = "black", "One organ type" = "honeydew3")) + 
  theme_minimal() +
  xlab("Magnitude of log fold change values") +
  geom_hline(yintercept = 1.5) + geom_hline(yintercept = 2.5) + geom_vline(xintercept = 0) +
  theme(panel.grid.major.y = element_blank(), axis.title.y = element_blank(), plot.background = element_rect(fill = "white"))

setwd(dir_out)

ggsave(sup_p1, filename = "Sup_fig1.png", width = 4320, height = 2430, units = "px", bg = "white")
ggsave(sup_p2, filename = "Sup_fig2.png", width = 4320, height = 2430, units = "px", bg = "white")
ggsave(sup_p3, filename = "Sup_fig3.png", width = 8640, height = 4320, units = "px", bg = "white")
ggsave(sup_p4, filename = "Sup_fig4.png", width = 2430, height = 2430, units = "px", bg = "white")



### Supplemental figure 5: Nematode expression in different tissue types as it relates to rhizobial identity
sup_p5 <- ComplexUpset::upset(t %>% 
                                filter(organism == "Nematode"),
                              c("exp_gall21", "exp_gall22"),
                              sort_intersections_by=c('degree', 'cardinality'),
                              base_annotation = list("Gene count" =
                                                       ComplexUpset::intersection_size() +
                                                       geom_hline(yintercept = nem_total_genes, color = "dark gray") +
                                                       annotate("text", y = nem_total_genes, x = 0.4, hjust = 1, vjust = 0.5, color = "dark gray", cex = 3,
                                                                label = paste0("Total genes in genome (", format(nem_total_genes, big.mark = ","), " genes)")) +
                                                       coord_cartesian(clip = "off")),
                              name = "Expression intersections",
                              height_ratio = 1:3,
                              set_sizes = FALSE,
                              labeller = ggplot2::as_labeller(
                                c("exp_gall21" = "Expressed in galls from hosts with Em1021",
                                  "exp_gall22" = "Expressed in galls from hosts with Em1022")),
                              wrap = TRUE) +
  ggtitle("Nematode gene expression across hosts with different rhizobia strains")

### Supplemental figure 6: Em1021 expression in different tissue types as it relates to rhizobial identity
sup_p6 <- ComplexUpset::upset(t %>% 
                                filter(organism == "Rhizob21"),
                              c("exp_nodp21", "exp_nodm21"),
                              sort_intersections_by=c('degree', 'cardinality'),
                              base_annotation = list("Gene count" =
                                                       ComplexUpset::intersection_size() +
                                                       geom_hline(yintercept = r21_total_genes, color = "dark gray") +
                                                       annotate("text", y = r21_total_genes, x = 0.4, hjust = 1, vjust = 0.5, color = "dark gray", cex = 3,
                                                                label = paste0("Total genes in genome (", format(r21_total_genes, big.mark = ","), " genes)")) +
                                                       coord_cartesian(clip = "off")),
                              name = "Expression intersections",
                              height_ratio = 1:3,
                              set_sizes = FALSE,
                              labeller = ggplot2::as_labeller(
                                c("exp_nodm21" = "Expressed in nodules from hosts infected with parasites",
                                  "exp_nodp21" = "Expressed in nodules from uninfected hosts")),
                              wrap = TRUE) +
  ggtitle("Rhizobe strain Em1021 gene expression across hosts with different parasite infection status")

### Supplemental figure 7: Em1021 expression in different tissue types as it relates to rhizobial identity
sup_p7 <- ComplexUpset::upset(t %>% 
                                filter(organism == "Rhizob22"),
                              c("exp_nodp22", "exp_nodm22"),
                              sort_intersections_by=c('degree', 'cardinality'),
                              base_annotation = list("Gene count" =
                                                       ComplexUpset::intersection_size() +
                                                       geom_hline(yintercept = r22_total_genes, color = "dark gray") +
                                                       annotate("text", y = r22_total_genes, x = 0.4, hjust = 1, vjust = 0.5, color = "dark gray", cex = 3,
                                                                label = paste0("Total genes in genome (", format(r22_total_genes, big.mark = ","), " genes)")) +
                                                       coord_cartesian(clip = "off")),
                              name = "Expression intersections",
                              height_ratio = 1:3,
                              set_sizes = FALSE,
                              labeller = ggplot2::as_labeller(
                                c("exp_nodm22" = "Expressed in nodules from hosts infected with parasites",
                                  "exp_nodp22" = "Expressed in nodules from uninfected hosts")),
                              wrap = TRUE) +
  ggtitle("Rhizobe strain Em1022 gene expression across hosts with different parasite infection status")


ggsave(sup_p5, filename = "Sup_fig5.png", width = 4320, height = 2430, units = "px", bg = "white")
ggsave(sup_p6, filename = "Sup_fig6.png", width = 4320, height = 2430, units = "px", bg = "white")
ggsave(sup_p7, filename = "Sup_fig7.png", width = 4320, height = 2430, units = "px", bg = "white")


sup_p8 <- ComplexUpset::upset(t %>%
                                filter(med_gano_r.sig | med_gano_i.sig | med_garo_r.sig | med_garo_i.sig | med_noro_r.sig | med_noro_i.sig | med_all_r_gall.sig) %>%
                                filter(organism == "Medicago"),
                              c("med_gano_r.sig", "med_gano_i.sig",
                                "med_garo_r.sig", "med_garo_i.sig",
                                "med_noro_r.sig", "med_noro_i.sig",
                                "med_all_r_gall.sig"),
                              sort_intersections_by=c('degree', 'cardinality'),
                              name = "Expression intersections",
                              labeller = ggplot2::as_labeller(
                                c("med_gano_r.sig" = "Similar response to rhizobia strains in galls and nodules",
                                  "med_gano_i.sig" = "Different response to rhizobia strains across galls and nodules",
                                  "med_garo_r.sig" = "Similar response to rhizobia strains in galls and roots",
                                  "med_garo_i.sig" = "Different response to rhizobia strains in galls and roots",
                                  "med_noro_r.sig" = "Similar response to rhizobia strains in nodules and roots",
                                  "med_noro_i.sig" = "Different response to rhizobia strains in nodules and roots",
                                  "med_all_r_gall.sig" = "Similar response across all tissues")),
                              wrap = TRUE,  
                              set_sizes = (
                                ComplexUpset::upset_set_size() + 
                                  geom_text(aes(label=..count..), hjust=1.1, stat='count') +
                                  expand_limits(y=1000))
)
sup_p8

ggsave(sup_p8, filename = "Sup_fig8.png", width = 3645, height = 2430, units = "px", bg = "white")
