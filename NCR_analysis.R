library(eulerr)

NCR_list <- read.csv(file.path(dir_medic, "MtrunA17r5.0-ANR-EGN-r1.9.b2g.gaf"), sep="\t", skip = 4, header = FALSE) %>%
  filter(str_detect(V10, "(NCR)")) %>% pull(V2)

NCR_list <- append(NCR_list, "MtrunA17_Chr7g0243321")
NCR_list <- append(NCR_list, "MtrunA17_Chr2g0306721")
NCR_list <- append(NCR_list, "MtrunA17_Chr1g0167541")

for (NCR in (t %>% filter(str_detect(legoo_acronym, "MtNCR")) %>% pull(gene))) {if (NCR %in% NCR_list) {} else {NCR_list <- append(NCR_list, NCR)}}



NCR_exp_nods     <- t %>% filter(gene %in% NCR_list, exp_nods) %>% pull(gene)
NCR_exp_nods21   <- t %>% filter(gene %in% NCR_list, exp_nods21) %>% pull(gene)
NCR_exp_nods22   <- t %>% filter(gene %in% NCR_list, exp_nods22) %>% pull(gene)
NCR_exp_nodp     <- t %>% filter(gene %in% NCR_list, exp_nodp) %>% pull(gene)
NCR_exp_nodp21   <- t %>% filter(gene %in% NCR_list, exp_nodp21) %>% pull(gene)
NCR_exp_nodp22   <- t %>% filter(gene %in% NCR_list, exp_nodp22) %>% pull(gene)
NCR_exp_nodm     <- t %>% filter(gene %in% NCR_list, exp_nodm) %>% pull(gene)
NCR_exp_nodm21   <- t %>% filter(gene %in% NCR_list, exp_nodm21) %>% pull(gene)
NCR_exp_nodm22   <- t %>% filter(gene %in% NCR_list, exp_nodm22) %>% pull(gene)
NCR_exp_gall21   <- t %>% filter(gene %in% NCR_list, exp_nodm22) %>% pull(gene)
NCR_exp_gall22   <- t %>% filter(gene %in% NCR_list, exp_nodm22) %>% pull(gene)

NCR_DE_nods_n     <- t %>% filter(gene %in% NCR_list, med_nods_n.sig) %>% pull(gene)
NCR_DE_nods_r     <- t %>% filter(gene %in% NCR_list, med_nods_r.sig) %>% pull(gene)
NCR_DE_nods_r21_n <- t %>% filter(gene %in% NCR_list, med_nods_r21_n.sig) %>% pull(gene)
NCR_DE_nods_r22_n <- t %>% filter(gene %in% NCR_list, med_nods_r22_n.sig) %>% pull(gene)
NCR_DE_nods_i     <- t %>% filter(gene %in% NCR_list, med_nods_i.sig) %>% pull(gene)
NCR_DE_nodp_r     <- t %>% filter(gene %in% NCR_list, med_nodp_r.sig) %>% pull(gene)
NCR_DE_nodm_r     <- t %>% filter(gene %in% NCR_list, med_nodm_r.sig) %>% pull(gene)
NCR_DE_gall_r     <- t %>% filter(gene %in% NCR_list, med_gall_r.sig) %>% pull(gene)


plot(euler(list(
  # "nods" = NCR_exp_nods, 
  # "nods21" = NCR_exp_nods21,
  # "nods22" = NCR_exp_nods22,
  # "nodp" = NCR_exp_nodp,
  "nodp21" = NCR_exp_nodp21,
  "nodp22" = NCR_exp_nodp22)), quantities = TRUE)

plot(euler(list(
  # "nods" = NCR_exp_nods, 
  # "nods21" = NCR_exp_nods21,
  # "nods22" = NCR_exp_nods22,
  # "nodp" = NCR_exp_nodp,
  "nodm21" = NCR_exp_nodm21,
  "nodm22" = NCR_exp_nodm22,
  "nodp21" = NCR_exp_nodp21,
   "nodp22" = NCR_exp_nodp22
  )), quantities = TRUE, fills = c("dark gray", "blue", "pink", "yellow"))

setdiff(NCR_exp_nodm21, NCR_exp_nodm22)
setdiff(NCR_exp_nodp21, NCR_exp_nodp22)
setdiff(NCR_exp_nodp21, NCR_exp_nodm21)
setdiff(NCR_exp_nodp22, NCR_exp_nodm22)

plot(euler(list(
  "nods_n" = NCR_DE_nods_n,
  "nods_r" = NCR_DE_nods_r,
  "nods_r21_n" = NCR_DE_nods_r21_n,
  "nods_r22_n" = NCR_DE_nods_r22_n,
  "nods_i" = NCR_DE_nods_i,
  "nodp_r" = NCR_DE_nodp_r)), quantities = TRUE)


ComplexUpset::upset(t %>% filter(gene %in% NCR_list) %>% select(gene, exp_nodm21, exp_nodm22, exp_nodp21, exp_nodp22),
                    intersect = c("exp_nodm21", "exp_nodm22", "exp_nodp21", "exp_nodp22"))
ComplexUpset::upset(t %>% filter(gene %in% NCR_list) %>% select(gene, med_nodp_r.sig, med_nods_r.sig, med_nods_n.sig, med_nods_i.sig,
                                                                med_nods_r21_n.sigu, med_nods_r21_n.sigd, med_nods_r22_n.sig),
                    intersect = c("med_nodp_r.sig", "med_nods_r.sig", "med_nods_n.sig", "med_nods_i.sig", "med_nods_r21_n.sigu", "med_nods_r21_n.sigd", "med_nods_r22_n.sig"))
ComplexUpset::upset(t %>% filter(gene %in% NCR_list) %>% select(gene, med_nodp_r.sig, med_nods_r.sig, med_nods_n.sig, med_nods_i.sig,
                                                                med_nods_r21_n.sigu, med_nods_r21_n.sigd, med_nods_r22_n.sig, exp_nods),
                    intersect = c("med_nodp_r.sig", "med_nods_r.sig", "med_nods_n.sig", "med_nods_i.sig", "med_nods_r21_n.sigu", "med_nods_r21_n.sigd", "med_nods_r22_n.sig", "exp_nods"))

ggplot(t %>% filter(gene %in% NCR_list)) + aes(x = med_nods_r21_n.ns_LFC, y = med_nods_r22_n.ns_LFC) +
  theme_classic() + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  geom_linerange(data = . %>% filter(med_nods_r21_n.sig | med_nods_i.sig | med_nods_r.sig | med_nodp_r.sig),
                 aes(xmin = med_nods_r21_n.ns_LFC - med_nods_r21_n.ns_SE,
                     xmax = med_nods_r21_n.ns_LFC + med_nods_r21_n.ns_SE,
                     color = interaction(med_nods_r21_n.sig, med_nods_i.sig)),
                 alpha = 0.8) +
  geom_linerange(data = . %>% filter(med_nods_r21_n.sig | med_nods_i.sig | med_nods_r.sig | med_nodp_r.sig),
                 aes(ymin = med_nods_r22_n.ns_LFC - med_nods_r22_n.ns_SE,
                     ymax = med_nods_r22_n.ns_LFC + med_nods_r22_n.ns_SE,
                     color = interaction(med_nods_r21_n.sig, med_nods_i.sig)),
                 alpha = 0.8) +
  geom_point(aes(fill = interaction(med_nods_r21_n.sig, med_nods_i.sig),
                 shape = interaction(med_nods_r.sig, med_nodp_r.sig),
                 size = interaction(med_nods_r.sig, med_nodp_r.sig)),
             alpha = 0.8) +
  scale_shape_manual(values = c("TRUE.TRUE" = 24,
                                "FALSE.TRUE" = 22,
                                "TRUE.FALSE" = 23,
                                "FALSE.FALSE" = 21)) +
  scale_fill_manual(values = c("FALSE.FALSE" = "gray",
                                "FALSE.TRUE" = "gold",
                                "TRUE.FALSE" = "red",
                                "TRUE.TRUE" = "green")) +
  scale_color_manual(values = c("FALSE.FALSE" = "gray",
                               "FALSE.TRUE" = "gold",
                               "TRUE.FALSE" = "red",
                               "TRUE.TRUE" = "green")) +
  scale_size_manual(values = c("TRUE.TRUE" = 5,
                               "FALSE.TRUE" = 5,
                               "TRUE.FALSE" = 5,
                               "FALSE.FALSE" = 2)) +
  guides(fill = guide_legend( 
    override.aes=list(shape = 21)))


deleted_genes_from_Shen_2023 <- c(
"MtrunA17_Chr3g0082791",
"MtrunA17_Chr3g0082801",
"MtrunA17_Chr3g0082811",
"MtrunA17_Chr3g0082821",
"MtrunA17_Chr3g0082831",
"MtrunA17_Chr3g0082841",
"MtrunA17_Chr3g0082851",
"MtrunA17_Chr3g0082861",
"MtrunA17_Chr3g0082871",
"MtrunA17_Chr3g0082881",
"MtrunA17_Chr3g0082891",
"MtrunA17_Chr3g0082901",
"MtrunA17_Chr3g0082911",
"MtrunA17_Chr3g0082921",
"MtrunA17_Chr3g0082931",
"MtrunA17_Chr3g0082941",
"MtrunA17_Chr3g0082951",
"MtrunA17_Chr3g1009661",
"MtrunA17_Chr3g0082961",
"MtrunA17_Chr3g0082971",
"MtrunA17_Chr3g0082981",
"MtrunA17_Chr3g1009664",
"MtrunA17_Chr3g0082991",
"MtrunA17_Chr3g1009666",
"MtrunA17_Chr3g0083001",
"MtrunA17_Chr3g0083011",
"MtrunA17_Chr3g1009670",
"MtrunA17_Chr3g0083021",
"MtrunA17_Chr3g0083031",
"MtrunA17_Chr3g0083041",
"MtrunA17_Chr3g0083051",
"MtrunA17_Chr3g0083061",
"MtrunA17_Chr3g0083071",
"MtrunA17_Chr3g0083081",
"MtrunA17_Chr3g0083091",
"MtrunA17_Chr3g0083101",
"MtrunA17_Chr3g0083111",
"MtrunA17_Chr3g0083121",
"MtrunA17_Chr3g0083131",
"MtrunA17_Chr3g1009692",
"MtrunA17_Chr3g0083141",
"MtrunA17_Chr3g0083151",
"MtrunA17_Chr3g0083161",
"MtrunA17_Chr3g0083171",
"MtrunA17_Chr3g0083181",
"MtrunA17_Chr3g0083191",
"MtrunA17_Chr3g0083201",
"MtrunA17_Chr3g0083211",
"MtrunA17_Chr3g0083221",
"MtrunA17_Chr3g0083231",
"MtrunA17_Chr3g0083241",
"MtrunA17_Chr3g0083251",
"MtrunA17_Chr3g0083261"
)
t %>% filter(gene %in% NCR_list) %>% filter(med_nods_r21_n.sig | med_nods_i.sig | med_nods_r.sig | med_nodm_r.sig | med_nodp_r.sig | med_nods_n.sig | med_nods_r22_n.sig | med_nods_r21_n.sig | med_nods_i.sig) %>% filter(gene %in% deleted_genes_from_Shen_2023) %>% select(gene, starts_with("med") & ends_with(".sig"), legoo_acronym, legoo_description, legoo_publication)
