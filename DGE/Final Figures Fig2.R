library(ggplot2)
library(ggpmisc)
library(tidyverse)
library(ggh4x)
library(eulerr)
library(ggpubr)
library(cowplot)
library(grid)
library(ggplotify)
library(gridExtra)
library(ComplexUpset)

setwd(dir_main)
source('Partially paired t test.R')
setwd(dir_out)






response_to_rhizobia_label <- expression(LFC~"( "["USDA1021"]^{underline("WSM1022")}~") in galls")
response_to_rhizobia_mag_label <- expression(Magnitude~of~LFC~"( "["USDA1021"]^{underline("WSM1022")}~") in galls")
response_to_rhizobia_gall_label <- expression(LFC~"( "["USDA1021"]^{underline("WSM1022")}~") in galls")
response_to_rhizobia_nods_label <- expression(LFC~"( "["USDA1021"]^{underline("WSM1022")}~") in nodules")
response_to_rhizobia_root_label <- expression(LFC~"( "["USDA1021"]^{underline("WSM1022")}~") in roots")

response_to_parasite_label <- expression(LFC~"( "["uninfected hosts"]^{underline("  infected hosts  ")}~") in nodules")
response_to_parasite_mag_label <- expression(Magnitude~of~LFC~"( "["uninfected hosts"]^{underline("  infected hosts  ")}~") in nodules")
response_to_parasite_em21_label <- expression(LFC~"( "["uninfected hosts"]^{underline("  infected hosts  ")}~") in nodules from hosts with USDA1021")
response_to_parasite_em22_label <- expression(LFC~"( "["uninfected hosts"]^{underline("  infected hosts  ")}~") in nodules from hosts with WSM1022")



gall_color <- "indianred1"
nodule_color <- "dodgerblue"
gall_and_nodule_color <- "mediumpurple1"
USDA1021_color <- "royalblue4"
WSM1022_color <- "steelblue1"
USDA1021_rhizo_color <- "royalblue4"
WSM1022_rhizo_color <- "darkslategray1"
USDA1021_and_WSM1022_color <- "seagreen3"
gano_ME_color <- "mediumpurple3"
gano_IE_color <- "olivedrab"
gano_ME_and_IE_color <- "ivory2"
nods_ME_color <- gano_ME_color
nods_IE_color <- gano_IE_color
nods_ME_and_IE_color <- gano_ME_and_IE_color
nematode_color <- "salmon4"



##################################


t_1v2_up <- t.test.partial(paired = t %>% filter(med_nods_r21_n.sigu | med_nods_r22_n.sigu) %>% select(med_nods_r21_n.ns_LFC, med_nods_r22_n.ns_LFC) %>% drop_na() %>% as.data.frame(),
                           unpaired.x = t %>% filter(med_nods_r21_n.sigu & !med_nods_r22_n.sigu) %>% select(med_nods_r21_n.ns_LFC) %>% drop_na() %>% pull(med_nods_r21_n.ns_LFC),
                           unpaired.y = t %>% filter(!med_nods_r21_n.sigu & med_nods_r22_n.sigu) %>% select(med_nods_r22_n.ns_LFC) %>% drop_na() %>% pull(med_nods_r22_n.ns_LFC))$p.value

t_mv2_up <- t.test.partial(paired = t %>% mutate(ave = (med_nods_i_em21_n.ns_LFC + med_nods_i_em22_n.ns_LFC)/2) %>% filter((med_nods_n.sig & ave > 0) | med_nods_r22_n.sigu) %>% select(ave, med_nods_r22_n.ns_LFC) %>% drop_na() %>% as.data.frame(),
                           unpaired.x = t %>% mutate(ave = (med_nods_i_em21_n.ns_LFC + med_nods_i_em22_n.ns_LFC)/2) %>% filter((med_nods_n.sig & ave > 0) & !med_nods_r22_n.sigu) %>% select(ave) %>% drop_na() %>% pull(ave),
                           unpaired.y = t %>% mutate(ave = (med_nods_i_em21_n.ns_LFC + med_nods_i_em22_n.ns_LFC)/2) %>% filter(!(med_nods_n.sig & ave > 0) & med_nods_r22_n.sigu) %>% select(med_nods_r22_n.ns_LFC) %>% drop_na() %>% pull(med_nods_r22_n.ns_LFC))$p.value

t_mv1_up <- t.test.partial(paired = t %>% mutate(ave = (med_nods_i_em21_n.ns_LFC + med_nods_i_em22_n.ns_LFC)/2) %>% filter((med_nods_n.sig & ave > 0) | med_nods_r21_n.sigu) %>% select(ave, med_nods_r21_n.ns_LFC) %>% drop_na() %>% as.data.frame(),
                           unpaired.x = t %>% mutate(ave = (med_nods_i_em21_n.ns_LFC + med_nods_i_em22_n.ns_LFC)/2) %>% filter((med_nods_n.sig & ave > 0) & !med_nods_r21_n.sigu) %>% select(ave) %>% drop_na() %>% pull(ave),
                           unpaired.y = t %>% mutate(ave = (med_nods_i_em21_n.ns_LFC + med_nods_i_em22_n.ns_LFC)/2) %>% filter(!(med_nods_n.sig & ave > 0) & med_nods_r21_n.sigu) %>% select(med_nods_r21_n.ns_LFC) %>% drop_na() %>% pull(med_nods_r21_n.ns_LFC))$p.value



t_1v2_down <- t.test.partial(paired = t %>% filter(med_nods_r21_n.sigd | med_nods_r22_n.sigd) %>% select(med_nods_r21_n.ns_LFC, med_nods_r22_n.ns_LFC) %>% drop_na() %>% as.data.frame(),
                             unpaired.x = t %>% filter(med_nods_r21_n.sigd & !med_nods_r22_n.sigd) %>% select(med_nods_r21_n.ns_LFC) %>% drop_na() %>% pull(med_nods_r21_n.ns_LFC),
                             unpaired.y = t %>% filter(!med_nods_r21_n.sigd & med_nods_r22_n.sigd) %>% select(med_nods_r22_n.ns_LFC) %>% drop_na() %>% pull(med_nods_r22_n.ns_LFC))$p.value

t_mv2_down <- t.test.partial(paired = t %>% mutate(ave = (med_nods_i_em21_n.ns_LFC + med_nods_i_em22_n.ns_LFC)/2) %>% filter((med_nods_n.sig & ave < 0) | med_nods_r22_n.sigd) %>% select(ave, med_nods_r22_n.ns_LFC) %>% drop_na() %>% as.data.frame(),
                             unpaired.x = t %>% mutate(ave = (med_nods_i_em21_n.ns_LFC + med_nods_i_em22_n.ns_LFC)/2) %>% filter((med_nods_n.sig & ave < 0) & !med_nods_r22_n.sigd) %>% select(ave) %>% drop_na() %>% pull(ave),
                             unpaired.y = t %>% mutate(ave = (med_nods_i_em21_n.ns_LFC + med_nods_i_em22_n.ns_LFC)/2) %>% filter(!(med_nods_n.sig & ave < 0) & med_nods_r22_n.sigd) %>% select(med_nods_r22_n.ns_LFC) %>% drop_na() %>% pull(med_nods_r22_n.ns_LFC))$p.value

t_mv1_down <- t.test.partial(paired = t %>% mutate(ave = (med_nods_i_em21_n.ns_LFC + med_nods_i_em22_n.ns_LFC)/2) %>% filter((med_nods_n.sig & ave < 0) | med_nods_r21_n.sigd) %>% select(ave, med_nods_r21_n.ns_LFC) %>% drop_na() %>% as.data.frame(),
                             unpaired.x = t %>% mutate(ave = (med_nods_i_em21_n.ns_LFC + med_nods_i_em22_n.ns_LFC)/2) %>% filter((med_nods_n.sig & ave < 0) & !med_nods_r21_n.sigd) %>% select(ave) %>% drop_na() %>% pull(ave),
                             unpaired.y = t %>% mutate(ave = (med_nods_i_em21_n.ns_LFC + med_nods_i_em22_n.ns_LFC)/2) %>% filter(!(med_nods_n.sig & ave < 0) & med_nods_r21_n.sigd) %>% select(med_nods_r21_n.ns_LFC) %>% drop_na() %>% pull(med_nods_r21_n.ns_LFC))$p.value




t_1v2_mag  <- t.test.partial(paired = t %>% filter(med_nods_r21_n.sigd | med_nods_r22_n.sigd) %>% select(med_nods_r21_n.ns_LFC, med_nods_r22_n.ns_LFC) %>% drop_na() %>% abs() %>% as.data.frame(),
                             unpaired.x = t %>% filter(med_nods_r21_n.sigd & !med_nods_r22_n.sigd) %>% select(med_nods_r21_n.ns_LFC) %>% drop_na() %>% abs() %>% pull(med_nods_r21_n.ns_LFC),
                             unpaired.y = t %>% filter(!med_nods_r21_n.sigd & med_nods_r22_n.sigd) %>% select(med_nods_r22_n.ns_LFC) %>% drop_na() %>% abs() %>% pull(med_nods_r22_n.ns_LFC))$p.value

t_mv2_mag <- t.test.partial(paired = t %>% mutate(ave = (med_nods_i_em21_n.ns_LFC + med_nods_i_em22_n.ns_LFC)/2) %>% filter(med_nods_n.sig | med_nods_r22_n.sig) %>% select(ave, med_nods_r22_n.ns_LFC) %>% drop_na() %>% abs() %>% as.data.frame(),
                            unpaired.x = t %>% mutate(ave = (med_nods_i_em21_n.ns_LFC + med_nods_i_em22_n.ns_LFC)/2) %>% filter(med_nods_n.sig & !med_nods_r22_n.sig) %>% select(ave) %>% drop_na() %>% pull(ave) %>% abs(),
                            unpaired.y = t %>% mutate(ave = (med_nods_i_em21_n.ns_LFC + med_nods_i_em22_n.ns_LFC)/2) %>% filter(!med_nods_n.sig & med_nods_r22_n.sig) %>% select(med_nods_r22_n.ns_LFC) %>% drop_na() %>% pull(med_nods_r22_n.ns_LFC) %>% abs())$p.value

t_mv1_mag <- t.test.partial(paired = t %>% mutate(ave = (med_nods_i_em21_n.ns_LFC + med_nods_i_em22_n.ns_LFC)/2) %>% filter(med_nods_n.sig | med_nods_r21_n.sig) %>% select(ave, med_nods_r21_n.ns_LFC) %>% drop_na() %>% abs() %>% as.data.frame(),
                            unpaired.x = t %>% mutate(ave = (med_nods_i_em21_n.ns_LFC + med_nods_i_em22_n.ns_LFC)/2) %>% filter(med_nods_n.sig & !med_nods_r21_n.sigd) %>% select(ave) %>% drop_na() %>% pull(ave) %>% abs(),
                            unpaired.y = t %>% mutate(ave = (med_nods_i_em21_n.ns_LFC + med_nods_i_em22_n.ns_LFC)/2) %>% filter(!med_nods_n.sig & med_nods_r21_n.sigd) %>% select(med_nods_r21_n.ns_LFC) %>% drop_na() %>% pull(med_nods_r21_n.ns_LFC) %>% abs())$p.value





change_notation <- function(x) {
  if (round(x, 3) == 0)
  {return_value <- "p <0.001"}
  else 
  {return_value <- paste("p =",round(x, 3))}
  return(return_value)
}



astrification <- function(input){
  return(case_when(input <= 0.0001 ~ "****",
                   input <= 0.001 ~ "***",
                   input <= 0.01 ~ "**",
                   input <= 0.05 ~ "*",
                   input > 0.05 ~ "ns"))
}

t_1v2_up   <- astrification(t_1v2_up)
t_1v2_down <- astrification(t_1v2_down)
t_1v2_mag  <- astrification(t_1v2_mag)
t_mv1_up   <- astrification(t_mv1_up)
t_mv1_down <- astrification(t_mv1_down)
t_mv1_mag  <- astrification(t_mv1_mag)
t_mv2_up   <- astrification(t_mv2_up)
t_mv2_down <- astrification(t_mv2_down)
t_mv2_mag  <- astrification(t_mv2_mag)










################
### FIGURE 1 ###
################



f2a_table <- rbind(t %>% filter(organism == "Medicago") %>% select(LFC = med_nods_r21_n.ns_LFC, sig = med_nods_r21_n.sig, mutualism_gene) %>% mutate(category = "Host genes in\nnodules with USDA1021") %>% mutate(mutualism_gene = if_else(mutualism_gene & sig, TRUE, FALSE)) %>% drop_na(),
                   t %>% filter(organism == "Medicago") %>% select(LFC = med_nods_r22_n.ns_LFC, sig = med_nods_r22_n.sig, mutualism_gene) %>% mutate(category = "Host genes in\nnodules with WSM1022") %>% mutate(mutualism_gene = if_else(mutualism_gene & sig, TRUE, FALSE)) %>% drop_na())

f2a_main <- ggplot(f2a_table %>% arrange(sig, mutualism_gene)) + 
  geom_vline(xintercept = 0) +
  aes(x = LFC, y = category, fill = interaction(sig, category), alpha = sig, color = mutualism_gene, stroke = mutualism_gene) +
  geom_point(shape = 21, position = position_jitter(height = 0.4, seed = 3)) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.1), guide = "none") +
  scale_fill_manual(name = "Differential expression",
                    values = c("FALSE.Host genes in\nnodules with USDA1021" = NA,
                               "TRUE.Host genes in\nnodules with USDA1021" = USDA1021_color,
                               "TRUE.Host genes in\nnodules with WSM1022" = WSM1022_color),
                    labels = c("FALSE.Host genes in\nnodules with USDA1021" = "Chage in expression\nnot statistically different\nfrom zero",
                               "TRUE.Host genes in\nnodules with USDA1021" = "Differential expression of\nhost genes in nodules\nwith USDA1021",
                               "TRUE.Host genes in\nnodules with WSM1022" = "Differential expression of\nhost genes in nodules\nwith WSM1022"),
                    drop = FALSE) +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5,
                             nrow = 2, byrow = TRUE,
                             override.aes = list(
                               alpha = c(0.25, 1, 1),
                               size = 3))) +
  theme_classic() + theme(axis.title.y = element_blank(), legend.position = "bottom", axis.ticks.y = element_blank()) +
  xlab(response_to_parasite_label) + ylim(c("Host genes in\nnodules with USDA1021","Host genes in\nnodules with WSM1022")) + theme(legend.position = "none") +
  annotate("text", x = -Inf, y = -Inf,
           vjust = 0.45, hjust = 0,
           color = WSM1022_color,
           label = paste(" n =", length(gl_med_nods_r22_n), "\n\n")) +
  annotate("text", x = -Inf, y = -Inf,
           vjust = 0.45, hjust = 0,
           color = USDA1021_color,
           label = paste(" n =", length(gl_med_nods_r21_n), "\n")) +
  scale_color_manual(values = c("TRUE" = "goldenrod", "FALSE" = "black")) +
  scale_discrete_manual(aesthetic = "stroke", 
                        values = c("TRUE" = 1, "FALSE" = 0.5))

f2a_main


f2a_inset <- as.ggplot(plot(euler(list(
  "Genes responsive to parasites in nodules with USDA1021" = unique(gl_med_nods_r21_n),
  "Genes responsive to parasites in nodules with WSM1022" = unique(gl_med_nods_r22_n))), 
  quantities = list(cex = .35, col = "white"), labels = FALSE,
  fill = c(USDA1021_color, WSM1022_color)))
f2a_inset


f2a <- ggdraw() + draw_plot(f2a_main) + draw_plot(f2a_inset, x = 0.8, y = 0.35, width = 0.29, height = 0.29)
f2a




f2b_table <- rbind(t %>% filter(organism == "Rhizob21") %>% select(LFC = r21_nods_n.ns_LFC, sig = r21_nods_n.sig) %>% mutate(category = "USDA1021\ngenes") %>% drop_na(),
                   t %>% filter(organism == "Rhizob22") %>% select(LFC = r22_nods_n.ns_LFC, sig = r22_nods_n.sig) %>% mutate(category = "WSM1022\ngenes") %>% drop_na())


f2b <- ggplot(f2b_table %>% arrange(sig)) + 
  geom_vline(xintercept = 0) +
  aes(x = LFC, y = category, fill = interaction(sig, category), alpha = sig) +
  geom_point(shape = 21, color = "black", position = position_jitter(height = 0.4)) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.1), guide = "none") +
  scale_fill_manual(name = "Differential expression",
                    values = c("FALSE.USDA1021\ngenes" = NA,
                               "FALSE.WSM1022\ngenes" = NA,
                               "TRUE.USDA1021\ngenes" = USDA1021_color,
                               "TRUE.WSM1022\ngenes" = WSM1022_color),
                    labels = c("FALSE.USDA1021\ngenes" = "Chage in expression not\nstatistically different from zero",
                               "FALSE.WSM1022\ngenes" = "Chage in expression not\nstatistically different from zero",
                               "TRUE.USDA1021\ngenes" = "Differential expression of\nUSDA1021 genes",
                               "TRUE.WSM1022\ngenes" = "Differential expression of\nWSM1022 genes"),
                    drop = FALSE) +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5,
                             nrow = 1, byrow = TRUE,
                             override.aes = list(
                               alpha = c(0.25, 0.25, 1, 1),
                               size = 3))) +
  theme_classic() + theme(axis.title.y = element_blank(), legend.position = "bottom", axis.ticks.y = element_blank()) +
  xlab(response_to_parasite_label) + ylim(c("USDA1021\ngenes","WSM1022\ngenes")) + theme(legend.position = "none")  +
  annotate("text", x = -Inf, y = -Inf,
           vjust = 0.45, hjust = 0,
           color = WSM1022_color,
           label = paste(" n =", length(gl_r22_nods_n), "\n\n")) +
  annotate("text", x = -Inf, y = -Inf,
           vjust = 0.45, hjust = 0,
           color = USDA1021_color,
           label = paste(" n =", length(gl_r21_nods_n), "\n"))
f2b




f2c_subset_all  <- "Differentially expressed\nhost genes"
f2c_subset_unin <- "Differentially expressed\nhost genes upregulated\nin uninfected hosts"
f2c_subset_in   <- "Differentially expressed\nhost genes upregulated\nin infected hosts"


f2c_table <- rbind(
  t %>% filter(med_nods_r21_n.sig) %>% select(med_nods_r21_n.ns_LFC) %>% abs() %>% dplyr::rename("LFC" = med_nods_r21_n.ns_LFC) %>% mutate(subset = f2c_subset_all, test = "USDA1021"),
  t %>% filter(med_nods_r22_n.sig) %>% select(med_nods_r22_n.ns_LFC) %>% abs() %>% dplyr::rename("LFC" = med_nods_r22_n.ns_LFC) %>% mutate(subset = f2c_subset_all, test = "WSM1022"),
  t %>% filter(med_nods_r21_n.sigu) %>% select(med_nods_r21_n.ns_LFC) %>% abs() %>% dplyr::rename("LFC" = med_nods_r21_n.ns_LFC) %>% mutate(subset = f2c_subset_in, test = "USDA1021"),
  t %>% filter(med_nods_r22_n.sigu) %>% select(med_nods_r22_n.ns_LFC) %>% abs() %>% dplyr::rename("LFC" = med_nods_r22_n.ns_LFC) %>% mutate(subset = f2c_subset_in, test = "WSM1022"),
  t %>% filter(med_nods_r21_n.sigd) %>% select(med_nods_r21_n.ns_LFC) %>% abs() %>% dplyr::rename("LFC" = med_nods_r21_n.ns_LFC) %>% mutate(subset = f2c_subset_unin, test = "USDA1021"),
  t %>% filter(med_nods_r22_n.sigd) %>% select(med_nods_r22_n.ns_LFC) %>% abs() %>% dplyr::rename("LFC" = med_nods_r22_n.ns_LFC) %>% mutate(subset = f2c_subset_unin, test = "WSM1022")) %>% as.data.frame()


f2c <- ggplot(f2c_table) + 
  aes(y = subset, x = LFC, fill = test) +
  geom_vline(xintercept = 0, color = "black") +
  geom_boxplot(outlier.shape = 21) + 
  ggplot2::annotate("rect", xmin = 0, xmax = 8, ymin = 1.5, ymax = 1.51, alpha = 0.4, color = "black") +
  ggplot2::annotate("rect", xmin = 0, xmax = 8, ymin = 2.5, ymax = 2.51, alpha = 0.4, color = "gray") +
  
  ggplot2::annotate("text", size = 3, x = 0, y = 0.8 ,    hjust = 1, label = paste0("n=", length(f2c_table %>% filter(subset == f2c_subset_all) %>% filter(test == "USDA1021") %>% pull(LFC)), " "), color = USDA1021_color) +
  ggplot2::annotate("text", size = 3, x = 0, y = 1.2, hjust = 1, label = paste0("n=", length(f2c_table %>% filter(subset == f2c_subset_all) %>% filter(test == "WSM1022") %>% pull(LFC)), " "), color = WSM1022_color) +
  ggplot2::annotate("text", size = 3, x = 0, y = 1.8,    hjust = 1, label = paste0("n=", length(f2c_table %>% filter(subset == f2c_subset_in) %>% filter(test == "USDA1021") %>% pull(LFC)), " "), color = USDA1021_color) +
  ggplot2::annotate("text", size = 3, x = 0, y = 2.2, hjust = 1, label = paste0("n=", length(f2c_table %>% filter(subset == f2c_subset_in) %>% filter(test == "WSM1022") %>% pull(LFC)), " "), color = WSM1022_color) +
  ggplot2::annotate("text", size = 3, x = 0, y = 2.8,    hjust = 1, label = paste0("n=", length(f2c_table %>% filter(subset == f2c_subset_unin) %>% filter(test == "USDA1021") %>% pull(LFC)), " "), color = USDA1021_color) +
  ggplot2::annotate("text", size = 3, x = 0, y = 3.2, hjust = 1, label = paste0("n=", length(f2c_table %>% filter(subset == f2c_subset_unin) %>% filter(test == "WSM1022") %>% pull(LFC)), " "), color = WSM1022_color) +
  
  geom_segment(aes(x = 0.3,   xend = 0.3,   y = 1.2, yend = 0.8), color = "black") +
  geom_segment(aes(x = 0.3,   xend = 0.3,   y = 2.2, yend = 1.8), color = "black") +
  geom_segment(aes(x = 0.3,   xend = 0.3,   y = 3.2, yend = 2.8), color = "black") +
  
  ggplot2::annotate("text", size = 4, x = 0.2, y = 3, hjust = 0.5, angle = 90, label = t_1v2_up) +
  ggplot2::annotate("text", size = 4, x = 0.2, y = 2, hjust = 0.5, angle = 90, label = t_1v2_down) +
  ggplot2::annotate("text", size = 4, x = 0.2, y = 1, hjust = 0.5, angle = 90, label = t_1v2_mag) +
  
  theme_classic() + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "bottom") +
  xlab(response_to_parasite_mag_label) + ylab("") +
  scale_fill_manual(name = "Subset",
                    values = c("USDA1021" = USDA1021_color,
                               "WSM1022" = WSM1022_color,
                               "Average main effect" = gano_ME_color),
                    labels = c("Average main effect" = "Average main effect in\ngalls and nodules")) +
  coord_cartesian(expand = TRUE, clip = "off") + expand_limits(x = -0.5) +
  guides(fill = guide_legend("Rhizobia", title.position = "top", title.hjust = 0.5)) +
  theme(legend.background = element_rect(fill = NA, color = NA),
        legend.position = c(.8,.9))
# f2c






f2d_table <- t %>% filter(exp_nods) %>% filter(organism == "Medicago") %>% filter(med_nods_i.baseMean > 50) %>%
  select(med_nods_i_em21_n.ns_LFC, med_nods_i_em21_n.ns_SE,
         med_nods_i_em22_n.ns_LFC, med_nods_i_em22_n.ns_SE,
         med_nods_n.sig, med_nods_i.sig, mutualism_gene) %>%
  mutate(em21LFC = ifelse(is.na(med_nods_i_em21_n.ns_LFC), 0, med_nods_i_em21_n.ns_LFC),
         em21SE = ifelse(is.na(med_nods_i_em21_n.ns_SE), 0, med_nods_i_em21_n.ns_SE),
         em22LFC = ifelse(is.na(med_nods_i_em22_n.ns_LFC), 0, med_nods_i_em22_n.ns_LFC),
         em22SE = ifelse(is.na(med_nods_i_em22_n.ns_SE), 0, med_nods_i_em22_n.ns_SE),
         RSIG = ifelse(is.na(med_nods_n.sig), FALSE, med_nods_n.sig),
         ISIG = ifelse(is.na(med_nods_i.sig), FALSE, med_nods_i.sig)) %>%
  mutate(mutualism_gene = if_else(mutualism_gene & (ISIG | RSIG), TRUE, FALSE)) %>%
  select(em21LFC, em21SE, em22LFC, em22SE, RSIG, ISIG, mutualism_gene)



f2d_main <- ggplot(f2d_table) + aes(x = em21LFC, y = em22LFC, 
                                    fill = interaction(RSIG, ISIG), 
                                    color = interaction(RSIG, ISIG), 
                                    alpha = (RSIG | ISIG)) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  geom_point(data = . %>% filter(!(RSIG | ISIG)), color = "black", shape = 21) +
  geom_linerange(data = . %>% filter(RSIG | ISIG),
                 aes(xmax = em21LFC + em21SE,
                     xmin = em21LFC - em21SE),
                 linewidth = 0.2, show.legend = TRUE) +
  geom_linerange(data = . %>% filter(RSIG | ISIG),
                 aes(ymax = em22LFC + em22SE,
                     ymin = em22LFC - em22SE),
                 linewidth = 0.2, show.legend = TRUE) +
  geom_point(data = . %>% filter(RSIG | ISIG), color = "black", shape = 21) +
  geom_point(data = . %>% filter(mutualism_gene), color = "goldenrod", shape = 21, fill = NA, stroke = 2) +
  scale_color_manual(name = "Main effect and interaction tests",
                     values = c("FALSE.FALSE" = "light gray",
                                "FALSE.TRUE" = nods_IE_color,
                                "TRUE.FALSE" = nods_ME_color,
                                "TRUE.TRUE" = "black"),
                     labels = c("FALSE.FALSE" = "No effect",
                                "FALSE.TRUE" = "Infection x strain interaction",
                                "TRUE.FALSE" = "Infection main effect",
                                "TRUE.TRUE" = "Main effect and interaction")) +
  scale_fill_manual(name = "Main effect and interaction tests",
                    labels = c("FALSE.FALSE" = "No effect",
                               "FALSE.TRUE" = "Infection x strain interaction",
                               "TRUE.FALSE" = "Infection main effect",
                               "TRUE.TRUE" = "Main effect and interaction"),
                    values = c("FALSE.FALSE" = "light gray",
                               "FALSE.TRUE" = nods_IE_color,
                               "TRUE.FALSE" = nods_ME_color,
                               "TRUE.TRUE" = nods_ME_and_IE_color)) +
  guides(fill = guide_legend("Main effect and interaction tests", 
                             title.position = "top", title.hjust = 0.5,
                             nrow = 4, byrow = TRUE,
                             override.aes = list(
                               linetype = c(0,1,1,1),
                               alpha = c(0.25, 1, 1, 1),
                               size = 3)),
         color = "none") +
  scale_alpha_manual(guide = "none", values = c("FALSE" = 0.15, "TRUE" = 0.8)) +
  theme_classic() + xlab(response_to_parasite_em21_label) + ylab(response_to_parasite_em22_label) + theme(legend.position = "bottom") +
  annotate("text", x = -Inf, y = -Inf,
           vjust = 0.45, hjust = 0,
           color = nods_ME_color,
           label = paste(" n =", length(gl_med_nods_n), "\n")) +
  annotate("text", x = -Inf, y = -Inf,
           vjust = 0.45, hjust = 0,
           color = nods_IE_color,
           label = paste(" n =", length(gl_med_nods_i), "\n\n")) +
  theme(legend.background = element_rect(fill = NA, color = NA),
        legend.position = c(.175, .925), 
        legend.spacing.y = unit(0.005, "in")
        ) +
  ggplot2::annotate("text", size = 3, x = 5, y = -Inf, hjust = 0.5, vjust = 0.45, color = gano_ME_color , label = "Infection\nmain effect\n\n") +
  ggplot2::annotate("text", size = 3, x = 6.75, y = -Inf, hjust = 0.5, vjust = 0.45, color = gano_IE_color , label = "Infection x\nstrain interaction\n\n") +
  coord_cartesian(clip = "off")

f2d_inset <- as.ggplot(plot(euler(list(
  "Genes with a similar response to parasites in nodules regardless of strain" = unique(gl_med_nods_r[gl_med_nods_r %in% setdiff(gl_med_nods_r, gl_med_nods_i)]),
  "Genes with strain dependent responses to parasites" = unique(gl_med_nods_i[gl_med_nods_i %in% setdiff(gl_med_nods_i, gl_med_nods_r)]))), 
  quantities = list(cex = .75, type = c("percent")), labels = FALSE,
  fill = c(nods_ME_color, nods_IE_color, nods_ME_and_IE_color)))


f2d <- ggdraw() + draw_plot(f2d_main) + draw_plot(f2d_inset, x = 0.8, y = 0.15, width = 0.2, height = 0.2)
f2d




f2a <- f2a + ggtitle("a") + 
  theme(plot.title = element_text(size = 18, hjust = 0), plot.title.position = "plot") +
  theme(plot.margin = margin(b=0, r = 7.5))
f2b <- f2b + ggtitle("b") + 
  theme(plot.title = element_text(size = 18, hjust = 0), plot.title.position = "plot") +
  theme(plot.margin = margin(b=0, t=0, r = 7.5))
f2c <- f2c + ggtitle("c") + 
  theme(plot.title = element_text(size = 18, hjust = 0), plot.title.position = "plot") +
  theme(plot.margin = margin(b=5, t=0, r = 7.5))
f2d <- f2d + ggtitle("d") + 
  theme(plot.title = element_text(size = 18, hjust = 0), plot.title.position = "plot") +
  theme(plot.margin = margin(t=5, r = 7.5))

f2 <- as.ggplot(
  grid.arrange(grobs = list(f2a, f2b, f2c, f2d),
               row = 4,
               heights = c(5,5,7,15))
) + theme(panel.background = element_rect(fill = 'white', color = 'white'))

f2_wide <- as.ggplot(
  grid.arrange(grobs = list(f2a, f2b, f2c, f2d),
               row = 1, col = 4,
               heights = c(5,5,7),
               widths = c(4,5),
              layout_matrix = rbind(c(1,4),
                                   c(2,4),
                                   c(3,4)))) + 
  theme(panel.background = element_rect(fill = 'white', color = 'white'))

f2_long <- as.ggplot(
  grid.arrange(grobs = list(f2a, f2b, f2c, f2d),
               row = 4, col = 1,
               heights = c(2,2,3,6))) + 
  theme(panel.background = element_rect(fill = 'white', color = 'white'))

setwd(dir_out)
ggsave(filename = "Figure2.png", plot = f2, device = "png", width = 2048, height = 4096, units = "px")
ggsave(filename = "Figure2_wide.png", plot = f2_wide, device = "png", height = 2560, width = 5120, units = "px", scale = 0.75, limitsize = FALSE)
ggsave(filename = "Figure2_long.png", plot = f2_long, device = "png", width = 2560, height = 5120, units = "px", scale = 0.75, limitsize = FALSE)
