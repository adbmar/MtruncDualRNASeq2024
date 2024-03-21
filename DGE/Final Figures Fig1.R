library(ggplot2)
# library(ggpmisc)
library(ggh4x)
library(eulerr)
library(ggpubr)
library(cowplot)
library(grid)
library(ggplotify)
library(gridExtra)
library(ComplexUpset)
library(tidyverse)

select <- dplyr::select

setwd(dir_main)
source('Partially paired t test.R')
setwd(dir_out)

# response_to_rhizobia_label <- expression(atop(Log[2]~fold~change~of~normalized~read~count,
#                                               across~rhizobia~strain~"( "["USDA1021"]^{underline("WSM1022")}~")"))
# 
# response_to_rhizobia_mag_label <- expression(atop(Magnitude~of~log[2]~fold~change~of~normalized~read~count,
#                                                   across~rhizobia~strain~"( "["USDA1021"]^{underline("WSM1022")}~")"))
# 
# response_to_rhizobia_gall_label <- expression(atop(Log[2]~fold~change~of~normalized~read~count,
#                                               across~rhizobia~strain~"in"~galls~"( "["USDA1021"]^{underline("WSM1022")}~")"))
# 
# response_to_rhizobia_nods_label <- expression(atop(Log[2]~fold~change~of~normalized~read~count,
#                                               across~rhizobia~strain~"in"~nodules~"( "["USDA1021"]^{underline("WSM1022")}~")"))
# 
# response_to_rhizobia_root_label <- expression(atop(Log[2]~fold~change~of~normalized~read~count,
#                                               across~rhizobia~strain~"in"~roots~"( "["USDA1021"]^{underline("WSM1022")}~")"))
# 
# 
# response_to_parasite_label <- expression(atop(Log[2]~fold~change~of~normalized~read~count,
#                                          across~parasite~status~"( "["uninfected hosts"]^{underline("infected hosts")}~")"))
# 
# response_to_parasite_em21_label <- expression(atop(Log[2]~fold~change~of~normalized~read~count,
#                                               across~parasite~status~"in"~hosts~with~USDA1021~"( "["uninfected hosts"]^{underline("infected hosts")}~")"))
# 
# response_to_parasite_em22_label <- expression(atop(Log[2]~fold~change~of~normalized~read~count,
#                                               across~parasite~status~"in"~hosts~with~WSM1022~"( "["uninfected hosts"]^{underline("infected hosts")}~")"))
# 




response_to_rhizobia_label <- expression(LFC~"( "["USDA1021"]^{underline("WSM1022")}~") in galls")
response_to_rhizobia_mag_label <- expression(Magnitude~of~LFC~"( "["USDA1021"]^{underline("WSM1022")}~")")
response_to_rhizobia_gall_label <- expression(LFC~"( "["USDA1021"]^{underline("WSM1022")}~") in galls")
response_to_rhizobia_nods_label <- expression(LFC~"( "["USDA1021"]^{underline("WSM1022")}~") in nodules")
response_to_rhizobia_root_label <- expression(LFC~"( "["USDA1021"]^{underline("WSM1022")}~") in roots")


response_to_parasite_label <- expression(LFC~"( "["uninfected hosts"]^{underline("infected hosts")}~")")
response_to_parasite_em21_label <- expression(LFC~"( "["uninfected hosts"]^{underline("infected hosts")}~") in nodules")
response_to_parasite_em22_label <- expression(LFC~"( "["uninfected hosts"]^{underline("infected hosts")}~") in nodules")




gall_color <- "indianred1"
nodule_color <- "dodgerblue"
gall_and_nodule_color <- "mediumpurple1"
USDA1021_color <- "cadetblue1"
WSM1022_color <- "gold1"
USDA1021_and_WSM1022_color <- "seagreen3"
gano_ME_color <- "mediumpurple3"
gano_IE_color <- "olivedrab"
gano_ME_and_IE_color <- "ivory2"
nods_ME_color <- "orange"
nods_IE_color <- "purple"
nods_ME_and_IE_color <- "floralwhite"
nematode_color <- "salmon4"

##################################

t_gvn_up <- t.test.partial(paired = t %>% filter(med_gall_r.sigu | med_nodp_r.sigu) %>% select(med_gall_r.ns_LFC, med_nodp_r.ns_LFC) %>% drop_na(),
                           unpaired.x = t %>% filter(med_gall_r.sigu & !med_nodp_r.sigu) %>% select(med_gall_r.ns_LFC) %>% drop_na() %>% pull(med_gall_r.ns_LFC),
                           unpaired.y = t %>% filter(!med_gall_r.sigu & med_nodp_r.sigu) %>% select(med_nodp_r.ns_LFC) %>% drop_na() %>% pull(med_nodp_r.ns_LFC))$p.value

t_mvn_up <- t.test.partial(paired = t %>% mutate(ave = (med_gano_i_gall_r.ns_LFC + med_gano_i_nodp_r.ns_LFC)/2) %>% filter((med_gano_r.sig & ave > 0) | med_nodp_r.sigu) %>% select(ave, med_nodp_r.ns_LFC) %>% drop_na(),
                           unpaired.x = t %>% mutate(ave = (med_gano_i_gall_r.ns_LFC + med_gano_i_nodp_r.ns_LFC)/2) %>% filter((med_gano_r.sig & ave > 0) & !med_nodp_r.sigu) %>% select(ave) %>% drop_na() %>% pull(ave),
                           unpaired.y = t %>% mutate(ave = (med_gano_i_gall_r.ns_LFC + med_gano_i_nodp_r.ns_LFC)/2) %>% filter(!(med_gano_r.sig & ave > 0) & med_nodp_r.sigu) %>% select(med_nodp_r.ns_LFC) %>% drop_na() %>% pull(med_nodp_r.ns_LFC))$p.value

t_mvg_up <- t.test.partial(paired = t %>% mutate(ave = (med_gano_i_gall_r.ns_LFC + med_gano_i_nodp_r.ns_LFC)/2) %>% filter((med_gano_r.sig & ave > 0) | med_gall_r.sigu) %>% select(ave, med_gall_r.ns_LFC) %>% drop_na(),
                           unpaired.x = t %>% mutate(ave = (med_gano_i_gall_r.ns_LFC + med_gano_i_nodp_r.ns_LFC)/2) %>% filter((med_gano_r.sig & ave > 0) & !med_gall_r.sigu) %>% select(ave) %>% drop_na() %>% pull(ave),
                           unpaired.y = t %>% mutate(ave = (med_gano_i_gall_r.ns_LFC + med_gano_i_nodp_r.ns_LFC)/2) %>% filter(!(med_gano_r.sig & ave > 0) & med_gall_r.sigu) %>% select(med_gall_r.ns_LFC) %>% drop_na() %>% pull(med_gall_r.ns_LFC))$p.value



t_gvn_down <- t.test.partial(paired = t %>% filter(med_gall_r.sigd | med_nodp_r.sigd) %>% select(med_gall_r.ns_LFC, med_nodp_r.ns_LFC) %>% drop_na(),
                             unpaired.x = t %>% filter(med_gall_r.sigd & !med_nodp_r.sigd) %>% select(med_gall_r.ns_LFC) %>% drop_na() %>% pull(med_gall_r.ns_LFC),
                             unpaired.y = t %>% filter(!med_gall_r.sigd & med_nodp_r.sigd) %>% select(med_nodp_r.ns_LFC) %>% drop_na() %>% pull(med_nodp_r.ns_LFC))$p.value

t_mvn_down <- t.test.partial(paired = t %>% mutate(ave = (med_gano_i_gall_r.ns_LFC + med_gano_i_nodp_r.ns_LFC)/2) %>% filter((med_gano_r.sig & ave < 0) | med_nodp_r.sigd) %>% select(ave, med_nodp_r.ns_LFC) %>% drop_na(),
                             unpaired.x = t %>% mutate(ave = (med_gano_i_gall_r.ns_LFC + med_gano_i_nodp_r.ns_LFC)/2) %>% filter((med_gano_r.sig & ave < 0) & !med_nodp_r.sigd) %>% select(ave) %>% drop_na() %>% pull(ave),
                             unpaired.y = t %>% mutate(ave = (med_gano_i_gall_r.ns_LFC + med_gano_i_nodp_r.ns_LFC)/2) %>% filter(!(med_gano_r.sig & ave < 0) & med_nodp_r.sigd) %>% select(med_nodp_r.ns_LFC) %>% drop_na() %>% pull(med_nodp_r.ns_LFC))$p.value

t_mvg_down <- t.test.partial(paired = t %>% mutate(ave = (med_gano_i_gall_r.ns_LFC + med_gano_i_nodp_r.ns_LFC)/2) %>% filter((med_gano_r.sig & ave < 0) | med_gall_r.sigd) %>% select(ave, med_gall_r.ns_LFC) %>% drop_na(),
                             unpaired.x = t %>% mutate(ave = (med_gano_i_gall_r.ns_LFC + med_gano_i_nodp_r.ns_LFC)/2) %>% filter((med_gano_r.sig & ave < 0) & !med_gall_r.sigd) %>% select(ave) %>% drop_na() %>% pull(ave),
                             unpaired.y = t %>% mutate(ave = (med_gano_i_gall_r.ns_LFC + med_gano_i_nodp_r.ns_LFC)/2) %>% filter(!(med_gano_r.sig & ave < 0) & med_gall_r.sigd) %>% select(med_gall_r.ns_LFC) %>% drop_na() %>% pull(med_gall_r.ns_LFC))$p.value




t_gvn_mag  <- t.test.partial(paired = t %>% filter(med_gall_r.sigd | med_nodp_r.sigd) %>% select(med_gall_r.ns_LFC, med_nodp_r.ns_LFC) %>% drop_na() %>% abs(),
                             unpaired.x = t %>% filter(med_gall_r.sigd & !med_nodp_r.sigd) %>% select(med_gall_r.ns_LFC) %>% drop_na() %>% abs() %>% pull(med_gall_r.ns_LFC),
                             unpaired.y = t %>% filter(!med_gall_r.sigd & med_nodp_r.sigd) %>% select(med_nodp_r.ns_LFC) %>% drop_na() %>% abs() %>% pull(med_nodp_r.ns_LFC))$p.value

t_mvn_mag <- t.test.partial(paired = t %>% mutate(ave = (med_gano_i_gall_r.ns_LFC + med_gano_i_nodp_r.ns_LFC)/2) %>% filter(med_gano_r.sig | med_nodp_r.sig) %>% select(ave, med_nodp_r.ns_LFC) %>% drop_na() %>% abs(),
                            unpaired.x = t %>% mutate(ave = (med_gano_i_gall_r.ns_LFC + med_gano_i_nodp_r.ns_LFC)/2) %>% filter(med_gano_r.sig & !med_nodp_r.sig) %>% select(ave) %>% drop_na() %>% pull(ave) %>% abs(),
                            unpaired.y = t %>% mutate(ave = (med_gano_i_gall_r.ns_LFC + med_gano_i_nodp_r.ns_LFC)/2) %>% filter(!med_gano_r.sig & med_nodp_r.sig) %>% select(med_nodp_r.ns_LFC) %>% drop_na() %>% pull(med_nodp_r.ns_LFC) %>% abs())$p.value

t_mvg_mag <- t.test.partial(paired = t %>% mutate(ave = (med_gano_i_gall_r.ns_LFC + med_gano_i_nodp_r.ns_LFC)/2) %>% filter(med_gano_r.sig | med_gall_r.sig) %>% select(ave, med_gall_r.ns_LFC) %>% drop_na() %>% abs(),
                            unpaired.x = t %>% mutate(ave = (med_gano_i_gall_r.ns_LFC + med_gano_i_nodp_r.ns_LFC)/2) %>% filter(med_gano_r.sig & !med_gall_r.sigd) %>% select(ave) %>% drop_na() %>% pull(ave) %>% abs(),
                            unpaired.y = t %>% mutate(ave = (med_gano_i_gall_r.ns_LFC + med_gano_i_nodp_r.ns_LFC)/2) %>% filter(!med_gano_r.sig & med_gall_r.sigd) %>% select(med_gall_r.ns_LFC) %>% drop_na() %>% pull(med_gall_r.ns_LFC) %>% abs())$p.value


astrification <- function(input){
  return(case_when(input <= 0.0001 ~ "****",
                   input <= 0.001 ~ "***",
                   input <= 0.01 ~ "**",
                   input <= 0.05 ~ "*",
                   input > 0.05 ~ "ns"))
  }

t_mvg_up   <- astrification(t_mvg_up)
t_mvg_down <- astrification(t_mvg_down)
t_mvg_mag  <- astrification(t_mvg_mag)
t_gvn_up   <- astrification(t_gvn_up)
t_gvn_down <- astrification(t_gvn_down)
t_gvn_mag  <- astrification(t_gvn_mag)
t_mvn_up   <- astrification(t_mvn_up)
t_mvn_down <- astrification(t_mvn_down)
t_mvn_mag  <- astrification(t_mvn_mag)


################
### FIGURE 1 ###
################

f1a_table <- rbind(t %>% filter(organism == "Medicago") %>% select(LFC = med_gall_r.ns_LFC, sig = med_gall_r.sig) %>% mutate(category = "Host\ngenes") %>% drop_na())

f1a <- ggplot(f1a_table %>% arrange(sig)) + 
  geom_vline(xintercept = 0, color = "gray") +
  aes(x = LFC, y = category, fill = interaction(sig, category), alpha = sig) +
  geom_point(shape = 21, color = "black", position = position_jitter(height = 0.4)) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.1), guide = "none") +
  scale_fill_manual(name = "Differential expression",
                    values = c("FALSE.Host\ngenes" = NA,
                               "TRUE.Host\ngenes" = gall_color),
                    labels = c("FALSE.Host\ngenes" = "Chage in expression not\nstatistically different from zero",
                               "TRUE.Host\ngenes" = "Differential expression of host genes")) +
  theme_classic() + theme(axis.title.y = element_blank(), legend.position = "bottom", axis.ticks.y = element_blank()) + 
  guides(fill = guide_legend( title.position = "top", 
                              title.hjust = 0.5,
                              nrow = 1, byrow = TRUE)) +
  xlab(response_to_rhizobia_label) + theme(legend.position = "none") +
  annotate("text", x = -Inf, y = -Inf,
           vjust = 0.45, hjust = 0,
           color = gall_color,
           label = paste(" n =", length(gl_med_gall_r), "\n"))
 # f1a


f1b_table <- rbind(t %>% filter(organism == "Nematode") %>% select(LFC = nem_gall_r.ns_LFC, sig = nem_gall_r.sig) %>% mutate(category = "Nematode\ngenes") %>% drop_na())

f1b <- ggplot(f1b_table %>% arrange(sig)) + 
  geom_vline(xintercept = 0) +
  aes(x = LFC, y = category, fill = interaction(sig, category), alpha = sig) +
  geom_point(shape = 21, color = "black", position = position_jitter(height = 0.4)) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.1), guide = "none") +
  scale_fill_manual(name = "Differential expression",
                    values = c("FALSE.Nematode\ngenes" = NA,
                               "TRUE.Nematode\ngenes" = nematode_color),
                    labels = c("FALSE.Nematode\ngenes" = "Chage in expression not\nstatistically different from zero",
                               "TRUE.Nematode\ngenes" = "Differential expression of nematode genes")) +
  theme_classic() + theme(axis.title.y = element_blank(), legend.position = "bottom", axis.ticks.y = element_blank()) + 
  guides(fill = guide_legend(title.position = "top", 
                             title.hjust = 0.5,
                             nrow = 1, byrow = TRUE)) +
  xlab(response_to_rhizobia_label) + theme(legend.position = "none") +
  annotate("text", x = -Inf, y = -Inf,
           vjust = 0.45, hjust = 0,
           color = nematode_color,
           label = paste(" n =", length(gl_nem_gall_r), "\n"))





f1c_subset_all<- "Differentially\nexpressed\nhost genes"
f1c_subset_21 <- "Host genes\nupregulated\nwith USDA1021"
f1c_subset_22 <- "Host genes\nupregulated\nwith WSM1022"

f1c_table <- rbind(
  t %>% filter(med_gall_r.sig) %>% select(med_gall_r.ns_LFC) %>% abs() %>% dplyr::rename("LFC" = med_gall_r.ns_LFC) %>% mutate(subset = f1c_subset_all, test = "Galls"),
  t %>% filter(med_nodp_r.sig) %>% select(med_nodp_r.ns_LFC) %>% abs() %>% dplyr::rename("LFC" = med_nodp_r.ns_LFC) %>% mutate(subset = f1c_subset_all, test = "Nodules"),
  t %>% filter(med_gall_r.sigu) %>% select(med_gall_r.ns_LFC) %>% abs() %>% dplyr::rename("LFC" = med_gall_r.ns_LFC) %>% mutate(subset = f1c_subset_22, test = "Galls"),
  t %>% filter(med_nodp_r.sigu) %>% select(med_nodp_r.ns_LFC) %>% abs() %>% dplyr::rename("LFC" = med_nodp_r.ns_LFC) %>% mutate(subset = f1c_subset_22, test = "Nodules"),
  t %>% filter(med_gall_r.sigd) %>% select(med_gall_r.ns_LFC) %>% abs() %>% dplyr::rename("LFC" = med_gall_r.ns_LFC) %>% mutate(subset = f1c_subset_21, test = "Galls"),
  t %>% filter(med_nodp_r.sigd) %>% select(med_nodp_r.ns_LFC) %>% abs() %>% dplyr::rename("LFC" = med_nodp_r.ns_LFC) %>% mutate(subset = f1c_subset_21, test = "Nodules")) %>% as.data.frame()

f1c <- ggplot(f1c_table) + 
  aes(y = subset, x = LFC, fill = test) +
  geom_vline(xintercept = 0, color = "black") +
  geom_boxplot(outlier.shape = 21) + 
  ggplot2::annotate("rect", xmin = 0, xmax = 9, ymin = 1.5, ymax = 1.51, alpha = 0.4, color = "black") +
  ggplot2::annotate("rect", xmin = 0, xmax = 9, ymin = 2.5, ymax = 2.51, alpha = 0.4, color = "gray") +
  
  ggplot2::annotate("text", size = 3, x = 0, y = 0.8, hjust = 1, label = paste0("n = ", length(f1c_table %>% filter(subset == f1c_subset_all) %>% filter(test == "Galls")   %>% pull(LFC)), " "), color = gall_color)   +
  ggplot2::annotate("text", size = 3, x = 0, y = 1.2, hjust = 1, label = paste0("n = ", length(f1c_table %>% filter(subset == f1c_subset_all) %>% filter(test == "Nodules") %>% pull(LFC)), " "), color = nodule_color) +
  ggplot2::annotate("text", size = 3, x = 0, y = 1.8, hjust = 1, label = paste0("n = ", length(f1c_table %>% filter(subset == f1c_subset_21)  %>% filter(test == "Galls")   %>% pull(LFC)), " "), color = gall_color)   +
  ggplot2::annotate("text", size = 3, x = 0, y = 2.2, hjust = 1, label = paste0("n = ", length(f1c_table %>% filter(subset == f1c_subset_21)  %>% filter(test == "Nodules") %>% pull(LFC)), " "), color = nodule_color) +
  ggplot2::annotate("text", size = 3, x = 0, y = 2.8, hjust = 1, label = paste0("n = ", length(f1c_table %>% filter(subset == f1c_subset_22)  %>% filter(test == "Galls")   %>% pull(LFC)), " "), color = gall_color)   +
  ggplot2::annotate("text", size = 3, x = 0, y = 3.2, hjust = 1, label = paste0("n = ", length(f1c_table %>% filter(subset == f1c_subset_22)  %>% filter(test == "Nodules") %>% pull(LFC)), " "), color = nodule_color) +

  geom_segment(aes(x = 0.7,   xend = 0.7,   y = 1.2, yend = 0.8), color = "black") +
  geom_segment(aes(x = 0.7,   xend = 0.7,   y = 2.2, yend = 1.8), color = "black") +
  geom_segment(aes(x = 0.7,   xend = 0.7,   y = 3.2, yend = 2.8), color = "black") +

  ggplot2::annotate("text", size = 4, x = 0.55, y = 3, hjust = 0.5, angle = 90, label = t_gvn_up) +
  ggplot2::annotate("text", size = 4, x = 0.55, y = 2, hjust = 0.5, angle = 90, label = t_gvn_down) +
  ggplot2::annotate("text", size = 4, x = 0.55, y = 1, hjust = 0.5, angle = 90, label = t_gvn_mag) +

  theme_classic() + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank()) +
  xlab(response_to_rhizobia_mag_label) + ylab("") +
  scale_fill_manual(name = "Subset",
                    values = c("Galls" = gall_color,
                               "Nodules" = nodule_color,
                               "Average main effect" = gano_ME_color),
                    labels = c("Average main effect" = "Average main effect in\ngalls and nodules")) +
  coord_cartesian(expand = TRUE, clip = "off") + expand_limits(x = -0.5) +
  guides(fill = guide_legend("Organ", title.position = "top", title.hjust = 0.5)) +
  theme(legend.background = element_rect(fill = NA, color = NA),
        legend.position = c(.9,.9))
 f1c





f1d_table <- t %>% filter(exp_gall | exp_nodp) %>% filter(organism == "Medicago") %>% filter(med_gano_i.baseMean > 50) %>%
  select(med_gano_i_gall_r.ns_LFC, med_gano_i_gall_r.ns_SE,
         med_gano_i_nodp_r.ns_LFC, med_gano_i_nodp_r.ns_SE,
         med_gano_r.sig, med_gano_i.sig) %>%
  mutate(nodpLFC = ifelse(is.na(med_gano_i_nodp_r.ns_LFC), 0, med_gano_i_nodp_r.ns_LFC),
         nodpSE = ifelse(is.na(med_gano_i_nodp_r.ns_SE), 0, med_gano_i_nodp_r.ns_SE),
         gallLFC = ifelse(is.na(med_gano_i_gall_r.ns_LFC), 0, med_gano_i_gall_r.ns_LFC),
         gallSE = ifelse(is.na(med_gano_i_gall_r.ns_SE), 0, med_gano_i_gall_r.ns_SE),
         RSIG = ifelse(is.na(med_gano_r.sig), FALSE, med_gano_r.sig),
         ISIG = ifelse(is.na(med_gano_i.sig), FALSE, med_gano_i.sig)) %>%
  select(nodpLFC, nodpSE, gallLFC, gallSE, RSIG, ISIG)

f1d_main <- ggplot(f1d_table) + aes(x = nodpLFC, y = gallLFC, 
                               fill = interaction(RSIG, ISIG), 
                               color = interaction(RSIG, ISIG), 
                               alpha = (RSIG | ISIG)) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  geom_point(data = . %>% filter(!(RSIG | ISIG)), color = "black", shape = 21) +
  geom_linerange(data = . %>% filter(RSIG | ISIG),
                 aes(xmax = nodpLFC + nodpSE,
                     xmin = nodpLFC - nodpSE),
                 linewidth = 0.2, show.legend = TRUE) +
  geom_linerange(data = . %>% filter(RSIG | ISIG),
                 aes(ymax = gallLFC + gallSE,
                     ymin = gallLFC - gallSE),
                 linewidth = 0.2, show.legend = TRUE) +
  geom_point(data = . %>% filter(RSIG | ISIG), color = "black", shape = 21) +
  scale_color_manual(name = "Main effect and interaction tests",
                     values = c("FALSE.FALSE" = "light gray",
                                "FALSE.TRUE" = gano_IE_color,
                                "TRUE.FALSE" = gano_ME_color,
                                "TRUE.TRUE" = "black"),
                     labels = c("FALSE.FALSE" = "No effect",
                                "FALSE.TRUE" = "Rhizobia x organ interaction",
                                "TRUE.FALSE" = "Rhizobia main effect",
                                "TRUE.TRUE" = "Main effect and interaction")) +
  scale_fill_manual(name = "Main effect and interaction tests",
                    labels = c("FALSE.FALSE" = "No effect",
                               "FALSE.TRUE" = "Rhizobia x organ interaction",
                               "TRUE.FALSE" = "Rhizobia main effect",
                               "TRUE.TRUE" = "Main effect and interaction"),
                    values = c("FALSE.FALSE" = "light gray",
                               "FALSE.TRUE" = gano_IE_color,
                               "TRUE.FALSE" = gano_ME_color,
                               "TRUE.TRUE" = gano_ME_and_IE_color)) +
  guides(fill = guide_legend("Main effect and interaction tests", 
                             title.position = "top", title.hjust = 0.5,
                             nrow = 4, byrow = TRUE,
                             override.aes = list(
                               linetype = c(0,1,1,1),
                               alpha = c(0.25, 1, 1, 1),
                               size = 3)),
         color = "none") +
  scale_alpha_manual(guide = "none", values = c("FALSE" = 0.15, "TRUE" = 0.8)) +
  theme_classic() + xlab(response_to_rhizobia_nods_label) + ylab(response_to_rhizobia_gall_label)

f1d_inset <- as.ggplot(plot(euler(list(
  "gano_r" = gl_med_garo_r[gl_med_garo_r %in% setdiff(gl_med_gano_r, gl_med_gano_i)],
  "gano_i" = gl_med_garo_i[gl_med_garo_i %in% setdiff(gl_med_gano_i, gl_med_gano_r)])), 
  quantities = list(type = c("percent"),cex = 0.75), labels = FALSE,
  fill = c(gano_ME_color, gano_IE_color, gano_ME_and_IE_color)))
f1d_inset

f1d_annotation <- data.frame(xpos = c(-Inf, -Inf),
                             ypos = c(-Inf, -Inf + 0.5),
                             annotateText = c(paste(" n =", length(gl_med_gano_r), "\n\n"), 
                                              paste(" n =", length(gl_med_gano_i), "\n")),
                             color_ann = c("TRUE.FALSE", "FALSE.TRUE"))

f1d_annotated <- f1d_main + 
  geom_text(data = f1d_annotation, 
            inherit.aes = FALSE, 
            aes(x = xpos, 
                y = ypos, 
                label = annotateText,
                color = color_ann),
            vjust = 0.45, hjust = 0) +
  theme(legend.background = element_rect(fill = NA, color = NA),
        legend.position = c(.225,.9),
        legend.spacing.y = unit(0.005, "in")
  ) +
  ggplot2::annotate("text", size = 3, x = 3.9, y = -Inf, hjust = 0.5, vjust = 0.45, color = gano_ME_color , label = "Rhizobia\nmain effect\n\n") +
  ggplot2::annotate("text", size = 3, x = 7,   y = -Inf, hjust = 0.5, vjust = 0.45, color = gano_IE_color , label = "Rhizobia x\norgan interaction\n\n")


  
f1d <- ggdraw() + draw_plot(f1d_annotated) + draw_plot(f1d_inset, x = 0.725, y = 0.2, width = 0.2, height = 0.2)
# f1d








f1e_table <- t %>% filter(exp_gall | exp_root) %>% filter(organism == "Medicago") %>% filter(med_garo_i.baseMean > 50) %>%
  select(med_garo_i_gall_r.ns_LFC, med_garo_i_gall_r.ns_SE,
         med_garo_i_root_r.ns_LFC, med_garo_i_root_r.ns_SE,
         med_garo_r.sig, med_garo_i.sig) %>%
  mutate(rootLFC = ifelse(is.na(med_garo_i_root_r.ns_LFC), 0, med_garo_i_root_r.ns_LFC),
         rootSE = ifelse(is.na(med_garo_i_root_r.ns_SE), 0, med_garo_i_root_r.ns_SE),
         gallLFC = ifelse(is.na(med_garo_i_gall_r.ns_LFC), 0, med_garo_i_gall_r.ns_LFC),
         gallSE = ifelse(is.na(med_garo_i_gall_r.ns_SE), 0, med_garo_i_gall_r.ns_SE),
         RSIG = ifelse(is.na(med_garo_r.sig), FALSE, med_garo_r.sig),
         ISIG = ifelse(is.na(med_garo_i.sig), FALSE, med_garo_i.sig)) %>%
  select(rootLFC, rootSE, gallLFC, gallSE, RSIG, ISIG)

f1e_main <- ggplot(f1e_table) + aes(x = rootLFC, y = gallLFC, 
                                    fill = interaction(RSIG, ISIG), 
                                    color = interaction(RSIG, ISIG), 
                                    alpha = (RSIG | ISIG)) +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
  geom_point(data = . %>% filter(!(RSIG | ISIG)), color = "black", shape = 21) +
  geom_linerange(data = . %>% filter(RSIG | ISIG),
                 aes(xmax = rootLFC + rootSE,
                     xmin = rootLFC - rootSE),
                 linewidth = 0.2, show.legend = TRUE) +
  geom_linerange(data = . %>% filter(RSIG | ISIG),
                 aes(ymax = gallLFC + gallSE,
                     ymin = gallLFC - gallSE),
                 linewidth = 0.2, show.legend = TRUE) +
  geom_point(data = . %>% filter(RSIG | ISIG), color = "black", shape = 21) +
  scale_color_manual(name = "Main effect and interaction tests",
                     values = c("FALSE.FALSE" = "light gray",
                                "FALSE.TRUE" = gano_IE_color,
                                "TRUE.FALSE" = gano_ME_color,
                                "TRUE.TRUE" = "black"),
                     labels = c("FALSE.FALSE" = "No effect",
                                "FALSE.TRUE" = "Rhizobia x organ interaction",
                                "TRUE.FALSE" = "Rhizobia main effect",
                                "TRUE.TRUE" = "Main effect and interaction")) +
  scale_fill_manual(name = "Main effect and interaction tests",
                    labels = c("FALSE.FALSE" = "No effect",
                               "FALSE.TRUE" = "Rhizobia x organ interaction",
                               "TRUE.FALSE" = "Rhizobia main effect",
                               "TRUE.TRUE" = "Main effect and interaction"),
                    values = c("FALSE.FALSE" = "light gray",
                               "FALSE.TRUE" = gano_IE_color,
                               "TRUE.FALSE" = gano_ME_color,
                               "TRUE.TRUE" = gano_ME_and_IE_color)) +
  guides(fill = guide_legend("Main effect and interaction tests", 
                             title.position = "top", title.hjust = 0.5,
                             nrow = 4, byrow = TRUE,
                             override.aes = list(
                               linetype = c(0,1,1,1),
                               alpha = c(0.25, 1, 1, 1),
                               size = 3)),
         color = "none") +
  scale_alpha_manual(guide = "none", values = c("FALSE" = 0.15, "TRUE" = 0.8)) +
  theme_classic() + xlab(response_to_rhizobia_root_label) + ylab(response_to_rhizobia_gall_label) + theme(legend.position = "bottom")
# f1e_main

f1e_inset <- as.ggplot(plot(euler(list(
  "garo_r" = gl_med_garo_r[gl_med_garo_r %in% setdiff(gl_med_garo_r, gl_med_garo_i)],
  "garo_i" = gl_med_garo_i[gl_med_garo_i %in% setdiff(gl_med_garo_i, gl_med_garo_r)])), 
  quantities = list(type = "percent",cex = 0.75), labels = FALSE,
  fill = c(gano_ME_color, gano_IE_color, gano_ME_and_IE_color)))
# f1e_inset

f1e_annotation <- data.frame(xpos = c(-Inf, -Inf),
                             ypos = c(-Inf, -Inf + 0.5),
                             annotateText = c(paste(" n =", length(gl_med_garo_r), "\n\n"), 
                                              paste(" n =", length(gl_med_garo_i), "\n")),
                             color_ann = c("TRUE.FALSE", "FALSE.TRUE"))
# f1e_annotation

f1e_annotated <- f1e_main + geom_text(data = f1e_annotation, 
                                 inherit.aes = FALSE, 
                                 aes(x = xpos, 
                                     y = ypos, 
                                     label = annotateText,
                                     color = color_ann),
                                 vjust = 0.45, hjust = 0) +
  theme(legend.background = element_rect(fill = NA, color = NA),
        legend.position = c(.225,.9),
        legend.spacing.y = unit(0.005, "in")) +
  ggplot2::annotate("text", size = 3, x = 3.9, y = -Inf, hjust = 0.5, vjust = 0.45, color = gano_ME_color , label = "Rhizobia\nmain effect\n\n") +
  ggplot2::annotate("text", size = 3, x = 6.3, y = -Inf, hjust = 0.5, vjust = 0.45, color = gano_IE_color , label = "Rhizobia x\norgan interaction\n\n")
  
# f1e_main

f1e <- ggdraw() + draw_plot(f1e_annotated) + draw_plot(f1e_inset, x = 0.725, y = 0.2, width = 0.2, height = 0.2)
# f1e







responsive_to_mutualist_strain <- t %>% filter(med_gano_r.sig |
                                                 med_garo_r.sig |
                                                 med_noro_r.sig |
                                                 med_gano_i.sig |
                                                 med_noro_i.sig |
                                                 med_garo_i.sig |
                                                 med_gall_r.sig |
                                                 med_nodp_r.sig |
                                                 med_nods_r.sig |
                                                 med_nods_i.sig |
                                                 med_root_r.sig )

cols <- c("med_noro_i.sig", "med_garo_i.sig", "med_gano_i.sig", "med_noro_r.sig", "med_garo_r.sig", "med_gano_r.sig")
f1f <- ComplexUpset::upset(responsive_to_mutualist_strain,
                           cols, sort_sets = FALSE,
                           name = "Response to mutualist strain\nacross tissue type combinations",
                           labeller = ggplot2::as_labeller(c(
  "med_gano_i.sig" = "Different response in galls and nodules",
  "med_garo_i.sig" = "Different response in galls and roots",
  "med_noro_i.sig" = "Different response in nodules and roots",
  "med_gano_r.sig" = "Similar response in galls and nodules",
  "med_garo_r.sig" = "Similar response in galls and roots",
  "med_noro_r.sig" = "Similar response in nodules and roots"
)),
stripes = c(gano_IE_color, gano_IE_color, gano_IE_color, gano_ME_color, gano_ME_color, gano_ME_color))




f1a <- f1a + ggtitle("a\n") + 
  theme(plot.title = element_text(size = 18, hjust = 0), plot.title.position = "plot") +
  theme(plot.margin = margin(b=5))
f1b <- f1b + ggtitle("b\n") + 
  theme(plot.title = element_text(size = 18, hjust = 0), plot.title.position = "plot") +
  theme(plot.margin = margin(b=5, t=5))
f1c <- f1c + ggtitle("c\n") + 
  theme(plot.title = element_text(size = 18, hjust = 0), plot.title.position = "plot") +
  theme(plot.margin = margin(b=5, t=5))
f1d <- f1d + ggtitle("d\n") + 
  theme(plot.title = element_text(size = 18, hjust = 0), plot.title.position = "plot") +
  theme(plot.margin = margin(b=5, t=5))
f1e <- f1e + ggtitle("e\n") + 
  theme(plot.title = element_text(size = 18, hjust = 0), plot.title.position = "plot") +
  theme(plot.margin = margin(b=5, t=5))
f1f <- as.ggplot(f1f) + ggtitle("f\n") + 
  theme(plot.title = element_text(size = 18, hjust = 0), plot.title.position = "plot")


f1 <- as.ggplot(
  grid.arrange(grobs = list(f1a, f1b, f1c, f1d, f1e, f1f),
               row = 6,
               heights = c(3,3,4,8,8,5))
  ) + theme(panel.background = element_rect(fill = 'white', color = 'white'))


f1_red <- as.ggplot(
  grid.arrange(grobs = list(f1a, f1b, f1c, f1d, f1e),
               row = 5,
               heights = c(3, 3, 5, 10, 10))
) + theme(panel.background = element_rect(fill = 'white', color = 'white'))

f1_red_v2 <- as.ggplot(
grid.arrange(grobs = list(f1a, f1b, f1c, f1d, f1e),
             row = 3,
             cols = 3,
             heights = c(3, 3, 5),
             layout_matrix = rbind(c(1,4,5),
                                   c(2,4,5),
                                   c(3,4,5))
             )) + theme(panel.background = element_rect(fill = 'white', color = 'white'))

f1_red_square <- as.ggplot(grid.arrange(grobs = list(f1a, f1b, f1c, f1d, f1e),
                                        row = 3,
                                        cols = 2,
                                        heights = c(1,1,3),
                                        layout_matrix = rbind(c(1,3),
                                                              c(2,3),
                                                              c(4,5))
)) + theme(panel.background = element_rect(fill = 'white', color = 'white'))

setwd(dir_out)
ggsave(filename = "Figure1.png", plot = f1, device = "png", width = 2048, height = 8192, units = "px")
ggsave(filename = "Figure1_reduced.png", plot = f1_red, device = "png", width = 2048, height = 5760, units = "px")
ggsave(filename = "Figure1_reduced_long.png", plot = f1_red_v2, device = "png", height = 2048, width = 5760, units = "px")
ggsave(filename = "Figure1_reduced_squ.png", plot = f1_red_square, device = "png", height = 3072, width = 4096, units = "px", scale = 0.8, limitsize = FALSE)
