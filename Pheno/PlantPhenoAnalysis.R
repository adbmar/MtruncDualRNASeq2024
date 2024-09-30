################################################################################
# Author: Addison Dale Buxton-Martin
# 
# Description: This script runs the exploratory phenotype analysis for ___
# This analysis will import phenotypic data included in its encasing directory
# and perform analyses as outlined by ___
#
# See other script for cleaned up and to the point analyses and plotting
#
################################################################################


###################################
### Loading requisite libraries ###
###################################

# loading requisite libraries
# Define a function to check, install, and load a library
check_and_load_libraries <- function(library_names) {
  # Ensure library_names is a vector
  if (!is.character(library_names)) {
    stop("The input must be a character vector of package names.")
  }
  
  for (library_name in library_names) {
    if (!requireNamespace(library_name, quietly = TRUE)) {
      # Try installing the package using base install.packages
      tryCatch({
        install.packages(library_name, dependencies = TRUE)
      }, error = function(e) {
        # If base install.packages fails, try BiocManager::install
        if (requireNamespace("BiocManager", quietly = TRUE)) {
          BiocManager::install(library_name, dependencies = TRUE)
        } else {
          stop("BiocManager is required but not installed.")
        }
      })
    }
    # Load the library
    library(library_name, character.only = TRUE)
  }
}

check_and_load_libraries(c(
  "ggplot2", "ggplotify", "tidyverse", "lme4", "lmerTest",
  "ggeffects", "gridExtra", "grid", "broom.mixed", "ggpubr",
  "DHARMa", "ggtext", "corrplot", "car", "cowplot"))



################################################################################
############################## DIRECTORY HANDLING ##############################
################################################################################
dir_phenomain <- if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  dirname(rstudioapi::getSourceEditorContext()$path)
} else {
  dirname(normalizePath(sys.frame(1)$ofile))
}

#manual back up for setting main directory)
################################################################################
# Edit the file path of dir_main to be the folder containing this script.
# The data should be in this main directory as specified in the above Data
# Requirements section in the heading of this script
# dir_phenomain <- file.path("~", "Desktop", "DGE_12")
################################################################################

dir_phenoexplore <- file.path(dir_phenomain, "Exploratory_analyses")
dir.create(dir_phenoexplore, showWarnings = FALSE)

dir_phenoout <- file.path(dir_phenomain, "out")
dir.create(dir_phenoout, showWarnings = FALSE)


####################
### Loading Data ###
####################
# note data should be in the same directory as this script

phenodata <- read.csv(file.path(dir_phenomain, "Pheno_data.csv"))
scandata <- read.csv(file.path(dir_phenomain, "RhizoVisionExplorer_data.csv"))

# Harmonizing formatting between files
scandata$File.Name <- gsub(".png", "", scandata$File.Name)
scandata$File.Name <- gsub("id", "id:", scandata$File.Name)
scandata <- scandata %>% drop_na() %>% rename("id" = "File.Name")

# Joining data and organizing table
data <- left_join(scandata, phenodata, by = "id") %>% rowwise() %>%
  mutate(Block = interaction(Block.Side, Block.Shelf, Block.Position)) %>%
  mutate(Total.Nods = as.numeric(Nods) + as.numeric(Nods.Removed)) %>%
  select(id, Nema, Rhizo, Block, Galls, Total.Nods,
         Below.ground.biomass, Above.ground.biomass, Volume.mm3,
         Number.of.Root.Tips, Total.Root.Length.mm,
         Network.Area.mm2, Average.Diameter.mm,
         Median.Diameter.mm, Maximum.Diameter.mm, Perimeter.mm)

# Setting columns as factor or numeric
factor_cols <- c("id", "Nema", "Rhizo", "Block")
numeric_cols <- setdiff(names(data), factor_cols)
data[factor_cols] <-  lapply(data[factor_cols], as.factor)
data[numeric_cols] <- lapply(data[numeric_cols], as.numeric)

# Creating alternative dataframes
data_scaled <- sapply(data[numeric_cols], scale)

####################
### Correlations ###
####################

#Plotting all correlations without the gall column
try(dev.off())
try(dev.off())
setwd(dir_phenoout)
jpeg(height = 800, width = 800, file = "Correlation_all.jpeg", type = "cairo")
cor_all <- 
  corrplot.mixed(cor(data[6:ncol(data)]), 
                 upper = 'ellipse', order = 'hclust', tl.pos = "lt",
                 bg = "gray", lower.col = COL2("RdBu", 6))
try(dev.off())
try(dev.off())




#NOTE: Network area, Total root length, and perimeter were all very highly correlated
#NOTE: And Volume, number of root tips were close behind. Together all five of these
# root traits had correlations with one another of 0.9 or higher.

# As expected, above ground biomass was highly correlated with root volume (0.8)

#Plotting all correlations for all nematode positive treatments
setwd(dir_phenoout)
jpeg(height = 800, width = 800, file = "Correlation_gall.jpeg", type = "cairo")
cor_gall <- 
  corrplot.mixed(cor((data %>% filter(Nema == "+"))[5:ncol(data)]), 
                 upper = 'ellipse', order = 'hclust', tl.pos = "lt",
                 bg = "gray", lower.col = COL2("RdBu", 6))
try(dev.off())
try(dev.off())

##########################
## Statistical Analyses ##
##########################
options(contrasts=c("contr.sum", "contr.poly"))

### Defining a function for compact and easy analysis
check_residuals <- function(model, return_residuals = FALSE, stop_if_nonnormal = TRUE){
  #This function will check that the residuals of a model are normally distributed
  #Intended to be used at the beginning of the anlaysis wrapper function
  #Default will stop analysis if there are any important assumptions of an ANOVA test that are not met
  residual_test <- testResiduals(model, plot = F)
  residual_test_uniformity <- residual_test$uniformity
  residual_test_dispersion <- residual_test$dispersion
  residual_test_outlier <- residual_test$outliers
  if(residual_test_uniformity$p.value < 0.05){
    print("! ! ! ! ! ! WARNING: Residuals after fitting to model did not pass the KS test for normality!")
    if (stop_if_nonnormal) {stop("Model assumptions not met")}
    }
  if(residual_test_dispersion$p.value < 0.05){
    print("! ! ! ! ! ! WARNING: Residuals after fitting to model are not homoscedastic!")
    if (stop_if_nonnormal) {stop("Model assumptions not met")}
  }
  if(residual_test_outlier$p.value < 0.05){
    print("! ! ! ! ! ! WARNING: Residuals after fitting to model had outliers!")
    if (stop_if_nonnormal) {stop("Model assumptions not met")}
    }
  if(return_residuals == TRUE){return(residual_test)}
}

analysis <- function(formula){
  #Function to run a recurring statistical analysis that I will use
  #Note that the function automatically adds in the random variable of the block design
  model_formula <- paste(as.character(formula), "+ (1|Block)")
  model <- lmer(data = data, formula = as.formula(model_formula))
  check_residuals(model)
  print("Summary:")
  print(summary(model))
  print("Anova")
  Anova(model, type = 3)
}

analysis_gallonly <- function(formula){
  #variation of the baoveabove function that runs the analysis on just the gall samples
  model_formula <- paste(as.character(formula), "+ (1|Block)")
  model <- lmer(data = data %>% filter(Nema == "+"), formula = as.formula(model_formula))
  check_residuals(model)
  print("Summary:")
  print(summary(model))
  print("Anova")
  Anova(model, type = 3)
}

analysis_gallaszero <- function(formula){
  #Variation of the above function that runs the analysis on data where uninfected hosts
  #are treated as zero as opposed to NA
  model_formula <- paste(as.character(formula), "+ (1|Block)")
  model <- lmer(data = data %>% mutate(Galls = if_else(is.na(Galls), 0, Galls)),
                formula = as.formula(model_formula))
  check_residuals(model)
  print("Summary:")
  print(summary(model))
  print("Anova")
  Anova(model, type = 3)
}

analysis("Median.Diameter.mm ~ Nema*Rhizo") #Nema but marginally rhizo, infected are higher
analysis("Average.Diameter.mm ~ Nema*Rhizo") #Nema but not rhizo, infected are higher
analysis("Perimeter.mm ~ Nema*Rhizo") #Nema but not rhizo, infeted are lower
analysis("Total.Root.Length.mm ~ Nema*Rhizo") #Nema but not rhizo, infected are lower
analysis("Number.of.Root.Tips ~ Nema*Rhizo") #Nema but not rhizo, infected are lower

analysis("Maximum.Diameter.mm ~ Nema*Rhizo") #No effect
analysis("Network.Area.mm2 ~ Nema*Rhizo") #No effect

analysis("Above.ground.biomass ~ Total.Nods") #No effect
analysis("Above.ground.biomass ~ Rhizo * Total.Nods") #No effect
analysis("Above.ground.biomass ~ Rhizo") #No effect
analysis("Above.ground.biomass ~ Nema") #No effect
analysis("Above.ground.biomass ~ Rhizo*Nema") #No effect
analysis("Above.ground.biomass ~ Nema + Below.ground.biomass") #No effect
analysis("Above.ground.biomass ~ Galls") #No effect
analysis("Above.ground.biomass ~ Galls + Below.ground.biomass") #No effect

analysis("Total.Nods ~ Below.ground.biomass")
analysis("Galls ~ Below.ground.biomass")
analysis("Total.Nods ~ Galls")
analysis("Total.Nods ~ Galls + Below.ground.biomass")
analysis("Galls ~ Total.Nods + Below.ground.biomass")
analysis("Galls ~ Rhizo")
analysis("Galls ~ Rhizo + Below.ground.biomass")

analysis("Total.Nods ~ Nema * Rhizo")
analysis("Total.Nods ~ Rhizo")
analysis("Total.Nods ~ Nema")
analysis("Total.Nods ~ Nema * Rhizo + Below.ground.biomass")

################
### Plotting ###
################

# Setting colors
WSM1022_color <- "steelblue1"
USDA1021_color <- "royalblue4"
Infected_color <- "salmon4"
Uninfected_color <- "olivedrab3"

# Organizing data into long format
data_long <- data %>% mutate(Nema = if_else(Nema == "-", "Uninfected", "Infected")) %>% 
  select(-all_of(c("id"))) %>%
  rename("Nodule count" = Total.Nods,
         "Gall count" = Galls,
         "Above ground biomass (mg)" = Above.ground.biomass,
         "Average diameter (mm)" = Average.Diameter.mm,
         "Maximum Diameter (mm)" = Maximum.Diameter.mm,
         "Median Diameter (mm)" = Median.Diameter.mm,
         "Network Area (mm^2)" = Network.Area.mm2,
         "Root Tips" = Number.of.Root.Tips,
         "Perimeter (mm)" = Perimeter.mm,
         "Total Root Length (mm)" = Total.Root.Length.mm,
         "Below Ground Biomass (mg)" = Below.ground.biomass,
         "Volume (mm^3)" = Volume.mm3) %>%
pivot_longer(cols = -c("Nema", "Rhizo", "Block"), names_to = "trait")
data_long$value <- as.numeric(data_long$value)


# Getting a legend to use for future plots
legend1 <- cowplot::get_plot_component( ggplot(data_long) + 
                                          aes(x = Rhizo, y = value, fill = Nema) + 
                                          geom_point(shape = 21) + 
                                          scale_fill_manual(name = "Nematode infection status",
                                                            values = c("Infected" = Infected_color,
                                                                       "Uninfected" = Uninfected_color)) + 
                                          theme(legend.position = "right"),
                                        "guide-box-right", return_all = TRUE)

legend2 <- cowplot::get_plot_component(ggplot(data_long) + 
                                         aes(x = Nema, y = value, fill = Rhizo) + 
                                         geom_point(shape = 21) + 
                                         scale_fill_manual(name = "Rhizobia strain",
                                                           values = c("USDA1021" = USDA1021_color,
                                                                      "WSM1022" = WSM1022_color)) + 
                                         theme(legend.position = "right"),
                                       "guide-box-right", return_all = TRUE)
legend <- as.ggplot(ggarrange(legend1, legend2, ncol = 1, nrow = 2)) + 
  theme(plot.margin = unit(c(0.75, 0.75, 0.75, 0.75), "cm"),panel.background = element_rect(color = "black"))
legend 
# Function to apply to graphs for consistent plotting settings
p_add <- function(p, ylabel){
  # Defining settings that are applied across both plots
  p <- p +
    geom_point(shape = 21, alpha = 0.35, size = 1.5,
               position = position_jitterdodge(
                 dodge.width = 0.5,
                 jitter.width = 0.1)) +
    stat_summary(aes(group = NA),
                 fun = mean,
                 alpha = 0.75,
                 geom = "point",
                 color = "white", stroke = 1.2,
                 shape = 21, fill = "lightgray",
                 size = 2.7) + 
    stat_summary(aes(group = NA),
                 fun = mean,
                 geom = "point",
                 color = "black",
                 size = 2.5) + 
    stat_summary(aes(group = NA),
                 fun.data = mean_cl_normal,
                 geom = "errorbar",
                 width = 0.25,
                 fun.args = list(mult = 1)) +
    theme_classic() + 
    guides(fill=guide_legend(nrow=2, 
                             byrow=TRUE,
                             override.aes=list(shape = 21))) +
    theme(legend.key.height = unit(0.1, "cm")) +
    ylab(ylabel) +
    scale_fill_manual(values = c("Infected" = Infected_color,
                                 "Uninfected" = Uninfected_color,
                                 "USDA1021" = USDA1021_color,
                                 "WSM1022" = WSM1022_color)) +
    coord_cartesian(clip = "off", ) + theme(plot.margin = margin(0.2,0.2,0.2,0.2, "cm"))
  return(p)
}

# Wrapper function for plotting phenotypic data
pheno_plotting <- function(plotted_trait, r_sig = NA, n_sig = NA, get_legend = FALSE, rhizo_only = FALSE){
  plotted_trait_string <- parse(text = plotted_trait %>% str_replace_all(" ", "~"))
  
  df <- data_long %>% filter(trait == plotted_trait)
  if(rhizo_only){df <- df %>% filter(Nema == "Infected")}
  sig_bar_y <- (max(df$value, na.rm = TRUE) - min(df$value, na.rm = TRUE))*0.05 + max(df$value, na.rm = TRUE)
  
  pr <- ggplot(df) + 
    aes(x = Rhizo, y = value,
        group = Nema, fill = Nema) +
    xlab("Rhizobia strain")
  
  pr <- p_add(pr, plotted_trait_string)
  
  if(rhizo_only == FALSE){
    pn <- ggplot(df) + 
      aes(x = Nema, y = value,
          group = Rhizo, fill = Rhizo) +
      xlab("Nematode infection status")
    pn <- p_add(pn, plotted_trait_string)
  } else {
    pn <- ggplot() + theme(panel.background = element_rect(fill = "white"))
  }
  
  
  if(!is.na(n_sig) & !rhizo_only){
    pn <- pn + geom_bracket(inherit.aes = FALSE,
                            xmin = unique(factor(df$Nema))[1],
                            xmax = unique(factor(df$Nema))[2],
                            y.position = sig_bar_y,
                            label = n_sig)
  }
  if(!is.na(r_sig)){
    pr <- pr + geom_bracket(inherit.aes = FALSE,
                            xmin = unique(factor(df$Rhizo))[1],
                            xmax = unique(factor(df$Rhizo))[2],
                            y.position = sig_bar_y,
                            label = r_sig)
  }
  
  if ((!is.na(n_sig) |  !is.na(r_sig) & !rhizo_only)){
    ymax <- ((max(df$value, na.rm = TRUE) - min(df$value, na.rm = TRUE)))*0.125 + max(df$value, na.rm = TRUE)
    ymin <- min(df$value, na.rm = TRUE)
    pr <- pr + ylim(c(ymin, ymax))
    pn <- pn + ylim(c(ymin, ymax))
  }
  title <- ggplot() +
    theme(panel.background = element_blank()) +
    labs(title = plotted_trait_string) + 
    theme(plot.title = element_text(size=16, hjust = 0.5, vjust = 1),
          plot.background = element_rect(fill = "lightgray"))
  if (get_legend == FALSE){
    p <- grid.arrange(title, 
                      pr + theme(legend.position = "none"),
                      pn + theme(legend.position = "none"),
                      nrow = 2, ncol = 2,
                      layout_matrix = rbind(c(1,1), c(2,3)),
                      heights = c(1, 10))
    return(as.ggplot(p) + theme(plot.margin = margin(0.2,0.2,0.2,0.2, "cm")))
  } else {return(
    cowplot::get_plot_component(pr + 
                                  theme(legend.position = "right"), "guide-box-right", return_all = TRUE)
  )}
}


# Root architecture traits
paved <- pheno_plotting("Average diameter (mm)", n_sig = "p = 0.008")
pbra <- pheno_plotting("Branch Points")
pmaxd <- pheno_plotting("Maximum Diameter (mm)")
pmedd <- pheno_plotting("Median Diameter (mm)", n_sig = "p = 0.039")
pneta <- pheno_plotting("Network Area (mm^2)")
pperi <- pheno_plotting("Perimeter (mm)", n_sig = "p = 0.018")
prtip <- pheno_plotting("Root Tips", n_sig = "p = 0.012")
prlen <- pheno_plotting("Total Root Length (mm)", n_sig = "p = 0.031")
pbgb <- pheno_plotting("Below Ground Biomass (mg)")
pvol <- pheno_plotting("Volume (mm^3)")

sup_plot <- ggarrange(paved, pbfreq, pbra, pmaxd, pmedd, pneta, pperi, prtip, prlen, pvol, legend,
                      ncol = 3, nrow = 4) 
ggsave(plot = sup_plot, 
       filename = file.path(dir_phenoout, "SupFig_RootArch.jpeg"),
       units = "px", 
       width = 8600, height = 6400, scale = 0.5)


bgb_vol_cor <- cor(data %>% select(`Below Ground Biomass (mg)`, `Volume (mm^3)`) %>% drop_na())
ggplot(data) + aes(x = `Volume (mm^3)`, y = `Below Ground Biomass (mg)`) +
  geom_point() + 
  annotation("text", x = Inf, y = Inf, label = paste("r =", round(bgb_vol_cor, 3)),
             hjust = 1.1, vjust = 1.1, size = 3) +
  theme_classic()


# Nodule/Gall count data and above ground biomass data
pabg <- pheno_plotting("Above ground biomass (mg)")
pnod <- pheno_plotting("Nodule count", r_sig = "p = 0.004", n_sig = "p = 0.002")
pgall <- pheno_plotting("Gall count", r_sig = "p = 0.024", rhizo_only = TRUE)

# Plot of total nodules and above ground biomass
pavn <- ggplot(data) +
  aes(x = Total.Nods,
      y = Above.ground.biomass,
      fill = Rhizo,
      color = Nema,
      shape = Nema) +
  geom_point(shape = 21, stroke = 0.8) + theme_classic() +
  geom_smooth(method = "lm", inherit.aes = FALSE, 
              aes(x = Total.Nods, y = Above.ground.biomass)) +
  scale_color_manual(guide = "none",
                     values = c("+" = Infected_color,
                                "-" = Uninfected_color)) +
  scale_fill_manual(guide = "none",
                    values = c("USDA1021" = USDA1021_color,
                               "WSM1022" = WSM1022_color)) +
  ylab("Above ground plant biomass (mg)") + xlab("Nodule count") +
  # annotate("text", x = 90, y = 0.025, hjust = 1, label = "r = 0.488") +
  # annotate("text", x = 90, y = 0.015, hjust = 1, label = "Effect of nodule count: p < 0.001") +
  # annotate("text", x = 90, y = 0.005, hjust = 1, label = "Effect of rhizobia strain: p = 0.666")
  annotate("text", hjust = 1, vjust = 0.2, size = 2.5,
           x = max(data$Total.Nods), y = min(data$Above.ground.biomass),
           label = "r = 0.488\nEffect of nodule count: p < 0.001\nEffect of rhizobia strain: p = 0.934\nEffect of nematode infection: p = 0.645")
pavn
# Plot of galls and above ground biomass
pgvn <- ggplot(data %>% filter(Nema == "+")) +
  aes(x = Galls,
      y = Above.ground.biomass,
      fill = Rhizo,
      color = Rhizo) +
  geom_point(color = "black", shape = 21) + theme_classic() +
  # geom_smooth(method = "lm", inherit.aes = FALSE, 
  #             aes(x = Galls, y = Above.ground.biomass)) +
  scale_fill_manual(name = "Rhizobia strain",
                    values = c("USDA1021" = USDA1021_color,
                               "WSM1022" = WSM1022_color)) +
  scale_fill_manual(guide = FALSE,
                    values = c("USDA1021" = USDA1021_color,
                               "WSM1022" = WSM1022_color)) +
  ylab("Above ground plant biomass (mg)") +
  annotate("text", hjust = 1, vjust = 0.2, size = 2.5,
           x = max(data$Galls), y = min(data$Above.ground.biomass),
           label = "r = 0.020\nEffect of gall count: p = 0.130\nEffect of rhizobia strain: p = 0.266")
pgvn

fig4 <- ggarrange(pabg, pnod, pgall, pavn, pgvn, legend, ncol = 3, nrow = 2)
fig4
ggsave(plot = fig4, 
       filename = file.path(dir_phenoout, "Fig4.jpeg"),
       units = "px", 
       width = 5400, height = 3600, scale = 0.666)








# Function for plotting residuals across rhizobia strains
residual_plotting_rhizo <- function(formula, ylabel){
  model <- lmer(data = data, formula = formula)
  residuals <- resid(model)
  residual_data <- data.frame(Rhizo = data$Rhizo, Nema = data$Nema, Residuals = residuals)
  residual_data <- residual_data %>% mutate(Nema = if_else(Nema == "-", "Uninfected", "Infected"))
  p <- ggplot(residual_data) + 
    aes(x = Rhizo,
        y = Residuals, 
        fill = Nema)
  
  test <- p_add(p, ylabel) + xlab("Rhizobia strain")
  return(test)
}

# Function for plotting residuals across nematode infection status
residual_plotting_nema <- function(formula, ylabel){
  model <- lmer(data = data, formula = formula)
  residuals <- resid(model)
  residual_data <- data.frame(Rhizo = data$Rhizo, Nema = data$Nema, Residuals = residuals)
  residual_data <- residual_data %>% mutate(Nema = if_else(Nema == "-", "Uninfected", "Infected"))
  p <- ggplot(residual_data) + 
    aes(x = Nema,
        y = Residuals, 
        fill = Rhizo)
  
  test <- p_add(p, ylabel) + xlab("Nematode infection status")
  return(test)
}

# Wrapper script for plotting residuals across both rhizobia and nematode
residual_plotting <- function(formula, trait, correction, model_as_string){
  pr <- residual_plotting_rhizo(formula, "Model residuals")
  pn <- residual_plotting_nema(formula, "Model residuals")
  
  title <- ggplot() +
    theme(panel.background = element_blank()) +
    labs(title = paste(trait, "corrected by", correction),
         subtitle = paste("Model formula:", model_as_string)) + 
    theme(plot.title = element_text(size=16, hjust = 0.5, vjust = 1),
          plot.subtitle = element_text(size = 8, hjust = 0.5, vjust = 1),
          plot.background = element_rect(fill = "lightgray"))
  
  p <- grid.arrange(title, 
                    pr + theme(legend.position = "none"),
                    pn + theme(legend.position = "none"),
                    nrow = 2, ncol = 2,
                    layout_matrix = rbind(c(1,1), c(2,3)),
                    heights = c(1, 10))
  return(as.ggplot(p) + theme(plot.margin = margin(0.2,0.2,0.2,0.2, "cm")))
}

residual_plotting("Above.ground.biomass ~ Rhizo*Nema + Volume.mm3 + (1|Block)",
                  "Above ground biomass", "root volume",
                  "Above ground biomass ~ Rhizobia strain × Nematode status + Root volume")


ggplot(data) +
  aes(x = Total.Nods,
      y = Above.ground.biomass,
      fill = Rhizo,
      color = Rhizo) +
  geom_point(color = "black", shape = 21) + theme_classic() +
  geom_smooth(method = "lm", inherit.aes = FALSE, 
              aes(x = Total.Nods, y = Above.ground.biomass)) +
  scale_fill_manual(name = "Rhizobia strain",
                    values = c("USDA1021" = USDA1021_color,
                               "WSM1022" = WSM1022_color)) +
  scale_fill_manual(guide = FALSE,
                    values = c("USDA1021" = USDA1021_color,
                               "WSM1022" = WSM1022_color)) +
  ylab("Above ground biomass") + xlab("Nodule count")

ggplot(data %>% filter(Nema == "+")) +
  aes(x = Volume.mm3,
      y = Galls,
      fill = Rhizo,
      color = Rhizo) +
  geom_point(color = "black", shape = 21) + theme_classic() +
  geom_smooth(method = "lm") +
  scale_color_manual(name = "Rhizobia strain",
                     values = c("USDA1021" = USDA1021_color,
                                "WSM1022" = WSM1022_color)) +
  scale_fill_manual(guide = FALSE,
                    values = c("USDA1021" = USDA1021_color,
                               "WSM1022" = WSM1022_color)) +
  ylab("Galls") + xlab(expression(Root~volume~(mm^3)))


dev.off()





