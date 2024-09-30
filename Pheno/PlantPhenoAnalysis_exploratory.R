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
  "DHARMa", "ggtext", "corrplot", "car"))



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
data <- left_join(scandata, phenodata, by = "id") %>%
  mutate(Block = interaction(Block.Side, Block.Shelf, Block.Position)) %>%
  mutate(Total.Nods = as.numeric(Nods) + as.numeric(Nods.Removed)) %>%
  select(id, Nema, Rhizo, Block, Galls, Total.Nods,
         Above.ground.biomass, Volume.mm3,
         Number.of.Root.Tips, Number.of.Branch.Points, Total.Root.Length.mm,
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
cor_all <- corrplot.mixed(cor(data[6:ncol(data)]), 
                          upper = 'ellipse', order = 'hclust', tl.pos = "lt",
                          bg = "gray", lower.col = COL2("RdBu", 6))

#NOTE: Network area, Total root length, and perimeter were all very highly correlated
#NOTE: And Volume, number of root tips were close behind. Together all five of these
# root traits had correlations with one another of 0.9 or higher.

# As expected, above ground biomass was highly correlated with root volume (0.8)

#Plotting all correlations for all nematode positive treatments
cor_gall <- corrplot.mixed(cor((data %>% filter(Nema == "+"))[5:ncol(data)]), 
                           upper = 'ellipse', order = 'hclust', tl.pos = "lt",
                           bg = "gray", lower.col = COL2("RdBu", 6))


##########################
## Exploratory analysis ##
##########################
options(contrasts=c("contr.sum", "contr.poly"))

setwd(dir_phenoexplore)

# Function to plot raw counts
analysis <- function(data, 
                     responding_variable_col, 
                     responding_variable_name,
                     correlate = FALSE,
                     cor_variable_col = NA,
                     cor_variable_name = NA,
                     nema = TRUE,
                     custom_formula = NA) {
  
  responding_variable_col_string <- deparse(substitute(responding_variable_col))
  cor_variable_col_string <- deparse(substitute(cor_variable_col))
  
  if (correlate == FALSE) {ylabel <- responding_variable_name}
  else {ylabel <- paste("Residual of", responding_variable_name,
                        "from correlation with", cor_variable_name)}
  
  plot0 <- ggplot(data) + aes(x = if(nema == TRUE){interaction(Rhizo, Nema)}
                              else{Rhizo},
                              y = {{ responding_variable_col }}) +
    geom_boxplot(outliers = FALSE) +
    geom_point(shape = 21, fill = "dark gray", alpha = 0.8) +
    ggtitle("Raw Data") +
    theme_classic()
  
  # Plotting raw data
  if (correlate == FALSE){
    plot2 <- ggplot(data) + aes(x = if(nema == TRUE){interaction(Rhizo, Nema)}
                                else{Rhizo},
                                y = {{ responding_variable_col }}) +
      geom_boxplot(outliers = FALSE) +
      geom_point(position = position_jitter(width = 0.2),
                 alpha = 0.5, shape = 21, fill = "dark gray") +
      xlab("Treatment Group") +
      ylab(ylabel) + 
      ggtitle(paste("Raw Data for", responding_variable_name)) +
      theme_classic()
  }
  else {
    res_model <- lmer(as.formula(paste(
      responding_variable_col_string, "~", cor_variable_col_string, "+ (1 | Block)")), data = data)
    plot2 <- ggplot(data %>%
                      mutate(
                        predicted = predict(res_model),
                        residuals = {{ responding_variable_col }} - predicted)) +
      aes(x = if(nema == TRUE){interaction(Rhizo, Nema)}
          else{Rhizo},
          y = residuals) +
      geom_hline(yintercept = 0) +
      geom_boxplot(outliers = FALSE) +
      geom_point(position = position_jitter(width = 0.2),
                 alpha = 0.5, shape = 21, fill = "dark gray") +
      xlab("Treatment Group") +
      ylab(ylabel) + 
      ggtitle(paste("Residuals after correlating with", cor_variable_name,"for", responding_variable_name)) +
      theme_classic()
  }
  
  # Making model
  if (is.na(custom_formula)) {
    if(nema == TRUE){
      model_formula <- as.formula(paste(responding_variable_col_string,
                                        "~ Rhizo * Nema + (1 | Block)"))
    }
    else{
      model_formula <- as.formula(paste(responding_variable_col_string,
                                        "~ Rhizo + (1 | Block)"))
    }
    if(correlate == TRUE){
      if(nema == TRUE) {
        model_formula <- as.formula(paste(responding_variable_col_string, 
                                          "~ Rhizo * Nema +",
                                          cor_variable_col_string, 
                                          "+ (1 | Block)"))}
      else {
        model_formula <- as.formula(paste(responding_variable_col_string,
                                          "~ Rhizo +",
                                          cor_variable_col_string,
                                          "+ (1 | Block)"))
      }
    }
  }
  else{
    model_formula <- as.formula(custom_formula)
  }
  
  
  model <- lmer(model_formula, data)
  
  # Testing assumptions of model
  plot1 <- as.ggplot(function() plot(simulateResiduals(model)))
  
  
  # Plotting predicted means
  ggp <- ggpredict(model, terms = if(nema == TRUE){c("Rhizo", "Nema")}
                   else{c("Rhizo")})
  ggp_df <- as.data.frame(ggp) %>% rename(c("Rhizo" = "x", "Nema" = "group"))
  
  # Plotting predicted group means
  plot3 <- ggplot() +
    geom_linerange(data = ggp_df, 
                   aes(
                     if (nema == TRUE) {x = interaction(Rhizo, Nema)} else {
                       x = Rhizo},
                     ymax = conf.high,
                     ymin = conf.low),
                   linetype = "dashed") +
    geom_point(data = data, aes(
      if (nema == TRUE) {x = interaction(Rhizo, Nema)} else {
        x = Rhizo},
      y = {{ responding_variable_col }}),
      shape = 21, fill = "dark gray", alpha = 0.6,
      position = position_jitter(width = 0.2)) +
    geom_point(data = ggp_df, aes(
      if (nema == TRUE) {x = interaction(Rhizo, Nema)} else {
        x = Rhizo},
      y = predicted),
      size = 3, fill = "blue", shape = 21) +
    
    geom_point(size = 3) +
    theme_classic() +
    xlab("Treatment Groups") +
    ggtitle(paste("Mixed Model Predicted Group Means for", responding_variable_name, "overtop of raw data")) +
    ylab(ylabel)
  
  # Testing model with Anova
  plot2 <- as.ggplot(function() 
    grid.table(as.data.frame(
      anova(model, type = 3)) %>%
        mutate(sig = if_else(`Pr(>F)` < 0.05, "*", " ")))) + 
    ggtitle(paste("ANOVA Results for", responding_variable_name))
  
  
  
  p <- as.ggplot(ggarrange(plotlist = list(plot1, plot2, plot3), nrow = 3, heights = c(5,2,3))) + 
    ggtitle(paste0("<span style='font-size: 22pt;'>", responding_variable_name, "</font>")) +
    theme(plot.title = element_markdown())
  return(p)
  
}

analysis(data, Total.Nods, "Total Nodules")
ggsave(height = 10.5, width = 8.5, filename = "TotalNodules_unfiletered.jpeg")

analysis(data, Above.ground.biomass, "Above Ground Biomass")
ggsave(height = 10.5, width = 8.5, filename = "AboveGroundBiomass_unfiletered.jpeg")

analysis(data, Volume.mm3, "Total Root Volume")
ggsave(height = 10.5, width = 8.5, filename = "TotalRootVolume_unfiletered.jpeg")

analysis(data, Number.of.Root.Tips, "Number of Root tips")
ggsave(height = 10.5, width = 8.5, filename = "RootTips_unfiletered.jpeg")

analysis(data, Number.of.Branch.Points, "Number of Branch Points")
ggsave(height = 10.5, width = 8.5, filename = "BranchPoints_unfiletered.jpeg")

analysis(data, Total.Root.Length.mm, "Total Root Length (mm)")
ggsave(height = 10.5, width = 8.5, filename = "RootLength_unfiletered.jpeg")

analysis(data, Network.Area.mm2, "Network Area (mm2)")
ggsave(height = 10.5, width = 8.5, filename = "NetworkArea_unfiletered.jpeg")

analysis(data, Average.Diameter.mm, "Average Diameter (mm)")
ggsave(height = 10.5, width = 8.5, filename = "AverageDaimeter_unfiletered.jpeg")

analysis(data, Median.Diameter.mm, "Median Diameter (mm)")
ggsave(height = 10.5, width = 8.5, filename = "MedianDiameter_unfiletered.jpeg")

analysis(data, Maximum.Diameter.mm, "Maximum Diameter (mm)")
ggsave(height = 10.5, width = 8.5, filename = "MaxDiameter_unfiletered.jpeg")

analysis(data, Perimeter.mm, "Perimeter (mm)")
ggsave(height = 10.5, width = 8.5, filename = "Perimeter_unfiletered.jpeg")

analysis(data %>% filter(Nema == "+"), Galls, "Gall Count", nema = FALSE)
ggsave(height = 10.5, width = 8.5, filename = "GallCount_unfiletered.jpeg")

analysis(data, Total.Nods, "Total Nodules", correlate = TRUE, Above.ground.biomass, "Above ground biomass", nema = TRUE)
ggsave(height = 10.5, width = 8.5, filename = "TotalNodules_AboveGroundBiomassCorrected_unfiletered.jpeg")

analysis(data, Total.Nods, "Total Nodules", correlate = TRUE, Volume.mm3, "Total Root Volume", nema = TRUE)
ggsave(height = 10.5, width = 8.5, filename = "TotalNodules_RootVolumeCorrected_unfiletered.jpeg")

analysis(data %>% filter(Nema == "+"), Galls, "Galls", correlate = TRUE, Above.ground.biomass, "Above ground biomass", nema = FALSE)
ggsave(height = 10.5, width = 8.5, filename = "Galls_AboveGroundBiomassCorrected_unfiletered.jpeg")

analysis(data %>% filter(Nema == "+"), Galls, "Galls", correlate = TRUE, Volume.mm3, "Total Root Volume", nema = FALSE)
ggsave(height = 10.5, width = 8.5, filename = "Galls_RootVolumeCorrected_unfiletered.jpeg")

analysis(data %>% filter(Nema == "+"), Above.ground.biomass, "Above ground biomass", correlate = TRUE, Galls, "Galls", nema = FALSE)
ggsave(height = 10.5, width = 8.5, filename = "AboveGroundBiomass_GallCorrected_unfiletered.jpeg")

analysis(data, Above.ground.biomass, "Above ground biomass", correlate = TRUE, Total.Nods, "Total Nodules", nema = TRUE)
ggsave(height = 10.5, width = 8.5, filename = "AboveGroundBiomass_NoduleCorrected_unfiletered.jpeg")

analysis(data, Above.ground.biomass, "Above ground biomass", custom_formula = "Above.ground.biomass ~ Total.Nods * Rhizo * Nema + (1 | Block)")
ggsave(height = 10.5, width = 8.5, filename = "Three_way_interaction_unfiletered.jpeg")

analysis(data, Above.ground.biomass, "Above ground biomass", custom_formula = "Above.ground.biomass ~ Total.Nods * Rhizo + Total.Nods * Nema  + (1 | Block)")
ggsave(height = 10.5, width = 8.5, filename = "Three_way_workaround_unfiletered.jpeg")

analysis(data, Above.ground.biomass, "Above ground biomass", custom_formula = "Above.ground.biomass ~ Total.Nods * Rhizo + Total.Nods * Nema  + (1 | Block)")
ggsave(height = 10.5, width = 8.5, filename = "Three_way_workaround_unfiletered.jpeg")


options(contrasts=c("contr.sum", "contr.poly"))

### Defining a function for compact and easy analysis
check_residuals <- function(model, return_residuals = FALSE){
  residual_test <- testResiduals(model, plot = F)
  residual_test_uniformity <- residual_test$uniformity
  residual_test_dispersion <- residual_test$dispersion
  residual_test_outlier <- residual_test$outliers
  if(residual_test_uniformity$p.value < 0.05){print("! ! ! ! ! ! WARNING: Residuals after fitting to model did not pass the KS test for normality!")}
  if(residual_test_dispersion$p.value < 0.05){print("! ! ! ! ! ! WARNING: Residuals after fitting to model are not homoscedastic!")}
  if(residual_test_outlier$p.value < 0.05){print("! ! ! ! ! ! WARNING: Residuals after fitting to model had outliers!")}
  if(residual_test_uniformity$p.value < 0.05 |
     residual_test_dispersion$p.value < 0.05 |
     residual_test_outlier$p.value < 0.05) {print("test")
    #print(residual_test)
  }
  if(return_residuals == TRUE){return(residual_test)}
}

analysis <- function(formula){
  model_formula <- paste(as.character(formula), "+ (1|Block)")
  model <- lmer(data = data, formula = as.formula(model_formula))
  check_residuals(model)
  print("Summary:")
  print(summary(model))
  print("Anova")
  Anova(model, type = 3)
}

analysis_gallonly <- function(formula){
  model_formula <- paste(as.character(formula), "+ (1|Block)")
  model <- lmer(data = data %>% filter(Nema == "+"), formula = as.formula(model_formula))
  check_residuals(model)
  print("Summary:")
  print(summary(model))
  print("Anova")
  Anova(model, type = 3)
}

analysis_gallaszero <- function(formula){
  model_formula <- paste(as.character(formula), "+ (1|Block)")
  model <- lmer(data = data %>% mutate(Galls = if_else(is.na(Galls), 0, Galls)),
                formula = as.formula(model_formula))
  check_residuals(model)
  print("Summary:")
  print(summary(model))
  print("Anova")
  Anova(model, type = 3)
}




analysis("Volume.mm3 ~ Nema*Rhizo") #No
analysis("Total.Root.Length.mm ~ Nema*Rhizo")
analysis("Network.Area.mm2 ~ Nema*Rhizo") #No
analysis("Perimeter.mm ~ Nema*Rhizo")
analysis("Number.of.Root.Tips ~ Nema*Rhizo")
analysis("Number.of.Branch.Points ~ Nema*Rhizo") #No
analysis("Average.Diameter.mm ~ Nema*Rhizo")
analysis("Median.Diameter.mm ~ Nema*Rhizo")
analysis("Maximum.Diameter.mm ~ Nema*Rhizo") #No


analysis_gallaszero("Volume.mm3 ~ Nema + Galls")
analysis_gallaszero("Total.Root.Length.mm ~ Nema + Galls")
analysis_gallaszero("Network.Area.mm2 ~ Nema + Galls")
analysis_gallaszero("Perimeter.mm ~ Nema + Galls")
analysis_gallaszero("Number.of.Root.Tips ~ Nema + Galls")
analysis_gallaszero("Number.of.Branch.Points ~ Nema + Galls") #No
analysis_gallaszero("Average.Diameter.mm ~ Nema + Galls") #No
analysis_gallaszero("Median.Diameter.mm ~ Nema + Galls") #No
analysis_gallaszero("Maximum.Diameter.mm ~ Nema + Galls")

#Severity of nematode infection seems to have a positive relationship with all these traits
#including the ones that were decreased by nematode infection
analysis("Volume.mm3 ~ Galls")
analysis("Total.Root.Length.mm ~ Galls")
analysis("Network.Area.mm2 ~ Galls")
analysis("Perimeter.mm ~ Galls")
analysis("Number.of.Root.Tips ~ Galls")
analysis("Number.of.Branch.Points ~ Galls") #No
analysis("Average.Diameter.mm ~ Galls") #No
analysis("Median.Diameter.mm ~ Galls") #No
analysis("Maximum.Diameter.mm ~ Galls")

#however this relationship appears to disapear when we correct for volume
#suggesting that perhaps that trend is simply driven by the effect
#of root volume
analysis("Volume.mm3 ~ Galls")
analysis("Total.Root.Length.mm ~ Galls + Volume.mm3") #No
analysis("Network.Area.mm2 ~ Galls + Volume.mm3") #No
analysis("Perimeter.mm ~ Galls + Volume.mm3") #No
analysis("Number.of.Root.Tips ~ Galls + Volume.mm3") #No
analysis("Number.of.Branch.Points ~ Galls + Volume.mm3") #No
analysis("Average.Diameter.mm ~ Galls + Volume.mm3") #No
analysis("Median.Diameter.mm ~ Galls + Volume.mm3") #No
analysis("Maximum.Diameter.mm ~ Galls + Volume.mm3") #No






# Does above ground biomass correlate with below ground volume?
analysis("Above.ground.biomass ~ Volume.mm3")
ggplot(data, aes(x = Volume.mm3, y = Above.ground.biomass)) + geom_point() + theme_classic() + 
  geom_smooth(method='lm')





# Do plants with bigger root systems have more galls?
analysis("Galls ~ Volume.mm3")
ggplot(data, aes(x = Volume.mm3, y = Galls)) + geom_point() + theme_classic() + 
  geom_smooth(method='lm')
# Hosts with larger root systems had more galls

# Do plants with bigger root systems have more nodules?
analysis("Total.Nods ~ Volume.mm3")
ggplot(data, aes(x = Volume.mm3, y = Total.Nods)) + geom_point() + theme_classic() + 
  geom_smooth(method='lm')
# Hosts with larger root systems had more nodules

# Do plants with bigger root systems have more total symbionts?
ggplot(data, aes(x = Volume.mm3, y = Total.Nods)) + geom_point() + theme_classic() + 
  geom_smooth(method='lm')


# Together This suggests that nematodes don't add more total symbionts
# and that symbionts are competing for space/resources on the root. This suggests that
# correcting for root volume is necessary to get at the symbiont counts. Additionally,
# because symbiont count is made up of both gall and nodules, and because non-nematode
# treatments will always have 0 galls, testing with full symbiont counts may not be
# appropriate. Better to treat galls and nodules separately.



# Does rhizobia strain impact gall count?
analysis("Galls ~ Rhizo") #Yes p=0.02377
analysis("Galls ~ Rhizo + Volume.mm3") #Yes p=0.01531

# Does the number of nodules impact number of galls?
analysis("Galls ~ Total.Nods + Volume.mm3") #marginally p=0.0508

# Does the type of rhizobia within a quantity of nodules impact number of galls?
analysis("Galls ~ Rhizo + Total.Nods + Volume.mm3") #Yes p=0.001003

# Is the effect of the different rhizobia strains on the number of galls different?
analysis("Galls ~ Rhizo*Total.Nods + Volume.mm3") #No

# Since we loose a good correlation with Volume as we add factors, let's go back and confirm
# the patterns above with a simpler model
analysis("Galls~Rhizo*Total.Nods") #Confirmed



# Does nematode status impact nodule count? (Note we want to keep Nema status in there since galls are potentially competing for space with nodules)
analysis("Total.Nods ~ Rhizo + Nema + Volume.mm3")
# Rhizobia strain: YES p=0.0012760
# Nematode status: YES p=0.0005157
ggplot(data) + aes(x = Volume.mm3, y = Total.Nods, fill = Rhizo) + geom_point(shape=21) + theme_classic()

# Is the effect of nematode status on nodule count different by strain?
analysis("Total.Nods ~ Rhizo*Nema + Volume.mm3") #No

# We loose a good correlation with Volume as we add factors, let's go back and confirm
# the patterns above with a simpler model
analysis("Total.Nods ~ Rhizo*Nema") #Confirmed

# Does severity of nematode infection matter for total nodules?
analysis("Total.Nods ~ Galls + Volume.mm3") #Yes p=0.04946

# Does severity of nematode infection and rhizobia strain matter for total nodules?
analysis("Total.Nods ~ Galls + Rhizo + Volume.mm3") #Yes
# Galls p=0.01442
# Rhizo p=0.00259
ggplot(data = data, aes(x = Galls, y = Total.Nods, fill = Rhizo)) + geom_point(shape = 21) + theme_classic()

# Does severity of nematode and strain matter for total nodules?
analysis("Total.Nods ~ Galls*Rhizo + Volume.mm3") #No


# Does infection with nematodes effect host fitness?
analysis("Above.ground.biomass ~ Nema") #No
analysis("Above.ground.biomass ~ Nema*Rhizo") #No
analysis("Above.ground.biomass ~ Nema*Rhizo + Volume.mm3") #No

# Does severity of either interaction effect host fitness?
analysis("Above.ground.biomass ~ Galls*Total.Nods + Volume.mm3") #No
analysis("Above.ground.biomass ~ Nema*Rhizo") #No
analysis("Above.ground.biomass ~ Nema*Rhizo + Volume.mm3") #No
analysis("Above.ground.biomass ~ Galls + Total.Nods") #No
analysis("Above.ground.biomass ~ Galls*Rhizo") #No
analysis("Above.ground.biomass ~ Galls*Rhizo + Volume.mm3") #No

analysis("Above.ground.biomass ~ Galls*Total.Nods") #No







########################
## Summative analysis ##
########################
# Nematode infection status but not rhizobia strain impacted several root architecture traits but not rhizobia strain
analysis("Volume.mm3 ~ Nema*Rhizo") #No
analysis("Total.Root.Length.mm ~ Nema*Rhizo")
analysis("Network.Area.mm2 ~ Nema*Rhizo") #No
analysis("Perimeter.mm ~ Nema*Rhizo")
analysis("Number.of.Root.Tips ~ Nema*Rhizo")
analysis("Number.of.Branch.Points ~ Nema*Rhizo") #No
analysis("Average.Diameter.mm ~ Nema*Rhizo")
analysis("Median.Diameter.mm ~ Nema*Rhizo")
analysis("Maximum.Diameter.mm ~ Nema*Rhizo") #No


# Hosts with larger root systems were more fit
# This is in line with expectations
analysis("Above.ground.biomass ~ Volume.mm3")
ggplot(data, aes(x = Volume.mm3, y = Above.ground.biomass)) + geom_point() + theme_classic() + 
  geom_smooth(method='lm')

# Hosts with larger root systems have more symbionts
# This is in line with expectations
analysis("Total.Nods ~ Volume.mm3")
analysis("Galls ~ Volume.mm3")
ggplot(data, aes(x = Volume.mm3, y = Total.Nods)) + geom_point() + theme_classic() + 
  geom_smooth(method='lm')
ggplot(data, aes(x = Volume.mm3, y = Galls)) + geom_point() + theme_classic() + 
  geom_smooth(method='lm')

# Nematode status does not impact total symbionts but does impact nodule count,
# suggesting that symbionts are competing for the same real estate on roots
analysis("Total.Nods ~ Nema")

# Plants with more nodules had more galls and plants with more galls had more nodules
analysis("Galls ~ Total.Nods") 
analysis("Total.Nods ~ Galls") 
analysis("Galls ~ Total.Nods + Volume.mm3") 
analysis("Total.Nods ~ Galls + Volume.mm3") 

#Rhizobia strain impacts the severity of nematode infection, with USDA1021 and
# uninfected plants having more galls. This holds true even when corrected for by root volume
analysis("Galls ~ Rhizo") 
analysis("Galls ~ Rhizo*Total.Nods + Volume.mm3") 
analysis("Galls ~ Rhizo*Total.Nods") 
ggplot(data, aes(x = Total.Nods, y = Galls, fill = Rhizo)) + geom_point(shape = 21) + theme_classic()
ggplot(data, aes(x = Rhizo, y = Galls)) + geom_point() + theme_classic()

# Number of nodules is impacted by both rhizobia strain and nematode status
# In addition, infection severity also impacts nodule counts
analysis("Total.Nods ~ Rhizo*Nema + Volume.mm3")
analysis("Total.Nods ~ Rhizo*Nema")
analysis("Total.Nods ~ Rhizo*Galls + Volume.mm3")
analysis("Total.Nods ~ Rhizo*Galls")
ggplot(data, aes(y = Total.Nods, x=Galls, fill = Rhizo)) + geom_point(shape = 21)

# Infection status and rhizobia strain do not impact above ground biomass
analysis("Above.ground.biomass ~ Nema*Rhizo")

# However, number of nodules does
analysis("Above.ground.biomass ~ Galls + Total.Nods")
analysis("Above.ground.biomass ~ Galls*Total.Nods")
analysis("Above.ground.biomass ~ Galls + Total.Nods + Volume.mm3")
analysis("Above.ground.biomass ~ Galls + Volume.mm3")

# But strain doesn't alter that effect of total nodules
analysis("Above.ground.biomass ~ Galls + Total.Nods*Rhizo")
analysis("Above.ground.biomass ~ Galls + Total.Nods*Rhizo + Volume.mm3*Rhizo")
analysis("Above.ground.biomass ~ Total.Nods*Rhizo")
analysis("Above.ground.biomass ~ Total.Nods*Rhizo + Volume.mm3*Rhizo")
analysis("Above.ground.biomass ~ Total.Nods*Rhizo*Nema")
analysis("Above.ground.biomass ~ Total.Nods*Rhizo + Nema*Rhizo")





