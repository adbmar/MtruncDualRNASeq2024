library("RColorBrewer")
library("pheatmap")
library("PoiClaClu")

### Dataset quality check ###
dataset_quality_check_plotting <- function(dds, pdf_title, nema_factor = FALSE){
  #Function to plot dataset quality metrics to a pdf in the current working directory
  #Meant for DESeq dataset objects and to assess their quality
  #File and directory handling
  setwd(dir_out)
  print(paste("Dataset quality check: Printing dataset quality check to", getwd(), pdf_title))
  pdf_title <- paste0(pdf_title, ".pdf")
  pdf(file=pdf_title)
  
  #Prepping colors and data
  color_scheme <- "Spectral"
  print("Dataset quality check: Transforming data into log data")
  rld <- rlog(dds)
  
  #Euclidian sample distances
  print("Dataset quality check: Preparing euclidian sample distance heat map")
  sampleDists <- dist( t( assay(rld) ) )
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(rld$tissue, rld$nema, rld$rhizo, sep="_")
  colors <- colorRampPalette( rev(brewer.pal(9, color_scheme)) )(255)
  df <- as.data.frame(colData(rld)[ , c("tissue","rhizo","nema")])
  pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors, annotation_col=df, main = "Euclidian sample distances")
  
  #Poisson sample distances
  print("Dataset quality check: Preparing poisson sample distance heat map")
  poisd <- PoissonDistance( t( counts(dds) ) )
  samplePoisDistMatrix <- as.matrix(poisd$dd)
  rownames(samplePoisDistMatrix) <- paste(rld$tissue, rld$nema, rld$rhizo, sep="_")
  colnames(samplePoisDistMatrix) <- NULL
  pheatmap(samplePoisDistMatrix, clustering_distance_rows=poisd$dd, clustering_distance_cols=poisd$dd, col=colors, main = "Poisson sample distances")
  
  #PCA plot
  print("Dataset quality check: Preparing PCA plots")
  PCAdata <- plotPCA(rld, intgroup = c("tissue", "rhizo", "nema"), returnData=TRUE)
  pcaplot <- ggplot(PCAdata) +
    geom_hline(yintercept = 0, color = "dark gray") + geom_vline(xintercept = 0, color = "dark gray") +
    aes(PC1, PC2, fill = nema, shape = rhizo, color = tissue, stroke = 2, alpha = 0.8) + 
    scale_color_manual(values = c("Root" = "tan", "Nodule" = "blue", "Gall" = "red")) +
    scale_shape_manual(values = c("Em1021" = 23, "Em1022" = 21)) + 
    geom_point(size=3) + 
    ggtitle("PCA Plot") +
    xlab(paste("PC1: ", round(attr(PCAdata, "percentVar")[1]*100, 3), "% variance")) + 
    ylab(paste("PC2: ", round(attr(PCAdata, "percentVar")[2]*100, 3), "% variance")) +
    theme_classic()
  if (nema_factor == TRUE){
    try(pcaplot <- pcaplot +
          scale_fill_manual(values = c("Nema" = "dark red", "No_Nema" = "white")) +
          guides(fill = guide_legend(override.aes=list(colour = c("Nema" = "dark red", "No_Nema" = "white")))))
  }
  nudge <- position_nudge(y = 1)
  print(pcaplot)
  pcaplot <- pcaplot + geom_text(aes(label = name), position = nudge) 
  print(pcaplot)
  
  #MDS plots
  print("Dataset quality check: Preparing MDS plots")
  sampleDists <- dist( t( assay(rld) ) )
  sampleDistMatrix <- as.matrix(sampleDists)
  eu_mds <- data.frame(cmdscale(sampleDistMatrix))
  eu_mds <- cbind(eu_mds, as.data.frame(colData(rld)))
  eu_mds$name <- rownames(eu_mds)
  eumds_plot <- ggplot(eu_mds) +
    geom_hline(yintercept = 0, color = "dark gray") + geom_vline(xintercept = 0, color = "dark gray") +
    aes(X1, X2, fill = nema, shape = rhizo, color = tissue, stroke = 3) + 
    scale_color_manual(values = c("Root" = "tan", "Nodule" = "blue", "Gall" = "red")) +
    scale_shape_manual(values = c("Em1021" = 23, "Em1022" = 21)) + 
    geom_point(size=3) + 
    ggtitle("Euclidian MDS plot")+
    theme_classic()
  if (nema_factor == TRUE){
    try(eumds_plot <- eumds_plot +
          scale_fill_manual(values = c("Nema" = "dark red", "No_Nema" = "white")) +
          guides(fill = guide_legend(override.aes=list(colour = c("Nema" = "dark red", "No_Nema" = "white")))))
  }
  nudge <- position_nudge(y = max(max(eu_mds$X1),max(eu_mds$X2))/15)
  print(eumds_plot)
  eumds_plot <- eumds_plot + geom_text(aes(label = name), position = nudge)
  print(eumds_plot)
  
  poisd <- PoissonDistance( t( counts(dds) ) )
  samplePoisDistMatrix <- as.matrix(poisd$dd)
  p_mds <- data.frame(cmdscale(samplePoisDistMatrix))
  p_mds <- cbind(p_mds, as.data.frame(colData(rld)))
  p_mds$names <- attr(dds, "colNames")
  pmds_plot <- ggplot(p_mds) +
    geom_hline(yintercept = 0, color = "dark gray") + geom_vline(xintercept = 0, color = "dark gray") +
    aes(X1, X2,  fill = nema, shape = rhizo, color = tissue, stroke = 3) + 
    scale_fill_manual(values = c("Nema" = "dark red", "No_Nema" = "white")) +
    scale_color_manual(values = c("Root" = "tan", "Nodule" = "blue", "Gall" = "red")) +
    scale_shape_manual(values = c("Em1021" = 23, "Em1022" = 21)) +
    geom_point(size=3) +  
    ggtitle("Poison MDS plot")+
    theme_classic()
  if (nema_factor == TRUE){try(pmds_plot <- pmds_plot + guides(fill = guide_legend(override.aes=list(colour = c("Nema" = "dark red", "No_Nema" = "white")))))}
  nudge <- position_nudge(y = max(max(p_mds$X1),max(p_mds$X2))/15)
  print(pmds_plot)
  eumds_plot <- pmds_plot + geom_text(aes(label = name), position = nudge)
  print(pmds_plot)
  
  #Gene clustering of most variable genes
  print("Dataset quality check: Preparing gene clustering diagrams")
  #top 20
  topVarGenes <- utils::head(order(-rowVars(assay(rld))),20)
  mat <- assay(rld)[ topVarGenes, ]
  df <- as.data.frame(colData(rld)[ , c("tissue", "rhizo", "nema")])
  pheatmap(mat, annotation_col=df, main = "Gene clustering of 20 most variable genes")
  
  #top50
  topVarGenes <- utils::head(order(-rowVars(assay(rld))),50)
  mat <- assay(rld)[ topVarGenes, ]
  df <- as.data.frame(colData(rld)[ , c("tissue", "rhizo", "nema")])
  pheatmap(mat, annotation_col=df, main = "Gene clustering of 50 most variable genes")
  
  #top100
  topVarGenes <- utils::head(order(-rowVars(assay(rld))),100)
  mat <- assay(rld)[ topVarGenes, ]
  df <- as.data.frame(colData(rld)[ , c("tissue", "rhizo", "nema")])
  pheatmap(mat, annotation_col=df, main = "Gene clustering of 100 most variable genes")
  
  #top500
  topVarGenes <- head(order(-rowVars(assay(rld))),500)
  mat <- assay(rld)[ topVarGenes, ]
  df <- as.data.frame(colData(rld)[ , c("tissue", "rhizo", "nema")])
  pheatmap(mat, annotation_col=df, main = "Gene clustering of 500 most variable genes")
  
  #top1000
  topVarGenes <- head(order(-rowVars(assay(rld))),1000)
  mat <- assay(rld)[ topVarGenes, ]
  df <- as.data.frame(colData(rld)[ , c("tissue", "rhizo", "nema")])
  pheatmap(mat, annotation_col=df, main = "Gene clustering of 1000 most variable genes")
  
  #Closing pdf file
  dev.off()
  print("Dataset quality check: done")
  setwd(dir_main)
}

  setwd(dir_main)
  dataset_quality_check_plotting(dds_medicago, "Medicago_dataset_quality_check", TRUE)
  try(dev.off(), silent = TRUE)
  try(dev.off(), silent = TRUE)
  dataset_quality_check_plotting(dds_medicago_nods_qc, "Medicago_nods_dataset_quality_check")
  try(dev.off(), silent = TRUE)
  try(dev.off(), silent = TRUE)
  dataset_quality_check_plotting(dds_medicago_gano_qc, "Medicago_gano_dataset_quality_check")
  try(dev.off(), silent = TRUE)
  try(dev.off(), silent = TRUE)
  dataset_quality_check_plotting(dds_medicago_garo_qc, "Medicago_garo_dataset_quality_check")
  try(dev.off(), silent = TRUE)
  try(dev.off(), silent = TRUE)
  dataset_quality_check_plotting(dds_medicago_noro_qc, "Medicago_noro_dataset_quality_check")
  try(dev.off(), silent = TRUE)
  try(dev.off(), silent = TRUE)
  dataset_quality_check_plotting(dds_medicago_gall_qc, "Medicago_gall_dataset_quality_check")
  try(dev.off(), silent = TRUE)
  try(dev.off(), silent = TRUE)
  dataset_quality_check_plotting(dds_nematode, "Nematode_dataset_quality_check")
  try(dev.off(), silent = TRUE)
  try(dev.off(), silent = TRUE)
  dataset_quality_check_plotting(dds_rhizob21, "RhizobialEm1021_dataset_quality_check")
  try(dev.off(), silent = TRUE)
  try(dev.off(), silent = TRUE)
  dataset_quality_check_plotting(dds_rhizob22, "RhizobialEm1022_dataset_quality_check")
  try(dev.off(), silent = TRUE)
  try(dev.off(), silent = TRUE)
