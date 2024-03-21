################################################################################
# Author: Addison Dale Buxton-Martin
# 
# Description: This script runs the primarily analysis for ___
# This analysis will take gene expression quantification data (output of
# HTSeq-union method) and perform various differential expression
# analyses as outlined in the article found at ___.
#
# Data Requirements: This script requires that it be placed in a directory that 
# also contains a medicago, nematode, and two rhizobe subdirectories where 
# genome annotations for these organisms are found. Directories must be entitled
# "Medicago", "Nema", "Rhizob21", and "Rhizob22". All annotations must be in the
# form of .gtf and entitled "annotat.gtf". Inside of each directories are
# subdirectories entitled "counts" containing the count data from HTSeq. Each of
# these count files is prefaced with the sample number and an underscore then
# "count". The main directory must also include a samples.csv file that outlines
# the sample types as specified by Michael Love's DESeq2 workflow.
################################################################################


# loading requisite libraries
library("DESeq2")
library("GenomicFeatures")
library("rtracklayer")
library("GenomicRanges")
library("tidyverse")




################################################################################
############################## DIRECTORY HANDLING ##############################
################################################################################

################################################################################
# Edit the file path of dir_main to be the folder containing this script.
# The data should be in this main directory as specified in the above Data
# Requirements section in the heading of this script
dir_main <- file.path("~", "Desktop", "DGE")
################################################################################

dir_medic <- file.path(dir_main, "Medicago")
dir_nema <- file.path(dir_main, "Nema")
dir_rhizo21 <- file.path(dir_main, "Rhizob21")
dir_rhizo22 <- file.path(dir_main, "Rhizob22")

dir_out <- file.path(dir_main, "out")
dir.create(dir_out, showWarnings = FALSE)

setwd(dir_main)



### LOADING DATA IF SCRIPT HAS PREVIOUSLY BEEN RUN ###

setwd(dir_out)
if (file.exists("Diff_express_results_table.RData"))
{
  load(file = "Diff_express_results_table.RData")
  t <- t %>% as.data.frame()
  setwd(dir_main)
  stop
}else{
  setwd(dir_main)
}

#### STOP HERE IF SCRIPT HAS PREVIOUSLY BEEN RUN ###







################################################################################
######################## SUBSET LISTS FOR DATA PROCESS #########################
################################################################################

sampleData <- read.csv(file.path(dir_main, "samples.csv"),
                       row.names = 1,
                       header = TRUE,
                       check.names = FALSE)

all_subset <- rownames(sampleData)

gall_subset <- sampleData %>% filter(tissue == "Gall") %>% rownames()
nodp_subset <- sampleData %>% filter(tissue == "Nodule" & nema == "Nema") %>% rownames()
nodm_subset <- sampleData %>% filter(tissue == "Nodule" & nema == "No_Nema") %>% rownames()
root_subset <- sampleData %>% filter(tissue == "Root") %>% rownames()

gall21_subset <- sampleData %>% filter(tissue == "Gall" & rhizo == "Em1021") %>% rownames()
nodp21_subset <- sampleData %>% filter(tissue == "Nodule" & nema == "Nema" & rhizo == "Em1021") %>% rownames()
nodm21_subset <- sampleData %>% filter(tissue == "Nodule" & nema == "No_Nema" & rhizo == "Em1021") %>% rownames()
root21_subset <- sampleData %>% filter(tissue == "Root" & rhizo == "Em1021") %>% rownames()

gall22_subset <- sampleData %>% filter(tissue == "Gall" & rhizo == "Em1022") %>% rownames()
nodp22_subset <- sampleData %>% filter(tissue == "Nodule" & nema == "Nema" & rhizo == "Em1022") %>% rownames()
nodm22_subset <- sampleData %>% filter(tissue == "Nodule" & nema == "No_Nema" & rhizo == "Em1022") %>% rownames()
root22_subset <- sampleData %>% filter(tissue == "Root" & rhizo == "Em1021") %>% rownames()

medicago_subset <- all_subset
rhizob21_subset <- sampleData %>% filter(rhizo == "Em1021" & tissue == "Nodule") %>% rownames()
rhizob22_subset <- sampleData %>% filter(rhizo == "Em1022" & tissue == "Nodule") %>% rownames()
nematode_subset <- sampleData %>% filter(tissue == "Gall") %>% rownames()

nods_subset <- sampleData %>% filter(tissue == "Nodule") %>% rownames()
nods21_subset <- sampleData %>% filter(tissue == "Nodule" & rhizo == "Em1021") %>% rownames()
nods22_subset <- sampleData %>% filter(tissue == "Nodule" & rhizo == "Em1022") %>% rownames()
garo_subset <- sampleData %>% filter(tissue == "Gall" | tissue == "Root") %>% rownames()
gano_subset <- sampleData %>% filter(tissue == "Gall" | (tissue == "Nodule" & nema == "Nema")) %>% rownames()
noro_subset <- sampleData %>% filter(tissue == "Root" | (tissue == "Nodule" & nema == "Nema")) %>% rownames()

alln_subset <- sampleData %>% filter((tissue == "Nodule" & nema == "Nema") | tissue == "Gall" | tissue == "Root") %>% rownames()




################################################################################
############################# COUNT DATA IMPORTING #############################
################################################################################

import_counts <- function(directory) {
  # Function that imports the count files into data frame
  # Will take all files whose names are prefixed by a number followed by a
  # subscript and then "count" (i.e. "1_count") and load them into a dataframe
  files <- list.files(file.path(directory, "counts"))
  if (length(files) == length(which(endsWith(files, "_count")))) {
    return_df <- data.frame(gene = character())
    for (file in files){
      sample_num <- as.numeric(substr(file, 1, nchar(file)-6))
      if (is.na(sample_num)) {
        print("ERROR: A file in the directory supplied to import counts does not match required naming scheme")
        exit()
      }
      file_path <- file.path(directory, "counts", file)
      print(paste("Importing count data from:", file_path))
      imported_file <- read.csv(file_path, sep = "\t", header = FALSE)
      colnames(imported_file) <- c("gene", paste0("s", sample_num))
      return_df <- full_join(imported_file, return_df, by = "gene")}
    rownames(return_df) <- return_df$gene
    return(return_df)
  } else {
    print("ERROR: a file in the directory supplied to import counts does not match required naming scheme")
  }}


#Importing counts and removing the read data for poor or ambiguous alignments
#and putting them into count tables
remove_rows <- c("__too_low_aQual", "__not_aligned",
                 "__no_feature", "__ambiguous", "__alignment_not_unique")
medicago_counts <- import_counts(dir_medic) %>% filter(!(gene %in% remove_rows))
nematode_counts <- import_counts(dir_nema) %>% filter(!(gene %in% remove_rows))
rhizob21_counts <- import_counts(dir_rhizo21) %>% filter(!(gene %in% remove_rows))
rhizob22_counts <- import_counts(dir_rhizo22) %>% filter(!(gene %in% remove_rows))



################################################################################
############################# GETTING GENE LENGTH ##############################
################################################################################

prok_gene_length_processing <- function(directory){
  # Function to pull gene length information from gtf file of a prokaryote
  GTF_file_path <- file.path(directory, "annotation.gtf")
  rangeddata <- import.gff(GTF_file_path)
  granges <- as(rangeddata, "GRanges")
  granges <- as_tibble(granges)
  gene_length <- granges %>%  
    dplyr::select(gene_id, type, start, end) %>% 
    filter(type == "gene") %>% 
    rowwise() %>% mutate(width = (max(start, end) - min(start, end))/1000) %>% 
    dplyr::select(gene_id, width)
  # gene_length <- as.data.table(gene_length)
  colnames(gene_length) <- c("gene", "widthInKB")
  return(gene_length)
}

euk_gene_length_processing <- function(directory){
  # Function to pull gene length information from gtf file of a eukaryote
  GTF_file_path <- file.path(directory, "annotation.gtf")
  txdb <- makeTxDbFromGFF(GTF_file_path, format = "gtf")
  edg <- exonsBy(txdb, by = "gene")
  gene_length <- as.data.frame(sum(width(edg)/1000))
  base::colnames(gene_length) <- c("widthInKB")
  gene_length$gene <- rownames(gene_length)
  rownames(gene_length) <- NULL
  return(gene_length)
}    

#building a table of gene information for each organism
medicago_genes <- euk_gene_length_processing(dir_medic) %>% mutate(organism = "Medicago")
nematode_genes <- euk_gene_length_processing(dir_nema) %>% mutate(organism = "Nematode")
rhizob21_genes <- prok_gene_length_processing(dir_rhizo21) %>% mutate(organism = "Rhizob21")
rhizob22_genes <- prok_gene_length_processing(dir_rhizo22) %>% mutate(organism = "Rhizob22")

#filtering out erroneous duplicate genes from gene length function
rhizob21_genes <- rhizob21_genes %>%
  group_by(gene) %>%
  summarize(widthInKB = max(widthInKB)) %>%
  ungroup() %>%
  mutate(organism = "Rhizob21")
rhizob22_genes <- rhizob22_genes %>% group_by(gene) %>% ungroup()

#joining count data tables and gene length data tables
medicago_pretable <- full_join(medicago_genes, medicago_counts, by = "gene") %>% mutate(organism = "Medicago")
nematode_pretable <- full_join(nematode_genes, nematode_counts, by = "gene") %>% mutate(organism = "Nematode")
rhizob21_pretable <- full_join(rhizob21_genes, rhizob21_counts, by = "gene") %>% mutate(organism = "Rhizob21")
rhizob22_pretable <- full_join(rhizob22_genes, rhizob22_counts, by = "gene") %>% mutate(organism = "Rhizob22")





################################################################################
#################### CALCULATING RPKM AND EXPRESSION STATUS ####################
################################################################################

rpkm_and_base_exp <- function(table, rpkm_threshold = 1, exp_threshold = 2){
  # Function to calculate RPKM for each gene in every sample and to determine
  # expression status from that RPKM value
  
  # getting list of all the samples in the passed table and making blank
  # dataframes of rpkm and expression status for later merging into final
  # returned table
  cols <- table %>% dplyr::select(-c("gene", "widthInKB", "organism")) %>% base::colnames()
  df_rpkm <- table %>% dplyr::select("gene")
  df_exp <- table %>% dplyr::select("gene")
  
  # iterating through samples in table to calculate rpkm and expression statuses
  for (col in cols) {
    t <- as_tibble(table)
    rpkm_col <- paste0(col, "_rpkm")
    t[[rpkm_col]] <- t[[col]] / sum(t[[col]], na.rm = TRUE) / t$widthInKB * 1000000
    exp_col <- paste0("exp_", col)
    t[[exp_col]] <- ifelse(t[[rpkm_col]] >= rpkm_threshold, TRUE, FALSE)  
    
    #joining with blank dataframes created earlier
    df_rpkm <- full_join(
      as_tibble(t %>% dplyr::select(c("gene", rpkm_col))),
      df_rpkm,
      by = "gene")
    df_exp <- full_join(
      as_tibble(t %>% dplyr::select(c("gene", exp_col))), 
      df_exp, 
      by = "gene")
  }
  
  #joining rpkm and expression tables
  df <- full_join(df_rpkm, df_exp, by = "gene")
  return_table <- full_join(table, df, by = "gene")
  
  # inner function to determine expression status in categories of tissues
  # this function will return the table with an added column indicating the 
  # expression status in the subset fed to it
  exp_row_sums <- function(table, subset, col_name){
    if (sum(subset %in% cols) == length(subset)){
      return_table <- table %>% 
        mutate("{col_name}" := ifelse(
          rowSums(
            dplyr::select(.,
                   paste0("exp_", subset)),
            na.rm = TRUE) >= exp_threshold, 
          TRUE, 
          FALSE))
    }
    return(return_table)
  }
  return_table <- exp_row_sums(return_table, gall21_subset, "exp_gall21")
  return_table <- exp_row_sums(return_table, gall22_subset, "exp_gall22")
  return_table <- exp_row_sums(return_table, nodp21_subset, "exp_nodp21")
  return_table <- exp_row_sums(return_table, nodp22_subset, "exp_nodp22")
  return_table <- exp_row_sums(return_table, nodm21_subset, "exp_nodm21")
  return_table <- exp_row_sums(return_table, nodm22_subset, "exp_nodm22")
  return_table <- exp_row_sums(return_table, root21_subset, "exp_root21")
  return_table <- exp_row_sums(return_table, root22_subset, "exp_root22")
  return(return_table)
}

# Adding in RPKM and expression status to data tables
medicago_table <- rpkm_and_base_exp(medicago_pretable)
rhizob21_table <- rpkm_and_base_exp(rhizob21_pretable)
rhizob22_table <- rpkm_and_base_exp(rhizob22_pretable)
nematode_table <- rpkm_and_base_exp(nematode_pretable)

expression_col <- function(table){
  # Function to determine the expression status of sample types that are
  # composite sample types (ex: gall sample type which is composed of galls
  # from hosts with Em1021 and galls from hosts with Em1022)
  
  expression <- function(table, subset, col_name){
    # Inner function to determine expression from given columns
    cols <- table %>% base::colnames()
    return_table <- table
    if (sum(subset %in% cols) == length(subset)){
      return_table <- return_table %>%
        mutate("{col_name}" := ifelse(rowSums(dplyr::select(., subset), na.rm = TRUE) > 0,
                                      TRUE,
                                      FALSE))
    }
    return(return_table)
  }
  
  rt <- table
  rt <- expression(rt, c("exp_gall21", "exp_gall22"), "exp_gall")
  rt <- expression(rt, c("exp_nodp21", "exp_nodp22"), "exp_nodp")
  rt <- expression(rt, c("exp_nodm21", "exp_nodm22"), "exp_nodm")
  rt <- expression(rt, c("exp_root21", "exp_root22"), "exp_root")
  rt <- expression(rt, c("exp_nodp21", "exp_nodm21"), "exp_nods21")
  rt <- expression(rt, c("exp_nodp22", "exp_nodm22"), "exp_nods22")
  rt <- expression(rt, c("exp_nodm", "exp_nodp", "exp_nods21", "exp_nods22"), "exp_nods")
  rt <- expression(rt, c("exp_gall", "exp_root"), "exp_garo")
  rt <- expression(rt, c("exp_gall", "exp_nodp"), "exp_gano")
  rt <- expression(rt, c("exp_nodp", "exp_root"), "exp_noro")
  rt <- expression(rt, c("exp_nods", "exp_root", "exp_gall"), "exp")
  return(rt)
}

#Adding expression status for sample types to table
medicago_table <- expression_col(medicago_table)
rhizob21_table <- expression_col(rhizob21_table)
rhizob22_table <- expression_col(rhizob22_table)
nematode_table <- expression_col(nematode_table)

rm(medicago_pretable, medicaog_genes, nematode_pretable, nematode_genes,
   rhizob21_pretable, rhizobe21_genes, rhizob22_pretable, rhizobe22_genes)




################################################################################
############################## CREATING GENE LISTS #############################
################################################################################

filter_reads <- function(counts, gene_list){
  # Function for pulling genes out by expression tatus
  returnCounts <- counts[which(counts$gene %in% gene_list), ]
  print(paste("Total in the list provided to keep:", length(gene_list)))
  print(paste("Total features with counts after filtering:", length(returnCounts$gene)))
  print(paste("Number of features removed:", length(counts$gene) - length(returnCounts$gene)))
  return(returnCounts)
}

gl_med_exp <- medicago_table %>% filter(exp_gall | exp_nods | exp_root) %>% filter(organism == "Medicago") %>% pull(gene)
gl_med_exp_gall <- medicago_table %>% filter(exp_gall) %>% filter(organism == "Medicago") %>% pull(gene)
gl_med_exp_nodp <- medicago_table %>% filter(exp_nodp) %>% filter(organism == "Medicago") %>% pull(gene)
gl_med_exp_nodm <- medicago_table %>% filter(exp_nodm) %>% filter(organism == "Medicago") %>% pull(gene)
gl_med_exp_root <- medicago_table %>% filter(exp_root) %>% filter(organism == "Medicago") %>% pull(gene)
gl_med_exp_nods_em21 <- medicago_table %>% filter(exp_nods21) %>% filter(organism == "Medicago") %>% pull(gene)
gl_med_exp_nods_em22 <- medicago_table %>% filter(exp_nods22) %>% filter(organism == "Medicago") %>% pull(gene)
gl_med_exp_gall_and_nodp <- medicago_table %>% filter(exp_gall & exp_nodp) %>% filter(organism == "Medicago") %>% pull(gene)
gl_med_exp_gall_and_root <- medicago_table %>% filter(exp_gall & exp_root) %>% filter(organism == "Medicago") %>% pull(gene)
gl_med_exp_nodp_and_root <- medicago_table %>% filter(exp_root & exp_nodp) %>% filter(organism == "Medicago") %>% pull(gene)
gl_med_exp_nodp_and_nodm <- medicago_table %>% filter(exp_nodp & exp_nodm) %>% filter(organism == "Medicago") %>% pull(gene)
gl_med_exp_gall_or_nodp <- medicago_table %>% filter(exp_gall | exp_nodp) %>% filter(organism == "Medicago") %>% pull(gene)
gl_med_exp_gall_or_root <- medicago_table %>% filter(exp_gall | exp_root) %>% filter(organism == "Medicago") %>% pull(gene)
gl_med_exp_nodp_or_root <- medicago_table %>% filter(exp_root | exp_nodp) %>% filter(organism == "Medicago") %>% pull(gene)
gl_med_exp_nodp_or_nodm <- medicago_table %>% filter(exp_nodp | exp_nodm) %>% filter(organism == "Medicago") %>% pull(gene)
gl_med_exp_all_and <- medicago_table %>% filter(exp_gall & exp_nodp & exp_root) %>% filter(organism == "Medicago") %>% pull(gene)
gl_med_exp_all_or <- medicago_table %>% filter(exp_gall | exp_nodp | exp_root) %>% filter(organism == "Medicago") %>% pull(gene)

gl_r21_exp_nods_em21 <- rhizob21_table %>% filter(exp_nods21) %>% filter(organism == "Rhizob21") %>% pull(gene)
gl_r22_exp_nods_em22 <- rhizob22_table %>% filter(exp_nods22) %>% filter(organism == "Rhizob22") %>% pull(gene)

gl_nem_exp_gall <- nematode_table %>% filter(exp_gall) %>% filter(organism == "Nematode") %>% pull(gene)




################################################################################
############################ DIFFERENTIAL EXPRESSION ###########################
################################################################################

DESeq_setup <- function(subset, counts, equation){
  # Function to set up a DESeq dataset
  # Multiple tests can be run from a single dataset and so I've separated the
  # creation of datset objects from the actual tests
  samples <- sampleData[subset, ]
  counts <- counts[, subset]
  dds <- DESeqDataSetFromMatrix(counts, samples, equation)
  dds <- DESeq(dds, quiet = TRUE)
  return(dds)
}

Run_DESeq <- function(dds, equation, reduced_eq, effect,
                      rhizo_ref = FALSE, nema_ref = FALSE, tissue_ref = FALSE,
                      return_df = TRUE) {
  # Function to actually run DE analyses from DESeq dataset objects
  if (rhizo_ref != FALSE) {dds$rhizo <- relevel(dds$rhizo, ref = rhizo_ref)}
  if (nema_ref != FALSE) {dds$nema <- relevel(dds$nema, ref = nema_ref)}
  if (tissue_ref != FALSE) {dds$tissue <- relevel(dds$tissue, ref = tissue_ref)}
  dds <- DESeq(dds, test = "LRT", full = equation, reduced = reduced_eq)
  results <- lfcShrink(dds, coef = effect, type = "ashr")
  # return_df parameter will tell function to return a dataframe of results
  # result columns are renamed slightly
  # column names with ns before them are tests without shrinkage applied
  # tests with shrinkage applied do not have the "ns_" prefix
  if (return_df == TRUE){
    results_ns <- as.data.frame(results(dds, name = effect))
    results_ns$gene <- rownames(results_ns)
    results_ns <- results_ns %>% dplyr::select("gene", "log2FoldChange", "lfcSE") %>%
      dplyr::rename("ns_LFC" = "log2FoldChange", "ns_SE" = "lfcSE")
    res <- as.data.frame(results)
    res$gene <- rownames(res)
    res <- res %>% dplyr::select(c("gene", "baseMean", "log2FoldChange", "lfcSE", "padj")) %>%
      dplyr::rename("LFC" = "log2FoldChange", "SE" = "lfcSE")
    final_res <- full_join(res, results_ns, by = "gene")
    final_res <- final_res %>%
      mutate(significance = ifelse(is.na(padj), "Inconclusive", ifelse(padj < 0.05, "Significant", "Not significant"))) %>%
      mutate(sig = ifelse(significance == "Significant", TRUE, FALSE)) %>%
      mutate(sigu = ifelse(sig & LFC > 0, TRUE, FALSE)) %>%
      mutate(sigd = ifelse(sig & LFC < 0, TRUE, FALSE)) %>%
      mutate(sigu_ns = ifelse(sig & ns_LFC > 0, TRUE, FALSE)) %>%
      mutate(sigd_ns = ifelse(sig & ns_LFC < 0, TRUE, FALSE))
    return(final_res)
  }
  else {
    return(results)
  }
}




################################################################################
##################### QUALITY CONTROL AND SIMPLE DATASETS ######################
################################################################################

# Datsets to draw normalized count data from
dds_medicago <- DESeq_setup(medicago_subset, medicago_counts, ~1)
dds_nematode <- DESeq_setup(nematode_subset, nematode_counts, ~1)
dds_rhizob21 <- DESeq_setup(rhizob21_subset, rhizob21_counts, ~1)
dds_rhizob22 <- DESeq_setup(rhizob22_subset, rhizob22_counts, ~1)

# Drawing normalized counts from datasets
medicago_norm_counts <- as.data.frame(DESeq2::counts(dds_medicago, normalized = TRUE))
nematode_norm_counts <- as.data.frame(DESeq2::counts(dds_nematode, normalized = TRUE))
rhizob21_norm_counts <- as.data.frame(DESeq2::counts(dds_rhizob21, normalized = TRUE))
rhizob22_norm_counts <- as.data.frame(DESeq2::counts(dds_rhizob22, normalized = TRUE))

colnames(medicago_norm_counts) <- paste0(colnames(medicago_norm_counts), "_norm")
colnames(nematode_norm_counts) <- paste0(colnames(nematode_norm_counts), "_norm")
colnames(rhizob21_norm_counts) <- paste0(colnames(rhizob21_norm_counts), "_norm")
colnames(rhizob22_norm_counts) <- paste0(colnames(rhizob22_norm_counts), "_norm")

ave_norm <- function(table){
  cols <- colnames(table)
  table$gene <- rownames(table)
  return_table <- table
  if ("s1_norm" %in% cols & "s2_norm" %in% cols & "s5_norm" %in% cols & "s6_norm" %in% cols & "s9_norm" %in% cols & "s10_norm" %in% cols){
    return_table <- return_table %>% rowwise() %>% mutate(ave_norm_root = mean(s1_norm, s2_norm, s5_norm, s6_norm, s9_norm, s10_norm))}
  if ("s13_norm" %in% cols & "s14_norm" %in% cols & "s17_norm" %in% cols & "s18_norm" %in% cols & "s21_norm" %in% cols & "s22_norm" %in% cols){
    return_table <- return_table %>% rowwise() %>% mutate(ave_norm_gall = mean(s13_norm, s14_norm, s17_norm, s18_norm, s21_norm, s22_norm))}
  if ("s25_norm" %in% cols & "s26_norm" %in% cols & "s29_norm" %in% cols & "s30_norm" %in% cols & "s33_norm" %in% cols & "s34_norm" %in% cols){
    return_table <- return_table %>% rowwise() %>% mutate(ave_norm_nodp = mean(s25_norm, s26_norm, s29_norm, s30_norm, s33_norm, s34_norm))}
  if ("s37_norm" %in% cols & "s38_norm" %in% cols & "s41_norm" %in% cols & "s42_norm" %in% cols & "s45_norm" %in% cols & "s46_norm" %in% cols){
    return_table <- return_table %>% rowwise() %>% mutate(ave_norm_nodm = mean(s37_norm, s38_norm, s41_norm, s42_norm, s45_norm, s46_norm))}
  if ("s1_norm" %in% cols & "s5_norm" %in% cols & "s9_norm" %in% cols){
    return_table <- return_table %>% rowwise() %>% mutate(ave_norm_root_rh21 = mean(s1_norm, s5_norm, s9_norm))}
  if ("s2_norm" %in% cols & "s6_norm" %in% cols & "s10_norm" %in% cols){
    return_table <- return_table %>% rowwise() %>% mutate(ave_norm_root_rh22 = mean(s2_norm, s6_norm, s10_norm))}
  if ("s13_norm" %in% cols & "s17_norm" %in% cols & "s21_norm" %in% cols){
    return_table <- return_table %>% rowwise() %>% mutate(ave_norm_gall_rh21 = mean(s13_norm, s17_norm, s21_norm))}
  if ("s14_norm" %in% cols & "s18_norm" %in% cols & "s22_norm" %in% cols){
    return_table <- return_table %>% rowwise() %>% mutate(ave_norm_gall_rh22 = mean(s14_norm, s18_norm, s22_norm))}
  if ("s25_norm" %in% cols & "s29_norm" %in% cols & "s33_norm" %in% cols){
    return_table <- return_table %>% rowwise() %>% mutate(ave_norm_nodp_rh21 = mean(s25_norm, s29_norm, s33_norm))}
  if ("s26_norm" %in% cols & "s30_norm" %in% cols & "s34_norm" %in% cols){
    return_table <- return_table %>% rowwise() %>% mutate(ave_norm_nodp_rh22 = mean(s26_norm, s30_norm, s34_norm))}
  if ("s37_norm" %in% cols & "s41_norm" %in% cols & "s45_norm" %in% cols){
    return_table <- return_table %>% rowwise() %>% mutate(ave_norm_nodm_rh21 = mean(s37_norm, s41_norm, s45_norm))}
  if ("s38_norm" %in% cols & "s42_norm" %in% cols & "s46_norm" %in% cols){
    return_table <- return_table %>% rowwise() %>% mutate(ave_norm_nodm_rh22 = mean(s38_norm, s42_norm, s46_norm))}
  if ("s25_norm" %in% cols & "s29_norm" %in% cols & "s33_norm" %in% cols & "s37_norm" %in% cols & "s41_norm" %in% cols & "s45_norm" %in% cols){
    return_table <- return_table %>% rowwise() %>% mutate(ave_norm_rh21 = mean(s25_norm, s29_norm, s33_norm, s37_norm, s41_norm, s45_norm))}
  if ("s26_norm" %in% cols & "s30_norm" %in% cols & "s34_norm" %in% cols & "s38_norm" %in% cols & "s42_norm" %in% cols & "s46_norm" %in% cols){
    return_table <- return_table %>% rowwise() %>% mutate(ave_norm_rh22 = mean(s26_norm, s30_norm, s34_norm, s38_norm, s42_norm, s46_norm))}
  return(return_table) 
}

medicago_norm_counts <- ave_norm(medicago_norm_counts)
nematode_norm_counts <- ave_norm(nematode_norm_counts)
rhizob21_norm_counts <- ave_norm(rhizob21_norm_counts)
rhizob22_norm_counts <- ave_norm(rhizob22_norm_counts)

medicago_final_table <- left_join(medicago_table, medicago_norm_counts, by = "gene")
nematode_final_table <- left_join(nematode_table, nematode_norm_counts, by = "gene")
rhizob21_final_table <- left_join(rhizob21_table, rhizob21_norm_counts, by = "gene")
rhizob22_final_table <- left_join(rhizob22_table, rhizob22_norm_counts, by = "gene")

# Datasets to draw quality control metrics from
dds_medicago_gano_qc <- DESeq_setup(base::intersect(medicago_subset, gano_subset), medicago_counts, ~1)
dds_medicago_noro_qc <- DESeq_setup(base::intersect(medicago_subset, noro_subset), medicago_counts, ~1)
dds_medicago_nods_qc <- DESeq_setup(base::intersect(medicago_subset, nods_subset), medicago_counts, ~1)
dds_medicago_garo_qc <- DESeq_setup(base::intersect(medicago_subset, garo_subset), medicago_counts, ~1)
dds_medicago_gall_qc <- DESeq_setup(base::intersect(medicago_subset, gall_subset), medicago_counts, ~1)





################################################################################
################################### DE SETUP ###################################
################################################################################

dds_medicago_nodp_r     <- DESeq_setup(nodp_subset,                             filter_reads(medicago_counts, gl_med_exp_nodp), ~batch+rhizo)
dds_medicago_nodm_r     <- DESeq_setup(nodm_subset,                             filter_reads(medicago_counts, gl_med_exp_nodm), ~rhizo)
dds_medicago_gall_r     <- DESeq_setup(gall_subset,                             filter_reads(medicago_counts, gl_med_exp_gall), ~batch+rhizo)
dds_medicago_root_r     <- DESeq_setup(root_subset,                             filter_reads(medicago_counts, gl_med_exp_root), ~rhizo)
dds_medicago_nods_r21_n <- DESeq_setup(intersect(nods_subset, rhizob21_subset), filter_reads(medicago_counts, gl_med_exp_nods_em21), ~batch+nema)
dds_medicago_nods_r22_n <- DESeq_setup(intersect(nods_subset, rhizob22_subset), filter_reads(medicago_counts, gl_med_exp_nods_em22), ~batch+nema)



dds_medicago_gano_r     <- DESeq_setup(c(gall_subset, nodp_subset),             filter_reads(medicago_counts, gl_med_exp_gall_and_nodp), ~batch+tissue+rhizo)
dds_medicago_gano_i     <- DESeq_setup(c(gall_subset, nodp_subset),             filter_reads(medicago_counts, gl_med_exp_gall_or_nodp), ~batch+tissue+rhizo+rhizo:tissue)

dds_medicago_garo_r     <- DESeq_setup(garo_subset,                             filter_reads(medicago_counts, gl_med_exp_gall_and_root), ~batch+tissue+rhizo)
dds_medicago_garo_i     <- DESeq_setup(garo_subset,                             filter_reads(medicago_counts, gl_med_exp_gall_or_root), ~batch+tissue+rhizo+rhizo:tissue)

dds_medicago_noro_r     <- DESeq_setup(noro_subset,                             filter_reads(medicago_counts, gl_med_exp_nodp_and_root), ~batch+tissue+rhizo)
dds_medicago_noro_i     <- DESeq_setup(noro_subset,                             filter_reads(medicago_counts, gl_med_exp_nodp_or_root), ~batch+tissue+rhizo+rhizo:tissue)

dds_medicago_nods_n     <- DESeq_setup(nods_subset,                             filter_reads(medicago_counts, gl_med_exp_nodp_and_nodm), ~batch+nema+rhizo)
dds_medicago_nods_r     <- DESeq_setup(nods_subset,                             filter_reads(medicago_counts, gl_med_exp_nodp_and_nodm), ~batch+nema+rhizo)
dds_medicago_nods_i     <- DESeq_setup(nods_subset,                             filter_reads(medicago_counts, gl_med_exp_nodp_or_nodm), ~batch+nema+rhizo+rhizo:nema)

dds_medicago_all_r      <- DESeq_setup(alln_subset,                             filter_reads(medicago_counts, gl_med_exp_all_and),      ~batch+rhizo+tissue)
dds_medicago_all_i      <- DESeq_setup(alln_subset,                             filter_reads(medicago_counts, gl_med_exp_all_or),      ~batch+rhizo+tissue+rhizo:tissue)

dds_rhizob21_nods_n     <- DESeq_setup(intersect(nods_subset, rhizob21_subset), filter_reads(rhizob21_counts, gl_r21_exp_nods_em21), ~batch+nema)
dds_rhizob22_nods_n     <- DESeq_setup(intersect(nods_subset, rhizob22_subset), filter_reads(rhizob22_counts, gl_r22_exp_nods_em22), ~batch+nema)
dds_nematode_gall_r     <- DESeq_setup(gall_subset,                             filter_reads(nematode_counts, gl_nem_exp_gall), ~batch+rhizo)



################################################################################
################################## RUNNING DE ##################################
################################################################################

med_nodp_r     <- Run_DESeq(dds_medicago_nodp_r,     ~batch+rhizo,                 ~batch,            "rhizo_Em1022_vs_Em1021", rhizo_ref = "Em1021")
med_nodm_r     <- Run_DESeq(dds_medicago_nodm_r,     ~rhizo,                       ~1,            "rhizo_Em1022_vs_Em1021", rhizo_ref = "Em1021")
med_gall_r     <- Run_DESeq(dds_medicago_gall_r,     ~batch+rhizo,                 ~batch,            "rhizo_Em1022_vs_Em1021", rhizo_ref = "Em1021")
med_root_r     <- Run_DESeq(dds_medicago_root_r,     ~rhizo,                       ~1,                "rhizo_Em1022_vs_Em1021", rhizo_ref = "Em1021")
med_nods_r21_n <- Run_DESeq(dds_medicago_nods_r21_n, ~batch+nema,                  ~batch,            "nema_Nema_vs_No_Nema",   nema_ref = "No_Nema")
med_nods_r22_n <- Run_DESeq(dds_medicago_nods_r22_n, ~batch+nema,                  ~batch,            "nema_Nema_vs_No_Nema",   nema_ref = "No_Nema")

med_gano_r        <- Run_DESeq(dds_medicago_gano_r, ~batch+tissue+rhizo,              ~batch+tissue,       "rhizo_Em1022_vs_Em1021", rhizo_ref = "Em1021", tissue_ref = "Nodule")
med_gano_i        <- Run_DESeq(dds_medicago_gano_i, ~batch+tissue+rhizo+rhizo:tissue, ~batch+rhizo+tissue, "tissueGall.rhizoEm1022", rhizo_ref = "Em1021", tissue_ref = "Nodule")
med_gano_i_nodp_r <- Run_DESeq(dds_medicago_gano_i, ~batch+tissue+rhizo+rhizo:tissue, ~batch+rhizo+tissue, "rhizo_Em1022_vs_Em1021", rhizo_ref = "Em1021", tissue_ref = "Nodule")
med_gano_i_gall_r <- Run_DESeq(dds_medicago_gano_i, ~batch+tissue+rhizo+rhizo:tissue, ~batch+rhizo+tissue, "rhizo_Em1022_vs_Em1021", rhizo_ref = "Em1021", tissue_ref = "Gall")

med_garo_r        <- Run_DESeq(dds_medicago_garo_r, ~batch+tissue+rhizo,              ~batch+tissue,       "rhizo_Em1022_vs_Em1021", rhizo_ref = "Em1021", tissue_ref = "Root")
med_garo_i        <- Run_DESeq(dds_medicago_garo_i, ~batch+tissue+rhizo+rhizo:tissue, ~batch+rhizo+tissue, "tissueGall.rhizoEm1022", rhizo_ref = "Em1021", tissue_ref = "Root")
med_garo_i_gall_r <- Run_DESeq(dds_medicago_garo_i, ~batch+tissue+rhizo+rhizo:tissue, ~batch+rhizo+tissue, "rhizo_Em1022_vs_Em1021", rhizo_ref = "Em1021", tissue_ref = "Gall")
med_garo_i_root_r <- Run_DESeq(dds_medicago_garo_i, ~batch+tissue+rhizo+rhizo:tissue, ~batch+rhizo+tissue, "rhizo_Em1022_vs_Em1021", rhizo_ref = "Em1021", tissue_ref = "Root")

med_noro_r        <- Run_DESeq(dds_medicago_noro_r, ~batch+tissue+rhizo,              ~batch+tissue,       "rhizo_Em1022_vs_Em1021", rhizo_ref = "Em1021", tissue_ref = "Root")
med_noro_i        <- Run_DESeq(dds_medicago_noro_i, ~batch+tissue+rhizo+rhizo:tissue, ~batch+rhizo+tissue, "tissueNodule.rhizoEm1022", rhizo_ref = "Em1021", tissue_ref = "Root")
med_noro_i_nodp_r <- Run_DESeq(dds_medicago_noro_i, ~batch+tissue+rhizo+rhizo:tissue, ~batch+rhizo+tissue, "rhizo_Em1022_vs_Em1021", rhizo_ref = "Em1021", tissue_ref = "Nodule")
med_noro_i_root_r <- Run_DESeq(dds_medicago_noro_i, ~batch+tissue+rhizo+rhizo:tissue, ~batch+rhizo+tissue, "rhizo_Em1022_vs_Em1021", rhizo_ref = "Em1021", tissue_ref = "Root")

med_nods_n         <- Run_DESeq(dds_medicago_nods_n,     ~batch+nema+rhizo,            ~batch+rhizo,      "nema_Nema_vs_No_Nema",   nema_ref = "No_Nema", rhizo_ref = "Em1021")
med_nods_r         <- Run_DESeq(dds_medicago_nods_r,     ~batch+nema+rhizo,            ~batch+nema,       "rhizo_Em1022_vs_Em1021", nema_ref = "No_Nema", rhizo_ref = "Em1021")
med_nods_i         <- Run_DESeq(dds_medicago_nods_i,     ~batch+nema+rhizo+rhizo:nema, ~batch+nema+rhizo, "nemaNema.rhizoEm1022",   rhizo_ref = "Em1021", nema_ref = "No_Nema")
med_nods_i_nodp_r  <- Run_DESeq(dds_medicago_nods_i,     ~batch+nema+rhizo+rhizo:nema, ~batch+nema+rhizo, "rhizo_Em1022_vs_Em1021",   rhizo_ref = "Em1021", nema_ref = "Nema")
med_nods_i_nodm_r  <- Run_DESeq(dds_medicago_nods_i,     ~batch+nema+rhizo+rhizo:nema, ~batch+nema+rhizo, "rhizo_Em1022_vs_Em1021",   rhizo_ref = "Em1021", nema_ref = "No_Nema")
med_nods_i_em21_n  <- Run_DESeq(dds_medicago_nods_i,     ~batch+nema+rhizo+rhizo:nema, ~batch+nema+rhizo, "nema_Nema_vs_No_Nema",   rhizo_ref = "Em1021", nema_ref = "No_Nema")
med_nods_i_em22_n  <- Run_DESeq(dds_medicago_nods_i,     ~batch+nema+rhizo+rhizo:nema, ~batch+nema+rhizo, "nema_Nema_vs_No_Nema",   rhizo_ref = "Em1022", nema_ref = "No_Nema")

med_all_r_root          <- Run_DESeq(dds_medicago_all_r,      ~batch+rhizo+tissue,              ~batch+tissue,       "rhizo_Em1022_vs_Em1021", rhizo_ref = "Em1021", tissue_ref = "Root")
med_all_r_gall          <- Run_DESeq(dds_medicago_all_r,      ~batch+rhizo+tissue,              ~batch+tissue,       "rhizo_Em1022_vs_Em1021", rhizo_ref = "Em1021", tissue_ref = "Nodule")
med_all_r_nodp          <- Run_DESeq(dds_medicago_all_r,      ~batch+rhizo+tissue,              ~batch+tissue,       "rhizo_Em1022_vs_Em1021", rhizo_ref = "Em1021", tissue_ref = "Gall")
med_all_i_gano     <- Run_DESeq(dds_medicago_all_i,      ~batch+rhizo+tissue+rhizo:tissue, ~batch+rhizo+tissue, "rhizoEm1022.tissueGall", rhizo_ref = "Em1021", tissue_ref = "Nodule")
med_all_i_gano_gall     <- Run_DESeq(dds_medicago_all_i,      ~batch+rhizo+tissue+rhizo:tissue, ~batch+rhizo+tissue, "rhizo_Em1022_vs_Em1021", rhizo_ref = "Em1021", tissue_ref = "Gall")
med_all_i_gano_nodp     <- Run_DESeq(dds_medicago_all_i,      ~batch+rhizo+tissue+rhizo:tissue, ~batch+rhizo+tissue, "rhizo_Em1022_vs_Em1021", rhizo_ref = "Em1021", tissue_ref = "Nodule")
med_all_i_noro     <- Run_DESeq(dds_medicago_all_i,      ~batch+rhizo+tissue+rhizo:tissue, ~batch+rhizo+tissue, "rhizoEm1022.tissueNodule", rhizo_ref = "Em1021", tissue_ref = "Root")
med_all_i_noro_nodp     <- Run_DESeq(dds_medicago_all_i,      ~batch+rhizo+tissue+rhizo:tissue, ~batch+rhizo+tissue, "rhizo_Em1022_vs_Em1021", rhizo_ref = "Em1021", tissue_ref = "Nodule")
med_all_i_noro_root     <- Run_DESeq(dds_medicago_all_i,      ~batch+rhizo+tissue+rhizo:tissue, ~batch+rhizo+tissue, "rhizo_Em1022_vs_Em1021", rhizo_ref = "Em1021", tissue_ref = "Root")
med_all_i_garo     <- Run_DESeq(dds_medicago_all_i,      ~batch+rhizo+tissue+rhizo:tissue, ~batch+rhizo+tissue, "rhizoEm1022.tissueGall", rhizo_ref = "Em1021", tissue_ref = "Root")
med_all_i_garo_gall     <- Run_DESeq(dds_medicago_all_i,      ~batch+rhizo+tissue+rhizo:tissue, ~batch+rhizo+tissue, "rhizo_Em1022_vs_Em1021", rhizo_ref = "Em1021", tissue_ref = "Gall")
med_all_i_garo_root     <- Run_DESeq(dds_medicago_all_i,      ~batch+rhizo+tissue+rhizo:tissue, ~batch+rhizo+tissue, "rhizo_Em1022_vs_Em1021", rhizo_ref = "Em1021", tissue_ref = "Root")

r21_nods_n         <- Run_DESeq(dds_rhizob21_nods_n,     ~batch+nema,                  ~batch,            "nema_Nema_vs_No_Nema",   nema_ref = "No_Nema")
r22_nods_n         <- Run_DESeq(dds_rhizob22_nods_n,     ~batch+nema,                  ~batch,            "nema_Nema_vs_No_Nema",   nema_ref = "No_Nema")
nem_gall_r         <- Run_DESeq(dds_nematode_gall_r,     ~batch+rhizo,                 ~batch,            "rhizo_Em1022_vs_Em1021", rhizo_ref = "Em1021")




################################################################################
################################ COMPILING DATA ################################
################################################################################

# Function to set up dataframes for merging into one table
# Function adds prefix to column names so that columns are decipherable 
# in the final merged table
process_cols <- function(table){
  pre <- substitute(table)
  rownames(table) <- table$gene
  table <- subset(table, select = -c(gene))
  colnames(table) <- paste0(pre, ".", colnames(table))
  table$gene <- rownames(table)
  rownames(table) <- NULL
  return(table)
}

# Processing columns
med_nodp_r_res     <- process_cols(med_nodp_r)
med_nodm_r_res     <- process_cols(med_nodm_r)
med_gall_r_res     <- process_cols(med_gall_r)
med_root_r_res     <- process_cols(med_root_r)
med_nods_r21_n_res <- process_cols(med_nods_r21_n)
med_nods_r22_n_res <- process_cols(med_nods_r22_n)

med_gano_r_res     <- process_cols(med_gano_r)
med_gano_i_res     <- process_cols(med_gano_i)
med_gano_i_gall_r_res     <- process_cols(med_gano_i_gall_r)
med_gano_i_nodp_r_res     <- process_cols(med_gano_i_nodp_r)

med_garo_r_res     <- process_cols(med_garo_r)
med_garo_i_res     <- process_cols(med_garo_i)
med_garo_i_gall_r_res     <- process_cols(med_garo_i_gall_r)
med_garo_i_root_r_res     <- process_cols(med_garo_i_root_r)

med_noro_r_res     <- process_cols(med_noro_r)
med_noro_i_res     <- process_cols(med_noro_i)
med_noro_i_nodp_r_res     <- process_cols(med_noro_i_nodp_r)
med_noro_i_root_r_res     <- process_cols(med_noro_i_root_r)


med_nods_n_res     <- process_cols(med_nods_n)
med_nods_r_res     <- process_cols(med_nods_r)
med_nods_i_res     <- process_cols(med_nods_i)
med_nods_i_nodm_r_res     <- process_cols(med_nods_i_nodm_r)
med_nods_i_nodp_r_res     <- process_cols(med_nods_i_nodp_r)
med_nods_i_em21_n_res     <- process_cols(med_nods_i_em21_n)
med_nods_i_em22_n_res     <- process_cols(med_nods_i_em22_n)

med_all_r_root_res <- process_cols(med_all_r_root)
med_all_r_gall_res <- process_cols(med_all_r_gall)
med_all_r_nodp_res <- process_cols(med_all_r_nodp)
med_all_i_gano_res <- process_cols(med_all_i_gano)
med_all_i_noro_res <- process_cols(med_all_i_noro)
med_all_i_garo_res <- process_cols(med_all_i_garo)
med_all_i_gano_gall_res <- process_cols(med_all_i_gano_gall)
med_all_i_noro_nodp_res <- process_cols(med_all_i_noro_nodp)
med_all_i_garo_gall_res <- process_cols(med_all_i_garo_gall)
med_all_i_gano_nodp_res <- process_cols(med_all_i_gano_nodp)
med_all_i_noro_root_res <- process_cols(med_all_i_noro_root)
med_all_i_garo_root_res <- process_cols(med_all_i_garo_root)

r21_nods_n_res     <- process_cols(r21_nods_n)
r22_nods_n_res     <- process_cols(r22_nods_n)
nem_gall_r_res     <- process_cols(nem_gall_r)


# merging tables to get final table 
tables <- plyr::rbind.fill(medicago_final_table, nematode_final_table, rhizob21_final_table, rhizob22_final_table)
results_list <- lst(tables, med_gall_r_res, med_nodp_r_res, med_nodm_r_res, med_root_r_res, med_nods_r21_n_res, med_nods_r22_n_res,
                    med_garo_r_res, med_garo_i_res, med_garo_i_gall_r_res, med_garo_i_root_r_res,
                    med_gano_r_res, med_gano_i_res, med_gano_i_gall_r_res, med_gano_i_nodp_r_res,
                    med_noro_r_res, med_noro_i_res, med_noro_i_root_r_res, med_noro_i_nodp_r_res,
                    med_nods_n_res, med_nods_r_res, med_nods_i_res,
                    med_nods_i_nodm_r_res, med_nods_i_nodp_r_res, med_nods_i_em21_n_res, med_nods_i_em22_n_res,
                    med_all_r_root_res, med_all_r_gall_res, med_all_r_nodp_res,
                    med_all_i_garo_res, med_all_i_gano_res, med_all_i_noro_res,
                    med_all_i_garo_gall_res, med_all_i_gano_gall_res, med_all_i_noro_nodp_res,
                    med_all_i_garo_root_res, med_all_i_gano_nodp_res, med_all_i_noro_root_res,
                    r21_nods_n_res, r22_nods_n_res, nem_gall_r_res)


t <- Reduce(function(x, y) merge(x, y, all=TRUE, by = "gene"), results_list)



################################################################################
###################### ADDING SELECT ANNOTATIONS TO TABLE ######################
################################################################################

#Adding LEGOO annotations
LeGOO <- read.csv(file.path(dir_medic, "LeGOO-Currated-GeneNames.csv"))
colnames(LeGOO)[1] <- "gene"
LeGOO <- LeGOO %>% dplyr::select(gene, ACRONYM, DESCRIPTION, Publication)
colnames(LeGOO) <- c("gene", "legoo_acronym", "legoo_description", "legoo_publication")
LeGOO <- LeGOO %>% filter(gene != "#N/A" & gene != "UNDEF")
LeGOO <- LeGOO %>% group_by(gene) %>% summarise(across(everything(), ~toString(.)))
t <- left_join(t, LeGOO, by="gene")

#Adding NCR status to table based on LEGOO annotations
NCR_list <- read.csv(file.path(dir_medic, "MtrunA17r5.0-ANR-EGN-r1.9.b2g.gaf"), sep="\t", skip = 4, header = FALSE) %>%
  filter(str_detect(V10, "(NCR)")) %>% pull(V2)

NCR_list <- append(NCR_list, "MtrunA17_Chr7g0243321")
NCR_list <- append(NCR_list, "MtrunA17_Chr2g0306721")
NCR_list <- append(NCR_list, "MtrunA17_Chr1g0167541")

for (NCR in (t %>% filter(str_detect(legoo_acronym, "MtNCR")) %>% pull(gene))) {if (NCR %in% NCR_list) {} else {NCR_list <- append(NCR_list, NCR)}}

t <- t %>% rowwise() %>% mutate(NCR = ifelse(gene %in% NCR_list, TRUE, FALSE))



################################################################################
##################### COEFFICIENT OF VARIATION CALCULATION #####################
################################################################################

t <- t %>% mutate(all21_cv = sd(pick(s1, s5, s9,  s13, s17, s21, s25, s29, s33, s37, s41, s45), na.rm = TRUE) / rowMeans(pick(s1, s5, s9,  s13, s17, s21, s25, s29, s33, s37, s41, s45), na.rm = TRUE))
t <- t %>% mutate(all22_cv = sd(pick(s2, s6, s10, s14, s18, s22, s26, s30, s34, s38, s42, s46), na.rm = TRUE) / rowMeans(pick(s2, s6, s10, s14, s18, s22, s26, s30, s34, s38, s42, s46), na.rm = TRUE))
t <- t %>% mutate(all_cv = sd(pick(s1, s2, s5, s6, s9, s10, s13, s14, s17, s18, s21, s22, s25, s26, s29, s30, s33, s34, s37, s38, s41, s42, s45, s46), na.rm = TRUE) / rowMeans(pick(s1, s2, s5, s6, s9, s10, s13, s14, s17, s18, s21, s22, s25, s26, s29, s30, s33, s34, s37, s38, s41, s42, s45, s46), na.rm = TRUE))

t <- t %>% mutate(G21_cv = sd(pick(s13, s17, s21), na.rm = TRUE) / rowMeans(pick(s13, s17, s21), na.rm = TRUE)) %>% mutate(G21_cv_greater_all = ifelse(G21_cv > all_cv, TRUE, FALSE))
t <- t %>% mutate(G22_cv = sd(pick(s14, s18, s22), na.rm = TRUE) / rowMeans(pick(s14, s18, s22), na.rm = TRUE)) %>% mutate(G22_cv_greater_all = ifelse(G22_cv > all_cv, TRUE, FALSE))
t <- t %>% mutate(G_cv = sd(pick(s13, s14, s17, s18, s21, s22), na.rm = TRUE) / rowMeans(pick(s13, s14, s17, s18, s21, s22), na.rm = TRUE))

t <- t %>% mutate(Np21_cv = sd(pick(s25, s29, s33), na.rm = TRUE) / rowMeans(pick(s25, s29, s33), na.rm = TRUE)) %>% mutate(Np21_cv_greater_all = ifelse(Np21_cv > all_cv, TRUE, FALSE))
t <- t %>% mutate(Np22_cv = sd(pick(s26, s30, s34), na.rm = TRUE) / rowMeans(pick(s26, s30, s34), na.rm = TRUE)) %>% mutate(Np22_cv_greater_all = ifelse(Np22_cv > all_cv, TRUE, FALSE))
t <- t %>% mutate(Np_cv = sd(pick(s25, s26, s29, s30, s33, s34), na.rm = TRUE) / rowMeans(pick(s25, s26, s29, s30, s33, s34), na.rm = TRUE))

t <- t %>% mutate(Nm21_cv = sd(pick(s37, s41, s45), na.rm = TRUE) / rowMeans(pick(s37, s41, s45), na.rm = TRUE)) %>% mutate(Nm21_cv_greater_all = ifelse(Nm21_cv > all_cv, TRUE, FALSE))
t <- t %>% mutate(Nm22_cv = sd(pick(s38, s42, s46), na.rm = TRUE) / rowMeans(pick(s38, s42, s46), na.rm = TRUE)) %>% mutate(Nm22_cv_greater_all = ifelse(Nm22_cv > all_cv, TRUE, FALSE))
t <- t %>% mutate(Nm_cv = sd(pick(s37, s38, s41, s42, s45, s46), na.rm = TRUE) / rowMeans(pick(s37, s38, s41, s42, s45, s46), na.rm = TRUE))

t <- t %>% mutate(R21_cv = sd(pick(s1, s5, s9), na.rm = TRUE) / rowMeans(pick(s1, s5, s9), na.rm = TRUE))  %>% mutate(R21_cv_greater_all = ifelse(R21_cv > all_cv, TRUE, FALSE))
t <- t %>% mutate(R22_cv = sd(pick(s2, s6, s10), na.rm = TRUE) / rowMeans(pick(s2, s6, s10), na.rm = TRUE))  %>% mutate(R22_cv_greater_all = ifelse(R22_cv > all_cv, TRUE, FALSE))
t <- t %>% mutate(R_cv = sd(pick(s1, s2, s5, s6, s9, s10), na.rm = TRUE) / rowMeans(pick(s1, s2, s5, s6, s9, s10), na.rm = TRUE))

t <- t %>% mutate(N21_cv = sd(pick(s25, s29, s33, s37, s41, s45), na.rm = TRUE) / rowMeans(pick(s25, s29, s33, s37, s41, s45), na.rm = TRUE))
t <- t %>% mutate(N22_cv = sd(pick(s26, s30, s34, s38, s42, s46), na.rm = TRUE) / rowMeans(pick(s26, s30, s34, s38, s42, s46), na.rm = TRUE))
t <- t %>% mutate(N_cv = sd(pick(s25, s26, s29, s30, s33, s34, s37, s38, s41, s42, s45, s46), na.rm = TRUE) / rowMeans(pick(s25, s26, s29, s30, s33, s34, s37, s38, s41, s42, s45, s46), na.rm = TRUE))

t <- t %>% mutate(cv_qc_greater = (G21_cv_greater_all | G22_cv_greater_all | Np21_cv_greater_all | Np22_cv_greater_all | Nm21_cv_greater_all | Nm22_cv_greater_all | R21_cv_greater_all | R22_cv_greater_all))

cv_percent_table <- rbind(
  "Rhizob22" = t %>% filter(organism == "Rhizob21") %>% select(ends_with("greater_all")) %>% colSums(na.rm = TRUE) / length(t %>% filter(organism == "Rhizob21") %>% pull(gene)),
  "Rhizob21" = t %>% filter(organism == "Rhizob22") %>% select(ends_with("greater_all")) %>% colSums(na.rm = TRUE) / length(t %>% filter(organism == "Rhizob22") %>% pull(gene)),
  "Medicago" = t %>% filter(organism == "Medicago") %>% select(ends_with("greater_all")) %>% colSums(na.rm = TRUE) / length(t %>% filter(organism == "Medicago") %>% pull(gene)),
  "Nematode" = t %>% filter(organism == "Nematode") %>% select(ends_with("greater_all")) %>% colSums(na.rm = TRUE) / length(t %>% filter(organism == "Nematode") %>% pull(gene)))
cv_percent_table




################################################################################
############################### WRITING DATA OUT ###############################
################################################################################

setwd(dir_out)
save(t, file = "Diff_express_results_table.RData")
write.csv(t, file = "Diff_express_results_table.csv")
write.csv(cv_percent_table, file = "CV_Table.csv")
setwd(dir_main)




################################################################################
############################ GENERATING GENE LISTS #############################
################################################################################

gl_med_gall_r <- t %>% filter(med_gall_r.sig) %>% pull(gene)
gl_med_nodp_r <- t %>% filter(med_nodp_r.sig) %>% pull(gene)
gl_med_root_r <- t %>% filter(med_root_r.sig) %>% pull(gene)
gl_med_nods_r21_n <- t %>% filter(med_nods_r21_n.sig) %>% pull(gene)
gl_med_nods_r22_n <- t %>% filter(med_nods_r22_n.sig) %>% pull(gene)
gl_med_gano_r <- t %>% filter(med_gano_r.sig) %>% pull(gene)
gl_med_gano_i <- t %>% filter(med_gano_i.sig) %>% pull(gene)
gl_med_garo_r <- t %>% filter(med_garo_r.sig) %>% pull(gene)
gl_med_garo_i <- t %>% filter(med_garo_i.sig) %>% pull(gene)
gl_med_noro_r <- t %>% filter(med_noro_r.sig) %>% pull(gene)
gl_med_noro_i <- t %>% filter(med_noro_i.sig) %>% pull(gene)
gl_med_nods_r <- t %>% filter(med_nods_r.sig) %>% pull(gene)
gl_med_nods_n <-t %>% filter(med_nods_n.sig) %>% pull(gene)
gl_med_nods_i <-t %>% filter(med_nods_i.sig) %>% pull(gene)
gl_med_all_r <- t %>% filter(med_all_r_gall.sig) %>% pull(gene)
gl_nem_gall_r <-t %>% filter(nem_gall_r.sig) %>% pull(gene)
gl_r21_nods_n <-t %>% filter(r21_nods_n.sig) %>% pull(gene)
gl_r22_nods_n <-t %>% filter(r22_nods_n.sig) %>% pull(gene)

gl_med_gall_r_up <- t %>% filter(med_gall_r.sig & med_gall_r.ns_LFC > 0) %>% pull(gene)
gl_med_nodp_r_up <- t %>% filter(med_nodp_r.sig & med_nodp_r.ns_LFC > 0) %>% pull(gene)
gl_med_root_r_up <- t %>% filter(med_root_r.sig & med_root_r.ns_LFC > 0) %>% pull(gene)
gl_med_nods_r21_n_up <- t %>% filter(med_nods_r21_n.sig & med_nods_r21_n.ns_LFC > 0) %>% pull(gene)
gl_med_nods_r22_n_up <- t %>% filter(med_nods_r22_n.sig & med_nods_r22_n.ns_LFC > 0) %>% pull(gene)
gl_med_gano_r_up <- t %>% filter(med_gano_r.sig & med_gano_r.ns_LFC > 0) %>% pull(gene)
gl_med_garo_r_up <- t %>% filter(med_garo_r.sig & med_garo_r.ns_LFC > 0) %>% pull(gene)
gl_med_noro_r_up <- t %>% filter(med_noro_r.sig & med_noro_r.ns_LFC > 0) %>% pull(gene)
gl_med_nods_r_up <- t %>% filter(med_nods_r.sig & med_nods_r.ns_LFC > 0) %>% pull(gene)
gl_med_nods_n_up <-t %>% filter(med_nods_n.sig & med_nods_n.ns_LFC > 0) %>% pull(gene)
gl_nem_gall_r_up <-t %>% filter(nem_gall_r.sig & nem_gall_r.ns_LFC > 0) %>% pull(gene)
gl_r21_nods_n_up <-t %>% filter(r21_nods_n.sig & r21_nods_n.ns_LFC > 0) %>% pull(gene)
gl_r22_nods_n_up <-t %>% filter(r22_nods_n.sig & r22_nods_n.ns_LFC > 0) %>% pull(gene)

gl_med_gall_r_down <- t %>% filter(med_gall_r.sig & med_gall_r.ns_LFC < 0) %>% pull(gene)
gl_med_nodp_r_down <- t %>% filter(med_nodp_r.sig & med_nodp_r.ns_LFC < 0) %>% pull(gene)
gl_med_root_r_down <- t %>% filter(med_root_r.sig & med_root_r.ns_LFC < 0) %>% pull(gene)
gl_med_nods_r21_n_down <- t %>% filter(med_nods_r21_n.sig & med_nods_r21_n.ns_LFC < 0) %>% pull(gene)
gl_med_nods_r22_n_down <- t %>% filter(med_nods_r22_n.sig & med_nods_r22_n.ns_LFC < 0) %>% pull(gene)
gl_med_garo_r_down <- t %>% filter(med_garo_r.sig & med_garo_r.ns_LFC < 0) %>% pull(gene)
gl_med_gano_r_down <- t %>% filter(med_gano_r.sig & med_gano_r.ns_LFC < 0) %>% pull(gene)
gl_med_noro_r_down <- t %>% filter(med_noro_r.sig & med_noro_r.ns_LFC < 0) %>% pull(gene)
gl_med_nods_r_down <- t %>% filter(med_nods_r.sig & med_nods_r.ns_LFC < 0) %>% pull(gene)
gl_med_nods_n_down <-t %>% filter(med_nods_n.sig & med_nods_n.ns_LFC < 0) %>% pull(gene)
gl_nem_gall_r_down <-t %>% filter(nem_gall_r.sig & nem_gall_r.ns_LFC < 0) %>% pull(gene)
gl_r21_nods_n_down <-t %>% filter(r21_nods_n.sig & r21_nods_n.ns_LFC < 0) %>% pull(gene)
gl_r22_nods_n_down <-t %>% filter(r22_nods_n.sig & r22_nods_n.ns_LFC < 0) %>% pull(gene)





################################################################################
################################## QC REPORTS ##################################
################################################################################
source("Quality Check.R")



