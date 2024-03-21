library("AnnotationForge")
library("stringr")
library("clusterProfiler")
library("GOfuncR")
library("tidyverse")





########################################
########################################
###!!!!RUN PRIMARY ANALYSIS FIRST!!!!###
########################################
########################################

# Directory handling
# Note: main directory must contain a directory for storing database files
enrich_dir_name <- "EnrichDatabase_withParents"
dir_enrichdb <- file.path(dir_main, enrich_dir_name)
dir.create(dir_enrichdb, showWarnings = FALSE)




# Importing gtf annotations to create databases for enrichment
setwd(dir_main)

gtf <- read.csv(file.path(dir_medic, "annotation.gtf"),
                col.names = c("SEQNAME", "SOURCE", "FEATURE", "START",
                              "END", "SCORE", "STRAND", "FRAME", "ATT"),
                sep = "\t")
gtf <- gtf %>% filter(FEATURE == "gene")

# following line basically uses a regular expression in the ATTribute column of
# gtf file to pull out the gene name that matches the gene names used in this
# analysis
gtf <- gtf %>%  mutate(GID = str_extract(ATT, '(?<=; gene_id )([^;]*)(?=)'))


# Importing EC annotations
b2g_annot <- read_tsv(file.path(dir_medic, "MtrunA17r5.0-ANR-EGN-r1.9.b2g.annot"), col_names = FALSE)
colnames(b2g_annot) <- c("GID", "ANNOT")
GO_b2g_annot <- b2g_annot[which(startsWith(b2g_annot$ANNOT, "GO:")), c(1,2)]
colnames(GO_b2g_annot) <- c("GID","GO")
EC_b2g_annot <- b2g_annot[which(startsWith(b2g_annot$ANNOT, "EC:")), c(1,2)]
colnames(EC_b2g_annot) <- c("GID","EC")
b2g_annot <- full_join(GO_b2g_annot, EC_b2g_annot, by = "GID")
rm(EC_b2g_annot)
rm(GO_b2g_annot)

#B2G Annotation Processing
gaf_annot <- read_tsv(file.path(dir_medic, "MtrunA17r5.0-ANR-EGN-r1.9.b2g.gaf"), col_names = FALSE, skip = 4)
colnames(gaf_annot) <- c("DB", "GID", "SYMBOL", "QUALIFIER", "GO", "DB:REFERENCE", "EVIDENCE", "WITH_OR_FROM", "ASPECT", "DB_OBJ_NAME", "DB_OBJ_SYN", "DB_OBJ_TYPE", "TAXON", "DATE", "ASSIGN_BY", "ANNOT_EXT", "GENE_PROD_ID")
gaf_annot <- gaf_annot %>% select(GID, everything())

## Makes an organism package
setwd(dir_enrichdb)

try(MtruncTolouseFile <- system.file("extdata", "MtruncTolouse_info.txt", package="AnnotationForge"))
try(MtruncTolouse <- read.table(MtruncTolouseFile, sep="\t"))

## Now prepare some data.frames
MTSym <- as.data.frame(gaf_annot %>% select("GID", "SYMBOL"))
MTChr <- as.data.frame(gtf %>% select("GID", "SEQNAME"))
colnames(MTChr) <- c("GID", "CHROMOSOME")

#Removing duplicated rows
MTSym <- dplyr::distinct(MTSym)
MTChr <- dplyr::distinct(MTChr)




if (file.exists("GO_term_database_with_parents_minus_duplicates.csv")){
  MTGO <- read.csv("GO_term_database_with_parents.csv") %>% select("GID", "GO", "EVIDENCE")
} else {
  MTGO <- as.data.frame(gaf_annot %>% select("GID", "GO", "EVIDENCE"))
  colnames(MTGO) <- c("GID", "GO", "EVIDENCE")
  
  MTGO <- dplyr::distinct(MTGO)
  
  # Populating parental terms in addition to child terms
  get_parents <- function(table){
    start_time <- Sys.time()
    print(paste(start_time,"Beginning to process input table"))
    out_list <- list()
    for (i in 1:nrow(table)){
      current_time <- Sys.time()
      if (i %% 1000 == 0) {print(paste(Sys.time(),"Processed",i,"rows of input table   ",current_time-start_time,"units have elapsed   rate is",i/as.numeric(current_time-start_time),"rows per unit"))}
      GO <- table[i,"GO"]
      GID <- table[i,"GID"]
      EVIDENCE <- table[i,"EVIDENCE"]
      parents <- get_parent_nodes(GO)
      i_table <- data.frame(GID = NA, GO = NA, EVIDENCE = NA)
      for (j in 1:nrow(parents)){
        go_term <- parents[j,2]
        i_table <- rbind(i_table, c(GO = go_term, GID = GID, EVIDENCE = EVIDENCE)) %>% drop_na()
      }
      out_list[[i]] <- i_table
    }
    out_table <- bind_rows(out_list)
    return(out_table)
  }
  
  MTGO_mod <- get_parents(MTGO)
  MTGO_mod2 <- MTGO_mod %>% rename(GO = GID, GID = GO)
  MTGO_mod2 <- MTGO_mod2 %>% relocate(GID)
  write.csv(MTGO_mod2, "GO_term_database_with_parents.csv", row.names = FALSE)
  MTGO_mod2 <- dplyr::distinct(MTGO_mod2)
  write.csv(MTGO_mod2, "GO_term_database_with_parents_minus_duplicates.csv", row.names = FALSE)
}

MTGO <- dplyr::distinct(MTGO)

if (!dir.exists(paste0("org.Mtruncatula",".eg.db"))){
makeOrgPackage(gene_info = MTSym, chromosome = MTChr, go = MTGO,
               version="0.0.1",
               maintainer="Some One <so@someplace.org>",
               author="Some One <so@someplace.org>",
               outputDir = ".",
               tax_id="3880",
               genus="Medicago",
               species="truncatula",
               goTable="go")

install.packages("./org.Mtruncatula.eg.db", type = "source", repos = NULL)} else {
  print("Database already created. Passing...")
  }

make_db_with_subset <- function(list, prefix){
  if (!dir.exists(paste0("org.Mtruncatula",prefix,".eg.db"))){
    subset_MTSym <- dplyr::distinct(MTSym[which(MTSym$GID %in% list),])
    subset_MTChr <- dplyr::distinct(MTChr[which(MTChr$GID %in% list),])
    subset_MTGO <- dplyr::distinct(MTGO[which(MTGO$GID %in% list),])
    species <- paste0("truncatula", prefix)
    db_name <- paste0("./org.Mtruncatula",prefix,".eg.db")
    makeOrgPackage(gene_info = subset_MTSym, chromosome = subset_MTChr, go = subset_MTGO,
                   version="0.0.1",
                   maintainer="Some One <so@someplace.org>",
                   author="Some One <so@someplace.org>",
                   outputDir = ".",
                   tax_id="3880",
                   genus="Medicago",
                   species = species,
                   goTable="go")
    install.packages(db_name, type = "source", repos = NULL)}
  else{print("Database already created. Passing...")}
}

make_db_with_subset(list = gl_med_exp_gall, prefix = "gall")
make_db_with_subset(list = gl_med_exp_nodp, prefix = "nodp")
make_db_with_subset(list = gl_med_exp_nods_em21, prefix = "nods21")
make_db_with_subset(list = gl_med_exp_nods_em22, prefix = "nods22")
make_db_with_subset(list = gl_med_exp_nodp_and_nodm, prefix = "nods")
make_db_with_subset(list = gl_med_exp_nodp_or_nodm, prefix = "nodsi")
make_db_with_subset(list = gl_med_exp_gall_and_nodp, prefix = "gano")
make_db_with_subset(list = gl_med_exp_gall_or_nodp, prefix = "ganoi")
make_db_with_subset(list = gl_med_exp_gall_and_root, prefix = "garo")
make_db_with_subset(list = gl_med_exp_gall_or_root, prefix = "garoi")
make_db_with_subset(list = gl_med_exp_nodp_and_root, prefix = "noro")
make_db_with_subset(list = gl_med_exp_nodp_or_root, prefix = "noroi")


library(org.Mtruncatula.eg.db)
library(org.Mtruncatulagall.eg.db)
library(org.Mtruncatulanodp.eg.db)
library(org.Mtruncatulanods21.eg.db)
library(org.Mtruncatulanods22.eg.db)
library(org.Mtruncatulanods.eg.db)
library(org.Mtruncatulanodsi.eg.db)
library(org.Mtruncatulagano.eg.db)
library(org.Mtruncatulaganoi.eg.db)
library(org.Mtruncatulagaro.eg.db)
library(org.Mtruncatulagaroi.eg.db)
library(org.Mtruncatulanoro.eg.db)
library(org.Mtruncatulanoroi.eg.db)

all_db <- org.Mtruncatula.eg.db
gall_db <- org.Mtruncatulagall.eg.db
nodp_db <- org.Mtruncatulanodp.eg.db
nods21_db <- org.Mtruncatulanods21.eg.db
nods22_db <- org.Mtruncatulanods22.eg.db
nods_db <- org.Mtruncatulanods.eg.db
nods_i_db <- org.Mtruncatulanodsi.eg.db
gano_db <- org.Mtruncatulagaro.eg.db
gano_i_db <- org.Mtruncatulagaroi.eg.db
garo_db <- org.Mtruncatulagaro.eg.db
garo_i_db <- org.Mtruncatulagaroi.eg.db

rm(gtf)
rm(gaf_annot)
rm(MTChr)
rm(MTGO)
rm(MTSym)
rm(b2g_annot)

########################
#GO ENRICHMENT ANALYSIS#
########################
dir_out <- file.path(dir_main, "out")
dir.create(dir_out, showWarnings = FALSE)
dir_out_enrich <- file.path(dir_out, "enrich_with_parents")
dir.create(dir_out_enrich, showWarnings = FALSE)
setwd(dir_out_enrich)


GO_results <- function(file_name, gene_list, orgdb, subdir = FALSE){
  GO_res <- function(file_name, gene_list, orgdb){
    GO_res <- enrichGO(gene = gene_list, OrgDb = orgdb, keyType = "GID", ont = "All", pvalueCutoff = 0.01, qvalueCutoff = 0.5, maxGSSize = 2000)
    df <- as.data.frame(GO_res)
    write.csv(df, file = paste(file_name, "_out.csv", sep = ""))
    write.table(df$ID, file = paste(file_name, "_out_list.txt", sep = ""), col.names = FALSE, row.names = FALSE, quote = FALSE)
    for (ontology in list("BP", "MF", "CC")){
      GO_res <- enrichGO(gene = gene_list, OrgDb = orgdb, keyType = "GID", ont = ontology, pvalueCutoff = 0.01, qvalueCutoff = 0.5, maxGSSize = 2000)
      df <- as.data.frame(GO_res)
      goterm_n <- length(rownames(as.data.frame(GO_res)))
      if (goterm_n > 0){
        if (goterm_n > 1){
          print(paste("Found",goterm_n,"enriched terms for", ontology))
          jpeg(file = paste(file_name, "_", ontology, "_tree", "_n", ".jpeg", sep = ""), 
               width = ifelse(goterm_n > 5, 
                              ifelse(goterm_n > 10,
                                     ifelse(goterm_n > 20,
                                            ifelse(goterm_n > 50, 
                                                   ifelse(goterm_n > 75,
                                                          ifelse(goterm_n > 100, 2400, 1800),
                                                          1200),
                                                   1080),
                                            720),
                                     480),
                              320),
               height = ifelse(goterm_n > 5,
                               ifelse(goterm_n > 10,
                                      ifelse(goterm_n > 20,
                                             ifelse(goterm_n > 50,
                                                    ifelse(goterm_n > 75,
                                                           ifelse(goterm_n > 100, 2400, 1800),
                                                           1800), 
                                                    1080), 
                                             720), 
                                      480), 
                               320))
          plot(goplot(GO_res, showCategory = goterm_n+10, color = "Count"))
          dev.off()
          jpeg(file = paste(file_name, "_", ontology, "_tree", "_padj", ".jpeg", sep = ""), 
               width = ifelse(goterm_n > 5, 
                              ifelse(goterm_n > 10,
                                     ifelse(goterm_n > 20,
                                            ifelse(goterm_n > 50, 
                                                   ifelse(goterm_n > 75,
                                                          ifelse(goterm_n > 100, 2400, 1800),
                                                          1200),
                                                   1080),
                                            720),
                                     480),
                              320),
               height = ifelse(goterm_n > 5,
                               ifelse(goterm_n > 10,
                                      ifelse(goterm_n > 20,
                                             ifelse(goterm_n > 50,
                                                    ifelse(goterm_n > 75,
                                                           ifelse(goterm_n > 100, 2400, 1800),
                                                           1800), 
                                                    1080), 
                                             720), 
                                      480), 
                               320))
          plot(goplot(GO_res, showCategory = goterm_n+10, color = "p.adjust"))
          dev.off()
          jpeg(file = paste(file_name, "_", ontology, "_dot", ".jpeg", sep = ""), width = 720, height = ifelse(goterm_n*32 < 640, 640, goterm_n*32))
          plot(dotplot(GO_res, show = goterm_n+10))
          dev.off()
        } else {
          print("1 enriched term found")}
      } else {
        print("No enrichment found")}
      try(dev.off(), silent = TRUE)
      try(dev.off(), silent = TRUE)
      return(df)
    }
  }
  starting_wd <- getwd()
  if (subdir != FALSE) {subdirectory <- subdir} else {subdirectory <- file_name}
  dir_GO_res <- file.path(starting_wd, subdirectory)
  dir.create(dir_GO_res, showWarnings = FALSE)
  setwd(dir_GO_res)
  df <- GO_res(file_name, gene_list, orgdb)
  return(df)
  setwd(starting_wd)
}




go_gall_r      <- GO_results("gall_r", gl_med_gall_r, gall_db)
go_gall_r_up   <- GO_results("gall_r_up", gl_med_gall_r_up, gall_db, "gall_r")
go_gall_r_down <- GO_results("gall_r_down", gl_med_gall_r_down, gall_db, "gall_r")

go_nodp_r      <- GO_results("nodp_r", gl_med_nodp_r, nods_i_db)
go_nodp_r_up   <- GO_results("nodp_r_up", gl_med_nodp_r_up, nods_i_db, "nodp_r")
go_nodp_r_down <- GO_results("nodp_r_down", gl_med_nodp_r_down, nods_i_db, "nodp_r")


go_nods_r21_n       <- GO_results("nods_r21_n", gl_med_nods_r21_n, nods_i_db)
go_nods_r21_n_up    <- GO_results("nods_r21_n_up", gl_med_nods_r21_n_up, nods_i_db, "nods_r21_n")
go_nods_r21_n_down  <- GO_results("nods_r21_n_down", gl_med_nods_r21_n_down, nods_i_db, "nods_r21_n")

go_nods_r22_n      <- GO_results("nods_r22_n", gl_med_nods_r22_n, nods_i_db)
go_nods_r22_n_up   <- GO_results("nods_r22_n_up", gl_med_nods_r22_n_up, nods_i_db, "nods_r22_n")
go_nods_r22_n_down <- GO_results("nods_r22_n_down", gl_med_nods_r22_n_down, nods_i_db, "nods_r22_n")

go_nods_r      <- GO_results("nods_r", gl_med_nods_r, nods_i_db)
go_nods_r_up   <- GO_results("nods_r_up", gl_med_nods_r_up, nods_i_db, "nods_r")
go_nods_r_down <- GO_results("nods_r_down", gl_med_nods_r_down, nods_i_db, "nods_r")

go_nods_n      <- GO_results("nods_n", gl_med_nods_n, nods_i_db, "nods_n")
go_nods_n_up   <- GO_results("nods_n_up", gl_med_nods_n_up, nods_i_db, "nods_n")
go_nods_n_down <- GO_results("nods_n_down", gl_med_nods_n_down, nods_i_db, "nods_n")

go_nods_i      <- GO_results("nods_i", gl_med_nods_i, nods_i_db, "nods_i")



go_gano_r      <- GO_results("gano_r", gl_med_gano_r, gano_i_db)
go_gano_r_up   <- GO_results("gano_r_up", gl_med_gano_r_up, gano_i_db, "gano_r")
go_gano_r_down <- GO_results("gano_r_down", gl_med_gano_r_down, gano_i_db, "gano_r")
go_gano_i      <- GO_results("gano_i", gl_med_gano_i, gano_i_db)

go_garo_r      <- GO_results("garo_r", gl_med_garo_r, garo_i_db)
go_garo_r_up   <- GO_results("garo_r_up", gl_med_garo_r_up, garo_i_db, "garo_r")
go_garo_r_down <- GO_results("garo_r_down", gl_med_garo_r_down, garo_i_db, "garo_r")
go_garo_i      <- GO_results("garo_i", gl_med_garo_i, garo_i_db)


