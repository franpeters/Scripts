library(TCGAbiolinks)
library(tidyverse)
library(DBI)
library(SummarizedExperiment)

#need to set the location of this to the files source pane location
datasets =read.csv("diseases.csv")
database <- "es_database.sqlite"
db <- dbConnect(RSQLite::SQLite(), database)
setwd("~/") # important for tcga biolinks

dbExecute(db,"CREATE TABLE `tcga_result` (
          id INT,
          sample TEXT,
          ensembl_gene_id TEXT,
          count FLOAT,
          type TEXT,
          tcga_id TEXT,
          PRIMARY KEY (id,sample,type),
          FOREIGN KEY (ensembl_gene_id) REFERENCES genes (ensembl_gene_id),
          FOREIGN KEY (tcga_id) REFERENCES tcga_meta (tcga_id)
);"
)

dbExecute(db,"CREATE TABLE `tcga_metadata` (
          tcga_id TEXT,
          patient TEXT,
          sample TEXT,
          code TEXT,
          definition TEXT,
          is_ffpe TEXT,
          primary_diagnosis TEXT,
          classification TEXT,
          progression TEXT,
          prior_treatment TEXT,
          PRIMARY KEY (tcga_id,patient,sample)
);"
)

i=1
datasets <- datasets$Disease.ID
for(i in 2:length(datasets)){
  query.exp <- GDCquery(
    project = paste("TCGA",datasets[i],sep="-"),
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification",
    workflow.type = "STAR - Counts"
  )
  
  GDCdownload(
    query = query.exp,
    files.per.chunk = 100
  )
  
  brca.exp <- GDCprepare(
    query = query.exp
  )
  
  # Which samples are Primary Tumor
  samples.primary.tumour <- brca.exp$barcode[brca.exp$shortLetterCode == "TP"]
  
  # which samples are solid tissue normal
  #
  samples.primary.tissue.normal <- brca.exp$barcode[brca.exp$shortLetterCode == "NT"]
  
  tpm <- assay(brca.exp,"tpm_unstrand")
  fpkm <- assay(brca.exp,"fpkm_unstrand")
  raw <- assay(brca.exp,"unstranded")
  
  tpm_tumour <- tpm[,samples.primary.tumour]
  tpm_tumour <- data.frame(tpm_tumour)
  tpm_tumour$ensembl_gene_id <- rownames(tpm_tumour)
  tpm_tumour$ensembl_gene_id<- gsub("\\..*","",tpm_tumour$ensembl_gene_id)
  tpm_tumour$type <- "tpm_count"
  fpkm_tumour <- fpkm[,samples.primary.tumour]
  fpkm_tumour <- data.frame(fpkm_tumour)
  fpkm_tumour$ensembl_gene_id <- rownames(fpkm_tumour)
  fpkm_tumour$ensembl_gene_id<- gsub("\\..*","",fpkm_tumour$ensembl_gene_id)
  fpkm_tumour$type <- "fpkm_count"
  raw_tumour <- raw[,samples.primary.tumour]
  raw_tumour <- data.frame(raw_tumour)
  raw_tumour$ensembl_gene_id <- rownames(raw_tumour)
  raw_tumour$ensembl_gene_id<- gsub("\\..*","",raw_tumour$ensembl_gene_id)
  raw_tumour$type <- "raw_count"
  if(ncol(raw_tumour)==2&ncol(tpm_tumour)==2&ncol(fpkm_tumour)==2){
    
  }else{
    # 
    # raw_tumour_gathered <- gather(raw_tumour,sample,count,colnames(raw_tumour)[1]:colnames(raw_tumour)[ncol(raw_tumour)-1],-ensembl_gene_id)
    # raw_tumour_gathered$type <- "raw"
    tumour_one <- rbind(tpm_tumour,fpkm_tumour,raw_tumour)
    tumour_one_gathered <- gather(tumour_one,sample,count,colnames(tumour_one)[1]:colnames(tumour_one)[ncol(tumour_one)-2],-ensembl_gene_id,-type)
    tumour_one_gathered$tcga_id <- datasets[i]
    tumour_one_gathered$id <- seq.int(nrow(tumour_one_gathered))
    
    dbAppendTable(db,"tcga_result",tumour_one_gathered)
    rm(tumour_one_gathered)
  }
  
 
  tpm_tissue.normal <- tpm[,samples.primary.tissue.normal]
  tpm_tissue.normal <- data.frame(tpm_tissue.normal)
  tpm_tissue.normal$ensembl_gene_id <- rownames(tpm_tissue.normal)
  tpm_tissue.normal$ensembl_gene_id<- gsub("\\..*","",tpm_tissue.normal$ensembl_gene_id)
  tpm_tissue.normal$type <- "tpm_count"
  fpkm_tissue.normal <- fpkm[,samples.primary.tissue.normal]
  fpkm_tissue.normal <- data.frame(fpkm_tissue.normal)
  fpkm_tissue.normal$ensembl_gene_id <- rownames(fpkm_tissue.normal)
  fpkm_tissue.normal$ensembl_gene_id<- gsub("\\..*","",fpkm_tissue.normal$ensembl_gene_id)
  fpkm_tissue.normal$type <- "fpkm_count"
  raw_tissue.normal <- raw[,samples.primary.tissue.normal]
  raw_tissue.normal <- data.frame(raw_tissue.normal)
  raw_tissue.normal$ensembl_gene_id <- rownames(raw_tissue.normal)
  raw_tissue.normal$ensembl_gene_id<- gsub("\\..*","",raw_tissue.normal$ensembl_gene_id)
  raw_tissue.normal$type <- "raw_count"
  
  if(ncol(raw_tissue.normal)==2&ncol(tpm_tissue.normal)==2&ncol(fpkm_tissue.normal)==2){
    
  }else{
  
  tissue.normal_one <- rbind(tpm_tissue.normal,fpkm_tissue.normal,raw_tissue.normal)
  tissue.normal_one_gathered <- gather(tissue.normal_one,sample,count,colnames(tissue.normal_one)[1]:colnames(tissue.normal_one)[ncol(tissue.normal_one)-2],-ensembl_gene_id,-type)
  tissue.normal_one_gathered$tcga_id <- datasets[i]
  
  tissue.normal_one_gathered$id <- seq.int(nrow(tissue.normal_one_gathered))
  
  dbAppendTable(db,"tcga_result",tissue.normal_one_gathered)
  rm(tissue.normal_one_gathered)
  }
  df <- data.frame(brca.exp$patient,
                   brca.exp$barcode,
                   brca.exp$shortLetterCode,
                   brca.exp$definition,
                   brca.exp$is_ffpe,
                   brca.exp$primary_diagnosis,
                   brca.exp$classification_of_tumor,
                   brca.exp$progression_or_recurrence,
                   brca.exp$prior_treatment)
  df$tcga_id <-datasets[i]
  
  colnames(df) <- c("patient","sample","code","definition","is_ffpe","primary_diagnosis","classification","progression","prior_treatment","tcga_id")
  df$sample <- gsub("-", "\\.",df$sample)
  dbAppendTable(db,"tcga_metadata",df)

  }
