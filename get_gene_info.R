#### get all genes' sequence and pdb file and uniprot

library(httr)
library(jsonlite)
library(xml2)
library(tidyverse)
library(biomaRt)
ensembl_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
atr <-listAttributes(ensembl_mart)

all_genes <- getBM( attributes = c("ensembl_gene_id","chembl","description","pdb"), mart = ensembl_mart)

all_genes_2 <- getBM( attributes = c("ensembl_gene_id","hgnc_symbol","uniprot_gn_id","entrezgene_id"), mart = ensembl_mart)

merge <- merge(all_genes,all_genes_2,"ensembl_gene_id",all.x=T,all.y=T)

# all_genes_3 <- getBM( attributes = c("ensembl_gene_id","peptide_location"), mart = ensembl_mart)
# 
# all_genes_final <- merge(merge,all_genes_3,"ensembl_gene_id",all.x=T,all.y=T)
# 
setwd("C:/Users/2634371/OneDrive - AlmacGroup/Documents/Projects/Everything_scape/Database")

merge$hgnc_symbol[merge$hgnc_symbol==""]  <- merge$ensembl_gene_id

merge <- merge %>% distinct(hgnc_symbol,ensembl_gene_id,.keep_all = T)
write.csv(merge,"merge_genes.csv")
drop_na <- merge %>% drop_na(uniprot_gn_id) 
uniprots <- drop_na$uniprot_gn_id
genes <- drop_na$ensembl_gene_id

res <- data.frame()
i=1
print(uniprots[i])
print(i)
requestURL <- paste("https://www.ebi.ac.uk/proteins/api/proteins/",uniprots[i],sep="")
r <- GET(requestURL, accept("application/json"))
stop_for_status(r)
content <- content(r)
test <- content$features
t <- test[[1]]
t <- t[c("type","category","begin","end","description")]
t <- data.frame(t)
features <- t

for(x in 2:length(test)){
  t <- test[[x]]
  t <- t[c("type","category","begin","end","description")]
  names(t) <- c("type","category","begin","end","description")
  if(length(is.na(t$description))==0){
    t <- data.frame(t[-5])
    t$description <- "None"
  }else{
    t <- data.frame(t)
  }
  
  colnames(t) <- c("type","category","begin","end","description")
  features <- rbind(features,t)
}
features$gene <- genes[i]
sequences <- data.frame(content$sequence)
sequences$gene <- genes[i]
i=2


for(i in 61855:length(uniprots)){
  print(uniprots[i])
  if(uniprots[i]!=""){
  print(genes[i])
  print(i)
  Sys.sleep(2) 
  requestURL <- paste("https://www.ebi.ac.uk/proteins/api/proteins/",uniprots[i],sep="")
  r <- GET(requestURL, accept("application/json"))
  content <- content(r)
  test <- content$features
  t <- test[[1]]
  t <- t[c("type","category","begin","end","description")]
  if (any(is.na(names(t)))==TRUE){

    new_list <- list(type = t[[1]], 
                    category = t[[2]],
                    begin = t[[3]],
                    end = t[[4]],
                    description = t[[5]])
    for(x in 1:5){
      if(length(new_list[[x]])==0){
        new_list[[x]] = "None"
      }
    }
    t <- data.frame(new_list)
  }else{
    t <- data.frame(t)
  }
  
  if (length(test)>1){
    for(x in 2:length(test)){
      t <- test[[x]]
      t <- t[c("type","category","begin","end","description")]
      
      if (any(is.na(names(t)))==TRUE){
        
        new_list <- list(type = t[[1]], 
                         category = t[[2]],
                         begin = t[[3]],
                         end = t[[4]],
                         description = t[[5]])
        for(x in 1:5){
          if(length(new_list[[x]])==0){
            new_list[[x]] = "None"
          }
        }
        t <- data.frame(new_list)
      }else{
        t <- data.frame(t)
      }
      t$gene <- genes[i]
      features <- rbind(features,t)
    }
  }
  if(length(content$sequence)>0){
    sequence <- data.frame(content$sequence)
    sequence$gene <- genes[i]
    sequences <- rbind(sequences,sequence)
  }

  }
}
sequences <- data.frame(sequences)

saveRDS(features,"genes_features.rds")
saveRDS(sequences,"genes_sequences.rds")
write.csv(features,"genes_features.csv")
write.csv(sequences,"genes_sequences.csv")


