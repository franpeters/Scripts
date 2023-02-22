## modify MaxQuant XML file quick script -- needs finalising

library("optparse") ## to make this script command line
library(tidyverse)
library(xml2)
## the inputs will be 
#annotation file -> same format as msstats requires 
#file type of raw files (.raw,.wiff etc.)
#folder path the raw files are in

raw_folder = "/work/fpeters/proteomics/esophageal/PXD020903/raw/"
annotation <- read.csv("annotation.csv")
file_type = ".raw"
bool = FALSE
fraction = 32767

annotation$filepaths <- paste0(raw_folder,annotation[,1],file_type)

x <- read_xml("new.xml")

for (i in 1:nrow(annotation)){
  x %>%
    xml_find_all(".//filePaths") %>%
    xml_add_child("string",annotation$filepaths[i]) 
  x %>%
    xml_find_all(".//experiments") %>%
    xml_add_child("string",annotation$Condition[i]) 
  x %>%
    xml_find_all(".//fractions") %>%
    xml_add_child("short",fraction) 
  x %>%
    xml_find_all(".//ptms") %>%
    xml_add_child("boolean",bool) 
  x %>%
    xml_find_all(".//paramGroupIndices") %>%
    xml_add_child("string",0) 
  x %>%
    xml_find_all(".//referenceChannel") %>%
    xml_add_child("string","") 
}

write_xml(x,"r_xml.xml")
