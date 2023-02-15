fetch_alphafold <-  function (uniprot,  path = NULL) {
  file_name <- paste("AF-",uniprot,"-F1-model_v4.pdb",sep="")
  url <- paste("https://alphafold.ebi.ac.uk/files/",file_name,sep="")
  download.file(url,file_name)
   if (!is.character(uniprot)) {
     stop("pdb must be character string. e.g.'1bna'")
   }
   if (is.null(path)) {
     path <- getwd()
   }

     readLines(file_name)
   
  }
