# Install and load

  install.packages("BiocManager")
  BiocManager::install("Biostring")
  
  library(Biostrings)
  library(stringr)
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(BiocManager)
  
#Read sequence file
  
  archivo <- "FASTA_SIN_SCAFFOLD"
  
  secuencias <- readDNAStringSet(archivo)

# ATG as initiation codon

  #Save first codon

    primeras_tres_letras <- subseq(secuencias, start = 1, end = 3)

    conteo <- table(primeras_tres_letras)
    
  #Select those that doesn't start with "ATG"
    
    primeras_tres_letras[!grepl(paste("ATG"), as.character(primeras_tres_letras))]
    
    filtrado_general <- primeras_tres_letras[!grepl(paste("ATG"), as.character(primeras_tres_letras))]
    
    myDNASet <- filtrado_general
    
    names(myDNASet) <- sub(".(XP|YP|NP)([0-9]+\\.[0-9]+)_.", "\\1_\\2", names(myDNASet))
    
  #Avoid factor conversion 
    
    dnaDataFrame <- data.frame(
      Name = names(myDNASet),
      Sequence = as.character(myDNASet),
      stringsAsFactors = FALSE
    )
  
  #Print a document with separated information
    
    dnaDataFrame <- dnaDataFrame %>%
      mutate(Content = str_extract_all(Name, "\\[.*?\\]") %>% map_chr(paste, collapse=" "))
    
    max_splits <- max(str_count(dnaDataFrame$Content, "\\["))
    
    dnaDataFrame <- dnaDataFrame %>%
      separate(Content, into = paste0("Col", 1:max_splits), sep = "\\] \\[", fill = "right", extra = "merge")
    
    dnaDataFrame <- dnaDataFrame %>%
      mutate(Name = str_replace_all(Name, "\\s*\\[.*?\\]", ""))
    
    write.table(dnaDataFrame, file = "NO_ATG.txt", sep = "\t", row.names= FALSE)

# STOP codon

  #Save last codon
    
    ultimo_codon <- subseq(secuencias, start = width(secuencias) - 2, end = width(secuencias))
    
    conteo3 <- table(ultimo_codon)

    #Select those that doesn't finish with "TAA", "TAG", "TGA"
    
    ultimo_codon[!grepl(paste("TAA", "TAG", "TGA"), as.character(ultimo_codon))]
    
    final_codon <- c("TAA", "TAG", "TGA")
    
        #(Need to do this because STOP has 3 combinations):
    
        STOP <- paste(final_codon, collapse = "|")
        
    #Continue with selection
        
    filtrado_general_STOP <- ultimo_codon[!grepl(pattern, as.character(ultimo_codon))]
        
    myDNASet_STOP <- filtrado_general_STOP
        
    names(myDNASet_STOP) <- sub(".(XP|YP|NP)([0-9]+\\.[0-9]+)_.", "\\1_\\2", names(myDNASet_STOP))
        
  #Avoid factor conversion
    
    dnaDataFrame_STOP <- data.frame(
      Name = names(myDNASet_STOP),
      Sequence = as.character(myDNASet_STOP),
      stringsAsFactors = FALSE
    )

    #Print a document with separated information

    dnaDataFrame_STOP <- dnaDataFrame_STOP %>%
      mutate(Content = str_extract_all(Name, "\\[.*?\\]") %>% map_chr(paste, collapse=" "))
    
    max_splits_STOP <- max(str_count(dnaDataFrame_STOP$Content, "\\["))
    
    dnaDataFrame_STOP <- dnaDataFrame_STOP %>%
      separate(Content, into = paste0("Col", 1:max_splits), sep = "\\] \\[", fill = "right", extra = "merge")
    
    dnaDataFrame_STOP <- dnaDataFrame_STOP %>%
      mutate(Name = str_replace_all(Name, "\\s*\\[.*?\\]", ""))
    
    write.table(dnaDataFrame_STOP, file = "NO_STOP.txt", sep = "\t", row.names= FALSE)

#Second codon
    
    #Save second codon
    
    segundo_codon <- subseq(secuencias, start = 4, end = 6)
    
    conteo2 <- table(segundo_codon)
    
    write.table(conteo2, file = "Segundo_codon.txt", sep = "\t", row.names= FALSE)
    