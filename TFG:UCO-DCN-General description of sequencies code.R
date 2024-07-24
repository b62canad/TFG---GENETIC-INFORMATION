# Install and load

install.packages("data.table")
library(data.table)

install.packages("seqinr")
library(seqinr)

library(dplyr)
library(refseqR)
library(rentrez)

# Base composition of nuclear DNA

  #Read sequence file

    xf <- readDNAStringSet("FASTA_SIN_SCAFFOLD")
  
  #Create DNAStringSet object to collect DNA sequencies
  
    sequences <- DNAStringSet(xf)
    
  #Calculate frequencies of nucleotids in sequences
    
    nucleotide_freq <- alphabetFrequency(sequences)
    
  #Transform into data frame
  
    n_f <- as.data.frame(nucleotide_freq)
  
  #See and print number of nucleotids
    
    total_adeninas <- sum(n_f["A"])
    total_timinas <- sum(n_f["T"])
    total_guaninas <- sum(n_f["G"])
    total_citocinas <- sum(n_f["C"])
    
    print(paste("Adeninas =", total_adeninas))
    print(paste("Timinas =", total_timinas))
    print(paste("Guaninas =", total_guaninas))
    print(paste("Citocinas =", total_citocinas))
    
  #Create document
  
    Base_composition <- data.frame(
      "Adeninas" = c("13202884"), "Timinas" = c("12610299"), "Guaninas" = c("9813416"), "Citocinas" = c("8246151")
    )
    
    write.table(Base_composition, file = "Base_composition.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Total and mean exon length
    
  #Read gff file
    
    archivo_gff <- "anotacion_filtrada.gff"
    dt <- fread(archivo_gff, header = FALSE, col.names = c('secuencia', 'fuente', 'tipo', 'inicio', 'fin', 'puntuacion', 'cadena', 'fase', 'atributos'))
    
  #Filter exons
    
    exones <- dt[tipo == 'exon']
    
  #Obtain length keeping in mind the chain direction
    
    exones[, longitud := abs(fin - inicio) + 1]
  
  #See the length of all exons
    
    longitud_total <- sum(exones$longitud)
    
  #Obtain exon mean length
    
    longitud_media <- mean(exones$longitud)
    
  #Print total lenght
    
    cat("Longitud total de exones:", longitud_total, "\n")
    
  #Print exon mean length
    
    cat("Longitud media de exones:", longitud_media, "\n")
    
  #Create a file with those data
    
    Longitud_exones <- data.frame(
      "Longitud_total" = c("68786334"), "Longitud_media" = c("290.4338")
    )
    
    write.table(Longitud_exones, file = "exon_length.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Mean number of exon and intron
    
    # Read FASTA file
    
    fasta_file <- read.fasta(file = "FASTA_SIN_SCAFFOLD" , seqtype = "DNA", as.string = TRUE)
    
    # Extract headline
    
    headers <- sapply(getAnnot(fasta_file), function(x) {
      if (grepl("(XP|YP|NP)_\\d+\\.\\d+", x)) {
        return(regmatches(x, regexpr("(XP|YP|NP)_\\d+\\.\\d+", x)))
      } else {
        return(NA)
      }
    })
    
    # Filter headline that aren't empty
    
    protein_ids <- headers[!is.na(headers)]
    
    # Transform results to dataframe
    
    encabezados <- data.frame(protein_ids)
    
    # Save results
    
    write.table(encabezados, "encabezados secuencias codificantes", quote = FALSE, row.names = FALSE)
    
    #Define function
    characterizeTable <- function(targets) {
      
      # Start vector for protein ID and number of exon
      protein_id <- vector("character")
      exon <- vector("integer") 
      
      # Extract protein ID and number of exon
      for(xi in as.character(targets)) {
        xplink <- entrez_link(dbfrom = "protein", id = xi, db = "gene")
        genesummary = entrez_summary(db = "gene", id = xplink$links[1])
        
        if(length(genesummary$genomicinfo$exoncount) != 0) {
          exon = c(exon, genesummary$genomicinfo$exoncount)  
        } else {
          exon = c(exon, 0)  
        }
        
        protein_id = c(protein_id, xi)
      }
      
      # Add column with number of intron 
      t1 <- data.frame(Protein_ID = protein_id, Exon_Count = exon)
      
      t1 <- t1 %>% mutate(Intron_count = exon - 1)
      
      t1
      
    }
    
    # Like "targets" is a vector we need to put a "$" symbol
    targets <- read.table("encabezados secuencias codificantes")
    
    # "n" is the number of exons that we have to study
    exones <- characterizeTable(targets$V1[1:n])
    
    write.table(exones, file = "exones_final.txt", sep = "\t", row.names= FALSE)
    
    #Obtain exon mean:
    
    #Create object:
    exon_final_count <- read.table ("exones_final.txt", header = TRUE, sep = "\t")
    
    # Show first dataframe lines and see column names 
    head(exon_final_count)
    
    # Transform the column to numeric
    exon_final_count$Exon_Count <- as.numeric(as.character(exon_final_count$Exon_Count))
    
    # Calculate mean of an especific column
    media_exon <- mean(exon_final_count$Exon_Count, na.rm = TRUE)      
    
    print(media_exon)     
    
    #Obtain intron mean:
    
    # Transform the column to numeric
    exon_final_count$Intron_count <- as.numeric(as.character(exon_final_count$Intron_count))
    
    # Calculate mean of an especific column
    media_intron <- mean(exon_final_count$Intron_count, na.rm = TRUE)      
    
    print(media_intron) 
    