#Install and load

  install.packages("BiocManager")
  BiocManager::install("Biostring")

  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(version = "3.18")

  library(Biostrings)

  install.packages("stringr")
  library(stringr)
  
  install.packages("dplyr")
  library(dplyr)
  
#Read FASTA and aminoacids
  
  xm <-  readBStringSet(filepath ="FASTA")
  Aminoacids <- read.table(file = "Aminoacids.txt", header = TRUE)
  
#Create an empty vector
  
  resultado_acumulado <- integer(0)
  
#Obtain nucleotid frequencies
  
  frecuencia_nucleotidos <- function(i){
    trinucleotideFrequency(DNAString(toString(xm[i])), step = 3)
  }
  
#Obtain all the nucleotids frequencies
  
  for (i in seq_along(xm)) {
    resultado_acumulado <- c(resultado_acumulado,frecuencia_nucleotidos(i)) 
  }
  
#Obtain result in dataframe
  
  resultado_dataframe <- data.frame(codon = names(resultado_acumulado), number = as.numeric(resultado_acumulado))
  
#Group codon and number of aparition
  
  resultado_dataframe %>% group_by(codon) %>% count(number)
  
#Obtain result
  
  resultado <- tapply(resultado_dataframe$number, resultado_dataframe$codon, sum)
  
  resultado_suma <- data.frame(codon = names(resultado), number = as.numeric(resultado))
  
  resultado_final <- merge(Aminoacids, resultado_suma, by.x = "codon", all.x = TRUE)
  
#Order result
  
  resultado_final_ordenado <- arrange(resultado_final, aa)
  
#Save as .rda and .txt format
  
  write.table(resultado_final_ordenado, file = "Aminoacids ssp.rda", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  write.table(resultado_final_ordenado, file = "Aminoacids ssp.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

  

