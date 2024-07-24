#Install and load

install.packages("BiocManager")
BiocManager::install("Biostring")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.18")

library(Biostrings)

install.packages("stringr")
library(stringr)

# Read document

readDNAStringSet("GCF_000331145.1_ASM33114v1_cds_from_genomic.fna")

# FASTA without scaffold

  #Filter sequences by "NC_" due to the fact that scaffold don't have "NC_"

patron_completo <- readDNAStringSet("GCF_000331145.1_ASM33114v1_cds_from_genomic.fna")

patron_completo [grep("NC_", names(patron_completo))]

NC <- patron_completo [grep("NC_", names(patron_completo))]

  #Create new FASTA without scaffolds

nueva_ruta_fasta <- "FASTA_SIN_SCAFFOLD"

writeXStringSet(NC, format = "fasta", filepath = nueva_ruta_fasta)

readDNAStringSet("FASTA_SIN_SCAFFOLD")

# FASTA only with chromosomes

FASTA_SIN_SCAFFOLD <- readDNAStringSet("FASTA_SIN_SCAFFOLD")

FASTA_Cr <- FASTA_SIN_SCAFFOLD[grep("NC_011163.1", names(FASTA_SIN_SCAFFOLD), invert = TRUE)]

  #Create new FASTA only with chromosomes sequences

nueva_ruta_fasta_Cr <- "FASTA_Cr"

writeXStringSet(FASTA_Cr, format = "fasta", filepath = nueva_ruta_fasta_Cr)

readDNAStringSet("FASTA_Cr")

# FASTA only with chloroplast

  # Filter sequences by "NC" where "NC_011163.1" is chloroplast

patron_completo_C <- readDNAStringSet("GCF_000331145.1_ASM33114v1_cds_from_genomic.fna")

patron_completo_C [grep("NC_011163.1", names(patron_completo_C))]

NC_C <- patron_completo_C [grep("NC_011163.1", names(patron_completo_C))]

  #Create new FASTA only with chloroplast sequences

nueva_ruta_fasta_C <- "FASTA_C"

writeXStringSet(NC_C, format = "fasta", filepath = nueva_ruta_fasta_C)

readDNAStringSet("FASTA_C")

#Chi square 

X_obs <- c(xi)
X_esp <- c(P)
chisq.test(X_obs, p=X_esp)

    #Example for Valine:

      Val_obs <- c(161574, 123664, 218089, 437911)
      Val_esp <- c(0.25, 0.25, 0.25, 0.25)
      chisq.test(Val_obs, p=Val_esp)
