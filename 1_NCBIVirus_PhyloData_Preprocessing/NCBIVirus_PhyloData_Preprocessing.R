# Loading Required Libraries
library(tidyverse); library(ape); library(seqinr); library(ggtree); library(anytime); library(phytools)

# Loading in Yellow Fever Virus Genome Metadata
yfv_metadata <- read.csv("1_NCBIVirus_PhyloData_Preprocessing/NCBIVirus_YFVAllSeqs_Metadata.csv")
yfv_metadata <- yfv_metadata %>%
  select(-Geo_Location) %>%
  mutate(Host = ifelse(Host == "Homo sapiens", "Human", ifelse(Host == "", "Unknown", "NHP"))) %>%
  mutate(Country = gsub(" ", "-", Country)) %>%
  mutate(Country = gsub("'", "-", Country))

# Creating Subsetting Index to Identify Sequences of Interest (LatAm origin in our case, and 
# precise date of collection)
countries <- c("Bolivia", "Brazil", "Ecuador", "Peru", "Trinidad-and-Tobago", "Venezuela")
country_subset_index <- yfv_metadata$Country %in% countries
date_subset_index <- nchar(yfv_metadata$Collection_Date) == 10
subset_index <- which(country_subset_index & date_subset_index)
yfv_metadata_subset <- yfv_metadata[subset_index, ]

# Creating Vector of Filenames
names <- c()
for (i in 1:nrow(yfv_metadata_subset)) {
  id <- yfv_metadata_subset[i, ]
  names <- c(names, paste(id$Accession, id$Country, id$Host, id$Collection_Date, sep = "_"))
}

# Splitting Up The Raw FASTA and Renaming, Then Concantenating the Seqs With Corrected File Names Back Together
all_seqs <- read.FASTA("1_NCBIVirus_PhyloData_Preprocessing/NCBIVirus_YFVAllSeqs.fasta")
for (i in 1:nrow(yfv_metadata_subset)) {
  index <- subset_index[i]
  temp_fasta <- all_seqs[index]
  new_name <- gsub("_", "|", names[i])
  names(temp_fasta) <- new_name
  write.FASTA(x = temp_fasta, file = paste0("1_NCBIVirus_PhyloData_Preprocessing/individual_fasta_seqs/", names[i], ".fasta"))
  if (i == 1) {
    all_seqs_proc <- temp_fasta
  } else {
    all_seqs_proc <- c(all_seqs_proc, temp_fasta)
  }
}
write.FASTA(x = all_seqs_proc, file = paste0("2_Multiple_Sequence_Alignment_Generation/NCBIVirus_YFVSelectSeqs.fasta"))

# Creating Tab-Delimited Version of Metadata for Annotation
names_vertical_bar <- gsub("_", "|", names)
yfv_metadata_subset$taxa <- names_vertical_bar 
yfv_metadata_subset <- yfv_metadata_subset %>%
  relocate(taxa)
readr::write_tsv(yfv_metadata_subset, "2_Multiple_Sequence_Alignment_Generation/NCBIVirus_YFVSelectSeqs_Metadata.tsv")


## Note that on windows PC, saving with "|" in file name doesn't work, so have to save with "_"
## but "|" used in fasta metadata line
# new_name <- sub("[^|]+", str_extract(names[i], "[^_]+"), names(temp_fasta)) # removing the version from accession id
