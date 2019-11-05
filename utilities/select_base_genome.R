# comparing_monolithic_UCE_fastas v0.2: select_base_genome.R

select_base_genome <- function(output_matrix_file) {
  cat("This script will take an output_matrix_file that results from monolithic.sh and\n")
  cat("use it to provide several comparisons of the best performing base genome\n")
  cat("To use:\n")
  cat("select_base_genome(output_matrix_file) e.g.\n")
  cat('select_base_genome("/Users/alanaalexander/Dropbox/beetles/grey/Ecol_evol/18Apr2018_blast_results/output_matrix_95.txt")')
  
  # Loading required libraries
  if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse')

  # Reading in output_matrix_file
  temp <- read_table2(output_matrix_file,col_names = TRUE)
  
  # Getting taxa_names
  taxa_names <- gsub("_longest_base_genome","",names(temp)[grepl("_longest_base_genome",names(temp))])

  # Filtering out loci that are problematic between/within base genomes
  temp <- filter(temp,is.na(between_taxa_problem))
  
  # Filtering out loci that have problems within taxa or that have no info for some taxa
  columns_to_check <- which(grepl("_longest_base_genome",names(temp)))
  
  for (i in columns_to_check) {
    temp <- temp %>% filter(temp[,i]!="problem_within" & !is.na(temp[,i]))
  }  
  
  # Base genome names
  base_genome_names <- names(temp)[grepl("ucelocus.txt$",names(temp))]
  
  # Number of good UCE loci that were found in all taxa where each base genome
  # had probes present
  print(paste(dim(temp)[1]," UCE loci were found in all taxa",sep=""))
  print("Of these, the following number of loci were present in each base genome (larger is better):")
  for (i in 1:length(base_genome_names)) {
    print(paste((sum(!is.na(temp[,i])))," loci are present in ",base_genome_names[i],sep=""))
  }
  print("")
  
  # Which base genomes recovered the greatest number of taxa
  recording_taxa_base <- names(temp)[grepl("ucelocus.txt",names(temp))][(!(names(temp)[grepl("ucelocus.txt",names(temp))] %in% base_genome_names))]
  print(paste(length(taxa_names)," taxa are present in the data set",sep=""))
  print("Of these, the following number of average taxa were matched to by each base genome over all UCE loci (larger is better):")
  for (i in 1:length(recording_taxa_base)) {
    av_count <- sum(str_count(as.matrix(temp[,which(names(temp)==recording_taxa_base[i])]),","),na.rm=TRUE)/dim(temp)[1]
    print(paste(av_count," average taxa are recovered over all UCE loci in ",base_genome_names[i],sep=""))
  }
  print("")
  
  # Which base genomes result in the longest recovered sequences across taxa
  print("The following proportion of loci had each base genome lead to the longest alignment (larger is better):")
  base_genome_names  <- gsub("/ucelocus.txt","",base_genome_names)
  longest_base_genomes <- which(grepl("_longest_base_genome$",names(temp)))
  for (i in 1:length(base_genome_names)) {
    av_count <- sum(str_count(as.matrix(temp[,longest_base_genomes]),fixed(base_genome_names[i])),na.rm=TRUE)/(dim(temp)[1]*length(taxa_names))
    print(paste(av_count," proportion of loci*taxa combinations where the following base genome gave longest alignment: ",base_genome_names[i],sep=""))
  }
  print("")
}  
