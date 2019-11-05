# comparing_monolithic_UCE_fastas v0.2: whittle_uce_probes.R

whittle_uce_probes <- function(uce_list_file,probe_fasta_file,file_type,taxa_filter,basename) {
  print("This code will whittle down your probe_fasta_file by loci that are present in your uce_list_file")
  print("(and that are not 'problematic' i.e. paralagous within or between taxa if you are using the output")
  print("from comparing_monolithic_UCE_fastas/monolithic.sh). To run this code:")
  print("whittle_uce_probes(uce_list_file,probe_fasta_file,file_type,basename)")
  print("where uce_list_file is an output file from comparing_monolithic_UCE_fastas/monolithic.sh or a list of uce loci,")
  print("probe_fasta_file is the output probes file from the phyluce pipeline,")
  print("file_type is the type of file used to whittle i.e. output from monolithic.sh = 'monolithic' or a list of uce loci = 'uce_list'")
  print('taxa_filter is an option that allows you to filter just to the UCE loci recovered across all taxa = "Y" (filter based on this) or "N" (do not filter)')
  print("basename is the name of the taxa that you designed your final probeset across if chosing the monolithic option (can be left blank if using uce_list") 
  print("e.g.")
  cat('whittle_uce_probes("/Users/alanaalexander/Dropbox/beetles/grey_whittled_probes/output_matrix_99.txt","/Users/alanaalexander/Dropbox/beetles/grey_whittled_probes/Adephaga_11Kv1.fasta","monolithic","Pterostichus.1")\n')
  print("e.g.")
  cat('whittle_uce_probes("C:\\Users\\Alana\\Dropbox\\beetles\\grey_whittled_probes\\uce_loci_from_baca_2017.txt","C:\\Users\\Alana\\Dropbox\\beetles\\grey_whittled_probes\\Coleoptera-UCE-1.1K-v1","uce_list")\n')
    
  if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse')
  
  # Reading in probe_fasta_file
  temp <- read_table(probe_fasta_file,col_names = FALSE)
  
  # "One-lining" sequence associated with each fasta header
  number_rows <- as.numeric(as.matrix(temp %>% filter(grepl(">",X1)) %>% count())[1,1])*2
  outputmatrix <- matrix(NA,ncol=1,nrow=number_rows)
  outputmatrix[seq(1,dim(outputmatrix)[1],2),1] <- as.matrix(temp %>% filter(grepl(">",X1)))
  
  outputmatrixpos <- 2
  i <- 2
  tempseq <- NULL
  while (i <= dim(temp)[1]) {
      if (!(grepl(">",temp[i,1]))) {
          tempseq <- paste(tempseq,temp[i,1],sep="")
      } else {
        outputmatrix[outputmatrixpos,1] <- tempseq
        tempseq <- NULL
        outputmatrixpos <- outputmatrixpos+2
      }
      i <- i+1
  }
  outputmatrix[(dim(outputmatrix)[1]),] <- tempseq
  
  if (file_type=="monolithic") {
    # Reading in the monolithic output
    temp <- read_table2(uce_list_file)
  
    # Filtering out loci that have between_taxa_problems
    temp <- filter(temp,is.na(between_taxa_problem))
  
    # Getting columns associated with "_longest_base_genome"
    columns_to_check <- which(grepl("_longest_base_genome",names(temp)))
    
    # Using these columns to eliminate rows where there are within taxa problems
    # and those that have NAs (taxa are missing) if taxa_filter selected
    for (i in columns_to_check) {
      if (taxa_filter=="N") {
        temp <- temp %>% filter(temp[,i]!="problem_within" | is.na(temp[,i]))
      } else {
        temp <- temp %>% filter(temp[,i]!="problem_within" & !is.na(temp[,i]))
      }  
    }  
  
    # Taking just the uce-locus names for the basename in our probe_fasta_file
    temp <- select(temp,grep(basename,names(temp),fixed=TRUE)[1])
  
    # Removing any is.na rows
    temp <- filter(temp,!is.na(temp[,1]))
  
 } else {
    temp <- read_table2(uce_list_file,col_names = FALSE)
 }   

  # Getting the lines of the probe file that are found in our list of "keeper" loci
  # First, extracting uce names from the probe file
  headerlines <- NULL
  headerlines <- outputmatrix[(grep(">",outputmatrix[,1])),1]
  headerlines <- gsub("_p.*","",gsub(">","",headerlines))
  
  # Then finding the headerlines that are in our list of "kept" uce loci
  keepheaderlines <- which(headerlines %in% as.matrix(temp))
  
  # Finding the number of our "kept" uce loci that are in the headerlines
  kept_loci <- length(which(as.matrix(temp) %in% headerlines)
    
  # Getting the position of the headerlines in the original output file
  keepheaderlines <- (keepheaderlines*2)-1
  # Adding in their associate sequences
  keepheaderlines <- c(keepheaderlines,(keepheaderlines+1))
  keepheaderlines <- sort(keepheaderlines)
  
  # subsetting these desired probes from the total probes file
  outputmatrix <- matrix(outputmatrix[keepheaderlines,],ncol=1)
  
  write.table(outputmatrix,"whittled_UCE_probes.fasta",col.names=FALSE,row.names=FALSE,quote=FALSE)
  
  print(paste((dim(outputmatrix)[1]/2)," probes targetting ",kept_loci," loci have been written out to ",getwd(),"/whittled_UCE_probes.fasta",sep=""))
}
