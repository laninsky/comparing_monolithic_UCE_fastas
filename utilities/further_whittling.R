# comparing_monolithic_UCE_fastas v0.2: further_whittling.R

further_whittling <- function(monolithic_file,probe_fasta_file,basename,no_of_loci) {
  print("This code will whittle down your probe_fasta_file by loci that are present in your monolithic file,")
  print("that are not 'problematic' i.e. paralagous within or between taxa, that are found in all taxa")
  print("characterized in the monolithic_file, and that are among the no_of_loci longest average loci in the")
  print("monolithic_file. To run this code:")
  print("further_whittling(monolithic_file,probe_fasta_file,basename,no_of_loci)")
  print("where monolithic_file is an output file from comparing_monolithic_UCE_fastas/monolithic.sh,")
  print("probe_fasta_file is the output probes file from the phyluce pipeline (or whittled_UCE_probes.fasta,")
  print("if you've already run whittle_uce_probes.R)")
  print("basename is the name of the taxa that you designed your final probeset across, and")
  print("no_of_loci is the desired number of loci that probes_fasta_file will be filtered to")
  print("e.g.")
  cat('further_whittling("/Users/alanaalexander/Dropbox/beetles/grey_whittled_probes/output_matrix_99.txt","/Users/alanaalexander/Dropbox/beetles/grey_whittled_probes/Adephaga_11Kv1.fasta","Pterostichus.1",6000)\n')
  print("e.g.")
  cat('further_whittling("output_matrix.txt","whittled_UCE_probes.fasta","liocanthydrus",6716)\n')
  
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
  
  # Reading in the monolithic output
  temp <- read_table2(monolithic_file)
    
  # Filtering out loci that have between_taxa_problems
  temp <- filter(temp,is.na(between_taxa_problem))
    
  # Getting columns associated with "_longest_base_genome"
  columns_to_check <- which(grepl("_longest_base_genome",names(temp)))
    
  # Using these columns to eliminate rows where there are within taxa problems
  # and those that have NAs (taxa are missing) 
  for (i in columns_to_check) {
    temp <- temp %>% filter(temp[,i]!="problem_within" & !is.na(temp[,i]))
  }  
  
  # Removing loci that were not characterized in our base genome of choice
  removerows <- which(is.na(temp[,(grep(basename,names(temp))[1])]))
  temp <- temp[-removerows,]
  
  # Sorting loci based on average loci length and taking the the top no_of_loci loci
  # First need to ensure the columns are numerical, then taking columns with length and
  # calculating the mean
  temp <- temp %>% mutate_at(vars(contains("_length")),funs(as.numeric)) %>% 
    mutate(av_length=select(.,contains("_length")) %>% rowMeans()) %>% 
    arrange(desc(av_length)) %>% slice(1:no_of_loci)

  # Taking just the uce-locus names for the basename in our probe_fasta_file
  temp <- select(temp,grep(basename,names(temp),fixed=TRUE)[1])
    
  # Getting the lines of the probe file that are found in our list of "keeper" loci
  # First, extracting uce names from the probe file
  headerlines <- NULL
  headerlines <- outputmatrix[(grep(">",outputmatrix[,1])),1]
  headerlines <- gsub("_p.*","",gsub(">","",headerlines))
  
  print(paste("Original probes file had ",(dim(outputmatrix)[1]/2)," probes targetting ",length(unique(headerlines))," loci",sep=""))
  
  # Then finding the headerlines that are in our list of "kept" uce loci
  keepheaderlines <- which(headerlines %in% as.matrix(temp))
  
  # Finding the number of our "kept" uce loci that are in the headerlines
  kept_loci <- length(which(as.matrix(temp) %in% headerlines[keepheaderlines]))
  
  # Getting the position of the headerlines in the original output file
  keepheaderlines <- (keepheaderlines*2)-1
  # Adding in their associate sequences
  keepheaderlines <- c(keepheaderlines,(keepheaderlines+1))
  keepheaderlines <- sort(keepheaderlines)
  
  # subsetting these desired probes from the total probes file
  outputmatrix <- matrix(outputmatrix[keepheaderlines,],ncol=1)
  
  write.table(outputmatrix,"further_whittled_UCE_probes.fasta",col.names=FALSE,row.names=FALSE,quote=FALSE)
  
  print(paste((dim(outputmatrix)[1]/2)," probes targetting ",kept_loci," loci have been written out to ",getwd(),"/further_whittled_UCE_probes.fasta. Information on these loci has been written out to filtered_output_matrix.txt",sep=""))
  
  # Generating a summary of the loci that have been kept
  temp <- read_table2(monolithic_file)
  headerlines <- NULL
  headerlines <- outputmatrix[(grep(">",outputmatrix[,1])),1]
  headerlines <- unique(gsub("_p.*","",gsub(">","",headerlines)))
  keeprows <- which(as.matrix(temp[,(grep(basename,names(temp))[1])]) %in% headerlines)
  temp <- temp[keeprows,]
  
  write.table(temp,"filtered_output_matrix.txt",col.names=TRUE,row.names=FALSE,quote=FALSE)
}
