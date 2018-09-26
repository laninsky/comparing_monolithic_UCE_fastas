whittle_uce_probes <- function(monolithic_output_name,probe_fasta_file,basename) {
  print("This code will whittle down your probe_fasta_file by loci that are present in your monolithic_output_name file")
  print("and that are not 'problematic' (i.e. paralagous) within or between taxa. To run this code:")
  print("whittle_uce_probes(monolithic_output_name,probe_fasta_file,basename)")
  print("where monolithic_output_name is an output file from comparing_monolithic_UCE_fastas/monolithic.sh")
  print("probe_fasta_file is the output probes file from the phyluce pipeline")
  print("and basename is the name of the taxa that you designed your final probeset across")
  print("e.g.")
  cat('whittle_uce_probes("/Users/alanaalexander/Dropbox/beetles/grey_whittled_probes/output_matrix_99.txt","/Users/alanaalexander/Dropbox/beetles/grey_whittled_probes/Adephaga_11Kv1.fasta","Pterostichus.1")\n')
  
  if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
  
  # Reading in probe_fasta_file
  temp <- readLines(probe_fasta_file)
  
  # "One-lining" sequence associated with each fasta header
  outputmatrix <- matrix(NA,ncol=1,nrow=length(grep(">",temp))*2)
  outputmatrix[seq(1,dim(outputmatrix)[1],2),1] <- temp[(grep(">",temp))]
  
  outputmatrixpos <- 2
  i <- 2
  tempseq <- NULL
  while (i <= length(temp)) {
      if (!(grepl(">",temp[i]))) {
          tempseq <- paste(tempseq,temp[i],sep="")
      } else {
        outputmatrix[outputmatrixpos,1] <- tempseq
        tempseq <- NULL
        outputmatrixpos <- outputmatrixpos+2
      }
      i <- i+1
  }
  outputmatrix[(dim(outputmatrix)[1]),] <- tempseq
  
  # Reading in the monolithic output
  temp <- read.table(monolithic_output_name,stringsAsFactors=FALSE,header=TRUE)
  
  # Filtering out loci that have between_taxa_problems
  temp <- filter(temp,is.na(between_taxa_problem))
  
  # Filtering out loci that have problems within taxa
  for (i in grep("_longest_base_genome",names(temp))) {
    temp <- filter(temp,temp[,i]!="problem_within")
  }  
  
  # Taking just the uce-locus names for the basename in our probe_fasta_file
  temp <- select(temp,grep(basename,names(temp))[1])
  
  # Removing any is.na rows
  temp <- filter(temp,!is.na(temp[,1]))
  
  # Getting the lines of the probe file that are found in our list of "keeper" loci
  headerlines <- NULL
  
  for (i in 1:dim(temp)[1]) {
    headerlines <- c(headerlines,which(grepl(temp[i,1],outputmatrix[,1])))
    if ((i %% 100)==0) {
     print(paste("We are ",round((i/dim(temp)[1])*100),"% through the file",sep=""))
    }
  }  
  
  headerlines <- unique(headerlines)
  headerlines <- sort(c(headerlines,(headerlines+1)))
  outputmatrix <- matrix(outputmatrix[headerlines,],ncol=1)
  
  write.table(outputmatrix,"whittled_UCE_probes.fasta",col.names=FALSE,row.names=FALSE,quote=FALSE)
  
  print(paste((dim(outputmatrix)[1]/2)," probes targetting ",(dim(temp)[1])," loci have been written out to ",getwd(),"/whittled_UCE_probes.fasta",sep=""))
}  
  
  
  
  
