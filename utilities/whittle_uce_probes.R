whittle_uce_probes <- function(uce_list_file,probe_fasta_file,basename,file_type) {
  print("This code will whittle down your probe_fasta_file by loci that are present in your uce_list_file")
  print("and that are not 'problematic' (i.e. paralagous) within or between taxa (if you are using the output")
  print("from comparing_monolithic_UCE_fastas/monolithic.sh. To run this code:")
  print("whittle_uce_probes(uce_list_file,probe_fasta_file,basename)")
  print("where uce_list_file is an output file from comparing_monolithic_UCE_fastas/monolithic.sh or a list of uce loci,")
  print("probe_fasta_file is the output probes file from the phyluce pipeline,")
  print("basename is the name of the taxa that you designed your final probeset across. You can use a place holder name") 
  cat('if you did not design the probes e.g. "whatever"\n')
  print("file_type is the type of file used to whittle i.e. output from monolithic.sh or a list of uce loci")
  print("e.g.")
  cat('whittle_uce_probes("/Users/alanaalexander/Dropbox/beetles/grey_whittled_probes/output_matrix_99.txt","/Users/alanaalexander/Dropbox/beetles/grey_whittled_probes/Adephaga_11Kv1.fasta","Pterostichus.1","monolithic")\n')
  print("e.g.")
  cat('whittle_uce_probes("C:\\Users\\Alana\\Dropbox\\beetles\\grey_whittled_probes\\uce_loci_from_baca_2017.txt","C:\\Users\\Alana\\Dropbox\\beetles\\grey_whittled_probes\\Coleoptera-UCE-1.1K-v1","whatever","file_list")\n')
    
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
  
  if (file_type=="monolithic") {
    # Reading in the monolithic output
    temp <- read.table(uce_list_file,stringsAsFactors=FALSE,header=TRUE)
  
    # Filtering out loci that have between_taxa_problems
    temp <- filter(temp,is.na(between_taxa_problem))
  
    # Filtering out loci that have problems within taxa
    for (i in grep("_longest_base_genome",names(temp))) {
      temp <- filter(temp,temp[,i]!="problem_within" | is.na(temp[,i]))
    }  
  
    # Taking just the uce-locus names for the basename in our probe_fasta_file
    temp <- select(temp,grep(basename,names(temp))[1])
  
    # Removing any is.na rows
    temp <- filter(temp,!is.na(temp[,1]))
  
 } else {
    temp <- read.table(uce_list_file,stringsAsFactors=FALSE,header=TRUE)
 }   

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
