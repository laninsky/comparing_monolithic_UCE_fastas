probe_fasta_file <- "/Users/alanaalexander/Dropbox/beetles/grey_whittled_probes/Adephaga_11Kv1.fasta"

whittle_uce_probes <- function(target_loci,probe_fasta_file) {
  
  
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
  
  
  
  seq_or_no[x
  
  
