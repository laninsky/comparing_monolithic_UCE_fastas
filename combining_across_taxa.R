file_list <- list.files(pattern=".summarized",recursive=TRUE)
taxa_list <- unlist(strsplit(file_list,"/",fixed=TRUE))
taxa_list <- taxa_list[seq((length(taxa_list)/length(file_list)),length(taxa_list),(length(taxa_list)/length(file_list)))]
taxa_list <- unlist(strsplit(taxa_list,"_blast.txt.summarized")) 

length_list <- paste(taxa_list,"_length",sep="")
base_list <- paste(taxa_list,"_longest_base",sep="")


temp_output <- as.matrix(read.table(file_list[1]))
problem_loci_coords <- which(temp_output[,(which(temp_output[1,]=="problem_locus"))]=="Y")
problem_loci <- temp_output[problem_loci_coords,]
temp_output <- temp_output[-problem_loci_coords,]

first_column <- c("rosetta_locus",rep(NA,(dim(temp_output)[1]-1)))
temp_output <- cbind(first_column,temp_output)
