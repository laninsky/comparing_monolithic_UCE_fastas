file_list <- list.files(pattern=".summarized",recursive=TRUE)

temp_output <- as.matrix(read.table(file_list[1]))
problem_loci_coords <- which(temp_output[,(which(temp_output[1,]=="problem_locus"))]=="Y")
problem_loci <- temp_output[problem_loci_coords,]
temp_output <- temp_output[-problem_loci_coords,]

first_column <- c("rosetta_locus",rep(NA,(dim(temp_output)[1]-1)))
temp_output <- cbind(first_column,temp_output)
