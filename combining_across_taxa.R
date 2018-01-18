file_list <- list.files(pattern=".summarized",recursive=TRUE)
taxa_list <- unlist(strsplit(file_list,"/",fixed=TRUE))
taxa_list <- taxa_list[seq((length(taxa_list)/length(file_list)),length(taxa_list),(length(taxa_list)/length(file_list)))]
taxa_list <- unlist(strsplit(taxa_list,"_blast.txt.summarized")) 

length_list <- paste(taxa_list,"_length",sep="")
longbase_list <- paste(taxa_list,"_longest_base_genome",sep="")
temp_output <- as.matrix(read.table(file_list[1]))
base_list <- temp_output[1,1:((which(temp_output[1,]=="max_length"))-1)]
pivot_col <- which(temp_output[1,]=="max_length")-1

output_matrix <- t(matrix(c(base_list,length_list,longbase_list,"between_taxa_problem")))
problem_taxa <- t(matrix(c("taxa","base_genome","uce_locus")))

for (i in file_list) { #1A
   taxa <- unlist(strsplit(i,"/",fixed=TRUE))
   taxa <- taxa[length(taxa)]
   taxa <- unlist(strsplit(taxa,"_blast.txt.summarized"))
   output_taxa <- which(grepl(taxa,output_matrix[1,])==TRUE)
   temp_output <- as.matrix(read.table(i)) 
      for (j in 2:(dim(temp_output)[1])) { #2A
         k <- 1
         while (k <= pivot_col) { #10A
            if(!(is.na(temp_output[j,k]))) { #11A
               if(!(temp_output[j,k] %in% output_matrix[,k])) { #3A
                  if(temp_output[j,which(temp_output[1,]=="problem_locus")]=="N") { #4A
                     temp_row <- c(temp_output[j,1:pivot_col],rep(NA,(length(output_matrix[1,])-pivot_col)))
                     temp_row[output_taxa[1]] <- temp_output[j,which(temp_output[1,]=="max_length")]
                     temp_row[output_taxa[2]] <- temp_output[j,which(temp_output[1,]=="which_base_gives_max" )]                                        
                     output_matrix <- rbind(output_matrix,temp_row)
                     k <- pivot_col+1
                  } else { #4AB this one is for problem loci - need to hold them somewhere until the end
                     break
                     for (m in 1:pivot_col) { #20A
                        if(!(is.na(temp_output[j,m]))) { #21A
                           if(length(which(problem_taxa[,1]==taxa & problem_taxa[,2]==output_matrix[1,m]))<1) { #22A
                              temp_row <- c(taxa,output_matrix[1,m],temp_output[j,m])
                              problem_taxa <- rbind(problem_taxa,temp_row)
                            } else { #22AB
                              break
                            } #22B
                         } #21B 
                      } #20A        
                     k <- pivot_col+1
                  } #4B
               } else { #3AB this is for when 
                  break
               } #3B
               break
           } #11B
           k <- k + 1 
        } #10B   
    } #2B       
} #1B       
    


problem_loci_coords <- which(temp_output[,(which(temp_output[1,]=="problem_locus"))]=="Y")
problem_loci <- temp_output[problem_loci_coords,]
temp_output <- temp_output[-problem_loci_coords,]

first_column <- c("rosetta_locus",rep(NA,(dim(temp_output)[1]-1)))
temp_output <- cbind(first_column,temp_output)
