file_list <- list.files(pattern=".summarized",recursive=TRUE)
taxa_list <- unlist(strsplit(file_list,"/",fixed=TRUE))
taxa_list <- taxa_list[seq((length(taxa_list)/length(file_list)),length(taxa_list),(length(taxa_list)/length(file_list)))]
taxa_list <- unlist(strsplit(taxa_list,"_blast.txt.summarized")) 

length_list <- paste(taxa_list,"_length",sep="")
longbase_list <- paste(taxa_list,"_longest_base_genome",sep="")
temp_output <- as.matrix(read.table(file_list[1]))
base_list <- temp_output[1,1:((which(temp_output[1,]=="max_length"))-1)]
pivot_col <- which(temp_output[1,]=="max_length")-1

output_matrix <- t(matrix(c(base_list,length_list,longbase_list,base_list,"between_taxa_problem")))
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
               if(temp_output[j,which(temp_output[1,]=="problem_locus")]=="N") { #4A
                  locus_present <- "No"
                  match_rows <- NULL
                  for (p in 1:pivot_col) {
                     if(!(is.na(temp_output[j,p]))) {
                        if(temp_output[j,p] %in% output_matrix[,p]) {
                           match_rows <- unique(sort(c(match_rows,which(output_matrix[,p] %in% temp_output[j,p])))) 
                           locus_present <- "Yes"
                        }   
                     }   
                  }   
                  if(locus_present=="No") { #3A
                     temp_row <- c(temp_output[j,1:pivot_col],rep(NA,(length(output_matrix[1,])-pivot_col)))
                     temp_row[output_taxa[1]] <- temp_output[j,which(temp_output[1,]=="max_length")]
                     temp_row[output_taxa[2]] <- temp_output[j,which(temp_output[1,]=="which_base_gives_max" )]
                     record_taxa <- which(output_matrix[1,] %in% output_matrix[1,which(!(is.na(temp_output[j,1:pivot_col])))])
                     record_taxa <- record_taxa[((length(record_taxa)/2)+1):length(record_taxa)]
                     temp_row[record_taxa] <- paste(taxa,",",sep="")
                     output_matrix <- rbind(output_matrix,temp_row)
                     k <- pivot_col+1
                  } else { #3AB 
                     temp_row <- output_matrix[match_rows[1],]
                     for (m in 1:length(temp_row)) {
                        for (p in match_rows) {
                           if (is.na(temp_row[m])) {
                              temp_row[m] <- output_matrix[p,m]
                           } else {
                              if (!(is.na(output_matrix[p,m]))) {
                                 if (!(temp_row[m]==output_matrix[p,m])) {
                                    temp_row[m] <- paste(temp_row[m],",",output_matrix[p,m],sep="")
                                    temp_row[m] <- gsub("NA","",temp_row[m])
                                    temp_row[length(temp_row)] <- "Y"
                                 }
                              }
                           }
                        }
                     }
                     for (m in 1:pivot_col) {                     
                        if (is.na(temp_row[m])) {
                           temp_row[m] <- temp_output[j,m]
                        } else {
                           if (!(is.na(temp_output[j,m]))) {
                              if (!(temp_row[m]==temp_output[j,m])) {
                                 temp_row[m] <- paste(temp_row[m],",",temp_output[j,m],sep="")
                                 temp_row[m] <- gsub("NA","",temp_row[m])
                                 temp_row[length(temp_row)] <- "Y"
                              }
                           }
                        }
                     }
                     temp_row[output_taxa[1]] <- temp_output[j,which(temp_output[1,]=="max_length")]
                     temp_row[output_taxa[2]] <- temp_output[j,which(temp_output[1,]=="which_base_gives_max" )]
                     record_taxa <- which(output_matrix[1,] %in% output_matrix[1,which(!(is.na(temp_output[j,1:pivot_col])))])
                     record_taxa <- record_taxa[((length(record_taxa)/2)+1):length(record_taxa)]
                     temp_row[record_taxa] <- paste(temp_row[record_taxa],taxa,",",sep="")
                     temp_row[record_taxa] <- gsub("NA","",temp_row[record_taxa])
                     output_matrix <- output_matrix[-match_rows,]
                     output_matrix <- rbind(output_matrix,temp_row)
                     k <- pivot_col+1         
                  } #3B
               } else { #4AB this one is for problem loci - need to hold them somewhere until the end
                  for (m in 1:pivot_col) { #20A
                     if(!(is.na(temp_output[j,m]))) { #21A
                        if(length(which(problem_taxa[,1]==taxa & problem_taxa[,2]==output_matrix[1,m]))<1) { #22A
                           temp_row <- c(taxa,output_matrix[1,m],temp_output[j,m])
                           problem_taxa <- rbind(problem_taxa,temp_row)
                         } else { #22AB - what to do if the length is more than 1
                           if (!(temp_output[j,m] %in% problem_taxa[(which(problem_taxa[,1]==taxa & problem_taxa[,2]==output_matrix[1,m])),3])) { #30A
                              temp_row <- c(taxa,output_matrix[1,m],temp_output[j,m])
                              problem_taxa <- rbind(problem_taxa,temp_row)
                           } #30B   
                         } #22B
                      } #21B 
                   } #20A        
                  k <- pivot_col+1
              } #4B
           } #11B
           k <- k + 1 
        } #10B   
    } #2B       
} #1B       
    
for (j in 2:(dim(problem_taxa)[1])) {
   problem_col <- (grep(problem_taxa[j,2],output_matrix[1,]))[1]
   problem_row <- which(output_matrix[,problem_col]==problem_taxa[j,3])
   taxa_col <- grep(problem_taxa[j,1],output_matrix[1,])
   output_matrix[problem_row,taxa_col] <- "problem_within"
   problem_col <- (grep(problem_taxa[j,2],output_matrix[1,]))[2]
   output_matrix[problem_row,problem_col] <- paste(output_matrix[problem_row,problem_col],problem_taxa[j,1],"_prob_w,",sep="")
   output_matrix[problem_row,problem_col] <- gsub("NA","",output_matrix[problem_row,problem_col])
}   
