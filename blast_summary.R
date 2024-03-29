# comparing_monolithic_UCE_fastas v0.2: blast_summary.R

for (i in list.files(pattern="_blast.txt",recursive=TRUE)) {
  temp <- as.matrix(read.csv(i, header=FALSE))
  base_genomes <- unique(sort(rbind(temp[,1],temp[,2])))
  output_matrix <- t(matrix(c(base_genomes,"max_length","which_base_gives_max","problem_locus")))
  for (j in 1:(dim(temp)[1])) { #6A
     # if temp[j,3] already in table
     if (temp[j,3] %in% output_matrix[1:(dim(output_matrix)[1]),which(output_matrix[1,]==temp[j,1])]) { #4A
        # what to do if both [j,3] and [j,4] are already there
        if (temp[j,4] %in% output_matrix[1:(dim(output_matrix)[1]),which(output_matrix[1,]==temp[j,2])]) { #5A
          # what to do if [j,3] or [j,4] are already in the matrix but on multiple rows (i.e. problematic)
          if (length((which(output_matrix[,which(output_matrix[1,]==temp[j,1])]==temp[j,3])))!=1 | length(which(output_matrix[,which(output_matrix[1,]==temp[j,2])]==temp[j,4]))!=1) { #IF1
            output_matrix[which(output_matrix[,which(output_matrix[1,]==temp[j,1])]==temp[j,3]),which(output_matrix[1,]=="problem_locus")] <- "Y"
            output_matrix[which(output_matrix[,which(output_matrix[1,]==temp[j,2])]==temp[j,4]),which(output_matrix[1,]=="problem_locus")] <- "Y"
            temp_row <- t(matrix(NA,nrow=dim(output_matrix)[2]))
            temp_row[1,which(output_matrix[1,]=="problem_locus")] <- "Y"
            temp_row[1,which(output_matrix[1,]==temp[j,1])] <- temp[j,3]
            temp_row[1,which(output_matrix[1,]==temp[j,2])] <- temp[j,4]
            # What to do if both base genomes have equally long loci
            if(as.numeric(temp[j,5])==as.numeric(temp[j,6])) { #2A
              temp_row[1,which(output_matrix[1,]=="max_length")] <- as.numeric(temp[j,5])
              temp_row[1,which(output_matrix[1,]=="which_base_gives_max")] <- paste("TIE_",gsub("/ucelocus.txt","",temp[j,1],fixed=TRUE),"_",gsub("/ucelocus.txt","",temp[j,2],fixed=TRUE),sep="")
              # What to do if genomes DON'T have the same length loci
            } else { #2AB
              # What to do if temp[j,5] is long than temp[j,6]
              if(as.numeric(temp[j,5])>as.numeric(temp[j,6])) { #1A
                temp_row[1,which(output_matrix[1,]=="max_length")] <- as.numeric(temp[j,5])
                temp_row[1,which(output_matrix[1,]=="which_base_gives_max")] <- gsub("/ucelocus.txt","",temp[j,1],fixed=TRUE)
                # What to do if temp[j,6] is long than temp[j,5]
              } else {  #1AB
                temp_row[1,which(output_matrix[1,]=="max_length")] <- as.numeric(temp[j,6])
                temp_row[1,which(output_matrix[1,]=="which_base_gives_max")] <- gsub("/ucelocus.txt","",temp[j,2],fixed=TRUE)
              } #1B
            } #2B
            output_matrix <- rbind(output_matrix,temp_row)
          } #IFB  
        # what to do if temp[j,3] is already in there, but temp[j,4] is not  
        } else { #5AB
          # if there is no value in the output_matrix for the [j,2] base genome
          if(all(is.na(output_matrix[(which(output_matrix[,which(output_matrix[1,]==temp[j,1])]==(temp[j,3]))),which(output_matrix[1,]==temp[j,2])]))) { #10A
            # for each of the rows in output_matrix that match this
            for (k in which(output_matrix[,which(output_matrix[1,]==temp[j,1])]==(temp[j,3]))) { #20A
              output_matrix[k,which(output_matrix[1,]==temp[j,2])] <- temp[j,4]
              # If the new base genome gives the same locus length
              if (as.numeric(output_matrix[k,which(output_matrix[1,]=="max_length")])==temp[j,6]) { #21A
                # If it already has tie in the title
                if (grepl("TIE_",output_matrix[k,which(output_matrix[1,]=="which_base_gives_max")])) { #22A
                  output_matrix[k,which(output_matrix[1,]=="which_base_gives_max")] <- paste(output_matrix[k,which(output_matrix[1,]=="which_base_gives_max")],"_",gsub("/ucelocus.txt","",temp[j,2],fixed=TRUE),sep="")
                } else { #22AB
                  output_matrix[k,which(output_matrix[1,]=="which_base_gives_max")] <- paste("TIE_",output_matrix[k,which(output_matrix[1,]=="which_base_gives_max")],"_",gsub("/ucelocus.txt","",temp[j,2],fixed=TRUE),sep="")
                } #22B
              } else { #21AB
                if (as.numeric(output_matrix[k,which(output_matrix[1,]=="max_length")])<as.numeric(temp[j,6])) { #23A
                  output_matrix[k,which(output_matrix[1,]=="which_base_gives_max")] <- gsub("/ucelocus.txt","",temp[j,2],fixed=TRUE)
                  output_matrix[k,which(output_matrix[1,]=="max_length")] <- as.numeric(temp[j,6])
                } #23B
              } #21B
           } #20B  
          # we have a problem because the locus from the [j,1] base genome is matching to multiple others at the blast threshold we used
          } else { #10AB
            output_matrix[which(output_matrix[,which(output_matrix[1,]==temp[j,1])]==temp[j,3]),which(output_matrix[1,]=="problem_locus")] <- "Y"
            temp_row <- t(matrix(NA,nrow=dim(output_matrix)[2]))
            temp_row[1,which(output_matrix[1,]=="problem_locus")] <- "Y"
            temp_row[1,which(output_matrix[1,]==temp[j,1])] <- temp[j,3]
            temp_row[1,which(output_matrix[1,]==temp[j,2])] <- temp[j,4]
            # What to do if both base genomes have equally long loci
            if(as.numeric(temp[j,5])==as.numeric(temp[j,6])) { #2A
              temp_row[1,which(output_matrix[1,]=="max_length")] <- as.numeric(temp[j,5])
              temp_row[1,which(output_matrix[1,]=="which_base_gives_max")] <- paste("TIE_",gsub("/ucelocus.txt","",temp[j,1],fixed=TRUE),"_",gsub("/ucelocus.txt","",temp[j,2],fixed=TRUE),sep="")
            # What to do if genomes DON'T have the same length loci
            } else { #2AB
              # What to do if temp[j,5] is long than temp[j,6]
              if(as.numeric(temp[j,5])>as.numeric(temp[j,6])) { #1A
                temp_row[1,which(output_matrix[1,]=="max_length")] <- as.numeric(temp[j,5])
                temp_row[1,which(output_matrix[1,]=="which_base_gives_max")] <- gsub("/ucelocus.txt","",temp[j,1],fixed=TRUE)
              # What to do if temp[j,6] is long than temp[j,5]
              } else {  #1AB
                temp_row[1,which(output_matrix[1,]=="max_length")] <- as.numeric(temp[j,6])
                temp_row[1,which(output_matrix[1,]=="which_base_gives_max")] <- gsub("/ucelocus.txt","",temp[j,2],fixed=TRUE)
              } #1B
           } #2B
           output_matrix <- rbind(output_matrix,temp_row)                
          }  #10B
        } #5B
     # if temp[j,3] NOT already in table   
     } else { #4AB
         # TO DO what to do if temp[j,4] is already in there, but temp[j,3] is not
        if (temp[j,4] %in% output_matrix[1:(dim(output_matrix)[1]),which(output_matrix[1,]==temp[j,2])]) { #3A
          # if there is no value in the output_matrix for the [j,1] base genome
          if(all(is.na(output_matrix[(which(output_matrix[,which(output_matrix[1,]==temp[j,2])]==(temp[j,4]))),which(output_matrix[1,]==temp[j,1])]))) { #10A
            # for each of the rows in output_matrix that match this
            for (k in which(output_matrix[,which(output_matrix[1,]==temp[j,2])]==(temp[j,4]))) { #20A
              output_matrix[k,which(output_matrix[1,]==temp[j,1])] <- temp[j,3]
              # If the new base genome gives the same locus length
              if (as.numeric(output_matrix[k,which(output_matrix[1,]=="max_length")])==temp[j,5]) { #21A
                # If it already has tie in the title
                if (grepl("TIE_",output_matrix[k,which(output_matrix[1,]=="which_base_gives_max")])) { #22A
                  output_matrix[k,which(output_matrix[1,]=="which_base_gives_max")] <- paste(output_matrix[k,which(output_matrix[1,]=="which_base_gives_max")],"_",gsub("/ucelocus.txt","",temp[j,1],fixed=TRUE),sep="")
                } else { #22AB
                  output_matrix[k,which(output_matrix[1,]=="which_base_gives_max")] <- paste("TIE_",output_matrix[k,which(output_matrix[1,]=="which_base_gives_max")],"_",gsub("/ucelocus.txt","",temp[j,1],fixed=TRUE),sep="")
                } #22B
              } else { #21AB
                if (as.numeric(output_matrix[k,which(output_matrix[1,]=="max_length")])<as.numeric(temp[j,5])) { #23A
                  output_matrix[k,which(output_matrix[1,]=="which_base_gives_max")] <- gsub("/ucelocus.txt","",temp[j,1],fixed=TRUE)
                  output_matrix[k,which(output_matrix[1,]=="max_length")] <- as.numeric(temp[j,5])
                } #23B
              } #21B 
            } #20B  
          # we have a problem because the locus from the [j,2] base genome is matching to multiple others at the blast threshold we used
          } else { #10AB
            output_matrix[which(output_matrix[,which(output_matrix[1,]==temp[j,2])]==temp[j,4]),which(output_matrix[1,]=="problem_locus")] <- "Y"
            temp_row <- t(matrix(NA,nrow=dim(output_matrix)[2]))
            temp_row[1,which(output_matrix[1,]=="problem_locus")] <- "Y"
            temp_row[1,which(output_matrix[1,]==temp[j,1])] <- temp[j,3]
            temp_row[1,which(output_matrix[1,]==temp[j,2])] <- temp[j,4]
            # What to do if both base genomes have equally long loci
            if(as.numeric(temp[j,5])==as.numeric(temp[j,6])) { #2A
              temp_row[1,which(output_matrix[1,]=="max_length")] <- as.numeric(temp[j,5])
              temp_row[1,which(output_matrix[1,]=="which_base_gives_max")] <- paste("TIE_",gsub("/ucelocus.txt","",temp[j,1],fixed=TRUE),"_",gsub("/ucelocus.txt","",temp[j,2],fixed=TRUE),sep="")
            # What to do if genomes DON'T have the same length loci
            } else { #2AB
              # What to do if temp[j,5] is long than temp[j,6]
              if(as.numeric(temp[j,5])>as.numeric(temp[j,6])) { #1A
                temp_row[1,which(output_matrix[1,]=="max_length")] <- as.numeric(temp[j,5])
                temp_row[1,which(output_matrix[1,]=="which_base_gives_max")] <- gsub("/ucelocus.txt","",temp[j,1],fixed=TRUE)
              # What to do if temp[j,6] is long than temp[j,5]
              } else {  #1AB
                temp_row[1,which(output_matrix[1,]=="max_length")] <- as.numeric(temp[j,6])
                temp_row[1,which(output_matrix[1,]=="which_base_gives_max")] <- gsub("/ucelocus.txt","",temp[j,2],fixed=TRUE)
              } #1B
           } #2B
           output_matrix <- rbind(output_matrix,temp_row)  
          } #10B
        # what to do if neither temp[j,3] or temp[j,4] exist     
        } else { #3AB    
            temp_row <- t(matrix(NA,nrow=dim(output_matrix)[2]))
            temp_row[1,which(output_matrix[1,]=="problem_locus")] <- "N"
            temp_row[1,which(output_matrix[1,]==temp[j,1])] <- temp[j,3]
            temp_row[1,which(output_matrix[1,]==temp[j,2])] <- temp[j,4]
            # What to do if both base genomes have equally long loci
            if(as.numeric(temp[j,5])==as.numeric(temp[j,6])) { #2A
              temp_row[1,which(output_matrix[1,]=="max_length")] <- as.numeric(temp[j,5])
              temp_row[1,which(output_matrix[1,]=="which_base_gives_max")] <- paste("TIE_",gsub("/ucelocus.txt","",temp[j,1],fixed=TRUE),"_",gsub("/ucelocus.txt","",temp[j,2],fixed=TRUE),sep="")
            # What to do if genomes DON'T have the same length loci
            } else { #2AB
              # What to do if temp[j,5] is long than temp[j,6]
              if(as.numeric(temp[j,5])>as.numeric(temp[j,6])) { #1A
                temp_row[1,which(output_matrix[1,]=="max_length")] <- as.numeric(temp[j,5])
                temp_row[1,which(output_matrix[1,]=="which_base_gives_max")] <- gsub("/ucelocus.txt","",temp[j,1],fixed=TRUE)
              # What to do if temp[j,6] is long than temp[j,5]
              } else {  #1AB
                temp_row[1,which(output_matrix[1,]=="max_length")] <- as.numeric(temp[j,6])
                temp_row[1,which(output_matrix[1,]=="which_base_gives_max")] <- gsub("/ucelocus.txt","",temp[j,2],fixed=TRUE)
              } #1B
           } #2B
           output_matrix <- rbind(output_matrix,temp_row)            
        } #3B
     }  #4B 
  }#6B   
write.table(output_matrix,paste(i,".summarized",sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)
}
       
    
