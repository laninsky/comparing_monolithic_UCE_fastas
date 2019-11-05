# comparing_monolithic_UCE_fastas v0.2: combining_across_taxa.R

# Loading required libraries
if (!require('tidyverse')) install.packages('tidyverse'); library('tidyverse')

# From the summarized blast results, getting a file list
file_list <- list.files(pattern=".summarized",recursive=TRUE)

# From the file_list, getting a list of the taxa names
taxa_list <- unlist(strsplit(file_list,"/",fixed=TRUE))
taxa_list <- taxa_list[seq((length(taxa_list)/length(file_list)),length(taxa_list),(length(taxa_list)/length(file_list)))]
taxa_list <- unlist(strsplit(taxa_list,"_blast.txt.summarized")) 

# Using these taxa names to create some more column variable names for our output
length_list <- paste(taxa_list,"_length",sep="")
longbase_list <- paste(taxa_list,"_longest_base_genome",sep="")

# Getting the base names for column list from the first of our summarized file list.
temp_output <- read_table2(file_list[1])
base_list <- names(temp_output[1,1:((which(names(temp_output)=="max_length"))-1)])

# Finding the maximum column number where uce loci are listed
pivot_col <- which(names(temp_output)=="max_length")-1

# Creating the names for our output tables
output_matrix <- NULL
output_matrix_names <- c(base_list,length_list,longbase_list,paste(base_list,"_2",sep=""),"between_taxa_problem")
problem_taxa <- tibble(taxa=character(),base_genome=character(),uce_locus=character())

# How far from the original UCE data the taxa presence/absence column is
how_many_cols <- length(output_matrix_names)-length(base_list)-1

# Stepping through each of the files
for (i in file_list) { #1A
  print("Up to:")
  print(i)
  print("Out of:")
  print(file_list)
  
  # Getting the current taxa name from the currrent file 
  taxa <- unlist(strsplit(i,"/",fixed=TRUE))
  taxa <- taxa[length(taxa)]
  taxa <- unlist(strsplit(taxa,"_blast.txt.summarized"))
  
  # Which columns the current taxa will inform in output_matrix.txt
  output_taxa <- which(grepl(paste(taxa,"_l",sep=""),output_matrix_names)==TRUE)
  
  # Reading in the summarized blast matches for the taxa we are up to
  temp_output <- read_table2(i) 
  
  # Extracting loci which are problematic within taxa and adding these
  # to our problem_taxa table to be dealt with later on
  temp_problematic_loci <- temp_output %>% filter(problem_locus=="Y") %>% select(c(-max_length,-which_base_gives_max,-problem_locus))
  temp_problematic_loci <- pivot_longer(temp_problematic_loci,names(temp_problematic_loci),names_to="base_genome",values_to="uce_locus") %>% filter(!(is.na(uce_locus)))
  temp_problematic_loci <- add_column(temp_problematic_loci,taxa=(rep(taxa,dim(temp_problematic_loci)[1])),.before = TRUE)
  problem_taxa <- bind_rows(problem_taxa,temp_problematic_loci) %>% distinct()
  
  # Removing problematic matches from the taxa we are up to
  temp_output <- read_table2(i) %>% filter(problem_locus=="N")
  
  # If output_matrix doesn't yet exist, creating it with the first sample
  if(is.null(output_matrix)) {
    output_matrix <- matrix(NA,nrow=dim(temp_output)[1],ncol=length(output_matrix_names))
    output_matrix[,1:pivot_col] <- as.matrix(temp_output[,1:pivot_col])
    output_matrix[,(output_taxa[1])] <- as.matrix(temp_output[,(pivot_col+1)])
    output_matrix[,(output_taxa[2])] <- as.matrix(temp_output[,(pivot_col+2)])
    output_matrix <- as_tibble(output_matrix)
    names(output_matrix) <- output_matrix_names
    for (j in 1:pivot_col) {
      output_matrix[(which(!is.na(output_matrix[,j]))),(j+how_many_cols)] <- taxa
    }
    # Or, if it does exist, moving into reconciling this taxon with the existing output_matrix file
  } else {
    # First creating a correctly formated taxa specific output matrix
    temp_output_matrix <- matrix(NA,nrow=dim(temp_output)[1],ncol=length(output_matrix_names))
    temp_output_matrix[,1:pivot_col] <- as.matrix(temp_output[,1:pivot_col])
    temp_output_matrix[,(output_taxa[1])] <- as.matrix(temp_output[,(pivot_col+1)])
    temp_output_matrix[,(output_taxa[2])] <- as.matrix(temp_output[,(pivot_col+2)])
    temp_output_matrix <- as_tibble(temp_output_matrix)
    names(temp_output_matrix) <- output_matrix_names
    for (j in 1:pivot_col) {
      temp_output_matrix[(which(!is.na(temp_output_matrix[,j]))),(j+how_many_cols)] <- taxa
    }
    
    # Then, obtaining rows that match between temp_output_matrix and output_matrix
    match_rows <- semi_join(temp_output_matrix,output_matrix,by=names(output_matrix)[1:pivot_col])
    
    # And rows that don't
    problem_rows <- anti_join(temp_output_matrix,output_matrix,by=names(output_matrix)[1:pivot_col])  
    
    # Then stepping through match_rows and adding data to output_matrix
    for (j in 1:dim(match_rows)[1]) {
      # Getting the row that matches in output_matrix
      for (k in 1:pivot_col) {
        row_index <- which(output_matrix[,k]==toString(match_rows[j,k]))
        if(length(row_index)>0) {
          break
        }
      }
      # Replacing any values that are NA in output_matrix with the equivalent
      # value from the match row we are up to (with the exception of which taxa
      # each base genome characterizes a locus in - we'll do this at the end)
      NAs_to_replace <- which(is.na(output_matrix[row_index,]))
      NAs_to_replace <- NAs_to_replace[which(NAs_to_replace < (dim(match_rows)[2]-pivot_col))]
      output_matrix[row_index,NAs_to_replace] <- match_rows[j,NAs_to_replace]
      
      # Adding the taxa to the base_genomes to show what taxa were found for each base_genome
      cols_to_paste <- which(!is.na(match_rows[j,(dim(match_rows)[2]-pivot_col):dim(match_rows)[2]]))+dim(match_rows)[2]-pivot_col-1
      output_matrix[row_index,cols_to_paste] <- gsub("NA,","",paste(output_matrix[row_index,cols_to_paste],match_rows[j,cols_to_paste],sep=","))
      
      # Progress meter
      if ((j%%1000)==0) {
        print(paste(round(((j/dim(match_rows)[1])*100),2),"% of the way through non-problematic loci for ",taxa,sep=""))
      }  
    }
    
    # Going through the problematic rows
    for (j in 1:dim(problem_rows)[1]) {
      
      # Progress meter
      if ((j%%100)==0) {
        print(paste(round(((j/dim(problem_rows)[1])*100),2),"% of the way through problematic loci for ",taxa,sep=""))
      }  
      
      row_index <- NULL
      # Getting the rows that match in output_matrix
      for (k in 1:pivot_col) {
        row_index <- c(row_index,grep(paste(paste(toString(problem_rows[j,k]),"$|",sep=""),paste(toString(problem_rows[j,k]),",",sep=""),sep=""),as.matrix(output_matrix[,k])))
      }
      
      row_index <- unique(row_index)
      
      # If no row_index found, the locus isn't problematic between taxa, it just has not
      # been present in the taxa previously assessed. It can therefore be added to the
      # bottom of output_matrix.txt
      
      if(length(row_index)==0) {
        output_matrix <- bind_rows(output_matrix,problem_rows[j,])
        # If the row_index is not 0 then we'll need to investigate why
      } else {
        # If multiple rows are present then loci is almost certainly problematic!
        # We'll deal with that below, but in the meantime...
        if(length(row_index)==1) {
          
          # Sometimes the rows don't match up because of NA values. In this case,
          # the rows can be combined as long as the rest of the information is consistent
          # First need to get the non-NAs found for up to pivot_col for both samples
          nonNAs <- c(which(!is.na(output_matrix[row_index,(1:pivot_col)])),which(!is.na(problem_rows[j,(1:pivot_col)])))[duplicated(c(which(!is.na(output_matrix[row_index,(1:pivot_col)])),which(!is.na(problem_rows[j,(1:pivot_col)]))))]
          
          # Then seeing if they match up over this
          if(all(output_matrix[row_index,nonNAs]==problem_rows[j,nonNAs])) {
            
            # Replacing any values that are NA in output_matrix with the equivalent
            # value from the match row we are up to (with the exception of which taxa
            # each base genome characterizes a locus in - we'll do this at the end)
            NAs_to_replace <- which(is.na(output_matrix[row_index,]))
            NAs_to_replace <- NAs_to_replace[which(NAs_to_replace < (dim(problem_rows)[2]-pivot_col))]
            output_matrix[row_index,NAs_to_replace] <- problem_rows[j,NAs_to_replace]
            
            # Adding the taxa to the base_genomes to show what taxa were found for each base_genome
            cols_to_paste <- which(!is.na(problem_rows[j,(dim(problem_rows)[2]-pivot_col):dim(problem_rows)[2]]))+dim(problem_rows)[2]-pivot_col-1
            output_matrix[row_index,cols_to_paste] <- gsub("NA,","",paste(output_matrix[row_index,cols_to_paste],problem_rows[j,cols_to_paste],sep=","))
            # Otherwise...there are some discrepancies that need to be addressed
            # First up is the situation where loci match up to different loci
            # depending on the base genome used
          } else {
            # For each column position
            for (k in 1:(dim(problem_rows)[2]-pivot_col-1)) {
              # If output_matrix is.na, then copying the value over from problem_rows
              if(is.na(output_matrix[row_index,k])) {
                output_matrix[row_index,k] <- problem_rows[j,k]
              } else {                
                # If problem row doesn't have an NA and values do not match up paste them together
                if (!is.na(problem_rows[j,k]) & output_matrix[row_index,k]!=problem_rows[j,k]) {
                  output_matrix[row_index,k] <- paste(output_matrix[row_index,k],problem_rows[j,k],sep=",")
                }
              }
            }
            # Adding the taxa to the base_genomes to show what taxa were found for each base_genome
            cols_to_paste <- which(!is.na(problem_rows[j,(dim(problem_rows)[2]-pivot_col):dim(problem_rows)[2]]))+dim(problem_rows)[2]-pivot_col-1
            output_matrix[row_index,cols_to_paste] <- gsub("NA,","",paste(output_matrix[row_index,cols_to_paste],problem_rows[j,cols_to_paste],sep=","))
            
            # Adding in that this locus is a between_taxa_problem
            output_matrix[row_index,dim(output_matrix)[2]] <- "Y"
          }
        } else {
          # If there are more than one records then there is almost certainly something squiffy 
          # going on! We'll use the problem_rows as the "row of record" so we can drop the 
          # other rows from the output_matrix. Generally this situation results when loci are
          # considered distinct for one genome, but synonymous for others.
          
          # For each matching row
          for (l in row_index) {
            # For each column position, up to where we list whether a base genome recovered 
            # a taxon or not....
            for (k in 1:(dim(problem_rows)[2]-pivot_col-1)) {
              # If output_matrix is.na, then copying the value over from problem_rows
              if(is.na(problem_rows[j,k])) {
                problem_rows[j,k] <- output_matrix[l,k]
              } else {                
                # If problem row doesn't have an NA and values do not match up paste them together
                if (!is.na(output_matrix[l,k]) & output_matrix[l,k]!=problem_rows[j,k]) {
                  problem_rows[j,k] <- paste(output_matrix[l,k],problem_rows[j,k],sep=",")
                }
              }
            }
            # For the taxa recording rows
            for (k in (dim(problem_rows)[2]-pivot_col):(dim(problem_rows)[2]-1)) {
              # If the output matrix row has info
              if (!is.na(output_matrix[l,k])) {
                # But the problem_rows entry doesn't...
                if (is.na(problem_rows[j,k])) {
                  # Then copying over the output_matrix entry
                  problem_rows[j,k] <- output_matrix[l,k]
                  # If both the output_matrix row and problem_rows already have entries then...
                } else {    
                  # If the output_matrix entry isn't already present, 
                  # paste it in separated by a comma
                  if (!grepl(as.matrix(output_matrix[l,k]),as.matrix(problem_rows[j,k]))) {
                    problem_rows[j,k] <- gsub("NA,","",paste(output_matrix[l,k],problem_rows[j,k],sep=","))    
                  }
                }
              }
            }
          }
          # Mentioning this is a problem row
          problem_rows[j,dim(problem_rows)[2]] <- "Y"
          # Dropping the rows recorded in problem_rows from output_matrix
          output_matrix <- output_matrix[-row_index,]
          # Adding the problem row to the bottom of the output_matrix
          output_matrix <- bind_rows(output_matrix,problem_rows[j,])
        }
      }
    }
  }
}  

# Now dealing with the loci that were problematic within taxa     
print("Now processing the loci that had problematic matches across base genomes within taxa")
for (j in 1:(dim(problem_taxa)[1])) {
  
  if ((j%%1000)==0) {
    print(paste(round(((j/dim(problem_taxa)[1])*100),2),"% of the way through within-problematic loci. This is the last step of the pipeline",sep=""))
  }  
  
  # which column records the uce locus for that base genome
  problem_col <- which(names(output_matrix) %in% problem_taxa[j,2])[1]
  # Which row has that uce locus 
  problem_row <- grep(paste(problem_taxa[j,3],"(,|$)",sep=""),as.matrix(output_matrix[,problem_col]))
  
  if(length(problem_row)>0) {
    # Which columns correspond to the taxa affected
    taxa_col <- grep(paste(problem_taxa[j,1],"_l",sep=""),names(output_matrix),fixed=TRUE)
    # Recording "problem_within" for all of those columns
    output_matrix[problem_row,taxa_col] <- "problem_within"
    
    # Getting location of the taxa recording column
    problem_col <- grep(toString(problem_taxa[j,2]),names(output_matrix),fixed=TRUE)[2]
    # And noting that this taxa was present, but with problems
    output_matrix[problem_row,problem_col] <- paste(output_matrix[problem_row,problem_col],",",problem_taxa[j,1],"_prob_w",sep="")
    output_matrix[problem_row,problem_col] <- gsub("NA","",output_matrix[problem_row,problem_col])
  }
}

write.table(output_matrix,"output_matrix.txt",quote=FALSE,row.names=FALSE,col.names=TRUE)
