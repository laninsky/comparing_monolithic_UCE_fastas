
intable <- read.table("temp",header=FALSE,stringsAsFactors=FALSE,sep="\t")

#"one-lining" the fasta file
rows <- dim(intable)[1]

to_write <- intable[1,1]

sequencepaste <- NULL

for (j in 2:rows) {
if ((length(grep(">",intable[j,1])))>0) {
to_write <- rbind(to_write,toupper(sequencepaste))
to_write <- rbind(to_write,intable[j,1])
sequencepaste <- NULL
} else {
sequencepaste <- paste(sequencepaste,intable[j,1],sep="")
}
}

to_write <- rbind(to_write,toupper(sequencepaste))

forsort <- matrix(NA,ncol=2,nrow=((dim(to_write)[1])/2))
forsort[,1] <- to_write[(seq(1,(dim(to_write)[1]),2)),1]
forsort[,2] <- to_write[(seq(2,(dim(to_write)[1]),2)),1]
forsort <- forsort[order(forsort[,1]),]
to_write[(seq(1,(dim(to_write)[1]),2)),1] <- forsort[,1]
to_write[(seq(2,(dim(to_write)[1]),2)),1] <- forsort[,2]

rm(forsort)
rm(intable)
rm(rows)
rm(sequencepaste)

# getting the unique taxa names in our file
nolines <- dim(to_write)[1]

title_lines <- to_write[(seq(1,nolines,2)),]
taxa <- unlist(strsplit(title_lines,"_"))[(seq(2,nolines,2))]
taxa <- unique(unlist(strsplit(taxa,"\\s+"))[(seq(1,nolines,2))])

# finding what base genome we are up to
basegenome <- gsub(".fasta","",as.matrix(read.table("current_name"))[1,1],fixed=TRUE)

# creating subfolders for each of those
for (i in taxa) {
dir.create(file.path(i), showWarnings = FALSE)
dir.create(file.path(i, basegenome), showWarnings = FALSE)
}

# to do, tick through the files, writing out loci with uce as the title to the appropriate folder

file.path(i,basegenome,"ucelocus.txt")




