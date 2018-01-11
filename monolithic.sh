ls *.fasta > fastalist.txt

no_fastas=`wc -l fastalist.txt | awk '{ print$1 }'`

# Separating out monolithic fastas by taxa the loci is characterized in (further grouped by basegenome used)
for i in `seq 1 $no_fastas`; do current_file=`head -n $i fastalist.txt | tail -n 1`; mv $current_file temp; echo $current_file > current_name; Rscript split_fasta_by_taxa.R; mv temp $current_file; done

# do some kind of list by folder (should give the individual taxa present)
basewd=`pwd`
for i in ls -d */; do cd $i; ls -d */ > base_genomes.txt; no_base_genomes=`wc -l base_genomes.txt | awk '{ print$1 }'`; for j in `seq 1 $(( no_base_genomes - 1 ))`; do blaster="`head -n $j base_genomes.txt | tail -n 1`"ucelocus.txt""; for k in `seq $(( j + 1 )) $no_base_genome`; do blastee="`head -n $k base_genomes.txt | tail -n 1`"ucelocus.txt""


# Then go into each folder and blast across genomes. Save this info per taxon



for i in taxa*fasta; do ... cross blast 

# Then take these files and combine them across taxa

