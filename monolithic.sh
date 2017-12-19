ls *.fasta > fastalist.txt

no_fastas=`wc -l fastalist.txt | awk '{ print$1 }'`

# will replace mv temp $current_file with rm temp after R code inserted before then
for i in `seq 1 $no_fastas`; do current_file=`head -n $i fastalist.txt | tail -n 1`; mv $current_file temp; 

# Rscript here

echo $current_file > current_name; mv temp $current_file; done




for i in taxa*fasta; do ... cross blast 

