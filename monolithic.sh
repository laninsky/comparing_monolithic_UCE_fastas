ls *.fasta > fastalist.txt

no_fastas=`wc -l fastalist.txt | awk '{ print$1 }'`

# Separating out monolithic fastas by taxa the loci is characterized in (further grouped by basegenome used)
for i in `seq 1 $no_fastas`; do current_file=`head -n $i fastalist.txt | tail -n 1`; mv $current_file temp; echo $current_file > current_name; Rscript split_fasta_by_taxa.R; mv temp $current_file; done

# This loop is going into each taxon's folder and doing a blast between the fastas recovered using different basegenomes
basewd=`pwd`
for i in `cat fastalist.txt | sed 's/.insilico.incomplete.fasta//g'`;
do cd $i;
ls -d */ > base_genomes.txt;
no_base_genomes=`wc -l base_genomes.txt | awk '{ print$1 }'`;
# For every base_genome except the last one create a blast database
# because it will be compared to the base_genomes that came
# before, and we don't want duplicated comparisons
for j in `seq 1 $(( no_base_genomes - 1 ))`;
do blaster="`head -n $j base_genomes.txt | tail -n 1`"ucelocus.txt"";
$BLASTPATH/makeblastdb -in $blaster -dbtype nucl;
outfile_name="`echo $i | sed 's/\///g'`"_blast.txt"";
# Take this blast database and compare it to every "subsequent"
# base genome (i.e. do not do the reciprocal blasts that have already occured)
for k in `seq $(( j + 1 )) $no_base_genomes`;
do blastee="`head -n $k base_genomes.txt | tail -n 1`"ucelocus.txt"";
$BLASTPATH/blastn -db $blaster -query $blastee -perc_identity $BLASTSIM -outfmt '10 sseqid qseqid slen qlen' > temp;
sed "s|^|$blaster,$blastee,|g" temp >> $outfile_name;
rm temp;
done;
done;
cd $basewd;
done

# Summarizing across the blast files for each taxa to determine presence/absence of each locus in each taxon 
# predicated on which base genome was used.

Rscript blast_summary.R

# Taking the taxa specific summaries and combining across them
Rscript combining_across_taxa.R
