Following the analysis of which UCE loci blasted to which UCE loci across base genomes within taxa, we wanted to extract "the good loci" for each sample in order to calculate coverage (using https://github.com/laninsky/reference_aligning_to_established_loci/tree/master/phase_everyone). To do this, we sorted the uce loci by "non-problematic status" (i.e one to one relationships of UCE loci found in different base genomes and different taxa) and by them being found in all 7 of the taxa we were evaluating. If there was more than one base genome that led to a locus being found across all 7 taxa, we used the base genome that resulted in the longest-alignment for the largest number of taxa. This gave us a list formatted like the following:
```
insilico-Omoglymmius-base/ucelocus.txt	uce-6195
insilico-Amphizoa-base/ucelocus.txt	uce-10032
insilico-Bembidion-base/ucelocus.txt	uce-12642
insilico-Bembidion-base/ucelocus.txt	uce-4506
```
Unimaginatively, we called this file file_locus.txt

Then in the same working directory we used for the previous step, we made a directory to hold all the fasta files:
```
mkdir fastas
```
And then carried out the following code to get a separate fasta file per locus, named by the taxa within it:
```
numberofloci=`wc -l file_locus.txt | awk '{print $1}'`
for i in amphizoa bemhap1 chlSer1 lioTuu1 omoHam1 pterMel1 traGib1;
do for j in seq 1 $numberofloci;
do currentline=`head -n $j file_locus.txt | tail -n 1`;
basefile=$i/`echo $currentline | awk '{print $1}'`;
locusname=`echo $currentline | awk '{print $2}'`;



grep -n pattern file.txt | cut -d : -f 1

```
