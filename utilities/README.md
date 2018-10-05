# Extracting the 'good' loci
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
for i in amphizoa bemHap1 chlSer1 lioTuu1 omoHam1 pterMel1 traGib1;
do for j in `seq 1 $numberofloci`;
do currentline=`head -n $j file_locus.txt | tail -n 1`;
basefile=$i/`echo $currentline | awk '{print $1}'`;
locusname=`echo $currentline | awk '{print $2}'`;
lineno=`grep -n ${locusname}$ $basefile | cut -d : -f 1`;
lineno=$((lineno+1));
echo ">""$i" >> fastas/uce-$j.fasta;
head -n $lineno $basefile | tail -n 1 >> fastas/uce-$j.fasta;
done;
done;
```

# Whittling the probes fasta file down to only a subset of loci
To use, paste the contents of whittle_uce_probes.R into your R console.

This code will whittle down your probe_fasta_file by loci that are present in your monolithic_output_name file and that are not 'problematic' (i.e. paralagous) within or between taxa. To run this code:"
```
whittle_uce_probes(uce_list_file,probe_fasta_file,basename,file_type)
```
where uce_list_file is either the output from comparing_monolithic_UCE_fastas/monolithic.sh or a list of uce loci (depending on what you put for file_type), probe_fasta_file is the output probes file from the phyluce pipeline, basename is the name of the taxa that you designed your final probeset across, and file_type is either "monolithic" for the output of comparing_monolithic_UCE_fastas/monolithic.sh or "file_list" for a flat file with a list of uce loci separated by each line. Make sure your loci in the list are named the same as in your probe file e.g.
```
uce-100
uce-1011
uce-1017
uce-1023
uce-1025

```
An example of running the code with monolithic.sh output e.g.
```
whittle_uce_probes("/Users/alanaalexander/Dropbox/beetles/grey_whittled_probes/output_matrix_99.txt","/Users/alanaalexander/Dropbox/beetles/grey_whittled_probes/Adephaga_11Kv1.fasta","Pterostichus.1","monolithic")
```

An example of running the code with a uce locus file e.g.
```
whittle_uce_probes("C:\\Users\\Alana\\Dropbox\\beetles\\grey_whittled_probes\\uce_loci_from_baca_2017.txt","C:\\Users\\Alana\\Dropbox\\beetles\\grey_whittled_probes\\Coleoptera-UCE-1.1K-v1","Pterostichus.1","file_list")
```
