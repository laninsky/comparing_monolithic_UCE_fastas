# Whittling the probes fasta file down to only a subset of loci
You may want to use the results of blasting between the different base genomes, to filter 'problematic' loci from your final data set. The code `whittle_uce_probes.R` will whittle down your probe_fasta_file based on uce loci in your monolithic_output_name file (i.e. the output from comparing_monolithic_UCE_fastas/monolithic.sh) that are not 'problematic' (i.e. paralagous) within or between taxa (file_type = "monolithic") or by loci that are present in a list of uce_loci (file_type = "uce_list": see below). To run this code, paste (or source) the contents of whittle_uce_probes.R into your R console and then call it by:
```
whittle_uce_probes(uce_list_file,probe_fasta_file,file_type,basename)
```
where: 
* `uce_list_file` either points to the output from comparing_monolithic_UCE_fastas/monolithic.sh ("output_matrix.txt" or whatever you renamed it) or a list of uce loci you've assembled (depending on what you put for file_type: "monolithic" or "uce_list") 
* `probe_fasta_file` is the output probes file from the phyluce pipeline (corresponding to the basename taxon and/or the taxon you made your "uce_list" across e.g. the naming of the uce loci should be the same between the list and your probes file)
* `file_type` is either "monolithic" for the output_matrix.txt file from comparing_monolithic_UCE_fastas/monolithic.sh (which you can feel free to rename as long as it is still formatted the same) or "uce_list" for a flat file with a list of uce loci separated by each line. Make sure your loci in the list are named the same as in your probe file e.g.
```
uce-100
uce-1011
uce-1017
uce-1023
uce-1025
```
* `basename` is the name of the taxa that you wish to select your final probeset across if you are inputting a "monolithic" file type


An example of running the code with monolithic.sh output e.g.
```
whittle_uce_probes("/Users/alanaalexander/Dropbox/beetles/grey_whittled_probes/output_matrix_99.txt","/Users/alanaalexander/Dropbox/beetles/grey_whittled_probes/Adephaga_11Kv1.fasta","monolithic","Pterostichus.1")
```

An example of running the code with a uce locus file e.g.
```
whittle_uce_probes("C:\\Users\\Alana\\Dropbox\\beetles\\grey_whittled_probes\\uce_loci_from_baca_2017.txt","C:\\Users\\Alana\\Dropbox\\beetles\\grey_whittled_probes\\Coleoptera-UCE-1.1K-v1","uce_list")
```

# Extracting the 'good' loci
You may then want to extract the sequences of each taxa for your whittled down list. For example, we wanted to extract "the good loci" for each sample in order to calculate coverage (using https://github.com/laninsky/reference_aligning_to_established_loci/tree/master/phase_everyone). To do this, in excel (yes, I should write some R code for this!) we sorted the uce loci by "non-problematic status" (i.e one to one relationships of UCE loci found in different base genomes and different taxa) and by them being found in all 7 of the taxa we were evaluating. If there was more than one base genome that led to a locus being found across all 7 taxa, we used the base genome that resulted in the longest-alignment for the largest number of taxa. This gave us a list formatted like the following:
```
insilico-Omoglymmius-base/ucelocus.txt	uce-6195
insilico-Amphizoa-base/ucelocus.txt	uce-10032
insilico-Bembidion-base/ucelocus.txt	uce-12642
insilico-Bembidion-base/ucelocus.txt	uce-4506
```
Unimaginatively, we called this file file_locus.txt

Then in the same working directory that contained the monolithic fasta files, we made a directory to hold all the resultant fasta files:
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
