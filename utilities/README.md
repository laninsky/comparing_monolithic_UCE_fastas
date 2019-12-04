# Selecting a base genome
Based on the output of monolithic.sh, for all 'good' loci (not problematic within or between taxa) select_base_genome.R will print off the number of these loci that were found with each base genome (bigger is better), which base genomes recovered the greatest number of taxa on average (bigger is better), and which base genomes resulted in the longest average recovered sequences across taxa (again, bigger is better). You can use this information to make a decision about which base genome you think is the best for developing your probes from. To run this code, paste (or source) the contents of select_base_genome.R into your R console and then call it by:
```
select_base_genome(output_matrix_file)
```
where output_matrix file is the output from comparing_monolithic_UCE_fastas/monolithic.sh ("output_matrix.txt" or whatever you renamed it) 

# Whittling the probes fasta file down to only a subset of loci
You may want to use the results of blasting between the different base genomes, to filter 'problematic' loci from your final data set. The code `whittle_uce_probes.R` will whittle down your probe_fasta_file based on uce loci in your monolithic_output_name file (i.e. the output from comparing_monolithic_UCE_fastas/monolithic.sh) that are not 'problematic' (i.e. paralagous) within or between taxa (file_type = "monolithic") or by loci that are present in a list of uce_loci (file_type = "uce_list": see below). It will further filter them to just loci that are found across all the taxa in your dataset. To run this code, paste (or source) the contents of whittle_uce_probes.R into your R console and then call it by:
```
whittle_uce_probes(uce_list_file,probe_fasta_file,file_type,taxa_filter,basename)
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
* `taxa_filter` is an option that allows you to filter just to the UCE loci recovered across all taxa if using the "monolithic" option = "Y" (filter based on this) or "N" (do not filter)
* `basename` is the name of the taxa that you wish to select your final probeset across if you are inputting a "monolithic" file type

An example of running the code with monolithic.sh output e.g.
```
whittle_uce_probes("/Users/alanaalexander/Dropbox/beetles/grey_whittled_probes/output_matrix_99.txt","/Users/alanaalexander/Dropbox/beetles/grey_whittled_probes/Adephaga_11Kv1.fasta","monolithic","Y","Pterostichus.1")
```

An example of running the code with a uce locus file e.g.
```
whittle_uce_probes("C:\\Users\\Alana\\Dropbox\\beetles\\grey_whittled_probes\\uce_loci_from_baca_2017.txt","C:\\Users\\Alana\\Dropbox\\beetles\\grey_whittled_probes\\Coleoptera-UCE-1.1K-v1","uce_list")
```

This script will output a file, `whittled_UCE_probes.fasta` into your working directory, containing the whittled down probe set, and also print off a summary of how many probes/loci are in that file.

# Further whittling your loci
whittle_uce_probes might not filter your loci down to a small enough number of loci/probes. That's where `further_whittling_length` and `further_whittling_random` come in. Both of these functions will take your `monolithic_file` (the output from comparing_monolithic_UCE_fastas/monolithic.sh ("output_matrix.txt" or whatever you renamed it)), `probe_fasta_file` (the output probes file from the phyluce pipeline), `basename` (the taxa that probe_fasta_file is designed across), and a `no_of_loci` (the number of loci you want to shoot for). `further_whittling_length` will take the longest average `no_of_loci` from the `monolithic_file` and output the probes for these in the file "further_whittled_UCE_probes_length.fasta", and their associated properties in "filtered_output_matrix_length.txt". `further_whittling_random` does something similar, but instead of taking the `no_of_loci` longest average loci, it just picks out `no_of_loci` random loci from the `monolithic_file` and outputs the probes for these in the file "further_whittled_UCE_probes_random.fasta", and their associated properties in "filtered_output_matrix_random.txt". 

An example of running the `further_whittling_length` code:
`further_whittling_length("output_matrix.txt","whittled_UCE_probes.fasta","liocanthydrus",3213)`  

An example of running the `further_whittling_random` code:
`further_whittling_random("output_matrix.txt","whittled_UCE_probes.fasta","liocanthydrus",3213)`  

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
