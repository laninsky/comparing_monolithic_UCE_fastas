# comparing_monolithic_UCE_fastas v0.1
### Blasting between monolithic UCE fasta files to find out which loci are shared between different base genomes

This repository came about because we had 7 different monolithic fasta files containing all UCE loci for all taxa, identified using different base genomes. The idea with this script was to create a "rosetta" stone via BLAST to match up loci identified using different base genomes, so that we could look at the effect of different base genomes on the ability for different UCE loci to be identified in different taxa (confused? The [step-by-step](https://github.com/laninsky/comparing_monolithic_UCE_fastas#what-is-it-doing) guide to what the script does might be helpful).

### What you need and how to run it
The code assumes you have a folder with your monolithic fasta files in it (and no other fasta files). It also assumes you have copied the necessary R-scripts and bash script (in this repository) into the same folder as well. It assumes the fasta header for each locus within these folders is similar to the following (underscore before taxa name, white space after, no underscore in the uce-locus name):
```
>uce-8084_lioTuu1 |uce-8084
```
If this isn't the case, you will need to modify lines 39-40 and 53-54 in split_fasta_by_taxa.R to match the regex needed to extract the taxa names. It also assumes that these taxa names match those in the "taxa name" part of the `.fasta` files. Please rename the .fasta files if this isn't the case.

You also need to point to the BLAST binaries folder by setting the following variable before starting
```
export BLASTPATH=path/to/bin
e.g. export BLASTPATH=/Users/alanaalexander/bin/ncbi-blast-2.7.1+/bin
```

Finally, you'll need to set the BLAST matching similarity percentage threshold e.g.
```
export BLASTSIM=95
```

After doing all of this, start it off by:
```
bash monolithic.sh
```

### What is it doing?
Let's say we have four monolithic fasta files corresponding to using the following four taxa as base genomes: `liocanthydrus`, `neohydrocoptus`, `sternocanthus`, and `suphisellus`
```
liocanthydrus_insilico-incomplete.fasta
neohydrocoptus-insilico-incomplete.fasta
sternocanthus-insilico-incomplete.fasta
suphisellus-insilico-incomplete.fasta
```
  
Let's peak inside `liocanthydrus_insilico-incomplete.fasta`:
```
>uce-50205_liocanthydrus |uce-50205
TTGATATTTTACAAAACCGTTAGAAATCTGAATCTCCGGCAACCGCTTACTGTCAAGTGT
...
>uce-18224_neohydrocoptus |uce-18224
AGTTTTATTCAAGTGAATTTTTTTCTCTCCTTAAACACAATGCCGAGGATAACAATGAAT
...
>uce-118947_sternocanthus |uce-118947
GTCCGCTCGTTCTCTGCGAAGAACGGGGCCCAGGCCATCCACGGTGAACAATGCGAGGAG
...
>uce-22825_suphisellus |uce-22825
ATTTTAACTCTCCATTATAAGAACTGAATTTCTACTGAACAGTGCCAACACAATTTACCT
```
Contained in this file is every single UCE locus identified in each of the four taxa, using liocanthydrus as the base genome.
  
The first major step ([Line 6](https://github.com/laninsky/comparing_monolithic_UCE_fastas/blob/master/monolithic.sh)) is to go into each of the monolithic fasta files, and pull out all the loci for each taxon and then dump them into separate directories. Inside these directories, the loci for that taxon are separated out into subdirectories corresponding to the base genome they are found in e.g. (for liocanthydrus again):
```
|-liocanthydrus
	 |- liocanthydrus_insilico-incomplete.fasta
		 |- ucelocus.txt
	 |- neohydrocoptus-insilico-incomplete.fasta
		 |- ucelocus.txt
	 |- sternocanthus-insilico-incomplete.fasta
		 |- ucelocus.txt
	 |- suphisellus-insilico-incomplete.fasta
		 |- ucelocus.txt
```

The next step ([Lines 8-35](https://github.com/laninsky/comparing_monolithic_UCE_fastas/blob/master/monolithic.sh)) is to blast between the UCE loci found using different base genomes, within a taxon. We are doing this because if, for example for liocanthydrus the locus `uce-10001` derived using neohydrocoptus as a base genome matches to two UCE loci derived using another base genome, say `uce-987` and `uce-118947` using suphisellus as the base genome, then that locus looks to be paralagous within the liocanthydrus genome (or at least similar enough it is likely to cause issues in downstream analyses).

The BLAST results from this are written out into the directory of each taxon as taxaname_blast.txt:
```
|-liocanthydrus
	 |- liocanthydrus_blast.txt
	 |- liocanthydrus_insilico-incomplete.fasta
		 |- ucelocus.txt.nsq
		 |- ucelocus.txt.nin
		 |- ucelocus.txt.nhr
		 |- ucelocus.txt
	 |- neohydrocoptus-insilico-incomplete.fasta
		 |- ucelocus.txt.nsq
		 |- ucelocus.txt.nin
		 |- ucelocus.txt.nhr
		 |- ucelocus.txt
	 |- sternocanthus-insilico-incomplete.fasta
		 |- ucelocus.txt.nsq
		 |- ucelocus.txt.nin
		 |- ucelocus.txt.nhr
		 |- ucelocus.txt
	 |- suphisellus-insilico-incomplete.fasta
		 |- ucelocus.txt     
|-neohydrocoptus
   	 |- neohydrocoptus_blast.txt
	 |- liocanthydrus_insilico-incomplete.fasta (collapsed)
	 |- neohydrocoptus-insilico-incomplete.fasta (collapsed)
	 |- sternocanthus-insilico-incomplete.fasta (collapsed)
	 |- suphisellus-insilico-incomplete.fasta  (collapsed) 
|-sternocanthus
   	 |- sternocanthus_blast.txt
	 |- liocanthydrus_insilico-incomplete.fasta (collapsed)
	 |- neohydrocoptus-insilico-incomplete.fasta (collapsed)
	 |- sternocanthus-insilico-incomplete.fasta (collapsed)
	 |- suphisellus-insilico-incomplete.fasta (collapsed)  
|-suphisellus
   	 |- suphisellus_blast.txt
	 |- liocanthydrus_insilico-incomplete.fasta (collapsed)
	 |- neohydrocoptus-insilico-incomplete.fasta (collapsed)
	 |- sternocanthus-insilico-incomplete.fasta (collapsed)
	 |- suphisellus-insilico-incomplete.fasta (collapsed) 
```
You might note that in the above file tree, tht within the liocanthydrus folder the subfolder suphisellus-insilico-incomplete.fasta does not have any \*.nsq, \*.nin, or \*.nhr files. This is because we only do one-way blast comparisons between the base genomes, and the suphisellus-derived UCE loci for liocanthydrus have already been compared to all the other genomes (who do have blast databases created for this comparison).

For each taxon, the blast results are then summarized to determine presence/absence of each locus in each taxon ([Line 40](https://github.com/laninsky/comparing_monolithic_UCE_fastas/blob/master/monolithic.sh)) e.g. for the liocanthydrus directory:
```
|-liocanthydrus
	 |- liocanthydrus_blast.txt
	 |- liocanthydrus_blast.txt.summarized
	 |- liocanthydrus_insilico-incomplete.fasta
	 |- neohydrocoptus-insilico-incomplete.fasta
	 |- sternocanthus-insilico-incomplete.fasta
	 |- suphisellus-insilico-incomplete.fasta
```

And then finally these results are summarized across all of the taxa ([Line 43](https://github.com/laninsky/comparing_monolithic_UCE_fastas/blob/master/monolithic.sh), including identification of "between_taxa_problem" loci that look paralagous when compared across taxa i.e. they appear to be single-copy when looking at each taxon separately, but between taxa comparisons do not show a one-to-one BLAST match) and written out as "output_matrix.txt":
```
|- output_matrix.txt
|- liocanthydrus (collapsed)
|- neohydrocoptus (collapsed)
|- sternocanthus (collapsed)
|- suphisellus (collapsed)
```
output_matrix.txt contains the following columns/groups of columns in the following order:
* Base genomes (first group of columns named with original fasta name suffixed by "/ucelocus.txt"). Contains the name of the UCE locus as found within that base genome. Generally will only be a "one to one" relationship, unless locus has a between_taxa_problem (loci different within a taxon then match to the same locus across taxa). If multiple loci are present for a base genome, they are separated by commas.	
* \_length (group of columns named taxaname_length). Contains maximum length of locus for that taxon across all base genomes. Problem_within indicates that within a taxon there was not a one-to-one match of loci across different base genomes. NA = this locus is not found in that taxon.
* \_longest_base_genome (group of columns named taxaname_longest_base_genome). The base genome that the locus of maximum length was recovered in for each taxon. Same deal for "problem_within" and NA as for the length columns.
* Base genomes (second group of columns named with original fasta name suffixed by "/ucelocus.txt"). Which taxa were found for each base genome for this locus. If taxa had a "problem_within", labeled taxa_prob_w.
* between_taxa_problem (single column). A locus with a "between_taxa_problem" == "Y" (see explanation above), otherwise NA.

### Downstream processing
The utilities folder in this repository contains scripts for extracting the 'good' UCE loci (loci that do not appear to be paralogous in any lineage, and that are found in every lineage) based on the output of monolithic.sh ("output_matrix.txt"), and also for whittling down the probe file to target just these loci.

### Programs/packages necessary for the pipeline:
```
BLAST:
Altschul, Stephen F., Gish, Warren, Miller, Webb, Myers, Eugene W., and Lipman, David J. (1990). Basic local alignment search tool. J. Mol. Biol. 215; 403-410. Gapped BLAST is described in Altschul, Stephen F., Madden, Thomas L., Schaffer, Alejandro A., Zhang, Jinghui, Zhang, Zheng, Miller, Webb, and Lipman, David J. (1997). Gapped BLAST and PSI-BLAST: a new generation of protein database search programs. Nucleic Acids Res. 25(17); 3389-3402.

R:
R: A Language and Environment for Statistical Computing. R Core Team. R Foundation for Statistical Computing. Vienna, Austria
```

Along with the programs above, to cite this pipeline:
```
Alexander, A. 2018. comparing_monolithic_UCE_fastas v0.0. Available from: https://github.com/laninsky/comparing_monolithic_UCE_fastas

Gustafson, G.T., Alexander, A., Sproul, J.S., Pflug, J.M., Maddison, D.R. and Short, A.E., 2019. Ultraconserved element (UCE) probe set design: Base genome and initial design parameters critical for optimization. Ecology and Evolution.
```

### Version history
0.1: Version used in Baca et al. TBD: made some modifications so that folders already existing in the directory were ignored, and to be more resistent to different naming schemes.  
0.0: Version used in Gustafson et al. (2019)
