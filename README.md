# comparing_monolithic_UCE_fastas
Blasting between monolithic UCE fasta files to find out which loci are shared between different base genomes

This repository came about because we had 8 different monolithic fasta files containing all UCE loci for all taxa, identified using different base genomes. The idea with this script was to create a "rosetta" stone via BLAST to match up loci identified using different base genomes, so that we could look at the effect of different base genomes on the ability for different UCE loci to be identified in different taxa.

The code assumes you have a folder with just your monolithic fasta files in it (and no other fasta files). It also assumes you have copied the necessary R-scripts (in this repository) into the same folder as well. It assumes the fasta header for each locus is similar to the following (underscore before taxa name, white space after):
```
>use-8084_lioTuu1 |uce-8084
```
If this isn't the case, you will need to modify lines 39-40 and 53-54 in split_fasta_by_taxa.R

You also need to point to the BLAST binaries folder by setting the following variable before starting
```
exoprt BLASTPATH=path/to/bin
e.g. export BLASTPATH=/Users/alanaalexander/bin/ncbi-blast-2.7.1+/bin
```
