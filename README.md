# comparing_monolithic_UCE_fastas
### Blasting between monolithic UCE fasta files to find out which loci are shared between different base genomes

This repository came about because we had 7 different monolithic fasta files containing all UCE loci for all taxa, identified using different base genomes. The idea with this script was to create a "rosetta" stone via BLAST to match up loci identified using different base genomes, so that we could look at the effect of different base genomes on the ability for different UCE loci to be identified in different taxa.

The code assumes you have a folder with just your monolithic fasta files in it (and no other fasta files). It also assumes you have copied the necessary R-scripts and bash script (in this repository) into the same folder as well. It assumes the fasta header for each locus is similar to the following (underscore before taxa name, white space after):
```
>use-8084_lioTuu1 |uce-8084
```
If this isn't the case, you will need to modify lines 39-40 and 53-54 in split_fasta_by_taxa.R

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
TBD
```
