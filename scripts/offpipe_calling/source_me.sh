#!/bin/bash

module purge

module load cgpcavemanwrapper/1.18.2
module load cgppindel/3.11.0


module load bcftools-1.9/python-3.11.6
module load ensembl_vep/103.1

#Agilent SureSelect V5 paded100 bp 
HUM_PAD_MERGED_BAITSET=/lustre/scratch125/casm/team113da/projects/DERMATLAS/metadata/references/baitset/DNA/GRCh38_WES5_canonical_pad100.merged.bed

