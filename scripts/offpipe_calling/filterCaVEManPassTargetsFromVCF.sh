#!/bin/bash
#This script needs in order the parameters:
#   1) sample.list : file with the sample.list file with the information of sample pairs in the format
#   2) EXTended_baits_bedfile: BED file with the coordinates of the baits used e.g. "/lustre/CASM_BAITS/SureSelect_Whole_Human_Exome_v5_GRCh38_liftover_160.sorted.w100.canonical.bed"
mkdir -p logs

samplelist=$1
target=$2


for f in `cat $samplelist`; do 
    echo $f;
    #This changes the names of the folders from Tumour-vs-Normal to Tumour_vs_Normal
    name=`echo $f | sed 's/-/_/g'`;
    bsub -e logs/filt.target.$f.e -o logs/filt.target.$f.o -M6000 -n 2 -R"select[mem>6000] rusage[mem=6000] span[hosts=1]" "module load bcftools/1.9 bgzip/1.9; bcftools view -f PASS $f/$name.muts.ids.mnv.flagged.v1.vcf.gz -Oz -o $f/$name.muts.ids.mnv.flagged.v1.target.pass.vcf.gz -T $target; tabix -p vcf $f/$name.muts.ids.mnv.flagged.v1.target.pass.vcf.gz" ;
done
