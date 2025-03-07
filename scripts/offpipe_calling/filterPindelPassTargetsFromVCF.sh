#!/bin/bash
#This script needs in order the parameters:
#   1) sample.list : file with the sample.list file with the information of sample pairs in the format
#   2) EXTended_baits_bedfile: BED file with the coordinates of the baits used e.g. "/lustre/CASM_BAITS/SureSelect_Whole_Human_Exome_v5_GRCh38_liftover_160.sorted.w100.canonical.bed"
mkdir -p logs

samplelist=$1
target=$2


for f in `cat $samplelist`; do 
    echo $f;
    if [[ ! -e $f ]]; then
		echo "File $f does not exist"
		exit
	fi
	outdir=`dirname $f`
	filename=`basename $f .vcf.gz`
	vcf_out=$outdir/$filename.target.pass.vcf.gz
	sample=`echo $filename | cut -f 1 -d "."`
    #This changes the names of the folders from Tumour-vs-Normal to Tumour_vs_Normal
    name=`echo $f | sed 's/-/_/g'`;
    #sleep 0.5s;
    bsub -e logs/filt.target.$sample.e -o logs/filt.target.$sample.o -M4000 -n 2 -R"select[mem>4000] rusage[mem=4000] span[hosts=1]" "module load bcftools-1.9/python-3.11.6; bcftools view -f PASS $f -Oz -o ${vcf_out} -T $target; tabix -p vcf ${vcf_out}" ;
done
