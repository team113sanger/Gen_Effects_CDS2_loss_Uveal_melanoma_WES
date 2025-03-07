#!/bin/bash

vcflist=$1

# Versions

ens_ver=103
species=homo_sapiens
assembly=GRCh38
vep=ensembl_vep/103.1

# Annotation resources

REFBASE=/lustre/scratch124/casm/team78pipelines/reference/human/GRCh38_full_analysis_set_plus_decoy_hla/

ref=$REFBASE/genome.fa
vep_cache=$REFBASE/vep/cache/103
gnomadfile=$REFBASE/vep/gnomad/v3.1.2/gnomad.genomes.v3.1.2.short.vcf.gz
clinvarfile=$REFBASE/vep/clinvar/20230121/clinvar_20230121.chr.canonical.vcf.gz
dbsnpfile=$REFBASE/vep/dbsnp/155/dbSNP155.GRCh38.GCF_000001405.39.mod.vcf.gz
cosmicfile=$REFBASE/vep/cosmic/v97/CosmicV97Coding_Noncoding.normal.counts.vcf.gz

# vep custom annotation command line

custom="--custom $cosmicfile,Cosmic,vcf,exact,0,CNT, --custom $clinvarfile,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT --custom $dbsnpfile,dbSNP,vcf,exact,0, --custom $gnomadfile,gnomAD,vcf,exact,0,FLAG,AF"

# Check directories

if [[ -z $vcflist ]]; then
	echo -e "\nUsage: $0 vcflist\n"
	exit 0
fi

mkdir -p logs

# Check input files and directories

for d in $REFBASE $vep_cache; do
	if [[ ! -d $d ]]; then
		echo "Directory $d does not exist"
		exit
	fi
done

for annotfile in $ref $gnomadfile $clinvarfile $dbsnpfile $cosmicfile; do
	if [[ ! -e $annotfile ]]; then
		echo "File $annotfile does not exist"
		exit
	fi
done

# Submit jobs

for f in `cat $vcflist`; do
	if [[ ! -e $f ]]; then
		echo "File $f does not exist"
		exit
	fi
	outdir=`dirname $f`
	filename=`basename $f .vcf.gz`
	vcf_out=$outdir/$filename.vep.vcf.gz
	sample=`echo $filename | cut -f 1 -d "."`
##	echo "$outdir $filename $vcf_out $f"
	bsub -n 4 -q normal -M4000 -R"select[mem>4000] rusage[mem=4000]" -e logs/vep.$sample.e -o logs/vep.$sample.o "module load bcftools-1.9/python-3.11.6 $vep; export SINGULARITY_BINDPATH=\"/lustre,/nfs/casm/team113da/users/kw10\"; vep -i $f --cache_version $ens_ver  -t SO --format vcf -o $vcf_out  --cache --dir $vep_cache --buffer 20000 --species $species --offline --symbol --biotype --vcf --sift b --no_stats --assembly $assembly  --flag_pick_allele_gene --canonical --hgvs --shift_hgvs 1 --fasta $ref $custom --compress_output bgzip --mane --numbers --protein --polyphen p --transcript_version --show_ref_allele --fork 4 --domains --force_overwrite && tabix -p vcf $vcf_out"
done

