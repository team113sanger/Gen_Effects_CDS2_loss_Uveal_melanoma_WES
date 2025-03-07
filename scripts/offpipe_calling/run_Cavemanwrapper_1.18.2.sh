#! /bin/bash

# THe script is aimed to be run as 
# run_Cavemanwrapper.sh tumour_PDID normal_PDID num_CPU_to_USE BAM_DIR
#Run the variant calling with CaVEMan  - for SNVs and MNVs
#INPUT parameters
export tumour=$1
export normal=$2
export CPU_TO_USE=$3
# BAM directory
export BAMDIR=$4

####################################################################################
# Special environment variables required for caveman.pl v1.18.2 cavemanwrapper v18.2
####################################################################################
# set this to your lustre area
export WORK=$PWD
# how many CPUs you want to use
export CPU_TO_USE
# #################  DERMATLAS ONLY special values for this variables
# this defines the root of the reference area for the relevant species/build used in DERMATLAS
export REF_BASE=/lustre/scratch124/casm/team78pipelines/reference/human/GRCh38_full_analysis_set_plus_decoy_hla
# the sequencing type being processed (pulldown|exome|genome|genomic|followup|targeted|rna_seq) - Whole Exome for DERMATLAS
export SEQ_TYPE=WXS
# the sequencing type being processed  [WGS|WXS|RNA|RNA-Seq|AMPLICON|TARGETED] - Whole Exome for the project
export SEQ_TYPE_TUM_PROT=WXS
# the sequencing type being processed  [WGS|WXS|RNA|RNA-Seq|AMPLICON|TARGETED] - Whole Exome for theproject
export SEQ_TYPE_NORM_PROT=WXS
# Ensembl version for annotation - DERMATLAS is  ENSv103
export ENSM_VER=103
# path to germline indels bedfile
export GERM_INDEL=/lustre/scratch124/casm/team113/projects/7688_7687_Gen_Effects_CDS2_loss_Uveal_melanoma/7688_3365_Gen_Effects_CDS2_loss_Uveal_melanoma_WES/resources/caveman/example/empty_germline_bed_for_caveman.bed

export EXCLUDE_CONTIGS_FILE=$REF_BASE/caveman/ignore_contigs_caveman.txt
ln -fs $REF_BASE/genome.fa $WORK/. && grep -vFwf $EXCLUDE_CONTIGS_FILE $REF_BASE/genome.fa.fai > $WORK/genome.fa.fai

# paths to the BAM/CRAM files on lustre
export MT_BAM=$BAMDIR/$tumour/Filtered_bams/$tumour.sample.dupmarked.mXfilt_Filtered.bam 
export WT_BAM=$BAMDIR/$normal/Filtered_bams/$normal.sample.dupmarked.mXfilt_Filtered.bam

# paths to the copy number bed  -Since ASCAT wasn't run on these samples I followed the recommended defaults in confluence 
# https://confluence.sanger.ac.uk/display/CIH/Manual+execution+-+cgpCaVEManWrapper
# export MT_CN=
# export WT_CN=$WORK/wildtype.cn.bed

#The name of the ouput folder
pair="$tumour"
mkdir -p $WORK/$pair
# LOGDIR= $WORK/$pair/logs
# mkdir -p $LOGDIR
printf "$1 \n $2 \n $3 \n BAMDIR: $BAMDIR \n MT_BAM: $MT_BAM \n Projdir:$WORK\n OUTPUT_DIR:$WORK/$pair \n "
echo "caveman.pl 
  -reference $WORK/genome.fa.fai 
  -outdir $WORK/$pair 
  -tumour-bam $MT_BAM 
  -normal-bam $WT_BAM 
  -ignore-file $REF_BASE/genome.gap.tab 
  -tum-cn-default 5 
  -norm-cn-default 2 
  -species Human 
  -species-assembly GRCh38 
  -flag-bed-files $REF_BASE/caveman/flagging 
  -germline-indel $GERM_INDEL 
  -unmatched-vcf $REF_BASE/caveman/unmatched_v5 
  -seqType $SEQ_TYPE 
  -tumour-protocol $SEQ_TYPE_TUM_PROT 
  -normal-protocol $SEQ_TYPE_NORM_PROT 
  -normal-contamination 0.1 
  -noflag
  -flagConfig $REF_BASE/caveman/flag.vcf.config.ini 
  -flagToVcfConfig $REF_BASE/caveman/flag.to.vcf.convert.ini 
  -annot-bed-files $REF_BASE/vagrent/e$ENSM_VER 
  -threads $CPU_TO_USE "

#Move to the Pair directory to run the caveman command 
# cd $WORK/$pair
#Run the variant calling with CaVEMan v1.18.2
module load cgpcavemanwrapper/1.18.2
caveman.pl \
  -reference $WORK/genome.fa.fai \
  -outdir $WORK/$pair \
  -tumour-bam $MT_BAM \
  -normal-bam $WT_BAM \
  -ignore-file $REF_BASE/genome.gap.tab \
  -tum-cn-default 5 \
  -norm-cn-default 2 \
  -species Human \
  -species-assembly GRCh38 \
  -flag-bed-files $REF_BASE/caveman/flagging \
  -germline-indel $GERM_INDEL \
  -unmatched-vcf $REF_BASE/caveman/unmatched_v5 \
  -seqType $SEQ_TYPE \
  -tumour-protocol $SEQ_TYPE_TUM_PROT \
  -normal-protocol $SEQ_TYPE_NORM_PROT \
  -normal-contamination 0.0  \
  -noflag \
  -flagConfig $REF_BASE/caveman/flag.vcf.config.ini \
  -flagToVcfConfig $REF_BASE/caveman/flag.to.vcf.convert.ini \
  -annot-bed-files $REF_BASE/vagrent/e$ENSM_VER \
  -threads $CPU_TO_USE 








