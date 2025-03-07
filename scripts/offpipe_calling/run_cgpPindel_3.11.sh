#! /bin/bash

# THe script is aimed to be run as 
# run_DERMATLAS_ffpe_pindel.sh tumour_PDID normal_PDID num_CPU_to_USE BAM_DIR REQUIRED_FLAG_DIR
#INPUT parameters
tumour=$1
normal=$2
CPU_TO_USE=$3
# BAM directory
export BAMDIR=$4
# Requiered Rules directory
export PINDELRULES=$5

####################################################################################
# Special environment variables required for cgpPindel v3.90
####################################################################################
# set this to your lustre area
export WORK=$PWD
# how many CPUs you want to use
export CPU_TO_USE
# #################  DERMATLAS ONLY special values for this variables
# this defines the root of the reference area for the relevant species/build used in DERMATLAS
export REF_BASE=/lustre/scratch124/casm/team78pipelines/canpipe/live/ref/human/GRCh38_full_analysis_set_plus_decoy_hla
# path to germline indels bedfile
export GERM_INDEL=/lustre/scratch124/casm/team113/projects/7688_7687_Gen_Effects_CDS2_loss_Uveal_melanoma/7688_3365_Gen_Effects_CDS2_loss_Uveal_melanoma_WES/resources/caveman/example/empty_germline_bed_for_caveman.bed

# the sequencing type being processed (WGS|WXS|TGS) WXS - Whole Exome for DERMATLAS
export SEQ_TYPE=WXS
# Ensembl version for annotation - DERMATLAS is  ENSv103
export ENSM_VER=103
# this indicates patterns of contigs to be IGNORED, different for every build
# this is cgppindel option -e chrUn%,HLA%,%_alt,%_random,chrM,chrEBV


# paths to the BAM/CRAM files on lustre
export MT_BAM=$BAMDIR/$tumour/Filtered_bams/$tumour.sample.dupmarked.mXfilt_Filtered.bam 
export WT_BAM=$BAMDIR/$normal/Filtered_bams/$normal.sample.dupmarked.mXfilt_Filtered.bam

#The name of the ouput folder
pair="$tumour"
mkdir -p $WORK/$pair

printf "$1 \n $2 \n $3 \n BAMDIR: $BAMDIR \n MT_BAM: $MT_BAM \n Projdir:$WORK\n OUTPUT_DIR:$WORK/$pair \n FLAGs_used: $PINDELRULES \n"

#Run the variant calling with Pindel
module load cgppindel/3.11.0
pindel.pl \
  -o $WORK/$pair \
  -r $REF_BASE/genome.fa \
  -t $MT_BAM \
  -n $WT_BAM \
  -s $REF_BASE/pindel/simpleRepeats.bed.gz \
  -f ${PINDELRULES} \
  -g $REF_BASE/vagrent/e${ENSM_VER}/codingexon_regions.indel.bed.gz \
  -u $REF_BASE/pindel/pindel_np.v5.gff3.gz \
  -st $SEQ_TYPE \
  -e chrUn%,HLA%,%_alt,%_random,chrM,chrEBV \
  -b $REF_BASE/shared/HiDepth_mrg1000_no_exon_coreChrs_v3.bed.gz \
  -sf $REF_BASE/pindel/softRulesFragment.lst \
  -c $CPU_TO_USE
