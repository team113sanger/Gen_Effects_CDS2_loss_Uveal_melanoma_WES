#! /bin/bash
# THe script is aimed to be runned as the first step in the MNV calling from Unflagged SNV VCF files  
# run_mnv_casmsmartphase0.1.8.sh tumour_PDID normal_PDID sample_name VCF_DIR BAMDIR
# run_mnv_casmsmartphase0.1.8.sh PD52540a PD52540b PD52540a_vs_D52540b
#INPUT parameters
tumour=$1
normal=$2
sample=$3
VCF_DIR=$4
BAMDIR=$5


#The name of the ouput folder
pair=$tumour"_vs_"$normal
samplefolder="$tumour"
# Create the path that the unflagged VCF is going to have
export VCF_IN=$VCF_DIR/$samplefolder/$pair.muts.ids.vcf.gz
export BAM=$BAMDIR/$tumour/Filtered_bams/$tumour.sample.dupmarked.mXfilt_Filtered.bam 
export OUTPHASED=$VCF_DIR/$samplefolder/$pair.phased
export BED=$VCF_DIR/$samplefolder/$pair.adjacent_snvs.bed

printf "Tumour:$1 \n Normal:$2 \n  Sample_name:$3 \n Pair: $pair \n VCF_DIR:$VCF_DIR \n BAMDIR: $BAMDIR \n VCF_IN: $VCF_IN \n BED_out: $BED \n OUTPUT: $OUTPHASED \n "

echo "module load casm-smart-phase/0.1.8;  smart-phase  -g $BED  -p TUMOUR -r $BAM -m0 -x -o $OUTPHASED -a $VCF_IN "

module load casm-smart-phase/0.1.8

smart-phase  -g $BED  -p TUMOUR -r $BAM -m0 -x -o $OUTPHASED -a $VCF_IN

