#! /bin/bash
# THe script is aimed to be runned as the first step in the MNV calling from Unflagged SNV VCF files  
# run_mnv_casmsmartphase0.1.8_mergemnvs.sh tumour_PDID normal_PDID sample_name VCF_DIR 
# run_mnv_casmsmartphase0.1.8_mergemnvs.sh PD52540a PD52540b PD52540a_vs_D52540b
#INPUT parameters
tumour=$1
normal=$2
sample=$3
VCF_DIR=$4


#The name of the ouput folder
pair=$tumour"_vs_"$normal
samplefolder="$tumour"
# Create the path that the unflagged VCF is going to have
export VCF_IN=$VCF_DIR/$samplefolder/$sample.muts.ids.vcf.gz
export OUTPUT_MNV_VCF=$VCF_DIR/$samplefolder/$sample.muts.ids.mnv.vcf.gz
export PHASED=$VCF_DIR/$samplefolder/$sample.phased
export BED=$VCF_DIR/$samplefolder/$sample.adjacent_snvs.bed

printf "Tumour:$1 \n Normal:$2 \n  Sample_name:$3 \n Pair: $pair \n VCF_DIR:$VCF_DIR \n OUTPUT_MNV_VCF: $OUTPUT_MNV_VCF \n VCF_IN: $VCF_IN \n BED_out: $BED \n Phased_input: $PHASED \n "

echo "module load casm-smart-phase/0.1.8;  casmsmartphase merge-mnvs -f $VCF_IN -o $MNV_VCF -p $PHASED -b $BED -c 0.1 -x 30 &>$sample.mergelog "

module load casm-smart-phase/0.1.8

casmsmartphase merge-mnvs -f $VCF_IN -o $OUTPUT_MNV_VCF -p $PHASED -b $BED -c 0.1 -x 30 &>$VCF_DIR/$samplefolder/$sample.mergelog 
