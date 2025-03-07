#! /bin/bash
# THe script is aimed to be runned as the first step in the MNV calling from Unflagged SNV VCF files  
# run_mnv_casmsmartphase_generate_bed.sh tumour_PDID normal_PDID sample_name VCF_DIR 
# run_mnv_casmsmartphase_generate_bed.sh PD52540a PD52540b PD52540a_vs_D52540b
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
export BED_OUT=$VCF_DIR/$samplefolder/$sample.adjacent_snvs.bed

printf "Tumour:$1 \n Normal:$2 \n  Sample_name:$3 \n Pair: $pair \n VCF_DIR:$VCF_DIR \n VCF_IN: $VCF_IN \n BED_out: $BED_OUT "

echo "module load casm-smart-phase/0.1.8;  casmsmartphase generate-bed --markhz -f $VCF_IN -o $BED_OUT "

module load casm-smart-phase/0.1.8

casmsmartphase generate-bed --markhz -f $VCF_IN -o $BED_OUT

