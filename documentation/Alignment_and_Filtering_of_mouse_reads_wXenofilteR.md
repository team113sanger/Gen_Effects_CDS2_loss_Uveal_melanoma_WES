# Alignment and Filtering of mouse reads from Human BAM file from Xenografted samples

## Overview

This document describes the steps taken to align and filter out mouse reads from Human BAM files from Xenografted samples. Mapping parameters used to align reads against the Human GRCh38 reference genome and the mouse reference used [NOD_ShiLtJ_V1_PDX](../reference/NOD_ShiLtJ_V1_PDX_ref/README.md) were the same.

All the scripts and code mentioned below can be found in the `scripts` directory.

## Alignment to the Human GRCh38 reference genome 

The WES sequencing data was aligned to the GRCh38 Human reference genome using `bwa-mem`. PCR duplicates were marked using `samtools markdup` function. The same process was applied for all the samples from both targeting experiments in this project. This process was preformed through an internal pipeline.

## Filtering of mouse reads with Xenofilter

This process was used to remove mouse reads from the Human BAM files with the whole exome sequencing data generated for the OMM2.5 xenografts.

### Required Environment variables and software

The following environment variables are required to be set before running the scripts:
- **PROJECTDIR**: The path to the project directory where this repo got cloned into
- **STUDY**: The study ID,  7688 for this analysis
- **PROJECTID**: The project ID, 3365s for this analysis

The following software is required to be installed and visible in the path before running the scripts:
- **R**: R `4.2.2`
- **samtools**: samtools `v1.14`
- **bwa-mem**: bew-mem `v0.7.17`
- **XenofilteR**: XenofilteR `v1.6`

- Load the following variables and software
```bash
PROJECTDIR=/lustre/7688_3365_Gen_Effects_CDS2_loss_Uveal_melanoma_WES
STUDY=7688
PROJECTID=3365

# Example of a source me file containing the required environment variables and software to load
# iRODs is used internally to access the sequencing data but this should be replaced for the human unfiltered BAM files 
source ${PROJECTDIR:?unset}/scripts/pdx_processing/source_me.sh
```
#### R environment
If you're interested in reproducing the R environment, for a the code used in R 4.2.2 run the following commands within `R v4.2.2`, change the path on `projectdir` to the path where the repository was cloned into:

```R
projectdir<-"/lustre/7688_3365_Gen_Effects_CDS2_loss_Uveal_melanoma_WES"
pdx_processing_dir<- file.path(projectdir,"scripts/pdx_processing")

setwd(pdx_processing_dir)
install.packages("renv")
library(renv)
renv::restore(lockfile=file.path(pdx_processing_dir, "renv.lock")) # To rebuild an environment from the renv.lockfile

```

### CRAM to FASTQ files

#### Manifest generation and import of sequencing data

To import the sequencing metadata from iRODS to generate manifests with the sequencing statistics and information
we ran the following: 

**IMPORTANT**: All the manifest generated can be found within the `metadata/manifests` directory.

```bash
PROJECTDIR=/lustre/7688_3365_Gen_Effects_CDS2_loss_Uveal_melanoma_WES
STUDY=7688
PROJECTID=3365

#load iRODS module 
module load IRODS/1.0
# initiate session 
iinit

#To call the script to build the manifests with project ID Name
mkdir -p ${PROJECTDIR:?unset}/metadata/manifests
/software/team113/dermatlas/R/R-4.2.2/bin/Rscript ${PROJECTDIR:?unset}/scripts/pdx_processing/Build_manifest_from_irods_cram_information.R --seqscape_proj_id ${STUDY} --outdir ${PROJECTDIR:?unset}/metadata/manifests
```

- After generating the manifest, we split the information to only contain the Xenografted samples (CDS2_Tumour) 

**OUTPUTS**:
- **7688_cram_manifest_INFO_from_iRODS_all.txt** : Contains the information of all the samples across both SW837 and OMM2.5 targeting experiments.
- **7688_cram_manifest_INFO_from_iRODS.txt** : Contains the information of the samples on the OMM2.5 targeting experiments.

```bash 
# Reformat the manifest and filter 
#Then keep only the  CDS2_Tumour samples and the header 
mv ${PROJECTDIR:?unset}/metadata/manifests/${STUDY}_cram_manifest_INFO_from_iRODS.txt ${PROJECTDIR:?unset}/metadata/manifests/${STUDY}_cram_manifest_INFO_from_iRODS_all.txt
grep -E 'CDS2_Tumour|sample' ${PROJECTDIR:?unset}/metadata/manifests/${STUDY}_cram_manifest_INFO_from_iRODS_all.txt >${PROJECTDIR:?unset}/metadata/manifests/${STUDY}_cram_manifest_INFO_from_iRODS.txt
```

- To generate the list of jobs to transform the CRAM files to fastq files the files we ran

Output: 
- `scripts/7688_cramtofastq_from_iRODs_jobs.sh` : Contains the list of jobs to transform the CRAM to fastq files

```bash
PROJECTDIR=/lustre/7688_3365_Gen_Effects_CDS2_loss_Uveal_melanoma_WES
STUDY=7688
cd ${PROJECTDIR:?unset}/scripts/pdx_processing/

# Load environment with requiring 
source ${PROJECTDIR:?unset}/scripts/pdx_processing/source_me.sh

#This script takes the cram manifest and generates the SH file with the jobs to import and transform to fastq all of the cram files from iRODs
Rscript ${PROJECTDIR:?unset}/scripts/pdx_processing/cramtofastq_from_iRODs_based_cram_manifest.R --manifest ${STUDY}_cram_manifest_INFO_from_iRODS.txt --projectdir ${PROJECTDIR:?unset} --studyID ${STUDY} --mem 16000

```

- Then we proceed to execute the jobs import the BAM files and transform them into fastqs using samtools

```bash
PROJECTDIR=/lustre/7688_3365_Gen_Effects_CDS2_loss_Uveal_melanoma_WES
STUDY=7688
cd ${PROJECTDIR:?unset}
# Load environment with requring 
source ${PROJECTDIR:?unset}/scripts/pdx_processing/source_me.sh
# Login into iRODs
#iinit 
sh ${PROJECTDIR:?unset}/scripts/${STUDY}_cramtofastq_from_iRODs_jobs.sh
```

#### Generate the mouse genome reference files and bwa index

To be able to generate the mouse referenced use the steps mentioned in the [**NOD_ShiLtJ_V1_PDX_ref**](../reference/NOD_ShiLtJ_V1_PDX_ref/README.md) README file. 


#### Mapping against the NOD_ShiLtJ_V1_PDX mouse reference genome with bwa-mem 

To generate the jobs to map the fastq files against the mouse reference genome, we used the `PDX_bwa_mem_mapping_jobs_from_master_manif.R` script:

- **INPUT**: Use the file : `metadata/manifests/7688_cram_manifest_INFO_from_iRODS_wbam_counts_qc.txt`

```bash
PROJECTDIR=/lustre/7688_3365_Gen_Effects_CDS2_loss_Uveal_melanoma_WES
STUDY=7688
cd ${PROJECTDIR:?unset}/scripts/pdx_processing

# Load environment with requiring 
source ${PROJECTDIR:?unset}/scripts/pdx_processing/source_me.sh

#This scrip
Rscript ${PROJECTDIR:?unset}/scripts/pdx_processing/PDX_bwa_mem_mapping_jobs_from_master_manif.R --manifest ${STUDY}_cram_manifest_INFO_from_iRODS_wbam_counts_qc.txt --projectdir ${PROJECTDIR:?unset} --referencedir ${PROJECTDIR:?unset}/reference/NOD_ShiLtJ_V1_PDX_ref
```
This will generate two outputs:
 1. `bwamem_mapping_perlanrun_to_NOD_PDXV1_tum_only_jobs.sh`: which contains with the list of jobs to perform the mapping against the mouse reference genome with **bwa-mem**.
 2. `samtools_psample_merge_nodv1_tum_only_jobs.sh`: which contains the list of jobs to merge and index the BAM files per sample.

Submit the remapping jobs with the mouse reference using **bwa-mem**

```bash
PROJECTDIR=/lustre/7688_3365_Gen_Effects_CDS2_loss_Uveal_melanoma_WES
STUDY=7688
cd ${PROJECTDIR:?unset}/scripts/pdx_processing

# Load environment with requring 
source ${PROJECTDIR:?unset}/scripts/pdx_processing/source_me.sh

/bin/sh bwamem_mapping_perlanrun_to_NOD_PDXV1_tum_only_jobs.sh
```

Submit the merging per sample

```bash
PROJECTDIR=/lustre/7688_3365_Gen_Effects_CDS2_loss_Uveal_melanoma_WES
STUDY=7688
cd ${PROJECTDIR:?unset}/scripts/pdx_processing

# Load environment with requring 
source ${PROJECTDIR:?unset}/scripts/pdx_processing/source_me.sh

/bin/sh samtools_psample_merge_nodv1_tum_only_jobs.sh
```

### Filter the mouse reads from the Xenofilter jobs for the samples

First, we generated the manifest with the final filtered file name and example of the Xenofilter jobs for each sample if a single per samples file could have been used using the `run_Xenofilter_from_WES_master_manif.R` script.  

However, given the amount of sequencing data generated for the grafted samples, we required to split both Human and mouse BAM files with the same matching read numbers and read names to be able to run XenofilteR successfully.  This is due to a known issue with XenofilteR when and one of it's dependencies, `Rsamtools`, see issue [here](https://github.com/NKI-GCF/XenofilteR/issues/7) 


We generate the manifest:

**INPUT**: `metadata/manifests/7688_cram_manifest_INFO_from_iRODS_wbam_counts_qc_psamp_mouse.txt`

```bash
PROJECTDIR=/lustre/7688_3365_Gen_Effects_CDS2_loss_Uveal_melanoma_WES
STUDY=7688

cd ${PROJECTDIR}/scripts/pdx_processing

# Load environment with requring 
source ${PROJECTDIR}/scripts/pdx_processing/source_me.sh

# Set the job to create the .sh XenofilteR job submissions 
Rscript ${PROJECTDIR}/scripts/pdx_processing/run_Xenofilter_from_WES_master_manif.R --manifest ${STUDY}_cram_manifest_INFO_from_iRODS_wbam_counts_qc_psamp_mouse.txt --projectdir ${PROJECTDIR} --outdir ${PROJECTDIR}/bams/WES_xfilt
```
**OUTPUTS**:
 `metadata/manifests/7688_cram_manifest_INFO_from_iRODS_wbam_counts_qc_psamp_mouse_xfb.txt`


#### Split input BAM files, split by read names from canapps bam file and run XenofilteR for PDX filtering from Manifest

To split the reads we used the script: **split_bam_files_and_get_xenofilter_jobs.R**. 

The script takes the a manifest with BAM file information and generates the jobs to take the unfiltered human BAM files, sort by read name, obtain a plain text file with the read names for all the reads present in the file, then it will split this file in the number of files required that have a maximum of `NREADS_SPLIT` per file. Subsequently it will split the Human and Mouse BAM files by read names and then generate the jobs to run XenofilteR to filter out the mouse reads for the matching Human & Mouse BAM files.

**INPUT**: `metadata/manifests/7688_cram_manifest_INFO_from_iRODS_wbam_counts_qc_psamp_mouse_xfb.txt`

```bash
PROJECTDIR=/lustre/7688_3365_Gen_Effects_CDS2_loss_Uveal_melanoma_WES
STUDY=7688

cd ${PROJECTDIR:?unset}

# Load environment with requiring 
source ${PROJECTDIR:?unset}/scripts/pdx_processing/source_me.sh

# Variables:
NCORES=4
NREADS_SPLIT=50000000

Rscript ${PROJECTDIR:?unset}/scripts/pdx_processing/split_bam_files_and_get_xenofilter_jobs.R --study_id ${STUDY:?unset} \
--manifest ${PROJECTDIR:?unset}/metadata/manifests/${STUDY}_cram_manifest_INFO_from_iRODS_wbam_counts_qc_psamp_mouse_xfb.txt \
--projectdir ${PROJECTDIR:?unset} \
--xfilter_outdir ${PROJECTDIR:?unset}/bams/WES_xfilt/NOD_PDXV1 \
--split_nreads ${NREADS_SPLIT:?unset} \
--nthreads ${NCORES} \

```

**OUTPUTS**:
The scripts will take the unfiltered human BAM files, sort by read name, split by read name  following output files:
- [`7688_cram_manifest_INFO_from_iRODS_wbam_counts_qc_psamp_mouse_xfb_part.txt`](../metadata/manifests/7688_cram_manifest_INFO_from_iRODS_wbam_counts_qc_psamp_mouse_xfb_part.txt): Manifest that contains the metadata and file name information of the split files for all the samples.

Runner files:
- **7688_bamsplittin_by_read_names_jobs.sh** : Contains the list of jobs to split the BAM files by read names using `samtools` 
- **7688_Xenofilter_by_read_names_parts_jobs.sh** : Contains the list of jobs used to run XenofilteR for the split Human and Mouse BAM files using the script [`scripts/pdx_processing/bam_Xenofilter_rg.R`](../scripts/pdx_processing/bam_Xenofilter_rg.R)
- **7688_Xenofilter_merged_filtered_parts_jobs.sh**: Contains the list of jobs to merge the filtered BAM files per sample and index them using `samtools`

To submit the BAM file splitting by read names:

```bash
PROJECTDIR=/lustre/7688_3365_Gen_Effects_CDS2_loss_Uveal_melanoma_WES
STUDY=7688

cd ${PROJECTDIR:?unset}

# Load environment with requring 
source ${PROJECTDIR:?unset}/scripts/pdx_processing/source_me.sh

bash ${PROJECTDIR:?unset}/scripts/pdx_processing/${STUDY}_bamsplittin_by_read_names_jobs.sh

```

Submit Xenofilter filtering of split bamfiles using `7688_Xenofilter_by_read_names_parts_jobs.sh`

```bash
PROJECTDIR=/lustre/7688_3365_Gen_Effects_CDS2_loss_Uveal_melanoma_WES
STUDY=7688

cd ${PROJECTDIR:?unset}

# Load environment with requring 
source ${PROJECTDIR:?unset}/scripts/pdx_processing/source_me.sh

bash ${PROJECTDIR:?unset}/scripts/pdx_processing/${STUDY}_Xenofilter_by_read_names_parts_jobs.sh
```

Submit filtering of BAM files with `XenofilteR`

```bash
PROJECTDIR=/lustre/7688_3365_Gen_Effects_CDS2_loss_Uveal_melanoma_WES
STUDY=7688

cd ${PROJECTDIR:?unset}

# Load environment with requring 
source ${PROJECTDIR:?unset}/scripts/pdx_processing/source_me.sh

bash ${PROJECTDIR:?unset}/scripts/pdx_processing/${STUDY}_Xenofilter_merged_filtered_parts_jobs.sh
```

### Collate filtered read stats and plots

Finally, once all jobs are complete we collate information on the proportion of filtered reads per file per sample and plot it using the script `collate_xfilter_stats_from_manif_and_plots.R `

INPUT: `metadata/manifests/7688_cram_manifest_INFO_from_iRODS_wbam_counts_qc_psamp_mouse_xfb_part.txt`

```bash
PROJECTDIR=/lustre/7688_3365_Gen_Effects_CDS2_loss_Uveal_melanoma_WES
STUDY=7688

cd ${PROJECTDIR:?unset}

# Load environment with requring 
source ${PROJECTDIR:?unset}/scripts/pdx_processing/source_me.sh

# Run the collation of the xenofiltered stats information and plots
Rscript ${PROJECTDIR:?unset}/scripts/pdx_processing/collate_xfilter_stats_from_manif_and_plots.R --study_id ${STUDY:?unset} \
--manifest ${PROJECTDIR:?unset}/metadata/manifests/${STUDY}_cram_manifest_INFO_from_iRODS_wbam_counts_qc_psamp_mouse_xfb_part.txt \
--projectdir ${PROJECTDIR:?unset} \
--xfilter_outdir ${PROJECTDIR:?unset}/analysis/bam_xfilterstats 
```

**OUTPUTS**:
- [`7688_cohort_per_part.xfiltstats`](../analysis/bam_xfilterstats/7688_cohort_per_part.xfiltstats): Contains the information of the proportion of filtered reads per file per sample.
- [`7688_cohort_per_sample.xfiltstats`](../analysis/bam_xfilterstats/7688_cohort_per_sample.xfiltstats): Contains the information of total of filtered reads per sample.
- [`7688_per_sample_pdx_xfiltstats_bpl.pdf`](../analysis/bam_xfilterstats/7688_per_sample_pdx_xfiltstats_bpl.pdf): Bar plot with the percentage of reads filtered per sample.
- [`7688_per_sample_pdx_xfiltstats_vpl.pdf`](../analysis/bam_xfilterstats/7688_per_sample_pdx_xfiltstats_vpl.pdf): Bar plot with the percentage of reads filtered per sample.

