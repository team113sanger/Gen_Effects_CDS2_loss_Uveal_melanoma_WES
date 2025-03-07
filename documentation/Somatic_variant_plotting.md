# Results of the somatic Variants plotting

## Overview

This document is provides the code used to filter the MAF file containing all the variants identified with the somatic callers that passed each callers standard filters and flagging by Variant Allele Frequency (VAF>0.04). Following this filter, it generates a MAF file per experiment, with the variants used in the paper for the results of _CDS2_ targeting experiments on SW837 cell-lines . Finally it also generates the plots for the VAF and depth of the variants found in the SW837 cell line. Results on the and OMM2.5 xenografted cell lines are included in this repository.

## Required dependencies

### R environment setup for somatic variant plotting
To reproduce the R environment used to generate the figures install the packages required to generate the final files and plots, the following steps are required run the following commands within `R v4.2.2`, change the path on `projectdir` to the path where the repository was cloned into:

Set the `projectdir` variable to the path where the repository was cloned into and run the following commands in R to install the required packages and restore the environment from the `renv.lock` file.

```R
# Run in R v4.2.2
projectdir<-"/lustre/7688_3365_Gen_Effects_CDS2_loss_Uveal_melanoma_WES"
plotsc_dir<- file.path(projectdir,"scripts/somatic_variant_plots")

setwd(plotsc_dir)
install.packages("renv")
library(renv)
# To rebuild an environment from the renv.lockfile
# use renv::restore(lockfile=file.path(plotsc_dir, "renv.lock"), prompt=FALSE) if you want R not to prompt you for confirmation
renv::restore(lockfile=file.path(plotsc_dir, "renv.lock"))

```

### MAF repository cloning 

Using `git` [(see here for how to install git)](https://github.com/git-guides/install-git) clone the `MAF` repository into the `scripts` directory of the project. The following commands will clone the repository `0.5.4` branch. This repository contains the script `plot_vaf_vs_depth_from_maf.R` used to generate the VAF plots.

Set the `PROJECTDIR` variable to the path where the repository was cloned into and run the following commands in the terminal:

```bash
PROJECTDIR=/lustre/7688_3365_Gen_Effects_CDS2_loss_Uveal_melanoma_WES
STUDY=7688
PROJECT=3365

cd ${PROJECTDIR}/scripts

# Clone the MAF repository
git clone -b '0.5.4' git@github.com:team113sanger/dermatlas_analysis_maf.git MAF

```

## Result reproduction

### Variant filtering and plotting

#### Required input files 

To generate the filtered MAF files, tile plots and VAF plots used in the paper, the following input files are required:
- [**7688_3365_sample_exp_metadata.tsv**](../metadata/7688_3365_sample_exp_metadata.tsv)
- [**keep_caveman_pindel_all.maf**](../analysis/variants_combined/version2/all/keep_caveman_pindel_all.maf): contains all the variants for all samples identified by the somatic callers; all variants in coding and splice sites regions including synonymous variant are included
- [**OMM2.5_mutated_cancer_gene_census.v97.genes_names.tsv**](../resources/COSMIC/OMM2.5_mutated_cancer_gene_census.v97.genes_names.tsv): contains the list of genes in the COSMIC v97 [Cancer Gene Census](https://cancer.sanger.ac.uk/census) list that were mutated in OMM2.5 cells. 

#### Execution
```bash
PROJECTDIR=/lustre/7688_3365_Gen_Effects_CDS2_loss_Uveal_melanoma_WES
STUDY=7688
PROJECT=3365

cd ${PROJECTDIR}/scripts/somatic_variant_plots/

Rscript ${PROJECTDIR:?unset}/scripts/somatic_variant_plots/somatic_mutations_oncoplots.R 

```
**OUTPUTS**:
The following files are generated in the `somatic_variant_plots` folder and are named according to the experiment and the cell line:
**SW837 cell line results**:
- **SW837_keep_matched_caveman_pindel_all_indelvaf_filt.maf**: MAF file containing the filtered variants found in the SW837 whole exome.
- **SW837_plots**: Folder containing tile plots generated for the SW837 cell line results
    - **SW837_lines_oncoplot.pdf**: OncoPlot showing all the mutations across the _CDS2_ targeted  SW837 cell lines - Extended Data Figure 12 
    - **SW837_lines_lolliplot_CDS2.pdf**: Lollipop plot showing the mutations in the _CDS2_ gene in the SW837 cell lines 

**OMM2.5 cell line results**:
- **OMM2.5_CDS2_Tumours_keep_matched_caveman_pindel_all_indelvaf_filt.maf**: MAF file containing the filtered variants found in the OMM2.5 whole exome.
- **OMM2.5_plots**: Folder containing tile plots generated for the OMM2.5 cell line results
- **CDS2_nontreated_vs_treated_plots**: Folder containing the lollipop plots generated with the comparison of mutations found between doxycycline treated and non treated OMM2.5 grafts whole exome. 

### VAF vs Depth plots for SW837 cell line

#### Required input files

To generate the VAF plot with the filtered variants found in SW837 cells used in the paper for Extended Data Figure 12: 

```bash
PROJECTDIR=/lustre/7688_3365_Gen_Effects_CDS2_loss_Uveal_melanoma_WES
STUDY=7688
PROJECT=3365

SOMAT_VAR_PDIR="${PROJECTDIR}/analysis/somatic_variant_plots"
SCRIPTDIR="${PROJECTDIR}/scripts/MAF"

cd ${SOMAT_VAR_PDIR}

#Get the list of samples 
grep "SW837_C9" ${PROJECTDIR}/metadata/${STUDY}_${PROJECT}-analysed_all.tsv >${SOMAT_VAR_PDIR}/SW837_sample_list.tsv

mkdir -p ${SOMAT_VAR_PDIR}/SW837_plots

cd ${SOMAT_VAR_PDIR}/SW837_plots

# Generate the VAF plots
Rscript ${SCRIPTDIR}/plot_vaf_vs_depth_from_maf.R --file ${PROJECTDIR}/analysis/somatic_variant_plots/SW837_keep_matched_caveman_pindel_all_indelvaf_filt.maf --genefile ${PROJECTDIR}/analysis/somatic_variant_plots/SW837_keep_matched_caveman_pindel_all_indelvaf_filt_gene_list.tsv --samplefile ${SOMAT_VAR_PDIR}/SW837_sample_list.tsv --colname "Hugo_Symbol" --suffix "genes" --width 10 --height 12 -ncol 5

```
**OUTPUTS**:
- **AF_vs_depth_recurrent_genes.pdf**: VAF vs depth plot showing all the mutations across the _CDS2_ targeted  SW837 cell lines - Extended Data Figure 12


## Results

The list of files generated present in the results [`somatic_variant_plots`](../analysis/somatic_variant_plots) folder are:


```bash
├── CDS2_nontreated_vs_treated_plots
│   ├── Mut_comparison_fisher_tests_min_nmut_samples_nontreated_vs_treated_res_tab.tsv
│   ├── OMM2.5_CDS2_Tumours_non_treated_vs_treated_fisher_test_between_both_Cohorts_0.1pval_forestplot.pdf
│   ├── OMM2.5_CDS2_Tumours_non_treated_vs_treated_fisher_test_between_both_Cohorts_0.1pval_lolliplots_OSGIN1.pdf
│   ├── OMM2.5_CDS2_Tumours_non_treated_vs_treated_fisher_test_between_both_Cohorts_0.1pval_lolliplots_SNX19.pdf
│   ├── OMM2.5_CDS2_Tumours_non_treated_vs_treated_fisher_test_between_both_Cohorts_lolliplots_CDS1.pdf
│   └── OMM2.5_CDS2_Tumours_non_treated_vs_treated_fisher_test_between_both_Cohorts_lolliplots_CDS2.pdf
├── OMM2.5_CDS2_Tumours_keep_matched_caveman_pindel_all_indelvaf_filt.maf
├── OMM2.5_plots
│   ├── OMM2.5_CDS2_Tumours_lines_oncoplot_allmuts.pdf
│   ├── OMM2.5_CDS2_Tumours_lines_oncoplot_cds_gRNA_exontar_and_exofftar_and_sig_mut_genes.pdf
│   ├── OMM2.5_CDS2_Tumours_lines_oncoplot_cds_gRNA_exontar_and_exofftar_compared_genes.pdf
│   ├── OMM2.5_CDS2_Tumours_lines_oncoplot_cds_gRNA_exontar_and_exofftar_genes_bplot_cgcv97.pdf
│   ├── OMM2.5_CDS2_Tumours_lines_oncoplot_cds_gRNA_exontar_and_exofftar_genes_bplot.pdf
│   ├── OMM2.5_CDS2_Tumours_lines_oncoplot_cds_gRNA_exontar_and_exofftar_genes_nobplot.pdf
│   └── OMM2.5_CDS2_Tumours_lines_oncoplot_cds_gRNA_PLC_Genes.pdf
├── SW837_keep_matched_caveman_pindel_all_indelvaf_filt_gene_list.tsv
├── SW837_keep_matched_caveman_pindel_all_indelvaf_filt.maf
├── SW837_plots
│   ├── AF_vs_depth_indel_samples.pdf
│   ├── AF_vs_depth_recurrent_genes.pdf
│   ├── AF_vs_depth_recurrent_sites.pdf
│   ├── AF_vs_depth_snv_samples.pdf
│   ├── mutation_types_barplot_proportion_samples.pdf
│   ├── mutation_types_barplot_samples.pdf
│   ├── SW837_lines_lolliplot_CDS2.pdf
│   ├── SW837_lines_oncoplot.pdf
│   ├── top_recurrently_mutated_genes.tsv
│   └── top_recurrently_mutated_sites.tsv
└── SW837_sample_list.tsv
```


