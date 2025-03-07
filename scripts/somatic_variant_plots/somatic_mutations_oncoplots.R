#!/usr/bin/env Rscript

progname<-"somatic_mutations_found_oncoplots.R"


#List of packages to use
library(here)
library(ggplot2)
library(maftools)
library(RColorBrewer)



projdir<-here()
setwd(projdir)
studyid<-7688
projectid<-3365

var_dir<-file.path(projdir, "analysis", "variants_combined", "version2","all" )
results_dir<-file.path(projdir, "analysis", "somatic_variant_plots")
swplot_dir<-file.path(results_dir, "SW837_plots")
omplot_dir<-file.path(results_dir, "OMM2.5_plots")
comp_plot_dir<-file.path(results_dir, "CDS2_nontreated_vs_treated_plots")
metadata_fname<-file.path(projdir, "metadata", paste(studyid, projectid, "sample_exp_metadata.tsv", sep = "_"))


# Create directories
dir.create(results_dir, recursive = TRUE)
dir.create(swplot_dir, recursive = TRUE)
dir.create(omplot_dir, recursive = TRUE)
dir.create(comp_plot_dir, recursive = TRUE)


#Get the total Frequency of reads mutated by variant type observe
get_mutated_freq_per_var_type_psample<-function(tmaf){
  if(is.data.frame(tmaf)!="MAF"){
    stop("tmaf is not a maftools object ")
  }
  tres<-tmaf@data |>
    dplyr::select(Tumor_Sample_Barcode, Hugo_Symbol,Variant_Classification, Variant_Type, VAF_tum) |>
    dplyr::group_by(Tumor_Sample_Barcode,Hugo_Symbol, Variant_Type) |> dplyr::summarise(Total_vaf_mut=sum(VAF_tum))
}


#########  Read the required input files - for SW figures
mdata<-read.csv(metadata_fname, header = TRUE, stringsAsFactors = FALSE, sep="\t")
sw_maf_mdata<-NULL
sw_maf_mdata<-mdata[grepl("SW", mdata$Sample_name), c("Sample_name", "CDS2_targeting", "CDS2_status", "Group")]
colnames(sw_maf_mdata)[1]<-"Tumor_Sample_Barcode"
# Read the maf files with all the mutations
temp_maf<-read.csv(file.path(var_dir, "keep_caveman_pindel_all.maf"), header = TRUE, stringsAsFactors = FALSE, sep="\t")
#Filter Because at least 80% of baits had +51X or more we will filter 
temp_maf<-temp_maf[temp_maf$VAF_tum>=0.04,]
temp_maf<-temp_maf[grepl("SW", temp_maf$Tumor_Sample_Barcode),]
write.table(temp_maf,file = file.path(results_dir, "SW837_keep_matched_caveman_pindel_all_indelvaf_filt.maf"), quote=FALSE, row.names=FALSE, col.names = TRUE, sep="\t" )

sw_all_maf<-read.maf(maf =temp_maf , clinicalData = sw_maf_mdata)

## PARSE the MAF file and get the sum of the VAFs of all the mutations by Variant type per gene per sample
sw_all_tot_vaf_persamp_perg<- sw_all_maf@data |>
  dplyr::select(Tumor_Sample_Barcode, Hugo_Symbol,Variant_Classification, Variant_Type, VAF_tum) |>
  dplyr::group_by(Tumor_Sample_Barcode,Hugo_Symbol, Variant_Type) |> dplyr::summarise(Total_vaf_mut=sum(VAF_tum))
#Set group colors 
group_col<-colorRampPalette(brewer.pal(n=12, "Paired")[1:2])(3)
group_col<- group_col[1:2]
names(group_col)<-unique(sw_all_maf@clinical.data$Group)

# SW837 lines
pdf(file.path(swplot_dir, "SW837_lines_oncoplot.pdf"), height = 9, width = 11)
oncoplot(sw_all_maf, 
         showTitle = FALSE,
         # All mutated genes plus CDS1 and exonic off targets of CDS2gRNA
         genes= c( unique(sw_all_maf@gene.summary$Hugo_Symbol), "CDS1","FRRS1L", "SMG1"), 
         showTumorSampleBarcodes=TRUE,
         drawRowBar=FALSE,
         drawColBar=FALSE, # that is the total number of mutations
         SampleNamefontSize = 1,
         barcode_mar= 15,
         additionalFeatureCex = 1,
         annotationFontSize = 1.5,
         legendFontSize = 1.5, 
         annotationColor = list(Group=group_col),
         clinicalFeatures = "Group",
         gene_mar=7,
         showPct=FALSE,
         draw_titv = F)
dev.off()

#Format aminoacid_changes to be used with lolliplot2 
sw_all_maf@data$amino_acid_change<- gsub("p.", "", sw_all_maf@data$HGVSp_Short, fixed = T)
#Replace "-" for NA 
sw_all_maf@data$amino_acid_change[sw_all_maf@data$amino_acid_change=="-"]<- NA
#Subset MAFs for each group
sw_cdtar_maf<-subsetMaf(sw_all_maf, tsb=c("SW837_C9_CDS2_gRNA_R1_DNA","SW837_C9_CDS2_gRNA_R2_DNA","SW837_C9_CDS2_gRNA_R3_DNA" ) )
sw_sftar_maf<-subsetMaf(sw_all_maf, tsb=c("SW837_C9_safe_gRNA_R1_DNA","SW837_C9_safe_gRNA_R2_DNA","SW837_C9_safe_gRNA_R3_DNA" ) )

#Lollipop plot for CDS2
pdf(file.path(swplot_dir, paste0("SW837_lines_lolliplot_CDS2.pdf")),
    height=5, width =9)
lollipopPlot2(m1 = sw_cdtar_maf, m2 = sw_sftar_maf, gene="CDS2", AACol1 = "amino_acid_change", AACol2 = "amino_acid_change",
              m1_name = 'SW837_C9 CDS2_gRNA', m2_name = 'SW837_C9 STG1_gRNA')
dev.off()

### Generate the gene list
write.table(unique(sw_all_maf@gene.summary$Hugo_Symbol), 
            file = file.path(results_dir, "SW837_keep_matched_caveman_pindel_all_indelvaf_filt_gene_list.tsv"),
            quote=FALSE, row.names=FALSE, col.names = FALSE, sep="\t" )



#########  Read the required input files - for OM Tumour figures
mdata<-read.csv(metadata_fname, header = TRUE, stringsAsFactors = FALSE, sep="\t")
om_maf_mdata<-NULL
om_maf_mdata<-mdata[grepl("CDS2_Tumour", mdata$Sample_name), c("Sample_name", "CDS2_targeting", "CDS2_status", "Group")]
om_maf_mdata<-om_maf_mdata[ om_maf_mdata$Group %in% c("CDS2_Tumour_Non-treated", "CDS2_Tumour_Treated"),]
colnames(om_maf_mdata)[1]<-"Tumor_Sample_Barcode"
# Read the maf files with all the mutations
temp_maf<-read.csv(file.path(var_dir, "keep_caveman_pindel_all.maf"), header = TRUE, stringsAsFactors = FALSE, sep="\t")
#Filter Because at least 80% of baits had +51X or more we will filter 
temp_maf<-temp_maf[temp_maf$VAF_tum>=0.04,]
temp_maf<-temp_maf[grepl("CDS2_Tumour", temp_maf$Tumor_Sample_Barcode),]
write.table(temp_maf,file = file.path(results_dir, "OMM2.5_CDS2_Tumours_keep_matched_caveman_pindel_all_indelvaf_filt.maf"), quote=FALSE, row.names=FALSE, col.names = TRUE, sep="\t" )

om_all_maf<-read.maf(maf =temp_maf , clinicalData = om_maf_mdata)

## PARSE the MAF file and get the sum of the VAFs of all the mutations by Variant type per gene per sample
om_all_tot_vaf_persamp_perg<- om_all_maf@data |>
  dplyr::select(Tumor_Sample_Barcode, Hugo_Symbol,Variant_Classification, Variant_Type, VAF_tum) |>
  dplyr::group_by(Tumor_Sample_Barcode,Hugo_Symbol, Variant_Type) |> dplyr::summarise(Total_vaf_mut=sum(VAF_tum))
#Set group colors 
group_col<-colorRampPalette(brewer.pal(n=12, "Paired")[3:4])(3)
group_col<- group_col[1:2]
names(group_col)<-unique(om_all_maf@clinical.data$Group)

#  OM Tumours - all muts
pdf(file.path(omplot_dir, "OMM2.5_CDS2_Tumours_lines_oncoplot_allmuts.pdf"), height = 10, width = 10)
oncoplot(om_all_maf, 
         showTitle = FALSE,
         showTumorSampleBarcodes=TRUE,
         drawRowBar=FALSE,
         drawColBar=FALSE, # that is the total number of mutations
         SampleNamefontSize = 0.9,
         barcode_mar= 15,
         additionalFeatureCex = 1,
         annotationFontSize = 1,
         legendFontSize = 1, 
         annotationColor = list(Group=group_col),
         clinicalFeatures = "Group",
         gene_mar=7,
         showPct=FALSE,
         draw_titv = F)
dev.off()


#  OM Tumours - CDS2 plus exonic off-targets
pdf(file.path(omplot_dir, "OMM2.5_CDS2_Tumours_lines_oncoplot_cds_gRNA_exontar_and_exofftar_genes_nobplot.pdf"), height = 10, width = 10)
oncoplot(om_all_maf, 
         showTitle = FALSE,
         genes = c("CDS2", "FRRS1L", "SMG1"),
         showTumorSampleBarcodes=TRUE,
         drawRowBar=FALSE,
         drawColBar=FALSE, # that is the total number of mutations
         SampleNamefontSize = 0.9,
         barcode_mar= 15,
         additionalFeatureCex = 1,
         annotationFontSize = 1,
         legendFontSize = 1, 
         annotationColor = list(Group=group_col),
         clinicalFeatures = "Group",
         gene_mar=7,
         showPct=FALSE,
         draw_titv = F)
dev.off()

#  OM Tumours - CDS2 plus exonic off-targets
pdf(file.path(omplot_dir, "OMM2.5_CDS2_Tumours_lines_oncoplot_cds_gRNA_exontar_and_exofftar_genes_bplot.pdf"), height = 10, width = 10)
oncoplot(om_all_maf, 
         showTitle = FALSE,
         genes = c("CDS2", "FRRS1L", "SMG1"),
         showTumorSampleBarcodes=TRUE,
         drawRowBar=TRUE,
         drawColBar=TRUE, # that is the total number of mutations
         SampleNamefontSize = 0.9,
         barcode_mar= 15,
         additionalFeatureCex = 1,
         annotationFontSize = 1,
         legendFontSize = 1, 
         annotationColor = list(Group=group_col),
         clinicalFeatures = "Group",
         gene_mar=7,
         showPct=FALSE,
         draw_titv = F)
dev.off()

# List of Mutated genes in OMM2.5 that are in the list of Cancer Genes reported on 
# the list of 736 genes on cancer Gene Census of COSMIC v97
cancer_gene_census_v97<-NULL
cancer_gene_census_v97<- read.csv(file=file.path(projdir, "resources/COSMIC/OMM2.5_mutated_cancer_gene_census.v97.genes_names.tsv"), header = F, stringsAsFactors = F, sep="\t")
colnames(cancer_gene_census_v97)<-c("Gene.Symbol")
cancer_gene_census_v97<-as.data.frame(cancer_gene_census_v97)

#  OM Tumours - CDS2 plus exonic off-targets
pdf(file.path(omplot_dir, "OMM2.5_CDS2_Tumours_lines_oncoplot_cds_gRNA_exontar_and_exofftar_genes_bplot_cgcv97.pdf"), height = 18, width = 11)
oncoplot(om_all_maf, 
         showTitle = FALSE,
         genes = cancer_gene_census_v97$Gene.Symbol[cancer_gene_census_v97$Gene.Symbol %in% om_all_maf@data$Hugo_Symbol],
         showTumorSampleBarcodes=TRUE,
         drawRowBar=TRUE,
         drawColBar=TRUE, # that is the total number of mutations
         SampleNamefontSize = 0.9,
         barcode_mar= 14,
         additionalFeatureCex = 1,
         annotationFontSize = 1.5,
         legendFontSize = 1.5, 
         annotationColor = list(Group=group_col),
         sortByAnnotation = TRUE, 
         clinicalFeatures = "Group",
         gene_mar=7,
         showPct=FALSE,
         draw_titv = F)
dev.off()


#Compare the two Groups 

#Format aminoacid_changes to be used with lolliplot2 
om_all_maf@data$amino_acid_change<- gsub("p.", "", om_all_maf@data$HGVSp_Short, fixed = T)
#Replace "-" for NA 
om_all_maf@data$amino_acid_change[om_all_maf@data$amino_acid_change=="-"]<- NA
#Subset MAFs for each group
om_nontr_maf<-subsetMaf(om_all_maf, tsb=c("CDS2_Tumour_Non-treated_R1_DNA","CDS2_Tumour_R1_DNA", "CDS2_Tumour_Non-treated_R2_DNA", "CDS2_Tumour_R2_DNA" ) )
om_tr_maf<-subsetMaf(om_all_maf, tsb=c("CDS2_Tumour_Treated_R1_DNA", "CDS2_Tumour_Treated_R2_DNA", "CDS2_Tumour_Treated_R3_DNA", "CDS2_Tumour_Treated_R4_DNA" ) )

#Compare the samples
npt.vs.rt <- mafCompare(m1 = om_nontr_maf, m2 = om_tr_maf, m1Name = 'CDS2_Non_treated', m2Name = 'CDS2_Treated', minMut = 4)
pdf(file.path(comp_plot_dir, "OMM2.5_CDS2_Tumours_non_treated_vs_treated_fisher_test_between_both_Cohorts_0.1pval_forestplot.pdf"), height=12, width = 10)
forestPlot(mafCompareRes = npt.vs.rt, pVal=0.05, geneFontSize =1.5, lineWidth = 2)
dev.off()

# Lolliplot2 - OSGIN1
pdf(file.path(comp_plot_dir, paste0("OMM2.5_CDS2_Tumours_non_treated_vs_treated_fisher_test_between_both_Cohorts_0.1pval_lolliplots_", npt.vs.rt$result$Hugo_Symbol[1], ".pdf")),
              height=5, width =9)
lollipopPlot2(m1 = om_nontr_maf, m2 = om_tr_maf, gene=npt.vs.rt$result$Hugo_Symbol[1], AACol1 = "amino_acid_change", AACol2 = "amino_acid_change",
              m1_name = 'CDS2_Non_treated', m2_name = 'CDS2_Treated')
dev.off()
pdf(file.path(comp_plot_dir, paste0("OMM2.5_CDS2_Tumours_non_treated_vs_treated_fisher_test_between_both_Cohorts_0.1pval_lolliplots_", npt.vs.rt$result$Hugo_Symbol[2], ".pdf")),
              height=5, width =9)
lollipopPlot2(m1 = om_nontr_maf, m2 = om_tr_maf, gene=npt.vs.rt$result$Hugo_Symbol[2], AACol1 = "amino_acid_change", AACol2 = "amino_acid_change",
              m1_name = 'CDS2_Non_treated', m2_name = 'CDS2_Treated')
dev.off()

pdf(file.path(comp_plot_dir, paste0("OMM2.5_CDS2_Tumours_non_treated_vs_treated_fisher_test_between_both_Cohorts_lolliplots_CDS2.pdf")),
    height=5, width =9)
lollipopPlot2(m1 = om_nontr_maf, m2 = om_tr_maf, gene="CDS2", AACol1 = "amino_acid_change", AACol2 = "amino_acid_change",
              m1_name = 'CDS2_Non_treated', m2_name = 'CDS2_Treated')
dev.off()

pdf(file.path(comp_plot_dir, paste0("OMM2.5_CDS2_Tumours_non_treated_vs_treated_fisher_test_between_both_Cohorts_lolliplots_CDS1.pdf")),
    height=5, width =9)
lollipopPlot2(m1 = om_nontr_maf, m2 = om_tr_maf, gene="CDS1", AACol1 = "amino_acid_change", AACol2 = "amino_acid_change",
              m1_name = 'CDS2_Non_treated', m2_name = 'CDS2_Treated')
dev.off()


write.table(npt.vs.rt$results, file = file.path(comp_plot_dir,"Mut_comparison_fisher_tests_min_nmut_samples_nontreated_vs_treated_res_tab.tsv" ), col.names = TRUE, row.names = FALSE,quote = FALSE, sep="\t" )

#  OM Tumours - CDS2 plus exonic off-targets
pdf(file.path(omplot_dir, "OMM2.5_CDS2_Tumours_lines_oncoplot_cds_gRNA_exontar_and_exofftar_compared_genes.pdf"), height = 8, width = 10)
oncoplot(om_all_maf, 
         showTitle = FALSE,
         genes= npt.vs.rt$result$Hugo_Symbol[1:2],
         showTumorSampleBarcodes=TRUE,
         drawRowBar=FALSE,
         removeNonMutated = FALSE, 
         drawColBar=FALSE, # that is the total number of mutations
         SampleNamefontSize = 1.05,
         barcode_mar= 17,
         additionalFeatureCex = 1,
         annotationFontSize = 1.5,
         legendFontSize = 1.5, 
         annotationColor = list(Group=group_col),
         sortByAnnotation = TRUE, 
         clinicalFeatures = "Group",
         gene_mar=7,
         showPct=FALSE,
         draw_titv = F)
dev.off()


#  OM Tumours - CDS2 plus exonic off-targets and genes that had pval below 0.05
pdf(file.path(omplot_dir, "OMM2.5_CDS2_Tumours_lines_oncoplot_cds_gRNA_exontar_and_exofftar_and_sig_mut_genes.pdf"), height = 8, width = 10)
oncoplot(om_all_maf, 
         showTitle = FALSE,
         genes= c("CDS2", "CDS1", "FRRS1L", "SMG1", npt.vs.rt$result$Hugo_Symbol[1:2] ),
         keepGeneOrder = TRUE,
         showTumorSampleBarcodes=TRUE,
         drawRowBar=FALSE,
         removeNonMutated = FALSE, 
         drawColBar=FALSE, # that is the total number of mutations
         SampleNamefontSize = 1.05,
         barcode_mar= 17,
         additionalFeatureCex = 1,
         annotationFontSize = 1.5,
         legendFontSize = 1.5, 
         annotationColor = list(Group=group_col),
         sortByAnnotation = TRUE, 
         clinicalFeatures = "Group",
         gene_mar=7,
         showPct=FALSE,
         draw_titv = F)
dev.off()

#### PLC GENES ONco plot
plc_genes<-c("PLCG1", "PLCGE1", "PLCD1","PLCB1", "PLCB2","PLCZ1")

#  OM Tumours - CDS2 plus exonic off-targets and genes that had pval below 0.05
pdf(file.path(omplot_dir, "OMM2.5_CDS2_Tumours_lines_oncoplot_cds_gRNA_PLC_Genes.pdf"), height = 8, width = 10)
oncoplot(om_all_maf, 
         showTitle = FALSE,
         genes= plc_genes,
         keepGeneOrder = TRUE,
         showTumorSampleBarcodes=TRUE,
         drawRowBar=FALSE,
         removeNonMutated = FALSE, 
         drawColBar=FALSE, # that is the total number of mutations
         SampleNamefontSize = 1.05,
         barcode_mar= 17,
         additionalFeatureCex = 1,
         annotationFontSize = 1.5,
         legendFontSize = 1.5, 
         annotationColor = list(Group=group_col),
         sortByAnnotation = TRUE, 
         clinicalFeatures = "Group",
         gene_mar=7,
         showPct=FALSE,
         draw_titv = F)
dev.off()




