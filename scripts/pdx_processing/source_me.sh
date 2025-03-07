#!/bin/bash

#Clean the environment
module purge 


# R and R libraries 
_R_BIN_PATH="/software/team113/dermatlas/R/R-4.2.2/bin"
export PATH="${_R_BIN_PATH:?empty-path-variable}:${PATH}"
export R_LIBS="/lustre/scratch124/casm/team113/projects/7688_7687_Gen_Effects_CDS2_loss_Uveal_melanoma/7688_3365_Gen_Effects_CDS2_loss_Uveal_melanoma_WES/renv/library/R-4.2/x86_64-pc-linux-gnu"

#load iRODS module 
module load IRODS/1.0

# Required samtools versions 
module load samtools-1.14/python-3.12.0
module load bwa-0.7.17
