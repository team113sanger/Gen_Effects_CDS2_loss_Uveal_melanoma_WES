# Instructions to generate the NOD_V1_PDX reference files

To be able to filter mouse reads on the xenografted samples, we generated a custom FASTA file referred as **NOD_V1_PDX**. This file was created using the NOD_ShiLtJ_V1 reference and GRCm39 Y Chromosome sequence. In addition, we collated a custom GTF file containing the gene annotations for the NOD_ShiLtJ_V1 reference and GRCm39 Y Chromosome, from ENSEMBL release 107. 

In this folder, you will find:

- The RMarkdown file used to generate the reference fasta file - see [PDX_from_LatAm_mouse_reference_generation_NOD_ShiLtJ_v1.Rmd](./PDX_from_LatAm_mouse_reference_generation_NOD_ShiLtJ_v1.Rmd)
- An HTML file with the output of the RMarkdown file run as well as the information and summary statistics of the reference metrics - download the file and open with a web browser. [PDX_from_LatAm_mouse_reference_generation_NOD_ShiLtJ_v1.html](./PDX_from_LatAm_mouse_reference_generation_NOD_ShiLtJ_v1.html)


## To generate the `NOD_V1_PDX` reference Fasta file

To regenerate the reference fasta file execute the RMarkdown file contained within the `reference/NOD_ShiLtJ_V1_PDX_ref` folder (it requires R and the path of the project directory where the repository was cloned into as the `projectdir` variable). 

### First ensure the required R packages are installed

```R
projectdir<-"/lustre/7688_3365_Gen_Effects_CDS2_loss_Uveal_melanoma_WES"
nod_ref_dir<- file.path(projectdir,"reference/NOD_ShiLtJ_V1_PDX_ref")

setwd(nod_ref_dir)
install.packages("renv")
library(renv)
renv::restore(lockfile=file.path(nod_ref_dir, "renv.lock")) # To rebuild an environment from the renv.lockfile

```

### Run the RMarkdown file to download and generate the FASTA file and the required annotation

```R
# Run R :
R
#Directory where git clone was performed
projectdir<-"/lustre/7688_3365_Gen_Effects_CDS2_loss_Uveal_melanoma_WES"

library(rmarkdown)
rmarkdown::render(input=file.path(
    projectdir, "reference/NOD_ShiLtJ_V1_PDX_ref/PDX_from_LatAm_mouse_reference_generation_NOD_ShiLtJ_v1.Rmd"
    ), params=list(generate=TRUE,refdir=file.path(projectdir,"reference")  )
)
```

### Generate the index file for `bwa`

This should be run in a machine with 4 threads and 16GB of RAM and requires to have `bwa` installed.
- `bwa-mem` version `0.7.17` [**here**](https://github.com/lh3/bwa)

```bash
PROJECTDIR=/lustre/7688_3365_Gen_Effects_CDS2_loss_Uveal_melanoma_WES
STUDY=7688

cd ${PROJECTDIR}/reference/NOD_ShiLtJ_V1_PDX_ref

bwa index ${PROJECTDIR:?unset}/reference/NOD_ShiLtJ_V1_PDX_ref/genome.fa 

```
