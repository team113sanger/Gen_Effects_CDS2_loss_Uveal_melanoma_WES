echo -e "#CHROM\tPOS\tREF\tALT\tdbSNP" > dbSNP155_common.tsv
zcat /lustre/scratch126/casm/team113da/projects/5534_Landscape_sebaceous_tumours_GRCh38_Remap/RESOURCES/dbSNP155.GRCh38.GCF_000001405.39.mod.vcf.gz |grep ";COMMON"|cut -f 1,2,3,4,5| awk '{print $0"\tCOMMON"}' >> dbSNP155_common.tsv
