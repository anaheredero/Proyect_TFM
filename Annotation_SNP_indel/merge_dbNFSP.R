args = commandArgs(trailingOnly=TRUE)

gene_vep_file <- args[1]
gene_dbnsfp_file <- args[2]
sample_name <- args[3]

gene_vep <- read.table(gene_vep_file,header=T,sep="\t", quote="")
gene_dbnsfp <- read.table(gene_dbnsfp_file,header=T,sep="\t", quote="")

table_comp <- merge(gene_vep,gene_dbnsfp,by="Gene",all.x=T,all.y=F)

final_table_name <- paste(sample_name,"vep_dbNSFP.txt", sep="_")
write.table(table_comp, file=final_table_name, sep="\t", na="-", row.names=F, quote=FALSE)
