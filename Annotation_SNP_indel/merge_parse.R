args = commandArgs(trailingOnly=TRUE)
variants_file <- args[1]
variants_annot_file <- args[2]
sample_name <- args[3]

variants <- read.table(variants_file,header=T,sep="\t",quote="")
variants_annot <- read.csv(variants_annot_file,header=T,sep="\t",quote="")

table_comp <- merge(variants,variants_annot,by="ID",all.x=F,all.y=F)

final_table_name <- paste(sample_name,"variants_annot_all.tab", sep="_")
write.table(table_comp, file=final_table_name, sep="\t", na="-", row.names=F, quote=FALSE)
