
# ANNOTATION OF SVs in Manta´s INVERSIONS converted vcfs using AnnotSV.

#conda activate svs-analysis
export ANNOTSV=/processing_Data/bioinformatics/pipelines/AnnotSV

cat ../../samples_id.txt | xargs -I % echo '$ANNOTSV/bin/AnnotSV -SVinputFile ../04-postprocessing_svs/%_svs.vcf -outputFile ./AnnotSV_full_%.tsv -genomeBuild GRCh38  -typeOfAnnotation full' > _01_AnnotSV_annot_full.sh

cat ../../samples_id.txt | xargs -I % echo '$ANNOTSV/bin/AnnotSV -SVinputFile ../04-postprocessing_svs/%_svs.vcf -outputFile ./AnnotSV_split_%.tsv -genomeBuild GRCh38 -typeOfAnnotation split' > _02_AnnotSV_annot_split.sh
