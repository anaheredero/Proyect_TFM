# Running Sarek in WGS and WES samples.

## GENOMES

ln -s 00-reads/ .
ln -s samples_id.txt .
cat samples_id.txt | xargs -I % echo -e "%\tXY\t0\t%\t1\t00-reads/%_R1.fastq.gz\t00-reads/%_R2.fastq.gz" > samples.tsv

#conda activate nextflowv20.07.1
echo "nextflow run /processing_Data/bioinformatics/pipelines/sarek/main.nf -bg -profile hpc_isciii --input 'samples.tsv' --skip_qc bamQC --tools 'HaplotypeCaller,Manta,TIDDIT' --no_gatk_spark --generate_gvcf --genome GRCh38 --save_reference --outdir 20210216_ANALYSIS01_GENOME/02-sarek -resume" > _00_sarek_genome.sh

## EXOMES
ln -s 00-reads/ .
ln -s samples_id_exomes.txt .

echo "nextflow run /processing_Data/bioinformatics/pipelines/sarek/main.nf -bg -profile hpc_isciii --input 'samples_exomes.tsv' --target_bed ../RAW/EXOMES/BED/S03723424_Padded_modif.bed --skip_qc bamQC --tools 'HaplotypeCaller' --no_gatk_spark --genome GRCh38 --save_reference --generate_gvcf --outdir 20210216_ANALYSIS02_EXOME/02-sarek -resume" > _01_sarek_exome.sh
