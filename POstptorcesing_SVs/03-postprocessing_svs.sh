#JOINT Diploid Sample Analysis ~ MANTA.

#conda activate nf-core-sarek-2.6-fix-gatk

# 1. Configure Manta
configManta.py --bam <path/to/01-sarek/Preprocessing/${sample1}/Recalibrated/${sample1}.recal.bam> --bam <path/to/01-sarek/Preprocessing/${sample2}/Recalibrated/${sample2}.recal.bam> --referenceFasta <path/to/Homo_sapiens_assembly38.fasta> --runDir ./00-joint_call_manta

# 2. Execute Manta workflow and it will output the results based on the configuration first step.
python ./00-joint_call_manta/runWorkflow.py

# 3. Converting BND reported by Manta into INVersions. 
gzip -d ./00-joint_call_manta/results/variants/diploidSV.vcf.gz
/processing_Data/bioinformatics/pipelines/manta-1.6.0.release_src/bin/libexec/convertInversion.py  /processing_Data/bioinformatics/pipelines/miniconda3/envs/nf-core-sarek-2.6-fix-gatk/bin/samtools <path/to/WholeGenomeFasta/Homo_sapiens_assembly38.fasta> ./00-joint_call_manta/results/variants/diploidSV.vcf > ./01-inversions_manta/Manta_inv_converted.vcf

# 4. Separate vcf into individuals samples
for file in ./01-inversions_manta/Manta_inv_converted.vcf; do for sample in \`bcftools query -l \$file\` ; do bcftools view -c1 -Ov -s \$sample -o \$sample.vcf \$file ; done ; done
