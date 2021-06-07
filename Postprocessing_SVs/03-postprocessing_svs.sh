#JOINT Diploid Sample Analysis ~ MANTA.

#conda activate nf-core-sarek-2.6-fix-gatk

# 1. Configure Manta
configManta.py --bam ../02-sarek/Preprocessing/HS233/Recalibrated/HS233.recal.bam --bam ../02-sarek/Preprocessing/HS239/Recalibrated/HS239.recal.bam --referenceFasta ../../../REFERENCES/WholeGenomeFasta/Homo_sapiens_assembly38.fasta --runDir ./00-joint_call_manta

# 2. Execute Manta workflow and it will output the results based on the configuration first step.
echo "python ./00-joint_call_manta/runWorkflow.py" > _00_joint_call_manta.sh

# 3. Converting BND reported by Manta into INVersions. 
gzip -d ./00-joint_call_manta/results/variants/diploidSV.vcf.gz
echo '/processing_Data/bioinformatics/pipelines/manta-1.6.0.release_src/bin/libexec/convertInversion.py  /processing_Data/bioinformatics/pipelines/miniconda3/envs/nf-core-sarek-2.6-fix-gatk/bin/samtools ../../../REFERENCES/WholeGenomeFasta/Homo_sapiens_assembly38.fasta ./00-joint_call_manta/results/variants/diploidSV.vcf  > ./01-inversions_manta/Manta_inv_converted.vcf' > _01_convertInversion.sh

# 4. Separate vcf into individuals samples
echo "for file in ./01-inversions_manta/Manta_inv_converted.vcf; do for sample in \`bcftools query -l \$file\` ; do bcftools view -c1 -Ov -s \$sample -o \$sample.vcf \$file ; done ; done" > _02_split_samples.sh
