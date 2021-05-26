<img src="assets/BU_ISCIII_logo.png" alt="logo" width="300" align="right"/>

# Introduction
An open-source pipeline to detect germline variants from whole genome and whole exome sequencing.

This pipeline 

# Pipeline summary
1. Sequencing quality control ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Map read to reference ([`BWA-mem`])
3. Mark duplicates ([`GATK MarkDuplicates`])
4. Base (Quality Score) Recalibration ([`GATK BaseRecalibrator`], [`GATK ApplyBQSR`])
5. Joint Variant Calling and post-processing
  5.1. SNPs/InDels ([`HaplotypeCaller --GVCF`])
    5.1.1.  Joint genotyping ([`GATK CombineGVCFs`], [`GATK GenotypeGVCFs`]) 
    5.1.2.  Quality hard-filtering ([`GATK SelectVariants`], [`GATK VariantFiltration`], [`GATK MergeVcfs`] (https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering))
  5.2. SVs ([`Manta`](https://github.com/Illumina/manta/blob/master/docs/userGuide/README.md))
    5.2.1.  Inversion conversion ([`Manta - convertInversions.py`](https://github.com/Illumina/manta/blob/master/docs/userGuide/README.md#inversions))
6. Annotation
  6.1. SNPs/InDels ([`VEP --plugin dbNSFP4.1a`](https://www.ensembl.org/info/docs/tools/vep/index.html))
  6.2. SVs ([`AnnotSV`](https://lbgi.fr/AnnotSV/))
7. Overall pipeline run summaries ([`MultiQC`](https://multiqc.info/))
