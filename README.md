<img src="assets/BU_ISCIII_logo.png" alt="logo" width="300" align="right"/>

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.10.0-brightgreen.svg)](https://www.nextflow.io/)

# Introduction
An open-source pipeline to detect germline variants from whole genome and whole exome sequencing.
This workflow was design using [`Sarek`](https://github.com/nf-core/sarek), a **nf-core** community pipeline, created to aim the same objective as ours, and the recommendations of the Best Practices guide from GATK and literature related to this topic.

Sarek pipeline is built using [`Nextflow`](https://www.nextflow.io/), a workflow tool to run tasks across multiple compute infrastructures in a very portable, scalable and reproducible manner. We used this workflow for the first steps, regarding: preprocessing and variant calling (just HaplotypeCaller in GVCF mode).
Further steps were carried out by our own, following GATK Best Practices for short variants discovery.

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
