<img src="/BU_ISCIII_logo.png" alt="logo" width="300" align="right"/>

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A519.10.0-brightgreen.svg)](https://www.nextflow.io/)

# Introduction
An open-source pipeline to detect germline variants from whole genome and whole exome sequencing.

This workflow was design using [`Sarek`](https://github.com/nf-core/sarek), a **nf-core** community pipeline, created to aim the same objective as ours, and the recommendations of the Best Practices guide from [`GATK`](https://gatk.broadinstitute.org/hc/en-us) and literature related to this topic.

Sarek pipeline is built using [`Nextflow`](https://www.nextflow.io/), a workflow tool to run tasks across multiple compute infrastructures in a very portable, scalable and reproducible manner. We used this workflow for the first steps, regarding: preprocessing and variant calling (just HaplotypeCaller in GVCF mode).
Further steps were carried out by our own, following GATK Best Practices for short variants discovery.

This pipeline consists of 3 main steps (preprocessing being tuned up):
* **01-sarek**, which involves quality test, mapping and variant calling (just HaplotypeCaller in GVCF mode).
* **02-postprocessing**, which involves joint genotyping and variant quality filtering for SNPs and in/dels and joint call and conversion of inversions, for SVs.
* **03-annot**, where SNPs and in/dels and SVs are annotated.

# Pipeline summary
1. Sequencing quality control ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Map read to reference ([`BWA-mem`])
3. Mark duplicates ([`GATK MarkDuplicates`])
4. Base (Quality Score) Recalibration ([`GATK BaseRecalibrator`], [`GATK ApplyBQSR`])
5. Joint Variant Calling and post-processing
    1. SNPs/InDels ([`HaplotypeCaller --GVCF`])
        1. Joint genotyping ([`GATK CombineGVCFs`], [`GATK GenotypeGVCFs`]) 
        2. Quality hard-filtering ([`GATK SelectVariants`](https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering), [`GATK VariantFiltration`](https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering), [`GATK MergeVcfs`](https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering))
    2. SVs ([`Manta`](https://github.com/Illumina/manta/blob/master/docs/userGuide/README.md))
        1.  Inversion conversion ([`Manta - convertInversions.py`](https://github.com/Illumina/manta/blob/master/docs/userGuide/README.md#inversions)
6. Annotation
    1. SNPs/InDels ([`VEP --plugin dbNSFP4.1a`](https://www.ensembl.org/info/docs/tools/vep/index.html))
    2. SVs ([`AnnotSV`](https://lbgi.fr/AnnotSV/))
7. Overall pipeline run summaries ([`MultiQC`](https://multiqc.info/))

# Pipeline execution guidelines
BU-ISCIII currently implements its pipelines using nextflow and singularity; however some of our pipelines are not migrated yet or some specific analysis don't have the enough entity to make a automatic pipeline (for example while we are testing something). In these cases we use a standardized organization for the analysis in order to be able to understand an replicate the execution and follow the results.

* **samples_id.txt file**: contains the sample identifiers for all the samples being analyzed in the project. Usually sample names are included in the fastq files with this format: `{sample_name}_S45_R1_001.fastq.gz`. One easy way for setting samples_id.txt file is locating yourself in ANALYSIS folder and executing:

    `find ../RAW_NC -name "*.fastq.gz" | cut -d "/" -f 3 | cut -d "_" -f 1 | sort -u > samples_id.txt`.
    
* **00-reads**: we use this folder as starting point for our analysis. It is useful for simplify posterior steps to rename fastq files to just `{sample_name}_R{1,2}.fastq.gz`. As we don't want to rename the original files we use this folder to make symbolic links with the new names. We can use this command for this, locate yourself inside 00-reads for easier use of ln command:

    `cat ../samples_id.txt | xargs -I % echo "ln -s ../../RAW/%_*R1*.fastq.gz" %_R1.fastq.gz" | bash`\
    `cat ../samples_id.txt | xargs -I % echo "ln -s ../../RAW/%_*R2*.fastq.gz" %_R2.fastq.gz" | bash`

* **01-preprocessing, 02-sarek, 03-postprocessing_variants, 04-postprocessing_svs, 05-annot_variants, 06-annot_svs**: next steps of the pipeline will follow having each one its own folder. Each folder must have a **lablog** file and the **scripts** needed for regenerate the results stored in that folder. Also a logs folder is created where all generated logs are stored.

* **lablog**: this file is located in each folder of the analysis, it indicates the commands that must be run for replicate the results found in its folder. Most of the times lablog does not contain the commands for running the analysis, instead it runs commands for generating the necessary scripts for performing the analysis. For example necessary orders for iterating among all the samples in the analysis will be set in this file. And of course it will describe all the steps needed for the analysis. Let's see an example, iamgine a lablog file contained in 01-fastQC folder will contain the fastqc command for all the samples in the project.

Imagine **samples_id.txt** looks like this:

`sample1`\
`sample2`\
`sample3`

lablog file inside 01-fastQC folder would look like this:

`## This line iterates using samples_id.txt file and executes the xargs part one time per line in the samples_id.txt file. In each iteration "%" works as a variable that contains each value in each line (sample1, sample2,...)`\
`cat ../samples_id.txt | xargs -I % echo "mkdir %;fastqc -o % ../00-reads/%_R1.fastq.gz ../00-reads/%_R2.fastq.gz" >> _00_rawfastqc.sh`\
`## If you feel more confortable with traditional loops the equivalent way using while would be:`\
`cat ../samples_id.txt | while read in; do echo "mkdir "$in";fastqc -o "$in" ../00-reads/"$in"_R1.fastq.gz ../00-reads/"$in"_R2.fastq.gz"; done >> _00_rawfastqc.sh`

When we execute lablog file: bash lablog we generate two scripts that will look like this: **_00_rawfastqc.sh**

`mkdir sample1;fastqc -o sample1 ../00-reads/sample1_R1.fastq.gz ../00-reads/sample1_R2.fastq.gz`\
`mkdir sample2;fastqc -o sample2 ../00-reads/sample2_R1.fastq.gz ../00-reads/sample2_R2.fastq.gz`\
`mkdir sample3;fastqc -o sample3 ../00-reads/sample3_R1.fastq.gz ../00-reads/sample3_R2.fastq.gz`


