#!/usr/bin/env nextflowv20.07.1

/*
================================================================================
                                  ---
================================================================================
*/

nextflow run ../../sarek/main.nf -bg --input '*.tsv' --tools HaplotypeCaller --no_gatk_spark --genome GRCh38 -profile hpc_isciii -resume

