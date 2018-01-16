# A Simple GATK Variant Call Workflow
The [GATK](https://software.broadinstitute.org/gatk/) (Genome Analysis Toolkit) has been developed in Broad Institute, identifying SNPs and indels in germline DNA and RNAseq data. Its scope is now expanding to include somatic variant calling tools and to tackle copy number (CNV) and structure variation (SV).  Refer to three original manuscripts published in early 2010's that cover both computational foundations underlying the GATK and "best practice" to tackle many practical issues.

## 1. Preparing GATK run
1. Define GATK prefix command
    1. Q27: Define a prefix command as we did for Picard 
    1. Q28: Print help command output to see which commands are available
 
## 2. Recalibrate Alignment ([BQSR](https://docs.google.com/file/d/0B2dK2q40HDWeLTFzNndsNDBuVms/preview))
1. Objectives
    1. The quality of base calls produced by the machines are subject to various sources of systematic technical error.
    1. It can lead to over- or under-estimated base quality scores in the data, affecting to all downstream analyses.
1. Method
    1. The program builds a model of covariation based on the data and a set of known variants
    1. It adjusts the base quality scores in the data based on the model.

1. Output    
    1. The method refines the base quality to achieve more accurate base qualities.
    1. The result improves the accuracy of our variant calls.
1. Reference
    1. [Best Practice](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_bqsr_BaseRecalibrator.php)
    1. [Algorithm](https://docs.google.com/file/d/0B2dK2q40HDWeLTFzNndsNDBuVms/preview)

## 3. [Haplotype Caller](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php)
1. Objectives
    1. This is the most essential algorithm of the pipeline to achieve the most accurate sensitive variant calls.
1. Method
    1. Given BAM file, the program builds a De Bruijn-like graph to reassemble genomic regions of interest, and identifies what are the possible haplotypes present in the data.
    1. The program then realigns each haplotype against the reference haplotype using the Smith-Waterman algorithm in order to identify potentially variant sites.
1. Example command line
    ```bash
    $gatk \
        -T HaplotypeCaller \
        -R ${REF_FILE} \
        -I ${ALN}.so.mdup.br.bam \
        --genotyping_mode DISCOVERY \
        --dbsnp $dbsnp_vcf \
        -stand_call_conf 30 \
        -nct ${NCPU} \
        -o ${OUTD}/gatkhc.vcf
    ```
## 4. [Variant Recalibrator](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_variantrecalibration_VariantRecalibrator.php)
1. Objective
    1. To build a recalibration model to score variant quality for filtering purposes
1. Method
    1. It creates a Gaussian mixture model by looking at the distribution of annotation values (QD, MQ, ReadPosRankSum, and etc.) over a high quality subset of the input call set,
    1. It scores all input variants according to the model,
    1. It filters variants based on score cutoffs identified in the first step.
1. Output
    1. FILTER field in the final VCF file is tagged by some useful classification information
1. Note: Since the FASTQ file we prepared is reduced for the purpose of demo, we will apply a hard filter instead of adaptive filter using training dataset and guassian mixture model.

## 5. Batch Script from BWA mem to GATK HC caller
1. `tmux`
1. `binf`
1. View `bwa_gatk_lite.sh` and see if you can understand the script
    ```bash
    3 min
    ```
1. `./bwa_gatk_lite.sh chr7a chr7`
1. `ctrl+b` and `d` # to detach the tmux window
