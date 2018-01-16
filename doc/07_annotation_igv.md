# 1. VCF file
1. What does it stand for?
    1. VCF stands for Variant Call Format
1. Variant Call Format (VCF) is a text file format (most likely stored in a compressed manner). It contains meta-information lines, a header line, and then data lines each containing information about a position in the genome.
1. Visit [here](http://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40/) to understand the file format
1. Visit [here](https://gatkforums.broadinstitute.org/gatk/discussion/1268/what-is-a-vcf-and-how-should-i-interpret-it) to interpret a VCF file

# Our first VCF file
## 2. Understand VCF file
1. View our GATK output VCF file
    1. `tmux attach`
    1. `less -S ${OUTD}/gatkhc_filt.vcf`
    1. Q29: Find a line something like 
        ```bash
        #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Sample..."
        ```
    1. CHROM, POS, ID, REF, ALT
        1. Q30: Can you find insertion type variant?
        1. Q31: Can you find deletion type variant?
    1. QUAL: The Phred-scaled probability that a REF/ALT polymorphism exists at this site given sequencing data.
    1. FILTER: Classificaiton label about the variant call. In general, it annotates a confidence(PASS or LowPass) about the variant.
    1. INFO: Annotation field
        1. What is FS?
        1. What is QD?
    1. Q32: FORMAT and sampleID1
        1. Find GT in FORMAT field.
        1. Find 1/1 in sampleID field
    1. Q33: Search for 21174798
        1. Check its GT
        1. 0/1 vs. 1/1
## Analyze VCF
### 3. Annotation
1. Objectives: The raw VCF file only contains the location of variant and how much sure it is correct variant. What do we do with this?   
1. Annotation:
    1. Bring all information/discovered knowledge knowing about the variants
    1. [ANNOVAR](http://annovar.openbioinformatics.org/en/latest/)
    1. open the annotated VCF table
        1. less -S ${ALN}.annovar.chr7.merged.hg19_multianno.txt
        1. Search for xxxxxx
        1. Look at the annotation in the line
1. View the corresponding BAM file
    1. Switch to the other window in tmux
    1. `binf`
    1. `igv.sh ${ALN}.so.mdup.br.hc.bam &`
    1. Zoom into the genomic region: 7:XXXXX-XXXXXX

## 4. De Novo Variants in Trio sample VCF
1. Consider, trio samples (proband, mom, and dad) are available
1. Refer to [Genotype Refinment in GATK best practice](https://software.broadinstitute.org/gatk/documentation/article.php?id=4727)
1. Finding De Novo Variants
    1. Clinical molecular diagnostic labs
    1. Familial samples
1. run GATK bash script
    1. `binf`, `cd examples`
    1. `./gatk_gt_refinement.sh trio trio.ped`
    1. `grep PossibleDenovo trio.cgp.filt.denovo.vcf | less -S`

# 5. ASSIGNMENT
1. Locate the second FASTQ file
    1. `$HOME/projects/bioinfo_2018/fastq/chr7b.*.fastq.gz`
1. Run FastQC program to check the quality of FASTQ files
1. Reuse `bwa_gatk_lite.sh`
1. Perform a structure variant analysis using a program called, `shear`
    1. Download [`shear`]((http://vk.cs.umn.edu/SHEAR/download.php?v=1.1.2))
    1. Check an installatin or README file to check the requirement
    1. Install
    1. Run `shear` on the sample to obtain a potential structure variant location
    1. Once you debug each command line, write it into a batch shell script.
1. Load two BAM files (the one we worked during the workshop session and the second BAM file)
1. Go to the region in chr7:xxxx-xxxx to see how two BAM read pileups are different each other. 
