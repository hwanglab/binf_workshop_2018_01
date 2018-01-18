## 1. Aligment quality evaluation and post processing
In this session, we understand SAM/BAM file format, mark duplicated reads, gather basic alignment stats to check any abnormality, getting familiar with `samtools` 
### 2. Sorting pileup reads
1. `fastq`; `cd chr7a`
1. Reads in FASTQ files are randomly scattered.
1. Reads aligned the target reference genome sequence are scattered randomly.
1. `samtools view -hS bwa.bam`
1. `samtools sort bwa.bam -o bwa.so.bam`
1. Q20: View the sorted BAM file. Can you see any difference visually?

### 3. Mark PCR duplicate
Duplicate reads are defined as originating from a single fragment of DNA, and it can arise during library construction using PCR. Here, we use Picard "MarkDuplicates" command that compares sequences in the 5' positions of both reads and read pairs in BAM file and marks duplicated reads to keep the only primary copy (with the sums of their base-quality scores).

1. Q21: Using Picard program, mark 5'end duplicated reads in the BAM file we worked last (`bwa.so.bam`). Assign `bwa.so.mdup.bam` to the output BAM file name.
    ```bash
    5 min
    ``` 
    (Hint:Check section 5.6 in alignment.md)

### 4. Index BAM file
1. Before we examine the qulaity of BAM file, we also need generate a BAM index file.
    1. `BAM=bwa.so.mdup.bam`
    1. `samtools index $BAM` #let's assume that you generated `bwa.so.mdup.bam` in the previous step
    1. Check if `${BAM}.bai` is generated
    
## 5. Mapping Stats
Before we start to a long journey with the BAM file, we want to make sure that both quantity and quality of read coverage are good enough before more downstream analyses.
 
### 6. Familiar with SAM/BAM format
1. `samtools view` for help
1. `samtools view -h $BAM | less -S`
1. Check [SAM/BAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf)
1. [SAM flag online lookup](https://broadinstitute.github.io/picard/explain-flags.html)

### 7. Visit UCSC Genome Browser
In general, we are interested in a certain genomic region to check a read depth there
1. Open a web browswer from your computer 
1. Go to UCSC genome browser hgTable website (https://genome.ucsc.edu/cgi-bin/hgTables)
    1. Select assembly: 2009 hg19/b37 
    1. Select group:"Gene and Gene Prediction"
    1. Select track:"RefSeq Gene"
    1. Select table:"refGene"
    1. Select region: "genome"
    1. Select output format:"gtf"
    1. Start to download
1. Open Microsoft Excel
    1. Drag the file you download and drop it into Microsoft Excel spreadsheet

### 8. Basic Alignment Statistics
1. Use QualiMap to generate an alignment QC metrics
    1. Come back to X11 terminal
    1. Run following commands
        ```bash
        fastq
        cd chr7a
        GTF=/Informatics_workshop/reference/annotations/refgene_b37.gtf
        qualimap bamqc -bam $BAM -gff $GTF -outfile qualimap.pdf -outformat PDF
        ```
    1. `cd bwa.so.mdup_stats`
    1. `xpdf qualimap.pdf &`
    1. Q25: Discuss with your group to understand the alignment report.
        ```bash
        5 min
        ```
    1. Explore more functions available in `qualimap`
        1. GUI supported
        1. The other functions to analyze RNASeq, BisulfitedSeq, or ChIPSeq alignment
### 9. Explore read depth and breadth in the annotated regions        

1. Generate BED file from the GTF used in qualimap
1. `fastq`; `cd chr7a` 
1. `grep 'exon' $GTF | ucscgtf_to_bed.py | sort -k1,1V -k2,2n -k3,3n > refgene.bed`
    1. `grep 'exon' $GTF`: Print out lines containing a word `exon` in `$GTF` file
    1. Use a simple python script (`ucscgtf_tobed.py`) to convert each line (output of the step a.) to BED file format
        1. Search "bed file format" in google
        1. Note that `start` is 0-based position but `end` is 1-based position.   
    1. `sort -k1,1V -k2,2n -k3,3n`: Sort the output (BED file) of the step b in an order of the 1st column (alphabet), 2nd (numeric), and 3rd column (numeric)
    1. `> refgene.bed`: Finally, pipe out (`>`) to a file instead of showing lines in the terminal screen

1. Q25a. View `refgene.bed` to see how it looks
1. Q25b. View `$GTF` to see how it looks
1. `samtools view -hb -q 5 -F 0x400 $BAM | bedtools coverage -hist -a refgene.bed -b - | gzip -fc > ${BAM}.cvg.hist.gz`
    1. `samtools view -hb -q 5 -F 0x400 $BAM` # to extract the ones of which mapping qulaity >= 5 and not marked "Duplicated" from the BAM file  
    1. `bedtools coverage -hist -b - -a refgene.bed` # run `bedtools coverage` and calculate both the depth and breadth of coverage of features (i.e., `refgene.bed`)
    1. `gzip -fc > ${BAM}.cvg.hist.gz` # this time, zip the output file (`-f`:overwrite if the output file already exists, `-c`: compress)
1. Q26: View the bedtool histogram output

### 10. Tips
1. Always, peek the other sources similar to what you want to do before you write your script
1. Rely on programs commonly used in community for many years
1. Start with a small sample set first
1. Update your script one by one
1. Modulize your script as much as possible so that you can reuse the intermediary result and glue them together with the other scripts 
