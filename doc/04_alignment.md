# 1. Alignments
Depending on the sequencing technology and applications, we can choose an appropriate alignment program. For example, several alignment software is available depending on read length, molecule type (DNA or RNA), transcriptome (gene expression) analysis, variant calls, alignment of a short read against a reference genome sequence, or alignment between two FASTA format sequence.

In this section, we introduce an NGS short aligner (BWA).

## 2. BWA for short DNA sequencing reads
Most short-read alignment algorithms follow the general procedures such as 1) building an index for a reference target genome sequence, 2) find multiple seed locations (lookup table), 3) extending seeds, 4) pairwise alignment, and 5) selecting the optimal solutions. The first two steps are the most critical regarding speed and accuracy. [BWA](http://bio-bwa.sourceforge.net/bwa.shtml) employes Burrows-Wheeler transform (BWT) algorithm that it uses suffix arrays that efficiently stores the dynamic size of k-mers for a large genome sequence. Refer to the [publications](https://arxiv.org/pdf/1303.3997.pdf) for more technical details.

## Indexing, indexing, and indexing a target reference genome sequence FASTA file 
### 3. Prepare BWT index
First, BWT index is required for a reference genome sequence before we align FASTQ files to the genome sequence.
1. `binf` and `cd ref`
1. `bwa index` #help for `bwa` subcommand `index`
1. `mkdir ./ref/bwa` #prepare bwa index directory
1. `bwa index -p bwa/chr7 chr7.fa`
1. Q18c: to check which files are created

### 4. Prepare FASTA index file containing the sequence length 
1. `samtools faidx chr7.fa` # we will revisit `samtools` again later
1. Q19: Check if an index file is successfully created

### 5. Prepare FASTA dictionary file also containing the sequence length but used by Picard
1. make sure that `java` is in your $PATH
    1. `which java`
    1. `java -v`
1. define `picard` java command prefix
    1. `picard="java -Xmx4g -Djava.io.tmpdir=${HOME}/tmp -jar ${HOME}/bin/picard.jar"`
    1. Picard is a java program. Many java program commands starts this format.
    1. heap memory size
    1. a temporary directory
    1. a jar file path
1. `$picard -h` # for help
1. `$picard CreateSequenceDictionary -h` # to look up run options for the subcommand `CreateSequenceDictionary`
1. `binf`
1. Define some environment variables and run the following lines.
    ```bash
    REFNAME=chr7
    REF_PREFIX=ref/${REFNAME}
    REF_FILE=${REF_PREFIX}.fa
    picard="java -Xmx4g -Djava.io.tmpdir=${HOME}/tmp -jar ${HOME}/bin/picard.jar"
    $picard CreateSequenceDictionary R=${REF_FILE} O=${REF_PREFIX}.dict
    ``` 

### 6. If statement in bash
1. We want to make sure that if a file necessary in a pipeline script exists. Otherwise, we exit the program run.
    ```bash
    # We assume that this script runs at $HOME/projects/bioinfo_2018
    if [ ! -f ${REF_FILE} ]; then
        echo "${REF_FILE} does not exist!"
        exit 1
    fi
    ```
1. We also want to create a directory if it does not exist.
    ```bash
    SAMPLE=chr7a
    OUTD=fastq/${SAMPLE}
        
    if [ ! -d $OUTD ]; then
        echo "Creating a working directory ..."
        mkdir $OUTD
        echo "Done."
    fi
    ```

## BWA mem
### 7. Explore BWA mem commandline options
1. `bwa mem` #help menu for subcommand `mem`
1. [*] means "optional"
1. <*> means "required"
```commandline
Usage: bwa mem [options] <idxbase> <in1.fq> [in2.fq]
Algorithm options:

       -t INT        number of threads [1]
       -k INT        minimum seed length [19]
       -w INT        band width for banded alignment [100]
       -d INT        off-diagonal X-dropoff [100]
       -r FLOAT      look for internal seeds inside a seed longer than {-k} * FLOAT [1.5]
       -y INT        seed occurrence for the 3rd round seeding [20]
       -c INT        skip seeds with more than INT occurrences [500]
       -D FLOAT      drop chains shorter than FLOAT fraction of the longest overlapping chain [0.50]
       -W INT        discard a chain if seeded bases shorter than INT [0]
       -m INT        perform at most INT rounds of mate rescues for each read [50]
       -S            skip mate rescue
       -P            skip pairing; mate rescue performed unless -S also in use

Scoring options:

       -A INT        score for a sequence match, which scales options -TdBOELU unless overridden [1]
       -B INT        penalty for a mismatch [4]
       -O INT[,INT]  gap open penalties for deletions and insertions [6,6]
       -E INT[,INT]  gap extension penalty; a gap of size k cost '{-O} + {-E}*k' [1,1]
       -L INT[,INT]  penalty for 5'- and 3'-end clipping [5,5]
       -U INT        penalty for an unpaired read pair [17]

       -x STR        read type. Setting -x changes multiple parameters unless overriden [null]
                     pacbio: -k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0  (PacBio reads to ref)
                     ont2d: -k14 -W20 -r10 -A1 -B1 -O1 -E1 -L0  (Oxford Nanopore 2D-reads to ref)
                     intractg: -B9 -O16 -L5  (intra-species contigs to ref)

Input/output options:

       -p            smart pairing (ignoring in2.fq)
       -R STR        read group header line such as '@RG\tID:foo\tSM:bar' [null]
       -H STR/FILE   insert STR to header if it starts with @; or insert lines in FILE [null]
       -j            treat ALT contigs as part of the primary assembly (i.e. ignore <idxbase>.alt file)

       -v INT        verbose level: 1=error, 2=warning, 3=message, 4+=debugging [3]
       -T INT        minimum score to output [30]
       -h INT[,INT]  if there are <INT hits with score >80% of the max score, output all in XA [5,200]
       -a            output all alignments for SE or unpaired PE
       -C            append FASTA/FASTQ comment to SAM output
       -V            output the reference FASTA header in the XR tag
       -Y            use soft clipping for supplementary alignments
       -M            mark shorter split hits as secondary

       -I FLOAT[,FLOAT[,INT[,INT]]]
                     specify the mean, standard deviation (10% of the mean if absent), max
                     (4 sigma from the mean if absent) and min of the insert size distribution.
                     FR orientation only. [inferred]
```

### 8. BWA mem command for paired-end DNA sequencing reads 
1. -R (read group): It is okay to run BWA without read group, but it is necessary when we run Picard and GATK SNP caller with the other samples.
    1. ID: Read group identifier containg a flowcell ID and a lane ID
    1. SM: Sample ID
    1. PL: Platform/technology used to produce the read
    1. LB: DNA preparation library identifier
    1. Refer [here](https://software.broadinstitute.org/gatk/documentation/article.php?id=6472) for more information about the group id option
1. -t: Number of threads to utilize
1. -M: mark shorter split hits as a secondary
    1. Considering an application to identify chimeric reads, we enable this option
1. bwa's output is SAM file and we need to convert BAM file to reduce the file size
1. All together, a command can be
```bash
ALN=${OUTD}/bwa
bwa mem \
    -M \
    -t 1 \
    -R '@RG\tID:flowcell.lane1\tSM:sampleID\tPL:ILLUMINA\tLB:libraryID' \
    ${BWA_REF_PREFIX} \
    fastq/${SAMPLE}.R1_001.fastq.gz \
    fastq/${SAMPLE}.R2_001.fastq.gz \
    | samtools view -bS - -o ${ALN}.bam
```
