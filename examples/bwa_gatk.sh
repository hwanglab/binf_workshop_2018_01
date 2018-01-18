#!/usr/bin/bash -l

SAMPLE=$1
REFNAME=$2

NCPU=1
REF_PREFIX=ref/${REFNAME}
REF_FILE=${REF_PREFIX}.fa

if [ ! -f ${REF_FILE} ]; then
    echo "${REF_FILE} does not exist!"
    exit 1
fi

if [ ! -f ${REF_FILE}.fai ]; then
	echo "Creating a reference index file ..."
    echo "samtools faidx ${REF_FILE}"
    samtools faidx ${REF_FILE}
    echo "Done."
fi

picard="java -Xmx4g -Djava.io.tmpdir=${HOME}/tmp -jar ${HOME}/bin/picard.jar"
if [ ! -f ${REF_PREFIX}.dict ]; then
	echo "Creating a reference dict file ..."
    echo "$picard CreateSequenceDictionary R=${REF_FILE} O=${REF_PREFIX}.dict"
    $picard CreateSequenceDictionary R=${REF_FILE} O=${REF_PREFIX}.dict
    echo "Done."
fi

if [ ! -f fastq/${SAMPLE}.R1_001.fastq.gz ]; then
    echo "fastq/${SAMPLE}.R1_001.fastq.gz does not exist!"
    exit 1
fi

if [ ! -f fastq/${SAMPLE}.R2_001.fastq.gz ]; then
    echo "fastq/${SAMPLE}.R2_001.fastq.gz does not exist!"
    exit 1
fi

BWA_REF_PREFIX=ref/bwa/${REFNAME}
if [ ! -f ${BWA_REF_PREFIX}.bwt ]; then
    echo "Creating BWT file for ${REF_FILE} ..."
    mkdir ref/bwa
    echo "bwa index -p ${BWA_REF_PREFIX} ${REF_FILE}"
    bwa index -p ${BWA_REF_PREFIX} ${REF_FILE}
    echo "Done."
fi

OUTD=fastq/${SAMPLE}
if [ ! -d $OUTD ]; then
	echo "Creating a working directory ..."
    mkdir $OUTD
    echo "Done."
fi

# -------------------------------------------------------------------
ALN=${OUTD}/bwa

echo "bwa mem -t $NCPU -R '@RG\tID:flowcell.lane1\tSM:sampleID\tPL:ILLUMINA\tLB:libraryID' ${BWA_REF_PREFIX} fastq/${SAMPLE}.R1_001.fastq.gz fastq/${SAMPLE}.R2_001.fastq.gz | samtools view -bS - -o ${ALN}.bam"
bwa mem -t $NCPU -R '@RG\tID:flowcell.lane1\tSM:sampleID\tPL:ILLUMINA\tLB:libraryID' ${BWA_REF_PREFIX} fastq/${SAMPLE}.R1_001.fastq.gz fastq/${SAMPLE}.R2_001.fastq.gz | samtools view -bS - -o ${ALN}.bam

echo "samtools sort ${ALN}.bam -o ${ALN}.so.bam"
samtools sort ${ALN}.bam -o ${ALN}.so.bam

echo "Deleting ${ALN}.bam"
# rm -rf ${ALN}.bam

echo "Dedeuping ..."
echo "$picard MarkDuplicates I=${ALN}.so.bam O=${ALN}.so.mdup.bam M=${ALN}.so.mdup.table"
$picard MarkDuplicates I=${ALN}.so.bam O=${ALN}.so.mdup.bam M=${ALN}.so.mdup.table
echo "Done."

echo "indexing bam file ..."
echo "samtools index ${ALN}.so.mdup.bam"
samtools index ${ALN}.so.mdup.bam
echo "Done."

GTF=/Informatics_workshop/reference/annotations/refgene_b37.gtf
echo "qualimap bamqc -bam ${ALN}.so.mdup.bam -gff $GTF -outfile qualimap.pdf -outformat PDF"
qualimap bamqc -bam ${ALN}.so.mdup.bam -gff $GTF -outfile qualimap.pdf -outformat PDF

echo "grep 'exon' $GTF | ucscgtf_to_bed.py | sort -k1,1V -k2,2n -k3,3n > fastq/${SAMPLE}/refgene.bed"
grep 'exon' $GTF | ucscgtf_to_bed.py | sort -k1,1V -k2,2n -k3,3n > fastq/${SAMPLE}/refgene.bed
BED=fastq/${SAMPLE}/refgene.bed

echo "samtools view -hb -q 5 -F 0x400 ${ALN}.so.mdup.bam | bedtools coverage -hist -b - -a $BED | gzip -fc > ${ALN}.so.mdup.bam.cvg.hist.gz"
samtools view -hb -q 5 -F 0x400 ${ALN}.so.mdup.bam | bedtools coverage -hist -b - -a $BED | gzip -fc > ${ALN}.so.mdup.bam.cvg.hist.gz

# -------------------------------------------------------------------

# echo "Collecting bam stats ..."
# echo "$picard CollectAlignmentSummaryMetrics R=${REF_FILE} I=${ALN}.so.mdup.bam O=${ALN}.so.mdup.bam.alnstats"
# $picard CollectAlignmentSummaryMetrics R=${REF_FILE} I=${ALN}.so.mdup.bam O=${ALN}.so.mdup.bam.alnstats

# target_bed=/home/hongc2/projects/bioinfo_2018/ref/S07604514_Covered_chr22_e200c.bed
# bait_bed=/home/hongc2/projects/bioinfo_2018/ref/S07604514_Covered_chr22.bed
#
# echo "$picard BedToIntervalList I=${target_bed} O=${target_bed}.ilist SD=${REF_PREFIX}.dict"
# $picard BedToIntervalList I=${target_bed} O=${target_bed}.ilist SD=${REF_PREFIX}.dict
#
# echo "$picard BedToIntervalList I=${bait_bed} O=${bait_bed}.ilist SD=${REF_PREFIX}.dict"
# $picard BedToIntervalList I=${bait_bed} O=${bait_bed}.ilist SD=${REF_PREFIX}.dict
#
# echo "$picard CollectHsMetrics NEAR_DISTANCE=0 R=${REF_FILE} I=${ALN}.so.mdup.bam TARGET_INTERVALS=${target_bed}.ilist BAIT_INTERVALS=${bait_bed}.ilist O=${ALN}.so.mdup.bam.hsstats"
# $picard CollectHsMetrics NEAR_DISTANCE=0 R=${REF_FILE} I=${ALN}.so.mdup.bam TARGET_INTERVALS=${target_bed}.ilist BAIT_INTERVALS=${bait_bed}.ilist O=${ALN}.so.mdup.bam.hsstats
#
# GAT=ref/annotation/refGene_e20_so_merged_chr22.bed
#
# echo "$picard BedToIntervalList I=${GAT} O=${GAT}.ilist SD=${REF_PREFIX}.dict"
# $picard BedToIntervalList I=${GAT} O=${GAT}.ilist SD=${REF_PREFIX}.dict
#
# echo "$picard CollectWgsMetricsWithNonZeroCoverage I=${ALN}.so.mdup.bam INTERVALS=${GAT}.ilist O=${ALN}.so.mdup.bam.ga.stats CHART=${ALN}.so.mdup.bam.ga.stats.pdf R=${REF_FILE}"
#
# $picard CollectWgsMetricsWithNonZeroCoverage I=${ALN}.so.mdup.bam INTERVALS=${GAT}.ilist O=${ALN}.so.mdup.bam.ga.stats CHART=${ALN}.so.mdup.bam.ga.stats.pdf R=${REF_FILE}

gatk="java -Djava.io.tmpdir=${HOME}/tmp -Xmx4G -jar ${HOME}/bin/GenomeAnalysisTK.jar"
GATKDB=/Informatics_workshop/reference/GATK_known_bundle
dbsnp_vcf=${GATKDB}/dbsnp_138.hg19_no_chr.vcf
hapmap_vcf=${GATKDB}/hapmap_3.3.hg19.sites_no_chr.vcf
omni_vcf=${GATKDB}/1000G_omni2.5.hg19.sites_no_chr.vcf
kG_vcf=${GATKDB}/1000G_phase1.snps.high_confidence.b37.vcf
mills_vcf=${GATKDB}/Mills_and_1000G_gold_standard.indels.hg19.sites_no_chr.vcf

echo "====================="
echo "GATK:Base Recalibrating ..."
$gatk \
    -T BaseRecalibrator \
    -R ${REF_FILE} \
    -I ${ALN}.so.mdup.bam \
    -knownSites $dbsnp_vcf \
    -knownSites $mills_vcf \
    -nct ${NCPU} \
    -o ${OUTD}/recal_data.table

$gatk \
    -T BaseRecalibrator \
    -R ${REF_FILE} \
    -I ${ALN}.so.mdup.bam \
    -knownSites $dbsnp_vcf \
    -knownSites $mills_vcf \
    -BQSR ${OUTD}/recal_data.table \
    -nct ${NCPU} \
    -o ${OUTD}/post_recal_data.table

$gatk \
    -T PrintReads \
    -R ${REF_FILE} \
    -I ${ALN}.so.mdup.bam \
    -BQSR ${OUTD}/post_recal_data.table \
    --defaultBaseQualities 25 \
    -nct ${NCPU} \
    -o ${ALN}.so.mdup.br.bam
    
echo "Done."


echo " ====================="
echo "GATK:HaplotypeCaller ..."    
$gatk \
    -T HaplotypeCaller \
    -R ${REF_FILE} \
    -I ${ALN}.so.mdup.br.bam \
    --genotyping_mode DISCOVERY \
    --dbsnp $dbsnp_vcf \
    -stand_call_conf 30 \
    -nct ${NCPU} \
    -o ${OUTD}/gatkhc.vcf
echo "Done."

echo "======================="
echo "GATK: applying adaptive filter ..."
#NOTE: this is an option for WES. For multi-sample or WGS, you can add DP and inbreeding_coeff
$gatk \
    -T VariantRecalibrator \
    -R ${REF_FILE} \
    -input ${OUTD}/gatkhc.vcf \
    -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap_vcf \
    -resource:omni,known=false,training=true,truth=true,prior=12.0 $omni_vcf \
    -resource:1000G,known=false,training=true,truth=false,prior=10.0 $kG_vcf \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp_vcf \
    -an QD \
    -an FS \
    -an SOR \
    -an MQ \
    -an MQRankSum \
    -an ReadPosRankSum \
    -mode SNP \
    -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
    -recalFile ${OUTD}/recalibrate_SNP.recal \
    -tranchesFile ${OUTD}/recalibrate_SNP.tranches \
    -rscriptFile ${OUTD}/recalibrate_SNP_plots.R \
    -nt $NCPU \

$gatk \
    -T ApplyRecalibration \
    -R ${REF_FILE} \
    -input ${OUTD}/gatkhc.vcf \
    -mode SNP \
    --ts_filter_level 99.5 \
    -recalFile ${OUTD}/recalibrate_SNP.recal \
    -tranchesFile ${OUTD}/recalibrate_SNP.tranches \
    -nt $ncpu \
    -o ${OUTD}/gatkhc_vrsnp.vcf

#NOTE: this is an option for WES. For multi-sample or WGS, you can add DP and inbreeding_coeff
$gatk \
    -T VariantRecalibrator \
    -R ${REF_FILE} \
    -input ${OUTD}/gatkhc_vrsnp.vcf \
    -resource:mills,known=true,training=true,truth=true,prior=12.0 $mills_vcf \
    -an QD \
    -an FS \
    -an SOR \
    -an MQRankSum \
    -an ReadPosRankSum \
    -mode INDEL \
    -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
#     --maxGaussians 4 \
    -recalFile ${OUTD}/recalibrate_INDEL.recal \
    -tranchesFile ${OUTD}/recalibrate_INDEL.tranches \
    -rscriptFile ${OUTD}/recalibrate_INDEL_plots.R \
    -nt $NCPU

$gatk \
    -T ApplyRecalibration \
    -R ${REF_FILE} \
    -input ${OUTD}/gatkhc_vrsnp.vcf \
    -mode INDEL \
    --ts_filter_level 99.0 \
    -recalFile ${OUTD}/recalibrate_INDEL.recal \
    -tranchesFile ${OUTD}/recalibrate_INDEL.tranches \
    -nt $NCPU \
    -o ${OUTD}/gatkhc_vr.vcf
    
# echo "===================="
# echo "GATK: Run Genotype Refinement ..."
# $gatk \
#     -R ${REF_FILE} \
#     -T CalculateGenotypePosteriors \
#     --supporting 1000G_phase3_v4_20130502.sites.vcf \
#     -ped trio.ped \
#     -V ${OUTD}/gatkhc_vr.vcf \
#     -o ${OUTD}/gatkhc_vr_cgp.vcf
# 
# $gatk \
#     -T VariantFiltration \
#     -R ${REF_FILE} \
#     -V ${OUTD}/gatkhc_vr_cgp.vcf \
#     -G_filter "GQ < 20.0" \
#     -G_filterName lowGQ \
#     -o ${OUTD}/gatkhc_vr_cgp_filt.vcf
# 
# $gatk \
#     -T VariantAnnotator \
#     -R ${REF_FILE} \
#     -V ${OUTD}/gatkhc_vr_cgp_filt.vcf \
#     -A PossibleDeNovo \
#     -ped trio.ped \
#     -o ${OUTD}/gatkhc_vr_cgp_filt_denovo.vcf

echo " ====================="
echo "ANNOVAR ..."
ANNOVAR_DIR=/Informatics_workshop/tools/annovar
# Annotate with Annovar
$ANNOVAR_DIR/table_annovar.pl ${OUTD}/gatkhc_vr.vcf $ANNOVAR_DIR/humandb/ -buildver hg19 -out ${OUTD}/gatkhc_vr.annovar -remove -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2015aug_all,1000g2015aug_eur,exac03,avsnp147,dbnsfp30a,cosmic80 -operation g,r,r,f,f,f,f,f,f,f -nastring . -vcfinput

echo "Check [${ALN}.so.mdup.br.hc_filt.annovar]"
echo "Done."

echo " =================="
echo "clean up ..."
rm -rf ${ALN}.bam
echo "Done."