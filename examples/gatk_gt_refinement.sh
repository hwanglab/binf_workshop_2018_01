#!/usr/bin/bash -l
VCF_PREF=$1
PED=$2

GATK_BUNDLE=/Informatics_workshop/reference/GATK_known_bundle
REF_FILE=${GATK_BUNDLE}/human_g1k_v37_decoy.fasta
kg34vcf=${GATK_BUNDLE}/1000G_phase3_v4_20130502.sites.vcf.gz

gatk="java -Djava.io.tmpdir=${HOME}/tmp -Xmx4G -jar ${HOME}/bin/GenomeAnalysisTK.jar"

echo "===================="
echo "GATK: Run Genotype Refinement (focusing on only chr2,6,22) ..."
$gatk \
 -R ${REF_FILE} \
 -T CalculateGenotypePosteriors \
 --supporting $kg34vcf \
 -ped ${PED} \
 -V ${VCF_PREF}.vcf \
 -L 2,6,22 \
 -o ${VCF_PREF}.cgp.vcf

$gatk \
 -T VariantFiltration \
 -R ${REF_FILE} \
 -V ${VCF_PREF}.cgp.vcf \
 -G_filter "GQ < 20.0" \
 -G_filterName lowGQ \
 -L 2,6,22 \
 -o ${VCF_PREF}.cgp.filt.vcf

$gatk \
 -T VariantAnnotator \
 -R ${REF_FILE} \
 -V ${VCF_PREF}.cgp.filt.vcf \
 -A PossibleDeNovo \
 -L 2,6,22 \
 -ped ${PED} \
 -o ${VCF_PREF}.cgp.filt.denovo.vcf
