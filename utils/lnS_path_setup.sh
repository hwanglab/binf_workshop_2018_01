#!/bin/bash -l

if [ ! -d ${HOME}/bin ]; then
    echo "Creating ${HOME}/bin directory ..."
    mkdir ${HOME}/bin
    echo "Done."
fi

ln -s /usr/local/tools/gatk/GenomeAnalysisTK.jar $HOME/bin/
ln -s /usr/local/tools/picard/picard.jar $HOME/bin/
ln -s /usr/local/tools/igv/igv.sh $HOME/bin/
ln -s /usr/local/tools/cutadapt/cutadapt $HOME/bin/
ln -s /usr/local/tools/bwa/bwa $HOME/bin/
ln -s /usr/local/tools/bedtools/bedtools $HOME/bin/
ln -s /usr/local/tools/fastqc/fastqc $HOME/bin/
ln -s /usr/local/tools/samtools/samtools $HOME/bin/
ln -s /Informatics_workshop/tools/annovar/table_annovar.pl $HOME/bin/
