#!/usr/bin/env bash

#SBATCH -p THIN
#SBATCH --mem-per-cpu=5gb
#SBATCH --cpus-per-task 10
#SBATCH --ntasks-per-node=1
#SBATCH --time=96:00:00

set -e


outfold=$1
sample=$2
mode=$3

ml load samtools
dorado=/orfeo/LTS/burlo/LT_storage/shared/tools/dorado-0.7.1-linux-x64/bin/dorado
mosdepth=/orfeo/LTS/burlo/LT_storage/shared/tools/mosdepth


preproc=${outfold}/1.PRE.PROCESSING/DEMULTIPLEXED
aligndir=${outfold}/2.ALIGNMENT

if [ ! -d ${aligndir} ];then
    mkdir -p ${aligndir}
fi


if [ ! -d ${aligndir}/${sample}/coverage ];then
    mkdir -p ${aligndir}/${sample}/coverage
fi


if [[ ${mode} == "MT" ]];then
    ref="/fast/burlo/nardone/mtDNA_resources/rcrs_mutserve.fasta"
    ${dorado} aligner ${ref} ${preproc}/${sample}/${sample}.bam | samtools addreplacerg -r "@RG\tID:${sample}\tSM:${sample}" - | samtools sort -O BAM -o ${aligndir}/${sample}/${sample}_MT.bam && samtools index ${aligndir}/${sample}/${sample}_MT.bam &&
    ${mosdepth} -n -x ${aligndir}/${sample}/coverage/${sample}_MT ${aligndir}/${sample}/${sample}_MT.bam
else
    ref="/fast/burlo/nardone/STRC_ref/STRC_hg38.fasta"
    ${dorado} aligner ${ref} ${preproc}/${sample}/${sample}.bam | samtools addreplacerg -r "@RG\tID:${sample}\tSM:${sample}" - | samtools sort -O BAM -o ${aligndir}/${sample}/${sample}_STRC.bam && samtools index ${aligndir}/${sample}/${sample}_STRC.bam &&
    ${mosdepth} -n -x -b /fast/burlo/nardone/STRC_ref/STRC_hg38.bed ${aligndir}/${sample}/coverage/${sample}_STRC ${aligndir}/${sample}/${sample}_STRC.bam &&
    awk -F "\t" 'BEGIN{FS=OFS}(NR == 1 || NR == 3){print $0}' ${aligndir}/${sample}/coverage/${sample}_STRC.mosdepth.summary.txt > ${aligndir}/${sample}/coverage/${sample}_STRC_only.mosdepth.summary.txt
fi





