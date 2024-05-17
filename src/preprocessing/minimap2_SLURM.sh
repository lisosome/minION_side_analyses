#!/usr/bin/env bash

#SBATCH -p THIN
#SBATCH --mem=150gb
#SBATCH --cpus-per-task 15
#SBATCH --ntasks-per-node=1
#SBATCH --time=96:00:00


fastq=$1
out=$2
ref=$3
sample=$4

module load samtools

source activate lr_process

if [[ ! -d ${out} ]];then
  mkdir -p ${out}
fi

if [[  ${ref} == "hg38"  ]];then
  fasta=/fast/burlo/fcrudele/resources/hgRef/GRCh38.p13/no_alt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
  minimap2 -y -t 15 -ax map-ont ${fasta} ${fastq} | samtools view -@ 15 -b | samtools addreplacerg -r "@RG\tID:${sample}\tSM:${sample}" | samtools sort -@ 15 --write-index -o ${out}/${sample}_${ref}.bam\#\#idx\#\#${out}/${sample}_${ref}.bai
  mkdir -p ${out}/coverage
  mosdepth -x -n -t 15 ${out}/coverage/${sample} ${out}/${sample}_${ref}.bam 
elif [[ ${ref} == "T2T" ]];then
  fasta=/fast/burlo/fcrudele/resources/hgRef/T2T-CHM13/chm13v2.0.fa
  minimap2 -y -t 15 -ax map-ont ${fasta} ${fastq} | samtools addreplacerg -r "@RG\tID:${sample}\tSM:${sample}" | samtools view -@ 15 -b | samtools sort -@ 15 --write-index -o ${out}/${sample}_${ref}.bam\#\#idx\#\#${out}/${sample}_${ref}.bai
  mkdir -p ${out}/coverage
  mosdepth -x -n -t 15 ${out}/coverage/${sample} ${out}/${sample}_${ref}.bam
else
  minimap2 -y -t 15 -ax map-ont ${ref} ${fastq} | samtools addreplacerg -r "@RG\tID:${sample}\tSM:${sample}" | samtools view -@ 15 -b | samtools sort -@ 15 --write-index -o ${out}/${sample}.minimap.bam\#\#idx\#\#${out}/${sample}.minimap.bai
  mkdir -p ${out}/coverage
  mosdepth -x -n -t 15 ${out}/coverage/${sample} ${out}/${sample}.minimap.bam
fi
