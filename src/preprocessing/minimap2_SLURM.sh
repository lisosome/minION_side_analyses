#!/usr/bin/env bash

#SBATCH -p THIN
#SBATCH --mem=150gb
#SBATCH --cpus-per-task 15
#SBATCH --ntasks-per-node=1
#SBATCH --time=96:00:00


fastq_dir=$1
out=$2
ref=$3
samplesheet=$4

REF_DIR=$(dirname ${ref})
REF_FASTA=$(basename ${ref})
INPUT_DIR=${fastq_dir}

module load singularity
module load samtools

if [[ ! -d ${out} ]];then
  mkdir -p ${out}
fi

#declare -A codes=()
#while read line;do
#  sample=$(echo "$line" | cut -d " " -f1)
#  bar=$(echo "$line" | cut -d " " -f2)
#  codes[${bar}]=${sample}
#done < ${samplesheet}

while read line;do

sample=$(echo "$line" | cut -d " " -f1)
INPUT_FAST="${sample}.trimmed_and_clean.fastq.gz"

singularity  exec --no-home --no-mount ${PWD} \
  -B ${REF_DIR} \
  -B ${INPUT_DIR} \
  -B ${out} \
  docker://lisosome/minion_side:latest \
  minimap2 -y -t 15 -ax map-ont ${REF_DIR}/${REF_FASTA} ${INPUT_DIR}/${sample}/${INPUT_FAST} | samtools view -@ 15 -b | samtools addreplacerg -r "@RG\tID:${sample}\tSM:${sample}" -  | samtools sort -@ 15 --write-index -o ${out}/${sample}.bam\#\#idx\#\#${out}/${sample}.bam.bai

if [[ ! -d ${out}/coverage ]];then
  mkdir -p ${out}/coverage
fi

singularity  exec --no-home --no-mount ${PWD} \
  -B ${out} \
  docker://lisosome/minion_side:latest \
  mosdepth -x -n -t 15 ${out}/coverage/${sample} ${out}/${sample}.bam

done < ${samplesheet}