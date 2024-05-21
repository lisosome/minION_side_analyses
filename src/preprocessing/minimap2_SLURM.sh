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
INPUT_DIR=$(dirname ${fastq_dir})

module load singularity


if [[ ! -d ${out} ]];then
  mkdir -p ${out}
fi

declare -A codes=()
while read line;do
  sample=$(echo "$line" | cut -d " " -f1)
  bar=$(echo "$line" | cut -d " " -f2)
  codes[${bar}]=${sample}
done < ${samplesheet}

for sam in ${codes[@]};do
sample=${codes[${sam}]}
INPUT_FAST="${sample}.trimmed_and_clean.fastq.gz"

singularity  exec --no-home --no-mount ${PWD} \
  -B ${REF_DIR}:/ref \
  -B ${INPUT_DIR}:/input \
  -B ${out}:/outdir \
  docker://lisosome/minion_side:latest \
minimap2 -y -t 15 -ax map-ont /ref/${REF_FASTA} /input/${sample}/${INPUT_FAST} \
| samtools addreplacerg -r "@RG\tID:${sample}\tSM:${sample}" \
| samtools view -@ 15 -b | samtools sort -@ 15 --write-index -o /outdir/${sample}.bam\#\#idx\#\#/outdir/${sample}.bam.bai

if [[ ! -d ${out} ]];then
  mkdir -p ${out}/coverage
fi

singularity  exec --no-home --no-mount ${PWD} \
  -B ${out}:/base \
  docker://lisosome/minion_side:latest \
mosdepth -x -n -t 15 /base/coverage/${sample} /base/${sample}.bam

done