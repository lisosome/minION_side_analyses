#!/usr/bin/env bash

# This script will perform the pre-processing step for ONT long reads from minION
# WARNING: At this day, Orfeo does not have some utilities such as tabix and gbzip installed in all the system, so all this operation will rely on a conda environment named "minimap2"

#SBATCH -p THIN
#SBATCH --mem=20gb
#SBATCH --cpus-per-task 5
#SBATCH --ntasks-per-node=1
#SBATCH --time=96:00:00

# This script assumes that the pre process is coming after the basecalling so will create a folder named "2.TRIMMING"
# This script assumes that the folder containing the basecalled fastqs will be called "1.BASECALLING/pass"
base=$1 # Basefolder
mode=$2 # Execution mode

case $mode in
  SINGLE_SAMPLE)
sample=$3
base_in=${base}/1.BASECALLING/pass
base_out=${base}/2.TRIMMING

if [[ ! -d ${base_out} ]];then
  mkdir -p ${base_out}
fi

module load conda
source activate lr_process

cat $(ls -rt ${base_in}/*.fastq) | NanoLyse | NanoFilt -q 10 -l 500 --headcrop 50 | bgzip -c > ${base_out}/${sample}.trimmed_and_clean.fastq.gz #&& bgzip ${base_out}/${sample}.trimmed_and_clean.fastq
;;
  BARCODE)
    samplesheet=$3 # A simple file with no header and the sample code and the corresponding barcode separated by a blank space. E.g. 20-2364 1

    base_in=${base}/1.BASECALLING/pass
    base_out=${base}/2.TRIMMING

    declare -A codes=()
    while read line;do
      sample=$(echo "$line" | cut -d " " -f1)
      bar=$(echo "$line" | cut -d " " -f2)
      codes[${bar}]=${sample}
    done < ${samplesheet}

 
  source activate lr_process

  for bar in ${!codes[@]};do
    sam=${codes[${bar}]}
    out_fold=${base_out}/${sam}
    mkdir -p ${out_fold}
    if [[ ${bar} -lt 10 ]];then
      workdir=${base_in}/barcode0${bar}
    else
      workdir=${base_in}/barcode${bar}
    fi

    zcat $(ls -rt ${workdir}/*.fastq.gz) | NanoLyse | NanoFilt -q 10 -l 500 --headcrop 50 > ${out_fold}/${sam}.trimmed_and_clean.fastq && bgzip ${out_fold}/${sam}.trimmed_and_clean.fastq
    done
    ;;
esac
