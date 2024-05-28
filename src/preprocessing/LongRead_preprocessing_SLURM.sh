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
samplesheet=$2 # A simple file with no header and the sample code and the corresponding barcode separated by a blank space. E.g. 20-2364 1

ml load singularity

base_in=${base}/1.BASECALLING/pass
base_out=${base}/2.TRIMMING

#declare -A codes=()
#while read line;do
#  sample=$(echo "$line" | cut -d " " -f1)
#  bar=$(echo "$line" | cut -d " " -f2)
#  codes[${bar}]=${sample}
#done < ${samplesheet}

 

#for bar in ${!codes[@]};do
#  sam=${codes[${bar}]}
#  out_fold=${base_out}/${sam}
#  mkdir -p ${out_fold}
#  if [[ ${bar} -lt 10 ]];then
#    workdir=${base_in}/barcode0${bar}
#  else
#    workdir=${base_in}/barcode${bar}
#  fi
while read line;do

  sam=$(echo "$line" | cut -d " " -f1)
  bar=$(echo "$line" | cut -d " " -f2)
  out_fold="${base_out}/${sam}"
  
  if [[ ! -d ${out_fold} ]];then
    mkdir -p ${out_fold}
  fi

  if [[ ${bar} -lt 10 ]];then
    workdir=${base_in}/barcode0${bar}
    fastqs=$(ls -rt ${workdir}/*.fastq.gz | wc -l)
    echo "Fastq number for sample ${sam}: ${fastqs}"
    #cd ${workdir};
    #fastqs=$(ls -rt *.fastq.gz | awk '{print "/workdir/"$0}')
    #echo "${fastqs}"
  else
    workdir=${base_in}/barcode${bar}
    fastqs=$(ls -rt ${workdir}/*.fastq.gz | wc -l)
    echo "Fastq number for sample ${sam}: ${fastqs}"
    #cd ${workdir};
    #fastqs=$(ls -rt *.fastq.gz | awk '{print "/workdir/"$0}')
    #echo "${fastqs}"
  fi
  singularity  exec --no-home --no-mount ${PWD} \
  -B ${workdir} \
  -B ${out_fold}:/out_fold \
  docker://lisosome/minion_side:latest \
  zcat $(ls -rt ${workdir}/*.fastq.gz) | NanoLyse | NanoFilt -q 10 -l 500 --headcrop 50 | gzip -c > ${out_fold}/${sam}.trimmed_and_clean.fastq.gz
done < ${samplesheet}
