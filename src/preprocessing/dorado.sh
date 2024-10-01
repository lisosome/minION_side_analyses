#!/usr/bin/env bash

#SBATCH -p GPU
#SBATCH --gpus=1
#SBATCH --mem-per-cpu=5gb
#SBATCH --cpus-per-task 20
#SBATCH --ntasks-per-node=1
#SBATCH --time=18:00:00

set -e
# Script to perform basecalling using dorado. This allow as to work directly with a BAM file and perform variant calling from there
ml load samtools
rawdata=$1  # Folder containing minION raw data
base_out=$2 # Folder that will store the output folders (e.g 1.ALIGNMENT 2.VARIANT_CALLING)
samplesheet=$3

dorado=/orfeo/LTS/burlo/LT_storage/shared/tools/dorado-0.7.1-linux-x64/bin/dorado

final_summary=$(find ${rawdata} -type f -name "final_summary*.txt")

if [[ -f ${final_summary} ]];then
    kit_name=$(grep "protocol=" ${final_summary} | cut -d ":" -f3)
else
    echo "Final summary file not found in ${rawdata}. Assuming SQK-NBD114-24 for demultiplexing"
    kit_name=SQK-NBD114-24
fi

args="-v --min-qscore 8 -x cuda:0 --kit-name ${kit_name} --trim all"
module load cuda/11.7
CUDA_VISIBLE_DEVICES=0
declare CUDA_VISIBLE_DEVICES

outfold=${base_out}/1.PRE.PROCESSING 

if [ ! -d ${base_out}/0.RAW.DATA ];then
    mkdir -p ${base_out}/0.RAW.DATA/
fi

if [ ! -d ${outfold}/DEMULTIPLEXED ];then
    mkdir -p ${outfold}/DEMULTIPLEXED
fi


find ${rawdata} -type f -name "*.pod5" -exec ln -s {} ${base_out}/0.RAW.DATA/ \;


${dorado} basecaller sup@v5.0.0 ${base_out}/0.RAW.DATA ${args} > ${outfold}/basecalled.bam &&
${dorado} demux --emit-summary --no-classify --no-trim -t 10 -o ${base_out}/1.PRE.PROCESSING/DEMULTIPLEXED ${outfold}/basecalled.bam &&
python3 /orfeo/LTS/burlo/LT_storage/shared/tools/parse_samplesheet.py -s ${samplesheet} -d ${base_out}/1.PRE.PROCESSING/DEMULTIPLEXED -k ${kit_name}

