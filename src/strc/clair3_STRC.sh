#!/usr/bin/env bash

#SBATCH -p THIN
#SBATCH --mem-per-cpu=5gb
#SBATCH --cpus-per-task 8
#SBATCH --ntasks-per-node=1
#SBATCH --time=96:00:00


ml load singularity

input_dir=$1
sample_code=$2
ref=$3
output_dir=$4

model_dir=/fast/burlo/nardone/clair3_model/
refdir=$(dirname ${ref})
refname=$(basename ${ref})
MODEL_NAME="r1041_e82_400bps_sup_v430"

singularity exec \
    -B ${input_dir}:/input \
    -B ${output_dir}:/output \
    -B ${refdir}:/ref \
    -B ${model_dir}:/model \
    docker://hkubal/clair3:latest \
    /opt/bin/run_clair3.sh \
    --bam_fn="/input/${sample_code}.bam" \    ## change your bam file name here
    --ref_fn=/ref/${refname} \       ## change your reference file name here
    --threads=8 \               ## maximum threads to be used
    --platform="ont" \                   ## options: {ont,hifi,ilmn}
    --model_path="/model/${MODEL_NAME}" \
    --output="/output"