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
outdir=$4

model_dir=/fast/burlo/nardone/clair3_model
refdir=$(dirname ${ref})
refname=$(basename ${ref})
MODEL_NAME="r1041_e82_400bps_sup_v430"

fai="${ref}.fai"

if [[ ! -f ${fai} ]];then
    ml load samtools
    samtools faidx ${ref}
fi

output_dir=${outdir}/${sample_code}

if [[ ! -d ${output_dir} ]];then mkdir -p ${output_dir};done

singularity exec \
    -B ${input_dir} \
    -B ${output_dir} \
    -B ${refdir} \
    -B ${model_dir} \
    docker://hkubal/clair3:latest \
    /opt/bin/run_clair3.sh \
    --bam_fn="${input_dir}/${sample_code}.bam" \
    --ref_fn=${refdir}/${refname} \
    --threads=8 \
    --platform="ont" \
    --model_path="${model_dir}/${MODEL_NAME}" \
    --output="${output_dir}" && mv ${output_dir}/merge_output.vcf.gz ${output_dir}/${sample_code}_STRC.vcf.gz &&
mv ${output_dir}/merge_output.vcf.gz.tbi ${output_dir}/${sample_code}_STRC.vcf.gz.tbi &&
rm -r ${output_dir}/log/ ${output_dir}/tmp/ ${output_dir}/full* ${output_dir}/pileup* ${output_dir}/*.log