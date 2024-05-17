#!/usr/bin/env bash

#SBATCH -p THIN
#SBATCH --mem-per-cpu=20gb
#SBATCH --cpus-per-task 5
#SBATCH --ntasks-per-node=1
#SBATCH --time=96:00:00

set -e

ml load singularity
ml load java
ml load bcftools

printHelp(){
    echo "mtDNA Variant Calling"
    echo 
    echo "Usage: mtDNA_variantCall.sh [options]"
    echo "Options:"
    echo "-h, --help - Show this help."
    echo "-b, --bam - Input BAM"
    echo "-s, --samplecode - Sample identifier (DNA code or other)"
    echo "-o, --outdir - Output Directory"
    echo "-t, --threads - Threads to use"
}


testNCreate(){
    fol=$1
    if [ ! -d $fol ];then
        mkdir -p $fol
    fi
}

testOutput(){
    file=$1
    step=$2
    sample=$3
    if [ -f $file ];then
        echo "$step completed for sample ${sample}"
    else
        echo "$step failed for sample $sample. Exiting..."
        exit 1
    fi
}



if [[ $# -eq 0 ]];then 
    printHelp
    exit 0
fi


while [[ $# -gt 0 ]];do
    arg=$1
    case ${arg} in
        
        -h|--help)
        printHelp
        exit 0
        ;;

        -b|--bam)
        bam=$2
        shift
        ;;

        -s|--samplecode)
        sample=$2
        shift
        ;;

        -o|--outdir)
        outd=$2
        shift
        ;;
        
        -t|--threads)
        threads=$2
        shift
        ;;
        
        *)
        echo "Unknown option: $1"
        printHelp
        exit 0
    esac
    shift

done

sif=/orfeo/cephfs/home/burlo/nardone/software/mtdna-server-2/singularity/quay.io-genepi-mtdna-server-2-v2.1.6.img
res_dir=/orfeo/cephfs/home/burlo/nardone/software/mtdna-server-2/files
annotation_file=${res_dir}/rCRS_annotation.txt
ref_mutserve=${res_dir}/rcrs_mutserve.fasta
ref_mutect2=${res_dir}/mt_contigs.fasta

SINGULARITY_TMPDIR=/fast/burlo/nardone/singularity/sing_temp
SINGULARITY_CACHEDIR=/fast/burlo/nardone/singularity/cache

export SINGULARITY_TMPDIR
export SINGULARITY_CACHEDIR

mutserve=/orfeo/cephfs/home/burlo/nardone/software/mutserve2/mutserve.jar

for fol in ${outd}/mutserve ${outd}/mutect2 ${outd}/filtering ${outd}/merging ${outd}/annotation;do
    if [ ! -d ${fol} ];then
        mkdir -p ${fol}
    fi

done

# Step 1: Mutserve2

echo "Mutserve2 calling for sample ${sample}"

wdir=${outd}/mutserve/work
testNCreate ${wdir}
bam_fol=$(dirname ${bam})
bam_file=$(basename ${bam})
ref_name=$(basename ${ref_mutserve})

singularity exec --no-home --no-mount ${PWD} \
    -B ${res_dir}:/reference \
    -B ${outd}/mutserve:/outdir \
    -B ${bam_fol}:/input \
    ${sif} \
    java -jar /opt/mutserve/mutserve.jar  \
        call \
        --level 0.01 \
        --reference /reference/${ref_name} \
        --mapQ 20 \
        --baseQ 20 \
        --output /outdir/work/${sample}.tmp.mutserve.vcf.gz \
        --no-ansi \
        --strand-bias 1.6 \
        --write-raw \
        /input/${bam_file} && \
        
        bcftools norm \
        -m-any \
        -f ${ref_mutserve} \
        -o ${wdir}/${sample}.tmp.norm.mutserve.vcf.gz -Oz \
        ${wdir}/${sample}.tmp.mutserve.vcf.gz && bcftools index -t ${wdir}/${sample}.tmp.norm.mutserve.vcf.gz &&

        #source activate /orfeo/cephfs/home/burlo/nardone/.conda/envs/lr_process

        echo "${bam_file} ${sample}" > ${wdir}/sample_renaming.tsv &&
        bcftools reheader -s ${wdir}/sample_renaming.tsv ${wdir}/${sample}.tmp.norm.mutserve.vcf.gz | bcftools view -Oz -o ${outd}/mutserve/${sample}.norm.mutserve.vcf.gz && \
        bcftools index -t ${outd}/mutserve/${sample}.norm.mutserve.vcf.gz


testOutput ${outd}/mutserve/${sample}.norm.mutserve.vcf.gz "Mutserve Calling" ${sample}


#if [ -f ${outd}/mutserve/${sample}.norm.mutserve.vcf.gz ];then
#    echo "Mutserve calling for ${sample} completed!"
#else
#    echo "Something went wrong during Mutserve calling for sample ${sample}."
#    exit 1
#fi
    
# Step 2: Mutect2. We will need singularity for this

wdir=${outd}/mutect2/work
testNCreate $wdir
bam_fol=$(dirname ${bam})
bam_file=$(basename ${bam})
ref_name=$(basename ${ref_mutect2})
echo "Mutect2 calling for sample ${sample}"

singularity exec --no-home --no-mount ${PWD} \
    -B ${res_dir}:/reference \
    -B ${outd}/mutect2:/outdir \
    -B ${bam_fol}:/input \
    ${sif} \
    gatk Mutect2 \
        -R /reference/${ref_name} \
        -L chrM \
        --min-base-quality-score 20 \
        -callable-depth 6 \
        --native-pair-hmm-threads ${threads} \
        --max-reads-per-alignment-start 0 \
        --tmp-dir /outdir/work \
        -I /input/${bam_file} \
        -O /outdir/work/${sample}.mutect2.raw.vcf.gz && \

singularity exec --no-home --no-mount ${PWD} \
    -B ${wdir}:/work \
    -B ${res_dir}:/reference \
    ${sif} \
    gatk FilterMutectCalls \
        -R /reference/${ref_name} \
        --min-reads-per-strand 2 \
        -V /work/${sample}.mutect2.raw.vcf.gz \
        --tmp-dir /work \
        -O /work/${sample}.mutect2.raw.filtered.vcf.gz && \
bcftools norm \
        -m-any \
        -f ${ref_mutect2} \
        -o ${outd}/mutect2/${sample}.norm.mutect2.vcf.gz -Oz \
        ${wdir}/${sample}.mutect2.raw.filtered.vcf.gz && bcftools index -t ${outd}/mutect2/${sample}.norm.mutect2.vcf.gz 

testOutput ${outd}/mutect2/${sample}.norm.mutect2.vcf.gz "Mutect2 calling" ${sample}

# Step 3: Filtering variants
filter_fold=${outd}/filtering
wdir=${filter_fold}/work
testNCreate ${wdir}

for method in mutserve mutect2;do

vcf_file=${outd}/${method}/${sample}.norm.${method}.vcf.gz

echo -e "ID\tFilter\tPos\tRef\tVariant\tVariantLevel\tMeanBaseQuality\tCoverage\tGT" \
        > ${wdir}/${sample}.${method}.txt

bcftools query -u \
        -f "${sample}\t%FILTER\t%POS\t%REF\t%ALT\t[%AF\t%BQ\t%DP\t%GT]\n" \
        ${vcf_file} >> ${wdir}/${sample}.${method}.txt
 
 if [[ ${method} == "mutserve" ]];then
    awk -F'\t' 'NR == 1 || (length($4) == 1 && length($5) == 1)' ${wdir}/${sample}.${method}.txt > ${wdir}/${sample}.${method}.tmp.txt
 else
    awk -F'\t' 'NR == 1 || ((length($4) > 1 || length($5) > 1) && length($4) != length($5))' ${wdir}/${sample}.${method}.txt > ${wdir}/${sample}.${method}.tmp.txt
 fi

awk 'BEGIN {OFS="\t"} {
        if (NR == 1) { print $0, "Type"; next }
        if ((length($4) > 1 || length($5) > 1) && length($4) != length($5)) { $10="3" }
        else if ($9 == "1") { $10="1" }
        else if ($9 == "0/1" || $9 == "1/0" || $9 == "0|1" || $9 == "1|0") { $10="2" }
        else { $10="UNKNOWN" }
        print
    }' ${wdir}/${sample}.${method}.tmp.txt > ${filter_fold}/${sample}.${method}.filtered.txt

done

for method in mutserve mutect2;do
testOutput ${filter_fold}/${sample}.${method}.filtered.txt "Filtering ${method} variants" ${sample} 
done

# Step 4 merging variants
source activate /orfeo/cephfs/home/burlo/nardone/miniconda3/envs/csvtk

wdir=${outd}/merging/work
testNCreate ${wdir}

filtered=$(ls ${filter_fold}/${sample}.*.filtered.txt | tr "\n" " ")

csvtk concat \
        -t ${filtered} \
        -T -o ${wdir}/${sample}.variants.concat.txt \
        --num-cpus ${threads} && \
csvtk sort \
        -t ${wdir}/${sample}.variants.concat.txt \
        -k ID:N -k Pos:n -k Ref:N -k Type:nr  -k Variant:N \
        -T -o ${wdir}/${sample}.variants.sorted.txt \
        --num-cpus ${threads} && \
singularity exec \
    -B ${outd}/merging:/base \
    ${sif} \
    java -jar /opt/VariantMerger.jar \
    /base/work/${sample}.variants.sorted.txt \
    --output /base/${sample}.mtDNA.variants.txt

testOutput ${outd}/merging/${sample}.mtDNA.variants.txt "Variant merging" ${sample}

# Step 5: Annotation


testNCreate ${outd}/annotation

java -jar ${mutserve} \
     annotate \
    --input ${outd}/merging/${sample}.mtDNA.variants.txt \
    --output ${outd}/annotation/${sample}.mtDNA.variants.annotated.txt \
    --annotation ${annotation_file}


testOutput ${outd}/annotation/${sample}.mtDNA.variants.annotated.txt "Variant Annotation" ${sample}

# Final Step - VCF merging

wdir=${outd}/vcf_merge/work
testNCreate ${wdir}

mut_vcf=${outd}/mutserve/${sample}.norm.mutserve.vcf.gz
mtec_vcf=${outd}/mutect2/${sample}.norm.mutect2.vcf.gz
merged_file=${outd}/merging/${sample}.mtDNA.variants.txt

bcftools view -v indels ${mtec_vcf} -Oz -o ${wdir}/${sample}.mutect2.indels.vcf.gz && \
    bcftools index -t  ${wdir}/${sample}.mutect2.indels.vcf.gz

while read line;do
    echo -e "chrM\t${line}" >> ${wdir}/pos_file.tsv
done < <(tail -n+2 ${merged_file} | cut -f3)

bcftools concat -a ${mut_vcf} ${wdir}/${sample}.mutect2.indels.vcf.gz | \
    bcftools sort -Oz -o ${wdir}/${sample}.SNPs_INDELs.vcf.gz && bcftools index -t ${wdir}/${sample}.SNPs_INDELs.vcf.gz &&
    bcftools view -R ${wdir}/pos_file.tsv ${wdir}/${sample}.SNPs_INDELs.vcf.gz -Oz -o ${outd}/vcf_merge/${sample}.filt.SNPs.INDELs.vcf.gz && \
    bcftools index -t ${outd}/vcf_merge/${sample}.filt.SNPs.INDELs.vcf.gz

testOutput ${outd}/vcf_merge/${sample}.filt.SNPs.INDELs.vcf.gz "VCF merging" ${sample}

