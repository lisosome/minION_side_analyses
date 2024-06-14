#!/usr/bin/env bash

#SBATCH -p THIN
#SBATCH --mem-per-cpu=10gb
#SBATCH --cpus-per-task 1
#SBATCH --ntasks-per-node=1
#SBATCH --time=96:00:00
set -e

printHelp(){
    echo "minION side analyses"
    echo 
    echo "Usage: minION_analyses.sh [options]"
    echo "Options:"
    echo "-h, --help - Show this help."
    echo "-r, --raw_data - Directory storing minION raw data"
    echo "-s, --samplesheet - Samplesheet specifying sample in the first column and barcode in the second. The file must be space separated"
    echo "-o, --outdir - Output Directory"
    echo "-m, --mode - Variant calling mode. Choose from [MT, STRC]"
    
}

if [[ $# -eq 0 ]];then 
    printHelp
    exit 0
fi

testNCreate(){
    fol=$1
    if [ ! -d $fol ];then
        mkdir -p $fol
    fi
}

while [[ $# -gt 0 ]];do
    arg=$1
    case ${arg} in
        -h|--help)
        printHelp
        exit 0
        ;;

        -r|--raw_data)
        raw_data=$2
        shift
        ;;

        -s|--samplesheet)
        samplesheet=$2
        shift
        ;;
        
        -o|--outdir)
        outdir=$2
        shift
        ;;
        
        -m|--mode)
        mode=$2
        shift
        ;;
        
        *)
        echo "Unknown option: $1"
        printHelp
        exit 0
    esac
    shift

done

logdir=${outdir}/Logs

testNCreate ${logdir}

SINGULARITY_TMPDIR=/fast/burlo/${USER}/singularity/sing_temp
SINGULARITY_CACHEDIR=/fast/burlo/${USER}/singularity/cache

testNCreate ${SINGULARITY_TMPDIR}
testNCreate ${SINGULARITY_CACHEDIR}
export SINGULARITY_TMPDIR
export SINGULARITY_CACHEDIR


BASE="$( cd "$( dirname "${BASH_SOURCE[0]}" 2>/dev/null)" && pwd 2>/dev/null)"
DIR=${BASE}
echo $DIR
preproc_dir=${DIR}/src/preprocessing
mtdna=${DIR}/src/mtdna
strc=${DIR}/src/strc


echo "1. Basecalling"

basecall=${preproc_dir}/guppy_cpu_SLURM.sh

basecalling=$(sbatch -A burlo --parsable -o ${logdir}/preprocessing.log -e ${logdir}/preprocessing.err ${basecall} ${raw_data} ${outdir} BARCODE) 
echo "Basecalling JOB ID: ${basecalling}"

echo "2. TRIMMING"

trim=${preproc_dir}/LongRead_preprocessing_SLURM.sh

trimming=$(sbatch -A burlo --parsable --dependency=afterok:${basecalling} -o ${logdir}/trimming.log -e ${logdir}/trimming.err ${trim} ${outdir} ${samplesheet})
echo "Trimming JOB ID: ${trimming}"

case $mode in
    MT )

    echo "3.ALIGNMENT for mtDNA analyses"
    align=${preproc_dir}/minimap2_SLURM.sh
    ref="/fast/burlo/nardone/mtDNA_resources/rcrs_mutserve.fasta"
    alignment=$(sbatch -A burlo --parsable --dependency=afterok:${trimming} -o ${logdir}/mt_DNA_alignment.log -e ${logdir}/mt_DNA_alignment.err ${align} ${outdir}/2.TRIMMING ${outdir}/3.ALIGNMENT ${ref} ${samplesheet})
    echo "mtDNA alignment JOB ID: ${alignment}"

    echo "4.VARIANT.CALLING for mtDNA"
    call=${mtdna}/mtDNA_variantCall.sh
    if [[ ! -d ${out} ]];then
        mkdir -p ${outdir}/4.VARIANT.CALLING
    fi

    #declare -A codes=()
    #while read line;do
    #  sample=$(echo "$line" | cut -d " " -f1)
    #  bar=$(echo "$line" | cut -d " " -f2)
    #  codes[${bar}]=${sample}
    #done < ${samplesheet}

    while read line;do
    
    sample=$(echo "$line" | cut -d " " -f1)
    
    calling=$(sbatch -A burlo --parsable --dependency=afterok:${alignment} -o ${logdir}/${sample}_mt_DNA_variant_calling.log -e ${logdir}/${sample}_mt_DNA__variant_calling.err ${call} -b ${outdir}/3.ALIGNMENT/${sample}.bam -s ${sample} -o ${outdir}/4.VARIANT.CALLING -t 5)
    echo "mtDNA variant calling for sample ${sample} JOB ID: ${calling}"
    
    done < ${samplesheet}
    ;;
    STRC )

    echo "3.ALIGNMENT for STRC analyses"
    align=${preproc_dir}/minimap2_SLURM.sh
    ref="/fast/burlo/nardone/STRC_ref/STRC_hg38.fasta"
    alignment=$(sbatch --parsable --dependency=afterok:${trimming} -o ${logdir}/STRC_alignment.log -e ${logdir}/STRC_alignment.err ${align} ${outdir}/2.TRIMMING ${outdir}/3.ALIGNMENT ${ref} ${samplesheet})
    echo "STRC alignment JOB ID: ${alignment}"

    testNCreate ${outdir}/4.VARIANT.CALLING

    echo "4.VARIANT.CALLING for STRC analyses"
    call=${strc}/clair3_STRC.sh
    
    #declare -A codes=()
    #while read line;do
    #  sample=$(echo "$line" | cut -d " " -f1)
    #  bar=$(echo "$line" | cut -d " " -f2)
    #  codes[${bar}]=${sample}
    #done < ${samplesheet}

    
    while read line;do

    sample=$(echo "$line" | cut -d " " -f1)

    
    calling=$(sbatch -A burlo --parsable --dependency=afterok:${alignment} -o ${logdir}/${sample}_STRC_variant_calling.log -e ${logdir}/${sample}_STRC_variant_calling.err ${call} ${outdir}/3.ALIGNMENT ${sample} ${ref} ${outdir}/4.VARIANT.CALLING)
    echo "STRC variant calling for sample ${sample} JOB ID: ${calling}"

    done < ${samplesheet}
    ;;
esac