#!/usr/bin/env bash

#SBATCH -p THIN
#SBATCH --mem-per-cpu=10gb
#SBATCH --cpus-per-task 1
#SBATCH --ntasks-per-node=1
#SBATCH --time=96:00:00
set -e

#Function to print the helper 
printHelp(){ 
    echo "minION side analyses"
    echo 
    echo "Usage: minION_analyses.sh [options]"
    echo "Options:"
    echo "-h, --help - Show this help."
    echo "-r, --raw_data - Directory storing minION raw data"
    echo "-s, --samplesheet - Samplesheet specifying sample in the first column and barcode in the second and mode in the third. Choose from [MT, STRC]. The file must be space separated"
    echo "-o, --outdir - Output Directory"
    
}

if [[ $# -eq 0 ]];then # If the number of arguments is equal to 0, print the helper and quit
    printHelp
    exit 0
fi

# Function to test the existance of a folder and eventually create it
testNCreate(){
    fol=$1
    if [ ! -d $fol ];then
        mkdir -p $fol
    fi
}


# Definition of the flags/arguments of the script. While the number of arguments is greather than 0, define the variable $args.
# Case the $args variable in several options. Every option is specified as a $2 positional argument due to the shift command
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
        
        *)
        echo "Unknown option: $1"
        printHelp
        exit 0
    esac
    shift

done

# Define some variables, create useful folders and export singularity env variables to specify the cache and temp folder

logdir=${outdir}/Logs

testNCreate ${logdir}

SINGULARITY_TMPDIR=/fast/burlo/${USER}/singularity/sing_temp
SINGULARITY_CACHEDIR=/fast/burlo/${USER}/singularity/cache

testNCreate ${SINGULARITY_TMPDIR}
testNCreate ${SINGULARITY_CACHEDIR}
export SINGULARITY_TMPDIR
export SINGULARITY_CACHEDIR

# Define the script directory to specify the absolute path of the other scripts
BASE="$( cd "$( dirname "${BASH_SOURCE[0]}" 2>/dev/null)" && pwd 2>/dev/null)"
DIR=${BASE}
preproc_dir=${DIR}/src/preprocessing
mtdna=${DIR}/src/mtdna
strc=${DIR}/src/strc

# First we will perform a check to see what kind of analyses will be performed based on the samplesheet
# The samplesheet structure is like:
# 23-1234 4 MT
# 24-5678 5 STRC

# We will just check using grep like this, counting the record resulting from grep

mt_analysis=$(grep "MT" ${samplesheet} | wc -l)
strc_analysis=$(grep "STRC" ${samplesheet} | wc -l)

testNCreate ${outdir}/samplesheets
# If the number stored in the variable is greater than 0 than create the specific samplesheet
if [[ ${mt_analysis} -gt 0 ]];then
    grep "MT" ${samplesheet} | cut -d " " -f-2 > ${outdir}/samplesheets/MT_samplesheet
fi
# If the number stored in the variable is greater than 0 than create the specific samplesheet
if [[ ${strc_analysis} -gt 0 ]];then
    grep "STRC" ${samplesheet} | cut -d " " -f-2 > ${outdir}/samplesheets/STRC_samplesheet
fi

# Basecalling, trimming and demultiplexing step. Dorado will perform all this.
# This step also uses a personalized python script to parse the samplesheet, rename the BAM files created with the demultiplexing or eventually merge bam belonging to the same sample but specified by different barcodes
basecalling=$(sbatch -J basecalling -A burlo --parsable -o ${logdir}/preprocessing.log -e ${logdir}/preprocessing.err ${preproc_dir}/dorado.sh ${raw_data} ${outdir} ${samplesheet})
echo "Preprocessing JOB ID: ${basecalling}"

align=${preproc_dir}/dorado_aligner.sh

# Central part of the script. Launching jobs for alignment and variant calling. Without translating literally the code, the following operations will be executed:
#   1. Checking for the existance of the MT and/or STRC samplesheet to launch the correct analyses with the right samples
#   2. Perform a check to see if there are samples with the same DNA code but with different barcode:
#       First we will select the DNA codes from the samplesheet and then we will leverage the sort and the uniq command to count the duplicated records. THe count will be stored as the firs column. And then we use awk to print out only the records that are duplicated and count them.
#       If the count of the duplicates is greather than 1, then another samplesheet will be created, containing the result of the `uniq -c` command, with DNA codes repeated only once
#   3. For every line (sample) in the samplesheet (being this the one created at the beginning or the one with de-duplicated samples) launch the alignment job and the appropriate variant calling job  


if [[ -f ${outdir}/samplesheets/MT_samplesheet ]];then
    mt_dups=$(cut -d " " -f1 ${outdir}/samplesheets/MT_samplesheet | sort | uniq -c | awk '($1 > 1){print $0}' | wc -l)
    if [[ ${mt_dups} -gt 0 ]];then
        cut -d " " -f1 ${outdir}/samplesheets/MT_samplesheet | sort | uniq -c > ${outdir}/samplesheets/MT_samplesheet_uniq
        while read line;do
            sample=$(echo "$line" | cut -d " " -f2)
            alignment=$(sbatch -J mt_align_${sample} -A burlo --parsable --dependency=afterok:${basecalling} -o ${logdir}/${sample}_mt_DNA_alignment.log -e ${logdir}/${sample}_mt_DNA_alignment.err ${align} ${outdir} ${sample} MT)
            echo "Mithocondrial alignment for sample ${sample} JOB ID: ${alignment}"
            call=${mtdna}/mtDNA_variantCall.sh
            if [[ ! -d ${outdir}/3.VARIANT.CALLING ]];then
                mkdir -p ${outdir}/3.VARIANT.CALLING
            fi
            calling=$(sbatch -A burlo --parsable --dependency=afterok:${alignment} -o ${logdir}/${sample}_mt_DNA_variant_calling.log -e ${logdir}/${sample}_mt_DNA_variant_calling.err ${call} -b ${outdir}/2.ALIGNMENT/${sample}/${sample}_MT.bam -s ${sample} -o ${outdir}/3.VARIANT.CALLING/${sample} -t 5)
            echo "mtDNA variant calling for sample ${sample} JOB ID: ${calling}"
        done < ${outdir}/samplesheets/MT_samplesheet_uniq
    else
        while read line;do
            sample=$(echo "$line" | cut -d " " -f1)
            alignment=$(sbatch -J mt_align_${sample} -A burlo --parsable --dependency=afterok:${basecalling} -o ${logdir}/${sample}_mt_DNA_alignment.log -e ${logdir}/${sample}_mt_DNA_alignment.err ${align} ${outdir} ${sample} MT)
            echo "Mithocondrial alignment for sample ${sample} JOB ID: ${alignment}"
            call=${mtdna}/mtDNA_variantCall.sh
            if [[ ! -d ${outdir}/3.VARIANT.CALLING ]];then
                mkdir -p ${outdir}/3.VARIANT.CALLING
            fi
            calling=$(sbatch -A burlo --parsable --dependency=afterok:${alignment} -o ${logdir}/${sample}_mt_DNA_variant_calling.log -e ${logdir}/${sample}_mt_DNA_variant_calling.err ${call} -b ${outdir}/2.ALIGNMENT/${sample}/${sample}_MT.bam -s ${sample} -o ${outdir}/3.VARIANT.CALLING/${sample} -t 5)
            echo "mtDNA variant calling for sample ${sample} JOB ID: ${calling}"
        done < ${outdir}/samplesheets/MT_samplesheet
    fi
fi


if [[ -f ${outdir}/samplesheets/STRC_samplesheet ]];then
    strc_dups=$(cut -d " " -f1 ${outdir}/samplesheets/STRC_samplesheet | sort | uniq -c | awk '($1 == 2){print $0}' | wc -l)
    if [[ ${strc_dups} -gt 0 ]];then
        cut -d " " -f1 ${outdir}/samplesheets/STRC_samplesheet | sort | uniq -c > ${outdir}/samplesheets/STRC_samplesheet_uniq
        while read line;do
            sample=$(echo "$line" | cut -d " " -f2)
            alignment=$(sbatch -J strc_align_${sample} -A burlo --parsable --dependency=afterok:${basecalling} -o ${logdir}/${sample}_STRC_alignment.log -e ${logdir}/${sample}_STRC_alignment.err ${align} ${outdir} ${sample} STRC)
            echo "STRC alignment for sample ${sample} JOB ID: ${alignment}"
            call=${strc}/clair3_STRC.sh
            if [[ ! -d ${outdir}/3.VARIANT.CALLING ]];then
                mkdir -p ${outdir}/3.VARIANT.CALLING
            fi
            ref="/fast/burlo/nardone/STRC_ref/STRC_hg38.fasta"
            calling=$(sbatch -A burlo --parsable --dependency=afterok:${alignment} -o ${logdir}/${sample}_STRC_variant_calling.log -e ${logdir}/${sample}_STRC_variant_calling.err ${call} ${outdir}/2.ALIGNMENT/${sample} ${sample} ${ref} ${outdir}/3.VARIANT.CALLING)
            echo "STRC variant calling for sample ${sample} JOB ID: ${calling}"

        done < ${outdir}/samplesheets/STRC_samplesheet_uniq
    else
        while read line;do
            sample=$(echo "$line" | cut -d " " -f1)
            alignment=$(sbatch -J strc_align_${sample} -A burlo --parsable --dependency=afterok:${basecalling} -o ${logdir}/${sample}_STRC_alignment.log -e ${logdir}/${sample}_STRC_alignment.err ${align} ${outdir} ${sample} STRC)
            echo "STRC alignment for sample ${sample} JOB ID: ${alignment}"
            call=${strc}/clair3_STRC.sh
            if [[ ! -d ${outdir}/3.VARIANT.CALLING ]];then
                mkdir -p ${outdir}/3.VARIANT.CALLING
            fi
            ref="/fast/burlo/nardone/STRC_ref/STRC_hg38.fasta"
            calling=$(sbatch -A burlo --parsable --dependency=afterok:${alignment} -o ${logdir}/${sample}_STRC_variant_calling.log -e ${logdir}/${sample}_STRC_variant_calling.err ${call} ${outdir}/2.ALIGNMENT/${sample} ${sample} ${ref} ${outdir}/3.VARIANT.CALLING)
            echo "STRC variant calling for sample ${sample} JOB ID: ${calling}"
        done < ${outdir}/samplesheets/STRC_samplesheet
    fi
fi
