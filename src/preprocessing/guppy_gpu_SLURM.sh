#!/usr/bin/env bash

# This script will perform the guppy basecalling using GPU nodes through the SLURM queue manager.
# This code is inspired by https://github.com/colindaven/guppy_on_slurm/blob/master/runbatch_gpu_guppy.sh

# set partition
#SBATCH -p GPU

# set run on x MB node only
#SBATCH --mem=192gb

# set run on bigmem node only
#SBATCH --cpus-per-task 48

# task per nodes
#SBATCH --ntasks-per-node=1

# set max wallclock time
#SBATCH --time=18:00:00

# set name of job
#SBATCH --job-name=guppy

# mail alert at start, end and abortion of execution
#SBATCH --mail-type=ALL
#SBATCH --mail-user=giuseppegiovanni.nardone@burlo.trieste.it


rawdata=$1 # Directory in which are located all the files produced by minION
baseout=$2 # Directory to use for the analyses
mode=$3

case $mode in
  SINGLE_SAMPLE)
module load ont-guppy-gpu/6.5.7
#guppy=/opt/area/shared/programs/x86_64/ont-guppy-gpu/6.2.1/bin/guppy_basecaller
config=/orfeo/opt/programs/x86/fedora36/ont-guppy-cpu/6.5.7/data/dna_r10.4.1_e8.2_400bps_5khz_sup.cfg

echo "Linking fast5 files in one directory..."
mkdir -p ${baseout}/0.RAW.DATA ${baseout}/1.BASECALLING
find ${rawdata} -type f -name "*.pod5" -exec ln -s {} ${baseout}/0.RAW.DATA/ \;

echo "Executing basecalling..."

guppy_params="--compress_fastq --recursive --num_callers 12 --gpu_runners_per_device 4 --chunks_per_runner 512 --chunk_size 3000"

guppy_basecaller -i ${baseout}/0.RAW.DATA -s ${baseout}/1.BASECALLING -c ${config} -x 'auto' ${guppy_params}
;;
BARCODE)
  kit=SQK-NBD114-24

 module load ont-guppy-gpu/6.5.7

config=/orfeo/opt/programs/x86/fedora36/ont-guppy-cpu/6.5.7/data/dna_r10.4.1_e8.2_400bps_5khz_sup.cfg

echo "Linking fast5 files in one directory..."
mkdir -p ${baseout}/0.RAW.DATA ${baseout}/1.BASECALLING
find ${rawdata} -type f -name "*.pod5" -exec ln -s {} ${baseout}/0.RAW.DATA/ \;

echo "Executing basecalling and demultiplexing..."

guppy_params="--compress_fastq --recursive --num_callers 12 --gpu_runners_per_device 12 --chunks_per_runner 208 --chunk_size 2000 --chunks_per_caller 10000 --barcode_kits ${kit} --enable_trim_barcodes --detect_barcodes --detect_mid_strand_barcodes"

guppy_basecaller -i ${baseout}/0.RAW.DATA -s ${baseout}/1.BASECALLING -c ${config} -x 'auto' ${guppy_params}
;;
esac

