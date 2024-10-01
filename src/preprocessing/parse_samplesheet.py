#!/usr/bin/env python3

# Script to parse the samplesheet and rename the demultiplexed bam files. Eventually, merge bam files from same samples but with different barcodes

import sys
import os
import argparse


#samplesheet = sys.argv[1]
#demultiplex_folder = sys.argv[2]
#kit_name = sys.argv[3]


def get_args():
        parser = argparse.ArgumentParser(
        description = "%(prog)s: Parse the samplesheet and performs all sorts of marvellous things!"
        )
        parser.add_argument('-s', '--samplesheet', type=str, help='Space separated samplesheet with three columns: DNA-code, barcode_number and type of analysis')
        parser.add_argument('-d', '--demul', type=str, help='Name of the folder in which the dorado created BAM has been demultiplexed')
        parser.add_argument('-k', '--kit_name', type=str, help='Name of the kit used for the minION sequencing')
        args = parser.parse_args()
        return args


def parse_file(filename):
    # Initialize an empty dictionary to store the data
    data_dict = {}

    # Open the file and read it line by line
    with open(filename, 'r') as file:
        for line in file:
            # Split each line into sample ID and barcode number
            sample_id, barcode, _ = line.strip().split()

            # If the sample ID is already in the dictionary, append the barcode to the list
            if sample_id in data_dict:
                data_dict[sample_id].append(barcode)
            # If the sample ID is not in the dictionary, create a new entry
            else:
                data_dict[sample_id] = [barcode]

    return data_dict

def barcode_constructor(barcode):
    for bar in barcode:
        if int(bar) < 10:
            bar = f"barcode0{bar}"
        else:
            bar = f"barcode{bar}"
    return bar


def main():
    args = get_args()
    samplesheet = args.samplesheet
    demultiplex_folder = args.demul
    kit_name = args.kit_name
    sample_dict = parse_file(samplesheet)
    for sample, barcode in sample_dict.items():
        outdir = os.path.join(demultiplex_folder, sample)
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        bar_name = barcode_constructor(barcode)
        if len(barcode) == 1:
            demultiplexed = f"{kit_name}_{bar_name}"
            demul = os.path.join(demultiplex_folder, demultiplexed)
            cmd = f"samtools index -@ 15 {demul}.bam && mv {demul}.bam {outdir}/{sample}.bam && mv {demul}.bam.bai {outdir}/{sample}.bam.bai"
            print(cmd)
            os.system(cmd)
        else:
            barcode_to_merge = [barcode_constructor(x) for x in barcode]
            demultiplexed = [os.path.join(demultiplex_folder, f"{kit_name}_{x}.bam") for x in barcode_to_merge]
            cmd = f"""samtools merge -@ 15 {" ".join(demultiplexed)} -o {outdir}/{sample}.bam && samtools index -@ 15 {outdir}/{sample}.bam"""
            print(cmd)
            os.system(cmd)
            if os.path.exists(os.path.join(outdir, f"{sample}.bam")):
                for f in demultiplexed:
                    cmd = f"rm {f}*"
                    os.system(cmd)
            else:
                print("Something went wrong with the merging. Exiting...")
                sys.exit(1)
    

if __name__ == '__main__':
    main()


    