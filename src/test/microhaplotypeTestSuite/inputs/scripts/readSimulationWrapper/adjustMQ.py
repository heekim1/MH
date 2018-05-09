#!/usr/bin/env python
__author__ = "pushkardakle"

"""This script will adjust the mapping quality values in a BAM file. This is required since we might have cutoffs
for mapping quality and using the currently generated fastq files with bwa/tmap is giving very low mapping quality
values. This script will give a random mapping quality between the given minimum and maximum
"""

import pysam
import argparse
import os
import logging
import sys
import random

def main():

    # Getting in the commandline arguments
    parser = argparse.ArgumentParser(description="Script to edit BAM file to optimize for different scenarios",
                                     epilog='EXAMPLE: python editBAM.py -b input.bam -i 60 -a 75 -s 400 -o output.bam')
    parser.add_argument('-b', '--bam', help='Input BAM file', dest='ibam', required=True, type=str)
    parser.add_argument('-i', '--min', help='Minimum Mapping Quality', dest='min', default=70, type=int)
    parser.add_argument('-a', '--max', help='Max Mapping Quality', dest='max', default=80, type=int)
    parser.add_argument('-s', '--seed', help='Seed for reproducibility', dest='seed', type=int)
    parser.add_argument('-o', '--out', help='Output BAM file', dest='obam', type=str)
    cmdArguments = parser.parse_args()

    # Setting up the logger
    logging.basicConfig(format='%(asctime)s %(levelname)s : %(message)s', level=logging.INFO)

    # Set input file path to os.path.abspath
    if os.path.isfile(cmdArguments.ibam):
        cmdArguments.ibam = os.path.abspath(cmdArguments.ibam)
    else:
        sys.exit(logging.error("The given bam file at %s does not exists. Exiting" % cmdArguments.ibam))

    # Set output file name if none given
    if cmdArguments.obam is None:
        cmdArguments.obam = cmdArguments.ibam.split('/')[-1].split('.')[0] + "_adjusted_mapq.bam"

    # Set the random seed if specified
    if cmdArguments.seed is not None:
        random.seed(cmdArguments.seed)

    # Opening the input file
    bamInput = pysam.AlignmentFile(cmdArguments.ibam, "rb")
    # Opening the output file
    bamOutput = pysam.AlignmentFile(cmdArguments.obam, 'wb', template=bamInput)

    # Iterate through the file and choose a random mapping quality
    for read in bamInput.fetch():
        if read.mapping_quality > 0:
            read.mapping_quality = random.randint(cmdArguments.min, cmdArguments.max)
        bamOutput.write(read)

    # Close the bam files
    bamInput.close()
    bamOutput.close()

if __name__ == '__main__':
    main()
