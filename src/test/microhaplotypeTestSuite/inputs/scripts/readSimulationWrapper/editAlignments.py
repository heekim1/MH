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

def main():

    # Getting in the commandline arguments
    parser = argparse.ArgumentParser(description="Script to edit alignments file to optimize for different scenarios",
                                     epilog='EXAMPLE: python editAlignments.py -b input.bam -i 60 -a 75 -s 400 -o output.bam')
    parser.add_argument('-b', '--bam', help='Input BAM file', dest='ibam', required=True, type=str)
    parser.add_argument('-p', '--pos', help='The position to be edited in chr:pos format', dest='pos')
    parser.add_argument('-o', '--out', help='Output BAM file', dest='obam', type=str)
    cmdArguments = parser.parse_args()

    # Setting up the logger
    logging.basicConfig(format='%(asctime)s %(levelname)s : %(message)s', level=logging.INFO)

    # Set input file path to os.path.abspath
    if os.path.isfile(cmdArguments.ibam):
        cmdArguments.ibam = os.path.abspath(cmdArguments.ibam)
    else:
        sys.exit(logging.error("The given bam file at %s does not exist. Exiting" % cmdArguments.ibam))

    # Set output file name if none given
    if cmdArguments.obam is None:
        cmdArguments.obam = cmdArguments.ibam.split('/')[-1].split('.')[0] + "_adjusted_align.bam"

    # Opening the input file
    bamInput = pysam.AlignmentFile(cmdArguments.ibam, "rb")
    # Opening the output file
    bamOutput = pysam.AlignmentFile(cmdArguments.obam, 'wb', template=bamInput)

    # Breaking down the specified positions
    alignmentRemovalPositionsList = cmdArguments.pos.split(',')

    rogueReads = []

    # Iterate through the file and find the rogue reads
    for alignmentRemovalPosition in alignmentRemovalPositionsList:
        chrom, pos = alignmentRemovalPosition.split(':')
        for pileupColumn in bamInput.pileup(chrom, int(pos)-1, int(pos), max_depth=800000):
            for pileupread in pileupColumn.pileups:
                rogueReads.append(pileupread.alignment.query_name)

    # Get only the unique rogue read names
    rogueReads = set(rogueReads)

    # Do not output the rogue reads
    for read in bamInput.fetch():
        if read.query_name not in rogueReads:
            bamOutput.write(read)

    # Close the bam files
    bamInput.close()
    bamOutput.close()

if __name__ == '__main__':
    main()
