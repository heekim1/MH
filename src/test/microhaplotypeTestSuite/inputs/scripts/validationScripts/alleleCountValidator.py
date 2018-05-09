#!/usr/bin/env python
from __future__ import division
__author__ = 'pushkardakle'

import pysam
import argparse
import pybedtools
import logging
import os
import sys
import math


def main():
    # Getting in the commandline arguments
    parser = argparse.ArgumentParser(description="Validate Read count")
    parser.add_argument('-b', '--bam', dest='bamFile', required=True)
    parser.add_argument('-e', '--bed', dest='bedFile', required=True)
    parser.add_argument('-i', '--info', dest='infoFile')
    parser.add_argument('-o', '--out', dest='outFileName', required=True)
    parser.add_argument('-m', '--minCov', dest='minCov', default=10)
    cmdArguments = parser.parse_args()

    # Settings for mapping quality and deletion filter
    mappingQualThresholdLower = 60
    mappingQualThresholdUpper = 254
    absoluteAnalyticalThreshold = 10
    deletionRemoval = True

    # Initializing the logger
    logging.basicConfig(format='%(asctime)s %(levelname)s : %(message)s', level=logging.INFO)

    # Validating the input files
    validateInputFiles(cmdArguments)

    # Process bed file
    bedIntervals = processBed(cmdArguments.bedFile)

    # Get microhaplotype counts at intervals
    shortReadsCount, readsCount, microhaplotypeCount = getMicrohaplotypeCounts(cmdArguments.bamFile, bedIntervals)

    # Process the Ae and probOfDetectMixture values if the info file in given in the arguments
    microhaplotypeInfoDict = {}
    if cmdArguments.infoFile is not None:
        microhaplotypeInfoDict = processInfoFile(cmdArguments.infoFile)

    # Now lets apply the filters coverage, average mapping quality and deletion filters and find the final alleles
    appendedMicroHapInfoDict, minimumNumberOfContributors = filterMicrohaplotypeCalls(microhaplotypeCount,
                                                                                      microhaplotypeInfoDict,
                                                                                      cmdArguments.minCov,
                                                                                      mappingQualThresholdLower,
                                                                                      mappingQualThresholdUpper,
                                                                                      deletionRemoval,
                                                                                      absoluteAnalyticalThreshold)

    # Creating output file and printing results
    outputFile = open(cmdArguments.outFileName, 'w')

    # Printing the headers
    outputFile.write('##Deletion Filter:%s, Mapping Quality Threshold:%d, Minimum Coverage:%f\n' % (str(deletionRemoval),
                     mappingQualThresholdLower, cmdArguments.minCov))
    outputFile.write('#The minimum number of contributor : %d\n' % minimumNumberOfContributors)
    outputFile.write('MarkerId	AvgAe\tProbOfDetMix\tNumOfAlleles\tNumOfContributors\tAlleleCoverage'
                     '\tAlleleCoverageDetail(-strand:+strand:avgMQ)\tSNPs\tMinimumCoverage\n')

    # Printing results for individual microhaplotypes
    for microhapName in microhaplotypeCount.keys():
        outputString = microhapName + '\t'

        # If we have values for Ae or probDetMixture from info file only then print them
        if 'Ae' in appendedMicroHapInfoDict[microhapName]:
            outputString += appendedMicroHapInfoDict[microhapName]['Ae'] + '\t'
        else:
            outputString += '\t'

        if 'probDetMixture' in appendedMicroHapInfoDict[microhapName]:
            outputString += appendedMicroHapInfoDict[microhapName]['probDetMixture'] + '\t'
        else:
            outputString += '\t'

        outputString += str(appendedMicroHapInfoDict[microhapName]['NumOfAlleles']) + '\t' \
                        + str(appendedMicroHapInfoDict[microhapName]['NumOfContributors']) \
                        + '\t' + appendedMicroHapInfoDict[microhapName]['AlleleCoverage']+'\t' + \
                        appendedMicroHapInfoDict[microhapName][r'AlleleCoverageDetail(-strand:+strand:avgMQ)'] + '\t'

        if 'rsidList' in appendedMicroHapInfoDict[microhapName]:
            outputString += appendedMicroHapInfoDict[microhapName]['rsidList'] + '\t'
        else:
            outputString += '\t'

        outputString += str(appendedMicroHapInfoDict[microhapName]['MinimumCoverage']) + '\n'

        outputFile.write(outputString)

    outputFile.close()

    # Printing the stdout results
    print "##readCount : %d, ShortReadCount : %d" % (readsCount, shortReadsCount)
    print "#MarkerId\tallele\tcoverage\tminus_cov\tplus_cov\tavg_MQ"
    for microhapName in microhaplotypeCount.keys():
        for indAlleles in microhaplotypeCount[microhapName]:
            print "%s\t%s\t%d\t%d\t%d\t%d" % (microhapName, indAlleles, microhaplotypeCount[microhapName][indAlleles]['count'],
                                              microhaplotypeCount[microhapName][indAlleles]['strandBiasNeg'],
                                              microhaplotypeCount[microhapName][indAlleles]['strandBiasPos'],
                                              microhaplotypeCount[microhapName][indAlleles]['avgMQ'])


def validateInputFiles(lCmdArguments):

    if not os.path.isfile(lCmdArguments.bamFile):
        logging.error("The bam file at %s was not found. Exiting" % lCmdArguments.bamFile)
        sys.exit(1)
    else:
        lCmdArguments.bam = os.path.abspath(lCmdArguments.bamFile)

    if not os.path.isfile(lCmdArguments.bedFile):
        logging.error("The bed file at %s was not found. Exiting" % lCmdArguments.bedFile)
        sys.exit(1)
    else:
        lCmdArguments.bam = os.path.abspath(lCmdArguments.bedFile)

    if lCmdArguments.infoFile is not None:
        if not os.path.isfile(lCmdArguments.infoFile):
            logging.error("The info file at %s was not found. Exiting" % lCmdArguments.infoFile)
            sys.exit(1)
        else:
            lCmdArguments.bam = os.path.abspath(lCmdArguments.infoFile)

    try:
        lCmdArguments.minCov = float(lCmdArguments.minCov)
    except ValueError:
        logging.error("The minimum coverage value %s is not a float. Exiting" % lCmdArguments.minCov)
        sys.exit(1)

    if lCmdArguments.minCov < 0.0 or lCmdArguments.minCov > 1.0:
        logging.error("The minimum coverage value %f isnt between 0 and 1" % lCmdArguments.minCov)
        sys.exit(1)


def processBed(lBedFile):
    """
    Process the bed file and return a dictionary of the main interval as key and actual snps as values
    :param lBedFile: The location of the bed file
    :return: A dictionary of main interval as key and actual snps as values
    """
    rBedDict = {}
    bedFileStringRep = ''

    # Handling the header in bam file.
    # Current Handling:- If it has track at the start the bed entry would be skipped.
    # @TODO - Make the bed header checking more generic rather than just track
    with open(lBedFile,'rU') as bedFileHandler:
        for singleBedLine in bedFileHandler:
            if not singleBedLine.startswith('track'):
                bedFileStringRep += singleBedLine

    bedFileObj = pybedtools.BedTool(bedFileStringRep, from_string=True)  # Giving string input to bedtools

    # Group the bed intervals by the microhaplotype name and create a new temporary grouped bed file
    groupedBedFile = bedFileObj.groupby(g=8, c=[1, 2], o=['distinct', 'distinct'])

    # Read in the grouped bed file and create a dict which has entries as follows
    # Key mh11KK-191,chr11,99880162,99880352': Value [99880162, 99880223, 99880281, 99880351]
    # The key gives the region which we need to look into
    # The values gives positions of snps in that region
    with open(groupedBedFile.fn) as f:
        for bedLine in f:
            bedLine = bedLine.rstrip('\n')
            lineContents = bedLine.split('\t')
            intervalPositions = map(int, lineContents[2].split(','))
            rBedDict[lineContents[0] + "," + lineContents[1] + "," + str(min(intervalPositions)) +
                     "," + str(max(intervalPositions) + 1)] = intervalPositions

    return rBedDict


def getMicrohaplotypeCounts(lBamFile, lBedIntervals):

    # @TODO: Add function documentation and flow comments

    rMicrohaplotypeCountDict = {}
    rshortReadsCount = 0
    rReadsCount = 0
    bamFile = pysam.AlignmentFile(lBamFile, 'rb')
    for indInterval in lBedIntervals.keys():
        microhapName, microhapChrom, microhapStart, microhapEnd = indInterval.split(',')
        microhapPositions = lBedIntervals[indInterval]
        alleleDict, strandDict, mapqDict = ({},{},{})
        for pileupColumn in bamFile.pileup(microhapChrom, int(microhapStart), int(microhapEnd), max_depth=800000):
            if pileupColumn.pos in microhapPositions:
                for pileupread in pileupColumn.pileups:
                    if pileupread.query_position is None:
                        alleleDict[pileupread.alignment.query_name] = alleleDict.get(
                            pileupread.alignment.query_name, '') + "_"
                    else:
                        alleleDict[pileupread.alignment.query_name] = alleleDict.get(
                            pileupread.alignment.query_name, '') + pileupread.alignment.query_sequence[
                            pileupread.query_position]
                    strandDict[pileupread.alignment.query_name] = pileupread.alignment.flag
                    mapqDict[pileupread.alignment.query_name] = pileupread.alignment.mapping_quality
        intervalDict = {}
        for alleleCall in alleleDict.keys():
            if len(alleleDict[alleleCall]) == len(microhapPositions):
                intervalDict.setdefault(alleleDict[alleleCall], {})
                intervalDict[alleleDict[alleleCall]]['count'] = intervalDict[alleleDict[alleleCall]].get('count', 0)+1
                intervalDict[alleleDict[alleleCall]].setdefault('strandBiasNeg', 0)
                intervalDict[alleleDict[alleleCall]].setdefault('strandBiasPos', 0)
                if strandDict[alleleCall] == 0:
                    intervalDict[alleleDict[alleleCall]]['strandBiasPos'] = \
                        intervalDict[alleleDict[alleleCall]].get('strandBiasPos', 0)+1
                elif strandDict[alleleCall] == 16:
                    intervalDict[alleleDict[alleleCall]]['strandBiasNeg'] = \
                        intervalDict[alleleDict[alleleCall]].get('strandBiasNeg', 0)+1
                intervalDict[alleleDict[alleleCall]].setdefault('avgMQ', [])
                intervalDict[alleleDict[alleleCall]]['avgMQ'].append(mapqDict[alleleCall])
                rReadsCount += 1
            else:
                rshortReadsCount += 1



        for intervalAllele in intervalDict.keys():
            intervalDict[intervalAllele]['avgMQ'] = int(sum(intervalDict[intervalAllele]['avgMQ'])/
                                                        len(intervalDict[intervalAllele]['avgMQ']))

        rMicrohaplotypeCountDict[microhapName] = intervalDict
    bamFile.close()

    return rshortReadsCount, rReadsCount, rMicrohaplotypeCountDict


def processInfoFile(lInfoFile):

    """
    This function will take in the microhaplotype info file and make dicts for ae and probablity of detection of mixture
    :param lInfoFile: The info file
    :return: Dicts for the Ae and probability of detecting mixture
    """
    # Return dicts initialization
    rMicrohaplotypeInfoDict = {}

    # Read in the info file
    with open(lInfoFile, 'rU') as infoFile:
        for infoLine in infoFile:
            infoLine = infoLine.rstrip('\n')
            infoList = infoLine.split('\t')
            microhaplotypeName = infoList[1]
            probDetMixture = infoList[2]
            aeValue = infoList[3]
            rsidList = infoList[8]
            rMicrohaplotypeInfoDict.setdefault(microhaplotypeName, {})
            rMicrohaplotypeInfoDict[microhaplotypeName]['Ae'] = aeValue
            rMicrohaplotypeInfoDict[microhaplotypeName]['probDetMixture'] = probDetMixture
            rMicrohaplotypeInfoDict[microhaplotypeName]['rsidList'] = rsidList

    return rMicrohaplotypeInfoDict


def filterMicrohaplotypeCalls(lMicrohaplotypeCountDict, lMicrohaplotypeInfoDict, lMinCov,
                              lMappingQualThresholdLower, lMappingQualThresholdUpper, lDeletionRemovalFlag,
                              lAbsoluteAnalyticalThreshold):
    """
    This function will apply the filters for mincov, mapping quality and deletion. Also linerize the allele calls
    :param lMicrohaplotypeCountDict: The dict for all the allele calls with their counts,
    strand values and mapping quality values
    :param lMicrohaplotypeInfoDict: The dict for all info like Ae, probDetMixture, Number of Alleles, Number of
    Contributors and final allele calls
    :param lMinCov: The coverage cutoff
    :param lMappingQualThresholdLower: The lower bound for mapping qual
    :param lMappingQualThresholdUpper: The upper bound for mapping qual
    :param lDeletionRemovalFlag: The deletion removal flag
    :param lAbsoluteAnalyticalThreshold: The absolute analytical threshold
    :return: The info dict with information added
    """

    # Now lets iterate over all these keys and see what we have in the counts dictionary
    # Adding support for minimum number of contributors
    minimum_number_of_contributors_list = [0]
    for microhaplotypeName in lMicrohaplotypeCountDict.keys():
        lMicrohaplotypeInfoDict.setdefault(microhaplotypeName, {})
        lMicrohaplotypeInfoDict[microhaplotypeName]['AlleleCoverage'], \
        lMicrohaplotypeInfoDict[microhaplotypeName][r'AlleleCoverageDetail(-strand:+strand:avgMQ)'], \
        lMicrohaplotypeInfoDict[microhaplotypeName]['MinimumCoverage'] = \
            _getAlleleCoverage(lMicrohaplotypeCountDict[microhaplotypeName], lMinCov, lMappingQualThresholdLower,
                               lMappingQualThresholdUpper, lDeletionRemovalFlag, lAbsoluteAnalyticalThreshold)

        if lMicrohaplotypeInfoDict[microhaplotypeName]['AlleleCoverage'] == '':
            lMicrohaplotypeInfoDict[microhaplotypeName]['NumOfAlleles'] = 0
        else:
            lMicrohaplotypeInfoDict[microhaplotypeName]['NumOfAlleles'] = len(lMicrohaplotypeInfoDict[microhaplotypeName]['AlleleCoverage'].split(','))

        if lMicrohaplotypeInfoDict[microhaplotypeName]['NumOfAlleles'] == 0:
            lMicrohaplotypeInfoDict[microhaplotypeName]['NumOfContributors'] = 0
        else:
            lMicrohaplotypeInfoDict[microhaplotypeName]['NumOfContributors'] = int(math.ceil(lMicrohaplotypeInfoDict[microhaplotypeName]['NumOfAlleles']/2))
            minimum_number_of_contributors_list.append(int(math.ceil(lMicrohaplotypeInfoDict[microhaplotypeName]['NumOfAlleles']/2)))

    return lMicrohaplotypeInfoDict, max(minimum_number_of_contributors_list)


def _getAlleleCoverage(lMicrohaplotypeCount, lMinCov, lMappingQualThresholdLower,
                       lMappingQualThresholdUpper, lDeletionRemovalFlag, lAbsoluteAnalyticalThreshold):
    """
    Look into a dict of all the called alleles for the particular microhaplotype. Depending on deletion flag include
    /remove it. Also put filters of avg mapping quality and coverage
    :param lMicrohaplotypeCount: A dict of all the called alleles for a particular microhaplotype
    :param lMinCov: The minimum coverage threshold
    :param lMappingQualThresholdLower: The lower bound for mapping qual
    :param lMappingQualThresholdUpper: The upper bound for mapping qual
    :param lDeletionRemovalFlag: The flag to indicate if deletions should be kept/removed
    :param lAbsoluteAnalyticalThreshold: The absolute analytical threshold
    :return: The allelic string representation
    """

    AlleleStringRepresentationList = []
    # Adding extended allelic representation
    AlleleStringRepresentationExtendedList = []
    # 05042016 - Pushkar - Adding analytical threshold
    rAnalyticalThreshold = int(lMinCov * sum([lMicrohaplotypeCount[allele]['count']
                                              for allele in lMicrohaplotypeCount.keys()]))

    for alleleCall in lMicrohaplotypeCount.keys():
        if lDeletionRemovalFlag and '_' in alleleCall:
            continue
        elif lMicrohaplotypeCount[alleleCall]['count'] < rAnalyticalThreshold or \
                        lMicrohaplotypeCount[alleleCall]['count'] < lAbsoluteAnalyticalThreshold:
            continue
        elif lMicrohaplotypeCount[alleleCall]['avgMQ'] < lMappingQualThresholdLower or \
                        lMicrohaplotypeCount[alleleCall]['avgMQ'] > lMappingQualThresholdUpper:
            continue
        else:
            AlleleStringRepresentationList.append('[' + alleleCall + ']' +
                                                  str(lMicrohaplotypeCount[alleleCall]['count']))
            AlleleStringRepresentationExtendedList.append('[' + alleleCall + ']' +
                                                           str(lMicrohaplotypeCount[alleleCall]['count']) +
                                                           "(" +
                                                           str(lMicrohaplotypeCount[alleleCall]['strandBiasNeg']) +
                                                           ":" +
                                                           str(lMicrohaplotypeCount[alleleCall]['strandBiasPos']) +
                                                           ":" +
                                                           str(lMicrohaplotypeCount[alleleCall]['avgMQ']) + ")")

    rAlleleStringRepresentation = ','.join(AlleleStringRepresentationList)
    rAlleleStringRepresentationExtended = ','.join(AlleleStringRepresentationExtendedList)

    return rAlleleStringRepresentation, rAlleleStringRepresentationExtended, rAnalyticalThreshold


if __name__ == '__main__':
    main()