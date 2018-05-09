#!/usr/bin/env python
__author__ = "pushkardakle"

"""This script will simulate NGS ion torrent reads according to given ratio for mixture analysis
It will do so by generating the number of samples specified in the ratio at varying coverages and then mixing them
It can run without any options by taking all the defaults from the config.ini file
For the options that you can configure at run time please run the script with the --help/-h argument
"""

import vcf
import ConfigParser
import argparse
import logging
import csv
import subprocess
import shlex
import sys
import os
# from datetime import datetime
from copy import deepcopy
import pybedtools
from collections import OrderedDict
import random

def main():

    # Getting in the commandline arguments
    parser = argparse.ArgumentParser(description="Reads simulator for mixture samples",
                                     epilog='EXAMPLE: python readsSimulator.py -f reference.fasta -b regions.bed '
                                            '-m mutations.csv -r 1:3:1 -c config.ini -s 1000 -j samplesim'
                                            ' -a postSimulationBWA.sh')
    parser.add_argument('-f', '--fasta', dest='referenceFile')
    parser.add_argument('-b', '--bed', dest='bedFile')
    parser.add_argument('-m', '--mutation', dest='mutationsFile')
    parser.add_argument('-i', '--hotspots', dest='hotspotsFile')
    parser.add_argument('-r', '--ratio', dest='mixtureRatio')
    parser.add_argument('-c', '--config', dest='configFile')
    parser.add_argument('-a', '--postSimulationScript', dest='postSimulationScript')
    parser.add_argument('-s', '--seed', dest='seed', default='')
    parser.add_argument('-j', '--jobname', dest='jobName')
    cmdArguments = parser.parse_args()

    # Setting up the logger
    logging.basicConfig(format='%(asctime)s %(levelname)s : %(message)s', level=logging.INFO)

    # Getting the directory where the script is located
    scriptDir = sys.path[0]

    # Getting the config file path
    # If config not given by user we will use the default config is script directory
    if cmdArguments.configFile is None:
        cmdArguments.configFile = scriptDir + '/config.ini'

    # Check if file exists
    if not os.path.isfile(cmdArguments.configFile):
        sys.exit(logging.error("Config file at %s not found. Exiting" % cmdArguments.configFile))

    # Reading in the config file
    config = ConfigParser.ConfigParser()
    config.read(cmdArguments.configFile)

    # Setting up the default options
    cmdArguments = set_defaults(cmdArguments, config)

    # Reading in the other options from the config(Sim options are added after checking with bed file)
    configVcfTemplate = os.path.abspath(config.get('defaultFiles', 'vcfTemplate'))
    configBaseCoverage = config.get('global_options', 'base_coverage')
    configSimulatorPath = os.path.abspath(config.get('global_options', 'simulator_path'))

    # Validating the arguments
    validate_params(cmdArguments, configVcfTemplate, configSimulatorPath)

    # Logging the final mutations file
    logging.info('Using the mutations file at %s' % cmdArguments.mutationsFile)

    # Adding in functionality to check if first read length is less than the smallest interval
    # If it is not we will set it to be less that the smallest interval
    if config.has_option('simulator', 'first_read_length'):
        config.set('simulator', 'first_read_length', checkFirstReadLength(cmdArguments.bedFile,
                                                                          config.get('simulator',
                                                                                     'first_read_length')))

    # After validating simulator options line read length the sim options are added in
    configSimOptions = config.get('simulator', 'sim_options')


    # Setting up the working directories
    # All the temporary files will go in a directory named as the jobname
    # Only the final fastq and bam processing files will be outside
    currentDir = os.getcwd()
    tempDir = currentDir + "/" + cmdArguments.jobName
    if not os.path.exists(tempDir):
        os.makedirs(tempDir)
    os.chdir(tempDir)

    # Making the vcf files required by simulator
    vcfFilesPathDict = makevcf(cmdArguments.mutationsFile, configVcfTemplate,
                               cmdArguments.mixtureRatio, cmdArguments.jobName)

    # Generating fastq files for individual contributors
    fastqFilesPathDict = makefastq(configSimulatorPath, configSimOptions, cmdArguments.mixtureRatio,
                                   cmdArguments.referenceFile, cmdArguments.bedFile, configBaseCoverage,
                                   vcfFilesPathDict, cmdArguments.seed, cmdArguments.jobName)

    # Changing the working dir back to original dir
    os.chdir(currentDir)

    # Merging the generated fastq files for individual contributors
    mergedFastqFile = mergefastq(fastqFilesPathDict, cmdArguments.jobName)
    logging.info("Find final output fastq at %s" % mergedFastqFile)

    # Print counts for each contributor
    printContributorCounts(cmdArguments.mixtureRatio, mergedFastqFile)

    # If mapping script specified execute it
    if cmdArguments.postSimulationScript is not None:
        logging.info('Calling post simulation script at %s' % cmdArguments.postSimulationScript)
        logging.info('**************** Post Simulation Script Output ****************')
        mapReads(cmdArguments.postSimulationScript, cmdArguments.referenceFile, mergedFastqFile,
                 cmdArguments.jobName, scriptDir)

    logging.info("Complete")


def set_defaults(argparseObj, configObj):
    """
    Check if values for parameters have been set in argparse else sets them according to config file
    :param argparseObj: An argparse object
    :param configObj: A config object
    :return: Modified argparse object
    """

    if argparseObj.mixtureRatio is None:
        argparseObj.mixtureRatio = configObj.get('global_options', 'default_mixture_ratio')

    if argparseObj.referenceFile is None:
        argparseObj.referenceFile = os.path.abspath(configObj.get('defaultFiles', 'defaultReference'))

    if argparseObj.bedFile is None:
        argparseObj.bedFile = os.path.abspath(configObj.get('defaultFiles', 'defaultBed'))


    if argparseObj.jobName is None:
        # argparseObj.jobname = "sim_reads_" + datetime.utcnow().strftime('%Y_%m_%d_%H_%M_%S_%f')[:-4]
        argparseObj.jobName = "sim_reads"

    if argparseObj.seed is not None:
        random.seed(argparseObj.seed)

    # Adding in a check for the hotspots file
    if configObj.getboolean('global_options', 'useMutationsFile') and argparseObj.mutationsFile is None:
        argparseObj.mutationsFile = os.path.abspath(configObj.get('defaultFiles', 'defaultMutations'))
    else:
        if argparseObj.hotspotsFile is None:
            argparseObj.hotspotsFile = os.path.abspath(configObj.get('defaultFiles', 'defaultHotspots'))

        if os.path.isfile(argparseObj.hotspotsFile):
            argparseObj.mutationsFile = _getMutationsFile(argparseObj.hotspotsFile, argparseObj.jobName)
        else:
            logging.error("The specified hotspots file at %s does not exist. Exiting" % argparseObj.hotspotsFile)
            sys.exit(0)

    return argparseObj


def _getMutationsFile(lHotSpotsFile, lJobName):
    """
    Create a mutations file from the given hotspots file
    :param lHotSpotsFile: Path to the hotspots file
    :param lJobName: Job Name
    :return: The path to created mutations file
    """

    # Number of contributors for which mutations should be generated
    numberOfContributor = 8

        # Getting in the bed file and a dictionary
    bedDict = _processBed(lHotSpotsFile)

    # Opening the output file
    outputFile = open(lJobName + '_generated_mutations.csv', 'w')

    # Writing the header
    outputFile.write('Chromosome,Position,Reference Allele,Alternate Alleles,Allele Frequency\n')

    # Initializing the allele and frequency possibility
    allelePossibilities = ['A', 'T', 'G', 'C', '-']
    frequencyPossibilities = ['0.5', '1']

    # Iterating over the keys to create the mutations csv file
    for chromosomePos in bedDict.keys():
        chromosome, position = chromosomePos.split(":")
        referenceBase = bedDict[chromosomePos].split(';')[0].split('=')[1]

        # Initialize the mutationsString and the frequenciesString
        mutationsString = ''
        frequenciesString = ''

        while len(mutationsString)<numberOfContributor:
            currentMutation = random.choice(allelePossibilities)
            if currentMutation != referenceBase:
                mutationsString += currentMutation

        for iterator in range(0, numberOfContributor):
            frequenciesString += random.choice(frequencyPossibilities) + '/'

        mutationsString = '/'.join(mutationsString)
        frequenciesString = frequenciesString.rstrip('/')

        outputFile.write('%s,%s,%s,%s,%s\n' % (chromosome, position, referenceBase, mutationsString, frequenciesString))

    outputFile.close()

    return os.path.abspath(lJobName + '_generated_mutations.csv')


def _processBed(lBedFile):
    """
    Process the bed file and return a dictionary of the main interval as key and actual snps as values
    :param lBedFile: The location of the bed file
    :return: A dictionary of main interval as key and actual snps as values
    """
    rBedDict = OrderedDict()
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
    groupedBedFile = bedFileObj.groupby(g=[8,1,3], c=7, o='first')

    # Read in the grouped bed file and create a dict which has entries as follows
    # Key mh11KK-191,chr11,99880162 : Value 'REF=A,ALT=C,ANCHOR=T'
    # The key gives the region which we need to look into
    # The values gives reference in that region
    with open(groupedBedFile.fn) as f:
        for bedLine in f:
            bedLine = bedLine.rstrip('\n')
            lineContents = bedLine.split('\t')
            rBedDict[lineContents[1]+ ':' + lineContents[2]] = lineContents[3]

    return rBedDict


def validate_params(argparseObj, lVcfTemplate, lSimulatorPath):
    """
    Validate the input parameters
    :param argparseObj: An argparse object
    :param lVcfTemplate: VCF template path
    :param lSimulatorPath: Simulator executable path
    :return: None
    """

    if not os.path.isfile(argparseObj.referenceFile):
        sys.exit(logging.error("Reference file at %s not found. Exiting" % argparseObj.referenceFile))

    if not os.path.isfile(argparseObj.bedFile):
        sys.exit(logging.error("Bed file at %s not found. Exiting" % argparseObj.bedFile))

    if not os.path.isfile(argparseObj.mutationsFile):
        sys.exit(logging.error("Mutations file at %s not found. Exiting" % argparseObj.mutationsFile))

    if not os.path.isfile(lVcfTemplate):
        sys.exit(logging.error("Vcf template file at %s not found. Exiting" % lVcfTemplate))

    if not os.path.isfile(lSimulatorPath):
        sys.exit(logging.error("Simulator file at %s not found. Exiting" % lSimulatorPath))


def checkFirstReadLength(lBedFilePath, lSpecifiedReadLength):
    """
    Check and find the lowest interval size in the bed file and set the read length to it after displaying a warning
    :param lBedFilePath: Path to the interval bed file
    :return: The read length
    """
    rRevisedReadLength = lSpecifiedReadLength
    intervalLengthsArr = []
    readLengthSubtractAmount = 2
    with open(lBedFilePath,'rU') as bedFileHandler:
        for singleBedLine in bedFileHandler:
            if not singleBedLine.startswith('track'):
                singleBedLineList = singleBedLine.split()
                try:
                    intervalLengthsArr.append(int(singleBedLineList[2]) - int(singleBedLineList[1]))
                except IndexError:
                    logging.warning('Bed file incorrectly formatted at %s. Skipping this entry' % singleBedLine)

    if int(lSpecifiedReadLength) > int(min(intervalLengthsArr) - readLengthSubtractAmount):
        rRevisedReadLength = min(intervalLengthsArr) - readLengthSubtractAmount
        logging.warning('The specified read length %s is larger than the size of smallest interval. '
                        'Resetting it to %d' % (lSpecifiedReadLength, rRevisedReadLength))

    return rRevisedReadLength


def makevcf(lMutatationsFile, lVcfTemplate, lMixtureRatio, lJobName):
    """
    Create vcf file according to values in the mutations file
    The header will be picked up from the template vcf file
    :param lMutatationsFile: Mutations csv file
    :param lVcfTemplate: Template vcf file for header
    :param lMixtureRatio: Mixture ratio
    :param lJobName: Jobname
    :return: A dict with paths of vcf files created
    """
    rVcfFilePathDict = {}
    contributorList = lMixtureRatio.split(":")
    currentContributor = 0
    NoAlleleSpecifiedChar = '-'

    # Reading in the vcf template file for header and a dummy record
    vcfTemplateObj = vcf.Reader(open(lVcfTemplate,'r'))
    dummyRecord = vcfTemplateObj.next()

    # Create a vcf file for each contributor
    for indContributor in contributorList:
        contributorStrRep = 'contributor' + str((currentContributor + 1))
        vcfWriteObj = vcf.VCFWriter(open(lJobName + "_" + contributorStrRep + ".vcf", 'w'), vcfTemplateObj)
        logging.info('Creating vcf file for %s' % contributorStrRep)

        with open(lMutatationsFile, 'rU') as csvFileName:
            csvFile = csv.DictReader(csvFileName, dialect=csv.excel)

            for indRow in csvFile:
                alternateAlleles = indRow['Alternate Alleles'].split('/')
                alleleFreq = indRow['Allele Frequency'].split('/')
                currentAlternateAllele, currentAlleleFreq = None, None
                try:
                    currentAlternateAllele = alternateAlleles[currentContributor]
                    currentAlleleFreq = alleleFreq[currentContributor]
                except IndexError:
                    logging.warning('No alternate allele or freq given for contributor %d at site %s:%s',
                                    contributorStrRep, indRow['Chromosome'], indRow['Position'])
                else:
                    currentRecord = deepcopy(dummyRecord)
                    if currentAlternateAllele != NoAlleleSpecifiedChar and currentAlleleFreq != NoAlleleSpecifiedChar:
                        currentRecord.CHROM = indRow['Chromosome']
                        currentRecord.POS = indRow['Position']
                        currentRecord.REF = indRow['Reference Allele']
                        currentRecord.ALT = currentAlternateAllele.split(',')  # pyvcf takes a list for ALT
                        currentRecord.INFO['AF'] = currentAlleleFreq
                        if float(currentAlleleFreq) >= 1.0:  # if HOM then pl=3 else pl=default 2 from template
                            currentRecord.INFO['pl'] = 3
                        # Adding support for mutation type i.e. SUBSTITUTE for SNP INSERT/DELETE for ins/del
                        # specifying currentRecord.ALT[0] as we have list as alternate allele while we want
                        # to check the length of the first alternate allele. Hard coding 0 in here as i dont think
                        # we will have scenario of specifying multiple alternate alleles for one contributor
                        if len(currentRecord.REF) == len(currentRecord.ALT[0]):
                            currentRecord.INFO['mt'] = 'SUBSTITUTE'
                        elif len(currentRecord.REF) > len(currentRecord.ALT[0]):
                            currentRecord.INFO['mt'] = 'DELETE'
                        elif len(currentRecord.REF) < len(currentRecord.ALT[0]):
                            currentRecord.INFO['mt'] = 'INSERT'
                        else:
                            logging.warning('Unknown mutation type observed for contributor %d at site %s:%s',
                                    contributorStrRep, indRow['Chromosome'], indRow['Position'])
                        vcfWriteObj.write_record(currentRecord)

        vcfWriteObj.close()
        rVcfFilePathDict[contributorStrRep] = os.path.abspath(lJobName + "_" + contributorStrRep + ".vcf")
        currentContributor += 1

    return rVcfFilePathDict


def makefastq(lSimulatorPath, lSimulatorOptions, lMixtureRatio,
              lReferenceFile, lBedFile, lBaseCoverage, lVcfFilesPathDict, lSeed, lJobName):
    """
    Run the simulator with given reference, bed and vcf files and create fastq files for each contributor
    :param lSimulatorPath: Path for simulator
    :param lSimulatorOptions: Base options to simulator
    :param lMixtureRatio: Mixture Ratio
    :param lReferenceFile: Reference File
    :param lBedFile: Bed file
    :param lBaseCoverage: Base coverage
    :param lVcfFilesPathDict: Dict from makevcf
    :param lSeed: Seed for random number generation
    :param lJobName: Job Name
    :return: A dict with paths of created fastq files
    """

    rFastqFilesPathDict = {}
    contributorList = lMixtureRatio.split(":")
    currentContributor = 0

    for indContributor in contributorList:
        indContributorCoverage = int(int(lBaseCoverage) * float(indContributor))
        simulatorSeed = ' -z ' + str(lSeed) if lSeed is not '' else lSeed
        contributorStrRep = 'contributor' + str((currentContributor + 1))

        logging.info('Coverage value for contributor %d is %d' % ((currentContributor+1), indContributorCoverage))
        logging.info('Creating fastq file for %s' % contributorStrRep)
        simulatorCommand = lSimulatorPath \
                            + " " + lSimulatorOptions \
                            + " -C " + str(indContributorCoverage) \
                            + " -P " + contributorStrRep \
                            + " -x " + lBedFile \
                            + simulatorSeed \
                            + " -v " + lVcfFilesPathDict[contributorStrRep] \
                            + " " + lReferenceFile \
                            + " " + lJobName + "_" + contributorStrRep
        logging.info("Executing \"%s\" command for %s" % (simulatorCommand, contributorStrRep))
        logging.info('**************** Simulator output for %s ****************' % contributorStrRep)
        simulatorCommandSplit = shlex.split(simulatorCommand)
        subprocess.check_output(simulatorCommandSplit)
        rFastqFilesPathDict[contributorStrRep] = os.path.abspath(lJobName + "_" + contributorStrRep + ".bfast.fastq")
        currentContributor += 1

    return rFastqFilesPathDict


def mergefastq(lFastqFilesPathDict, lJobName):
    """
    Merge the contributor fastq files into a single fastq file
    :param lFastqFilesPathDict: A dict of fastq file paths to each contributor
    :param lJobName: Job Name
    :return: Path for final merged fastq file
    """

    with open(lJobName + '_merged.fastq', 'w') as outFile:
        for inFileName in lFastqFilesPathDict.values():
            with open(inFileName) as inFile:
                outFile.write(inFile.read())
    return os.path.abspath(lJobName + '_merged.fastq')


def mapReads(lPostSimulationScript, lReferenceFile, lFinalFastqFile, lJobName, lScriptDir):
    """
    Call mapping script with the reference file and generated fastq file
    :param lMappingScript: Mapping script to call
    :param lReferenceFile: The reference fasta file
    :param lFinalFastqFile: The final generated fastq file
    :return: None
    """
    postSimulationCommand = "sh " + lPostSimulationScript + " " + lReferenceFile + " " \
                     + lFinalFastqFile + " " + lJobName + " " + lScriptDir
    postSimulationCommandSplit = shlex.split(postSimulationCommand)
    subprocess.check_output(postSimulationCommandSplit)


def printContributorCounts(lMixtureRatio, lFinalFastqFile):
    """
    Print the generated read counts for each contributor
    :param lMixtureRatio: Mixture Ratio
    :param lFinalFastqFile: Path for final merged fastq file
    :return: None
    """
    contributorList = lMixtureRatio.split(":")
    currentContributor = 0

    for indContributor in contributorList:
        contributorStrRep = 'contributor' + str((currentContributor + 1))
        contributorCount = subprocess.check_output('grep -c ' + contributorStrRep + " " + lFinalFastqFile, shell=True)
        contributorCount = contributorCount.rstrip('\n')
        logging.info("Number of total generated reads for %s is %s" % (contributorStrRep, contributorCount))
        currentContributor += 1


if __name__ == '__main__':
    main()


