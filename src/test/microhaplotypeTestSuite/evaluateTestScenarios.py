#!/usr/bin/env python
__author__ = "pushkardakle"

"""This script will compare the outputs from the Microhaplotyper.jar file and the validation scripts
"""

# Importing system modules
import argparse
import logging
import os
import sys
import shlex
import subprocess
from datetime import datetime



def main():

    # Getting in the commandline arguments
    parser = argparse.ArgumentParser(description="Script to compare microhaplotyper and validation results",
                                     epilog='EXAMPLE: python evaluateTestScenarios.py -d /path/to/scenarioDir'
                                            ' -j /path/to/Microhaplotyper.jar -o output.csv')
    parser.add_argument('-d', '--dir', help='The dir with scenarios', dest='dir', type=str, required=True)
    parser.add_argument('-j', '--jar', help='Path to Microhaplotyper.jar. Default: Jar if exists in current dir',
                        dest='jar', type=str)
    parser.add_argument('-o', '--out', help='Output file name. Default comparison_output_timestamp.csv',
                        dest='out', type=str)
    cmdArguments = parser.parse_args()

    # Initializing the logger
    logging.basicConfig(format='%(asctime)s %(levelname)s : %(message)s', level=logging.INFO)

    # Remove any trailing / from dir
    cmdArguments.dir = cmdArguments.dir.rstrip('/')

    # Validate the cmd arguments
    cmdArguments = validateCmdArgumentsAndSetDefaults(cmdArguments)

    # Creating the output directory
    if not os.path.exists("./" + os.path.splitext(cmdArguments.out)[0] + '_jar'):
        os.mkdir("./" + os.path.splitext(cmdArguments.out)[0] + '_jar')

    # Create the output file
    comparisonOutputFileHandler = open(cmdArguments.out, 'w')
    comparisonOutputFileHandler.write('Scenario,Scenario Description,Scenario Result,Result Details\n')

    # Navigate the test scenario directory
    for scenarioDirName in next(os.walk(cmdArguments.dir))[1]:
        if os.path.isfile(cmdArguments.dir + '/' + scenarioDirName + '/' + scenarioDirName + "_readme.txt"):
            currentScenarioDir = cmdArguments.dir + '/' + scenarioDirName
            readmeDict = {}
            readmeFile = open(currentScenarioDir + '/' + scenarioDirName + "_readme.txt", 'rU')
            for readmeLine in readmeFile:
                readmeLine = readmeLine.rstrip('\n')
                readmeField, readmeValue = readmeLine.split(':')
                readmeDict[readmeField] = readmeValue

            # Set the output file name
            jarOutputName = "./" + os.path.splitext(cmdArguments.out)[0] + '_jar' + "/" + scenarioDirName \
                            + '_microhaplotyper.txt'

            truthFileOutputName = currentScenarioDir + '/' + readmeDict['Truth file']

            # Execute the microhaplotyper.jar file on current scenario files
            microhaplotyperJarCommand = 'java -jar ' + cmdArguments.jar \
                                        + ' -bam ' + currentScenarioDir + '/' + readmeDict['Bam File'] \
                                        + ' -bed ' + currentScenarioDir + '/' + readmeDict['Bed File'] \
                                        + ' -minCov ' + readmeDict['Minimum Coverage'] \
                                        + ' -o ' + jarOutputName

            # If the info file is given add it too into the executed command
            if readmeDict['Info File'] != "-":
                microhaplotyperJarCommand += ' -info ' + currentScenarioDir + '/' + readmeDict['Info File']

            # Execute the created commandline
            logging.info('Generating jar output for scenario %s' % scenarioDirName)
            logging.info('Commandline to generate the jar file %s' % microhaplotyperJarCommand)
            microhaplotyperJarCommandSplit = shlex.split(microhaplotyperJarCommand)
            subprocess.check_output(microhaplotyperJarCommandSplit)

            # Compare truth files and microhaplotyper jar output
            comparisonCheckFlag, comparisonCheckDict = compareMicrohaplotyperWithTruth(jarOutputName,
                                                                                       truthFileOutputName)
            comparisonCheckString  = ''
            for failCategories in comparisonCheckDict:
                comparisonCheckString += failCategories + ":" + comparisonCheckDict[failCategories].rstrip(',') + "|"

            comparisonOutputFileHandler.write('%s,%s,%s,\"%s\"\n' % (readmeDict['Scenario Name'], readmeDict['Scenario Description'], comparisonCheckFlag, comparisonCheckString))

        else:
            logging.warning("No readme file found for scenario %s. Skipping it" % scenarioDirName)
    comparisonOutputFileHandler.close()


def compareMicrohaplotyperWithTruth(lJarOutputFilePath, lTruthFilePath):
    """
    Compare the outputs from the validation script and the microhaplotyper jar file
    :param lJarOutputFilePath: The jar output file
    :param lTruthFilePath: The truth file
    :return: Comparison Flag(True/False), Category wise dict of all failures if any
    """

    rComparisonCheckFlag = 'PASS'
    rComparisonCheckDict = {}
    jarDict, jarHeadersList = resultsFileToDict(lJarOutputFilePath)
    truthDict, truthHeadersList = resultsFileToDict(lTruthFilePath)

    for markerName in truthDict.keys():
        if markerName in jarDict.keys():
            for headerElement in truthHeadersList:
                if headerElement == 'AlleleCoverage' and \
                                set(jarDict[markerName][headerElement].split(',')) != \
                                set(truthDict[markerName][headerElement].split(',')):

                    rComparisonCheckFlag = 'FAIL'
                    rComparisonCheckDict[headerElement] = rComparisonCheckDict.get(headerElement, '') \
                                                          + 'Non matching %s/%s/%s for marker %s,' % \
                                                            (headerElement, str(jarDict[markerName][headerElement]),
                                                          str(truthDict[markerName][headerElement]), markerName)
                elif headerElement != 'AlleleCoverage' and abs(float(jarDict[markerName][headerElement]) -
                                                                       float(truthDict[markerName][headerElement])) > 0:
                    rComparisonCheckFlag = 'FAIL'
                    rComparisonCheckDict[headerElement] = rComparisonCheckDict.get(headerElement, '') \
                                                          + 'Non matching %s/%s/%s for marker %s,' % \
                                                            (headerElement, str(jarDict[markerName][headerElement]),
                                                          str(truthDict[markerName][headerElement]), markerName)
        else:
            rComparisonCheckFlag = 'FAIL'
            rComparisonCheckDict['Markers'] = rComparisonCheckDict.get('Markers', '') + 'No entry for marker %s,' % markerName

    return rComparisonCheckFlag, rComparisonCheckDict


def resultsFileToDict(lInputResultsFile):
    rResultsDict = {}
    resultFileHandler = open(lInputResultsFile, 'rU')
    # Ignoring the line for filters
    resultsFilterLine = resultFileHandler.readline()

    # Getting the line for number of contributors
    numberOfContributorsLine = resultFileHandler.readline()

    # Getting in the header line
    resultsHeaderList = resultFileHandler.readline().rstrip('\n').split('\t')
    resultsHeaderList.pop(0)

    # Now iterating over the rest of the lines
    for resultLine in resultFileHandler:
        resultLine = resultLine.rstrip('\n')
        resultLineList = resultLine.split('\t')
        markerName = resultLineList.pop(0)
        rResultsDict[markerName] = dict(zip(resultsHeaderList, resultLineList))

    return rResultsDict, resultsHeaderList


def validateCmdArgumentsAndSetDefaults(lCmdArguments):
    """
    Validate the input arguments and set default values
    :param lCmdArguments: The argparse object
    :return: The edited argparse object
    """

    if not os.path.exists(lCmdArguments.dir):
        logging.error("The given scenarios dir at %s does not exist. Exiting" % lCmdArguments.dir)
        sys.exit(1)

    if lCmdArguments.jar is None:
        lCmdArguments.jar = os.getcwd() + '/Microhaplotyper.jar'

    if not os.path.isfile(lCmdArguments.jar):
        logging.error('The jar file at %s does not exist. Exiting' % lCmdArguments.jar)
        sys.exit(1)

    if lCmdArguments.out is None:
        # lCmdArguments.out = "comparison_output_" + datetime.utcnow().strftime('%Y_%m_%d_%H_%M_%S_%f')[:-4]
        lCmdArguments.out = "comparison_output.csv"
    return lCmdArguments


if __name__ == '__main__':
    main()