#!/usr/bin/env python
__author__ = "pushkardakle"

"""This suite generation script will generate all the test scenarios for microhaplotype and also generate the expected
truth outputs for it. It will take a tsv file as input along with a base config file and generate a directory with
all the scenarios.
"""


# Import system modules
import argparse
import logging
import os
import sys
import shlex
import subprocess
import csv
import ConfigParser
import glob
import shutil
from datetime import datetime


def main():

    # Getting in the commandline arguments
    parser = argparse.ArgumentParser(description="Script to generate the test suite for microhaplotypes",
                                     epilog='EXAMPLE: python generateTestScenarios.py -s scenarios.csv'
                                            ' -c default_config.ini -d /home/input -o /home/output')
    parser.add_argument('-s', '--scenarios', help='The test scenarios tsv file. Default: baseScenarios.tsv in scenarios directory', dest='scenarios',
                        type=str)
    parser.add_argument('-c', '--config', help='Config file. Default: baseConfig.ini in configs directory',
                        dest='config', type=str)
    parser.add_argument('-d', '--dir', help='Base directory for input files. Default: input in Script Directory',
                        dest='dir', type=str)
    parser.add_argument('-o', '--outdir', help='Output directory. Default: Script Directory',
                        dest='outdir', type=str)
    cmdArguments = parser.parse_args()

    # Initializing the logger
    logging.basicConfig(format='%(asctime)s %(levelname)s : %(message)s', level=logging.INFO)

    # Configuring some variable names
    generatedTestScenariosScratchDirName = 'microhaplotypeTestScenariosScratch'
    generatedTestScenariosDirName = 'microhaplotypeTestScenarios'
    defaultConfigName = 'baseConfig.ini'
    defaultInputDirectoryName = 'inputs'
    defaultTestScenariosName = 'baseScenarios.tsv'

    # Getting the directory where the script is located
    scriptDir = sys.path[0]

    # Validating the input files
    cmdArguments = validateInputsAndSetDefaults(cmdArguments, scriptDir, defaultConfigName, defaultTestScenariosName,
                                                defaultInputDirectoryName)

    # Getting in the base working directory
    if cmdArguments.outdir is not None and os.path.exists(cmdArguments.outdir):
        baseWorkingDir = cmdArguments.outdir.rstrip('/')
    else:
        baseWorkingDir = os.getcwd()


    # Creating a directory for the generated test scenarios
    # @TODO - Check for directory permissions that need to be set up
    if not os.path.exists(baseWorkingDir + '/' + generatedTestScenariosScratchDirName):
        os.mkdir(baseWorkingDir + '/' + generatedTestScenariosScratchDirName)
    logging.info("Scratch directory for generated scenarios: %s" %
                 baseWorkingDir + '/' + generatedTestScenariosScratchDirName)

    # Also create a directory to store the final bam, bai, microhaplotype and readme file
    if not os.path.exists(baseWorkingDir + '/' + generatedTestScenariosDirName):
        os.mkdir(baseWorkingDir + '/' + generatedTestScenariosDirName)
    logging.info("Directory for generated scenarios: %s" %
                 baseWorkingDir + '/' + generatedTestScenariosDirName)

    # Set up the current working directory to the test scenarios directory created in earlier instruction
    os.chdir(baseWorkingDir + '/' + generatedTestScenariosScratchDirName)

    # Reading in the scenarios file row by row for each scenario
    # U added for portability
    with open(cmdArguments.scenarios, 'rU') as scenariosFileHandle:
        scenariosFile = csv.DictReader(scenariosFileHandle, delimiter='\t')  # This is a tsv file
        # Processing for each scenario
        for indvScenario in scenariosFile:

            if not indvScenario['scenarioName'].startswith('#'):

                # Create a directory for the scenario and navigate to it
                if not os.path.exists("./" + indvScenario['scenarioName']):
                    os.mkdir("./" + indvScenario['scenarioName'])

                os.chdir("./" + indvScenario['scenarioName'])

                # Creating an updated config file for current scenario
                updatedConfigPath = generateScenarioConfig(indvScenario, cmdArguments.config, cmdArguments.dir)

                logging.info("Running simulation or scenario: %s" % indvScenario['scenarioName'])
                logging.info("Scenario Description: %s" % indvScenario['scenarioDescription'])
                # Run the simulation wrapper for current scenario
                simulatedBamFilePath = executeReadSimulationWrapper(indvScenario, updatedConfigPath)

                # Copy the generated bam, bai to new dir and also write up a readme
                copyToFinalDir(baseWorkingDir + '/' + generatedTestScenariosDirName, simulatedBamFilePath,
                               indvScenario, updatedConfigPath,
                               baseWorkingDir, generatedTestScenariosScratchDirName)

                # Navigate back to parent directory
                os.chdir(baseWorkingDir + '/' + generatedTestScenariosScratchDirName)
            else:
                logging.info("Skipping data generation for scenario: %s" % indvScenario['scenarioName'])


def validateInputsAndSetDefaults(lCmdArguments, lScriptDir, lConfigName, lTestScenariosName, lInputDirName):

    """
    This is an function to do validation and setting up of default values incase they are not specified for config and
    input dir
    :param lCmdArguments: An argparse object
    :param lScriptDir: Path to directory where script is located
    :param lConfigName: The default config file name
    :param lTestScenariosName: The default test scenarios file name
    :param lInputDirName: The default input directory name
    :return: Updated argparse object
    """

    # If the config file is not specified then take one in script dir. Check finally if it exists
    if lCmdArguments.config is None:
        lCmdArguments.config = lScriptDir + "/configs/" + lConfigName

    lCmdArguments.config = os.path.abspath(lCmdArguments.config)

    if not os.path.isfile(lCmdArguments.config):
        sys.exit(logging.error("The config file at %s was not found. Exiting" % lCmdArguments.config))

    # If the scenarios file is not specified then take one in script dir. Check finally if it exists
    if lCmdArguments.scenarios is None:
        lCmdArguments.scenarios = lScriptDir + "/scenarios/" + lTestScenariosName

    lCmdArguments.scenarios = os.path.abspath(lCmdArguments.scenarios)

    if not os.path.isfile(lCmdArguments.scenarios):
        sys.exit(logging.error("The scenarios file at %s was not found. Exiting" % lCmdArguments.scenarios))

    # If the directory with /inputs is not specified set it to inputs in script dir. Check finally if it exists
    if lCmdArguments.dir is None:
        lCmdArguments.dir = lScriptDir + "/" + lInputDirName

    if not os.path.exists(lCmdArguments.dir):
        sys.exit(logging.error("The inputs directory at %s was not found. Exiting" % lCmdArguments.dir))

    return lCmdArguments


def generateScenarioConfig(lIndvScenario, lConfigFilePath, lInputsDir):

    """
    Compare the scenarios file and the default config file and generate a config file for the current scenario
    :param lIndvScenario: Dict from csv for an individual scenario row
    :param lConfigFilePath: The path to default config file
    :param lInputsDir: The path to directory with inputs
    :return: Path to generated scenario file
    """

    # Reading in default config file
    configObj = ConfigParser.RawConfigParser()
    configObj.read(lConfigFilePath)

    # Setting in the base directory to inputs dir
    configObj.set('defaultFiles', 'baseDir', lInputsDir)

    # Setting the simulator parameters
    if lIndvScenario['simulatorOptions'] != "":
        for indvSimulatorEditedParams in lIndvScenario['simulatorOptions'].split(';'):
            indvSimulatorEditedParamsList = indvSimulatorEditedParams.split('=')
            try:
                configObj.set('simulator', indvSimulatorEditedParamsList[0], indvSimulatorEditedParamsList[1])
            except IndexError:
                logging.warning("The simulator parameter set %s was specified incorrectly. Ignoring it" %
                                indvSimulatorEditedParams)

    # Setting the values for all other parameters specified in the config file
    for sectionName in configObj.sections():
        for (parameterName, parameterValue) in configObj.items(sectionName):
            if sectionName + ":" + parameterName in lIndvScenario and lIndvScenario[
                                sectionName + ":" + parameterName] != '':
                configObj.set(sectionName, parameterName, lIndvScenario[sectionName + ":" + parameterName])

    # Writing in this edited config file to scenario directory
    with open('./' + lIndvScenario['scenarioName'] + "_config.ini", 'wb') as updatedConfigFile:
        configObj.write(updatedConfigFile)

    # Return the absolute path of created config file
    return os.path.abspath('./' + lIndvScenario['scenarioName'] + "_config.ini")


def copyToFinalDir(lFinalDirPath, lSimulatedBamFilePath, lIndvScenario, lConfigPath, lBaseWorkingDir,
                   lTestScenarioScratchDirName):

    """
    Copy the final bam, bai and microhaplotype bed files to scenario dir and also create a readme
    :param lFinalDirPath: The final test scenario directory
    :param lSimulatedBamFilePath: The path of simulated bam file for current scenario
    :param lIndvScenario: A dict of details for the current scenario
    :param lConfigPath: The path to config file used to generate the scenario
    :param lBaseWorkingDir: The path to base directory where scratch outputs are stored
    :param lTestScenarioScratchDirName: The name of the directory within test scenarios scratch which
    is storing scratch results
    :return: None
    """

    # Create output directory for the current scenario
    if not os.path.exists(lFinalDirPath + "/" + lIndvScenario['scenarioName']):
        os.mkdir(lFinalDirPath + "/" + lIndvScenario['scenarioName'])

    currentScenarioOutputDir = lFinalDirPath + "/" + lIndvScenario['scenarioName']

    # Copy the bam, bai and microhaplotype bed and config files.
    if os.path.isfile(lSimulatedBamFilePath):
        shutil.copy(lSimulatedBamFilePath, currentScenarioOutputDir)

    if os.path.isfile(lSimulatedBamFilePath + ".bai"):
        shutil.copy(lSimulatedBamFilePath + ".bai", currentScenarioOutputDir)

    if os.path.isfile(lConfigPath):
        shutil.copy(lConfigPath, currentScenarioOutputDir)
    else:
        logging.error('Updated Scenario config file at %s does not exist. Exiting' % lConfigPath)
        sys.exit(1)

    # Reading in the config file to get values for validation script path, microhaplotype bed file and info file
    configObj = ConfigParser.ConfigParser()
    configObj.read(lConfigPath)

    lMicrohapValidationScriptPath = os.path.abspath(configObj.get('defaultFiles', 'microhaplotypeValidationScript'))
    lMicrohapBedFilePath = os.path.abspath(configObj.get('defaultFiles', 'defaultHotspots'))
    lMicrohapInfoFilePath = os.path.abspath(configObj.get('defaultFiles', 'microhaplotypeInfoFile'))
    lMicrohapMinCov=configObj.get('defaultFiles', 'microhaplotypeMinCov')

    # Copying in the Microhaplotype bed file and info file
    if os.path.isfile(lMicrohapBedFilePath):
        shutil.copy(lMicrohapBedFilePath, currentScenarioOutputDir)
    else:
        logging.error('Microhaplotype bed file at %s does not exist. Exiting' % lMicrohapBedFilePath)
        sys.exit(1)

    if os.path.isfile(lMicrohapInfoFilePath):
        shutil.copy(lMicrohapInfoFilePath, currentScenarioOutputDir)

    # Create the readme files with details of the generated scenarios
    readmeFile = open(currentScenarioOutputDir + '/' + lIndvScenario['scenarioName'] + "_readme.txt", 'w')
    readmeFile.write('Scenario Name:' + lIndvScenario['scenarioName'] + '\n')
    readmeFile.write('Scenario Description:' + lIndvScenario['scenarioDescription'] + '\n')
    readmeFile.write('Scenario Scratch Dir:' + lBaseWorkingDir + '/' +
                     lTestScenarioScratchDirName +
                     '/' + lIndvScenario['scenarioName'] + '\n')
    readmeFile.write('Run Date:' + datetime.utcnow().strftime('%Y_%m_%d_%H_%M_%S_%f')[:-4] + '\n')
    readmeFile.write('Bam File:' + os.path.basename(lSimulatedBamFilePath) + '\n')
    readmeFile.write('Bed File:' + os.path.basename(lMicrohapBedFilePath) + '\n')
    readmeFile.write('Info File:' + os.path.basename(lMicrohapInfoFilePath) + '\n')
    readmeFile.write('Minimum Coverage:' + lMicrohapMinCov + '\n')
    readmeFile.write('Truth file:' + lIndvScenario['scenarioName'] + "_truth.txt" + '\n')
    readmeFile.close()

    # Generate the truth.txt
    validationScriptCommand = 'python ' + lMicrohapValidationScriptPath + \
                              ' --bam ' + lSimulatedBamFilePath + \
                              ' --bed ' + lMicrohapBedFilePath + \
                              ' --o ' + currentScenarioOutputDir + \
                              '/' + lIndvScenario['scenarioName'] + "_truth.txt" + \
                              ' --minCov ' + lMicrohapMinCov

    if os.path.basename(lMicrohapInfoFilePath) != "-":
        validationScriptCommand += ' --info ' + lMicrohapInfoFilePath

    validationScriptCommandSplit = shlex.split(validationScriptCommand)

    subprocess.check_output(validationScriptCommandSplit)


def executeReadSimulationWrapper(lIndvScenario, lConfigPath):

    """
    Execute the reads simulation wrapper on the created file and with the additional options if any given
    in the scenario file.
    :param lIndvScenario: Dict from csv for an individual scenario row
    :param lConfigPath: The path of default config file
    :return: Path to generated bam file
    """

    # Read in the config file to get the path for readSimulationWrapper.py file
    configObj = ConfigParser.ConfigParser()
    configObj.read(lConfigPath)

    # Creating the base command for calling read simulation wrapper
    readSimulationWrapperCommand = "python " + configObj.get('defaultFiles', 'readSimulationWrapper') + \
                                   " -c " + lConfigPath + \
                                   " -j " + lIndvScenario['scenarioName'] + "_reads" + \
                                   " -a " + configObj.get('defaultFiles', 'postSimulationScript')

    # Adding the random generation seed if specified in scenario
    if 'simulationSeed' in lIndvScenario and lIndvScenario['simulationSeed'] != "":
        readSimulationWrapperCommand += " -s " + lIndvScenario['simulationSeed']

    # Split generated command
    readSimulationWrapperCommandSplit = shlex.split(readSimulationWrapperCommand)

    # Execute the generated command
    subprocess.check_output(readSimulationWrapperCommandSplit)

    # Get the most appended bam file
    rLatestBamFile = os.path.abspath(max(glob.iglob('*.bam'), key=os.path.abspath))

    # Get the microhaplotype bed file path
    # rMicrohaplotypeBedFile = os.path.abspath(configObj.get('defaultFiles', 'microhaplotypeBedFile'))

    return rLatestBamFile


if __name__ == '__main__':
    main()