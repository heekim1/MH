### Microhaplotype test suite

This is the test suite to generate and test scenarios for microhaplotype mixture analysis. The suite consists of two main components a test scenarios generator **(generateTestScenarios.py)** and an evaluator **(evaluateTestScenarios.py)**. The test suite generation script will also generate an expected truth file for every scenario. This truth file will be evaluated against the output from the Microhaplotyper.jar file and a report in csv format will be created indicating pass/fail status for all the scenarios. In case of a failure the differences will be highlighted in result details.

## Installing Prerequisites
### Installing dwgsim
The simulation tool [dwgsim](https://github.com/nh13/DWGSIM) can be installed on Ubuntu and Mac by the following steps
#### Mac
```bash
git clone https://github.com/nh13/DWGSIM.git
cd DWGSIM
git submodule init
git submodule update
make
sudo cp dwgsim /usr/bin/
sudo chmod +x /usr/bin/dwgsim
```
#### Ubuntu 14.04
```bash
git clone https://github.com/nh13/DWGSIM.git
cd DWGSIM
git submodule init
git submodule update
sudo apt-get install libncurses5-dev
make
sudo cp dwgsim /usr/bin/
sudo chmod +x /usr/bin/dwgsim
```
The commands are also included in the installation scripts for dwgsim at *scripts/dwgsim_install_script_mac.sh* and *scripts/dwgsim_install_script_ubuntu.sh*.
If you have installed dwgsim at any other path than /usr/bin/dwgsim please update the path in the config.ini file for parameter *simulator_path*.
It would also be helpful to go through the dwgsim wiki at https://github.com/nh13/DWGSIM/wiki/Simulating-Reads-with-DWGSIM

### Installing bedtools
The pybedtools module required bedtools. To install it follow the below given steps
```bash
# Mac
sudo brew install bedtools
# Ubuntu
sudo apt-get install bedtools
```

### Installing the mapping tool(BWA/TMAP)
The script can optionally use BWA/TMAP for mapping. Please install either mapping tool and add it to path. For more details on installation visit for bwa: https://github.com/lh3/bwa, tmap: https://github.com/iontorrent/TMAP
IF you choose BWA as the mapping tool change the following line in the config file at /configs/baseConfig.ini
```bash
# from
postSimulationScript=%(baseDir)s/scripts/postSimulationScripts/postSimulationTMAP.sh
# to
postSimulationScript=%(baseDir)s/scripts/postSimulationScripts/postSimulationBWA.sh
```

### Installing python modules
The script will require one python module i.e. pyvcf as a prerequisite. If you have pip installed on you system then installation should be straightforward as

```bash
[sudo] pip install numpy pandas pyvcf pysam pybedtools 
```

If you don't have pip on you system you can install it using the following commands
#### Mac
```bash
wget https://bootstrap.pypa.io/get-pip.py
sudo get-pip.py
```

#### Ubuntu 14.04
```bash
sudo apt-get install python-pip
```

## Installing the Microhaplotype Test Suite
Once the prerequisites are installed you can use the microhaplotype test generation and evaluation scripts. You need to do a few changes in the configuration files so that the scripts can easily find the required files. You can do these changes to the config file at */configs/baseConfig.ini*
```bash
# change
defaultReference=.
# to
defaultReference=<path where the reference human fasta file is located>
# e.g. scriptDir=/Users/user1/Documents/hg19.fasta
# Please use fasta file from ion community at - http://s3.amazonaws.com/IonTorrent/reference/hg19.zip

# If you have dwgsim installed at a different location then change path for dwgsim too
# change
simulator_path=/usr/bin/dwgsim
# to
simulator_path=<path/to/dwgsim/executable>
```

## Running the Microhaplotype Test Suite
The microhaplotypeTestSuite has two components generateTestScenarios.py and evaluateTestScenarios.py. The generateTestScenarios.py script will take a file and a base config file as input. Optionally if the user has changed the location of /inputs/ directory he can also specify that. The user can specify an output directory to store all the generated scenarios(default is microhaplotypeTestScenarios dir in run directory). The script will create two directories microhaplotypeTestScenarios for the scenarios and microhaplotypeTestScenariosScratch for the background files in case required for debugging.

After generating the scenarios the user can run the evaluateTestScenarios.py script to compare the results to expected outputs and give a report highlighting the pass/fail status. It will also create a directory to store all the output files from the microhaplotyper.jar file in case required for comparison. The evaluation script mainly takes two inputs - the location of generated scenarios from earlier script and location of the microhaplotyper.jar. The user can optionally change the name of output file

An example run of the suite would be as follows:
```bash
python generateTestScenarios.py -s ./scenarios/trialScenarios.tsv -c ./configs/trialConfig.ini
python evaluateTestScenarios.py -j /home/user1/Microhaplotyper.jar -d /home/microhaplotypeTestScenarios
```

## generateTestScenarios.py command line options
* **'-s'/'--scenarios'** - [STRING-FILEPATH,OPTIONAl,DEFAULT=/scenarios/baseScenarios.ini]
The path to file giving details of the scenarios to be executed. Please refer to section on scenarios file below for more details

* **'-c', '--config'** - [STRING-FILEPATH,OPTIONAl,DEFAULT=/configs/baseConfig.ini]
The path to file given the base configuration options. Please refer to section on config file for more details

* **'-d', '--dir'** - [STRING-DIRPATH,OPTIONAl,DEFAULT=/inputs]
The script takes in all the input files from an included folder(/inputs/). If you change to location of this folder or are running the script from a different location then specfiy the new path here

* **'-o', '--outdir'** -  [STRING-DIRPATH,OPTIONAl,DEFAULT=<current directory>]
The script will by default create two folders microhaplotypeTestScenarios and microhaplotypeTestScenariosScratch for the scenarios and other background files respectively in the current working directory. If you want to change the location where these folders shall be create specify the new location/dir here

## evaluateTestScenarios.py command line options
* **'-d', '--dir'** - [STRING-DIRPATH,REQUIRED]
This is the pointer to the location where the scenarios generated by generateTestScenarios.py are located. Please give path ending with /microhaplotypeTestScenarios.

* **'-j', '--jar'** - [STRING-FILEPATH,REQUIRED]
Path to the Microhaplotyper.jar file for which the comparison with truth is required

* **'-o', '--out'** - [STRING-FILEPATH,OPTIONAl,DEFAULT=./comparison_output.csv]
Optionally change the output file name created by the evaluation script

## Extended configuration with baseConfig.ini
The test scenario generation script takes a .ini file as input. This file specfied all the basic configuration options that will be applied to all the generated scenarios. If the user wants to change these base configuration settings for any scenario he can specfiy the changes in the scearios.tsv file. Changes specfied in the scenarios file will override the ones in the .ini file.

The config file used variable substitution as in if a variable baseDir is specfied earlier we can re use its value in the rest of the configuration options by including it as %(baseDir)s

Details of the configuration options available in the .ini file are as follows:-

- **[defaultFiles]**
  - **baseDir** = This will be autofilled by the wrapper as the location of the /inputs dir
  - **defaultReference** = The reference .fasta/fa file
  - **defaultBed** = The bed file to restrict the generated reads to
  - **defaultMutations** = A csv file listing all the mutations required. This file will be only considered if the flag global_options:useMutationsFile is set to true. This is a simple file giving a list on mutations that we want each contributor to have. Refer to the mutations file section below for details on it.
  - **defaultHotspots** = This is the microhaplotyper hotspots file. This file is used by the generation scripts to generate a set of mutations for each contributor on run time
  - **vcfTemplate** = This is the template used for vcf header. Unless required dont change this path

  - **readSimulationWrapper** = This is the location of the read simulation wrapper script.
  - **postSimulationScript** = This is the location of the script which will be run once the fastq files are generated by dwgsim. You can do mapping and any post processing if required on the bam files in these scripts.
  - **microhaplotypeValidationScript** = The location of the validation script which will generate the expected output
  - **microhaplotypeInfoFile** = The location of the info file used by the Microhaplotyper.jar file to get information of Ae and ProbofDetMixture values.
  - **microhaplotypeMinCov** = This is the values --minCov parameter for Microhaplotyper.jar


- **[simulator]**

This section is for dwgsim options. You can include any other options for dwgsim that you want or change the values of current options.
The options are used the generate the dwgsim option commandline given by the parameter **sim_options**

If you want to add a parameter lets say rate of mutations you could add it as
```bash
rateOfMutations = 0.1
sim_options = -e %(base_error_rate)s -1 %(first_read_length)s ... -r %(rateOfMutations)s
```
Thus the general syntax for additions is %(paramName)s

Some options like the coverage, reference file, bed file, vcf file are added in by the wrapper later.


- **[global_options]**
  - **base_coverage** = The base coverage for mixtures. If a contributor is 1x it will have this coverage. If 2x then coverage for it is 2*(base_coverage) value
  - **simulator_path** = Path for the dwgsim executable. Change it if dwgsim is not at /usr/bin/dwgsim
  - **default_mixture_ratio** = The default mixture ratio
  - **useMutationsFile** = Flag(True/False) which will tell the generator to ignore the hotspots file and use the Mutations file instead for the mutations in contributors.

## The scenarios file
The scenarios file is a tab separated tsv file. You can edit it in a text editor for minor edits by enabling the option of view whitespace if available. In order to do major edits you can open it in Microsoft Excel by specifying the delimiter as tab and all columns as text.

This is a dynamic file in a way that you can override any of the base config options by including a column for it in this file. IF you want to override the vcfTemplate option from the base config for a few scenarios you can insert a column by the name defaultFiles:vcftemplate. You only need to specify new values for the scenarios where you want to make a change. If left blank the generator will pick values from the config file.

You also have the option of disabling the execution of some of the scenarios by adding a # to the start of the scenario name. So if you want to skip the scenario s2_r11_c2 while generation change its name to #s2_r11_c2 

*Note: Incase you are adding a new parameter to the scenarios file please use the write parameter name in lowercase. For example if you want to add parameter postSimulationScript to the scenarios file add it as defaultFiles:postsimulationscript and not defaultFiles:postSimulationScript

## Mutations csv file
This is a csv file specifying the mutations to be incorporated in the generated reads.
The csv file has the following 5 columns

1. **Chromosome** - The chromosome where the mutation should be placed
2. **Position** - The position of the specified chromosome
3. **Reference Allele** - The allele in the reference file at specified position
4. **Alternate Allele** - The alternate alleles to be included in place of the reference. For multiple contributors separate the alternate alleles by '/'. For keeping the allele same as reference keep it as '-'. For e.g. T/-/C indicates you want contributor 1 to have T, contributor 2 to have allele same as reference and contributor3 to have allele C.
5. **Allele frequency** - The frequency of the alternate alleles. The possible values are [1,0.5,-]. 1 will indicate a Homozygous mutation, 0.5 will indicate heterozygous mutation and - will be ignored as alternate allele will be same as reference.

You can specify values for multiple contributors in the csv file and use only part of them. i.e. you can specify alternate alleles as T/A/C but have the mixture ratio as 1:2 in which case values for one first two contributors will be considered.

## The comparison output file
The comparison output file is a csv file which will give PASS/FAIL status  for all the executed scenarios. Incase of a failure it will give the details of the differences

## Maintainence notes
The main script is the one under scripts/readSimulationWrapper.py. It is the one that will genenrate the required bam files. All other scripts will just give it the required output for other scenarios. You can use the generateTestScenarios script as is or use the wrapper to generate your scenarios. Instead of using the evalulateTestScenarios you can prefer to use the execute_nightly_test script for evalulation purposes once you have generated the test scenarios.
