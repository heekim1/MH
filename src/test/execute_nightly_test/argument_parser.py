import argparse

import utils


class ArgumentParser(object):
    DEFAULT_SCENARIOS_DIR = '/Users/daklep/Documents/Projects/8_HID_NGS/2_Data_Simulation/13_microhaplotypeTestSuite/microhaplotypeTestScenarios'
    DEFAULT_JAR_LOCATION = '/Users/daklep/Documents/Projects/3_Converge_Build/converge/microhaplotype/target/' \
                           'Microhaplotyper-jar-with-dependencies.jar'
    DEFAULT_OUTPUT_DIR = './output'

    def __init__(self):
        self.args = self.get_commandline_arguments()

    def get_commandline_arguments(self):
        description = "Script to execute nightly tests for microhaplotype validation"
        epilog = "EXAMPLE: python main.py --scenariosdir /path/to/scenarios/dir " \
                 "--mhjarlocation /path/to/microhaplotyper.jar --outdir /path/to/output/dir"

        parser = argparse.ArgumentParser(description=description, epilog=epilog)

        parser.add_argument('-i', '--scenariosdir', dest='scenariosdir', default=self.DEFAULT_SCENARIOS_DIR,
                            type=lambda x: utils.is_valid_dir(parser, x))
        parser.add_argument('-j', '--mhjarlocation', dest='mhjarlocation', default=self.DEFAULT_JAR_LOCATION,
                            type=lambda x: utils.is_valid_file(parser, x))
        parser.add_argument('-o', '--outdir', dest='outdir', default=self.DEFAULT_OUTPUT_DIR,
                            type=lambda x: utils.is_valid_dir(parser, x))
        return parser.parse_args()
