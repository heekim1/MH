import logging
import os
import shlex
import subprocess
from datetime import datetime

from argument_parser import ArgumentParser
from compare_microhaplotype_outputs.compare_mh_results import CompareMHResults
from microhaplotype_truth_generator.main import Microhaplotyper
from microhaplotype_truth_generator.output_microhaplotypes import OutputMicrohaplotypes


class ExecuteNightlyTest(object):
    """
    Execute nightly testing for HID Genotyper
    This will execute the following functionality:-
    1. List all the scenarios in the given scenarios directory
    2. Read in the readme.txt file in the scenario directory and get all the paramters
    3. Create an output directory
    4. Generate the expected truth and the jar output
    5. Compare both and generate a results output file
    """

    def __init__(self, scenarios_dir, microhaplotype_jar_location, output_dir):

        """
        :param scenarios_dir: Path to the directory where scenarios are located
        :param microhaplotype_jar_location: Location of the microhaplotype jar file
        :param output_dir: The path where output results will be stored
        """

        self.scenarios_dir = scenarios_dir
        self.microhaplotype_jar_location = microhaplotype_jar_location
        self.parent_output_dir = output_dir
        self.scenarios_output_dir = self._get_scenarios_dir(self.parent_output_dir)
        self.comparison_results = {}

    def test_scenarios(self):
        for sub_dirname in os.listdir(self.scenarios_dir):
            logging.info("Evaluation scenario {}".format(sub_dirname))
            readme_file_location = self.scenarios_dir + '/' + sub_dirname + '/' + sub_dirname + '_readme.txt'
            if os.path.isfile(readme_file_location):
                run_params_dict = self._get_run_params(readme_file_location)
                current_scenario_input_dir = self.scenarios_dir + '/' + sub_dirname
                current_scenario_output_dir = self._create_current_scenario_dir(sub_dirname)
                truth_location = self._generate_truth(run_params_dict, current_scenario_input_dir,
                                                      current_scenario_output_dir)
                jar_output_location = self._generate_jar_output(run_params_dict, current_scenario_input_dir,
                                                                current_scenario_output_dir)
                comparison_obj = CompareMHResults(truth_location, jar_output_location)
                comparison_obj.compare()
                self.comparison_results[sub_dirname] = comparison_obj

    def output_results(self):
        output_file_timestamp = self.scenarios_output_dir.split(os.path.sep)[-1].split('-')[-1]
        with open(self.scenarios_output_dir + '/mh_simulated_' +
                          output_file_timestamp + '_validation_output.csv', 'w') as output:
            header = ["Scenario Name", "Number of Contributors Check", "Number of Contributors Check Detail",
                      "Results Check", "Results Details"]
            output.write(','.join(header) + '\n')
            for scenario_name in sorted(self.comparison_results.keys()):
                output.write(scenario_name + ',' + str(self.comparison_results[scenario_name]) + '\n')

    @staticmethod
    def _generate_truth(run_params, input_dir, output_dir):
        truth_obj = Microhaplotyper(bam=input_dir + '/' + run_params['Bam File'],
                                    bed=input_dir + '/' + run_params['Bed File'],
                                    info=input_dir + '/' + run_params['Info File'],
                                    analytical_threshold=float(run_params['Minimum Coverage']),
                                    out=output_dir + '/truth.txt')
        truth_obj.get_microhaplotypes()
        output_obj = OutputMicrohaplotypes(truth_obj)
        output_obj.output_microhaplotypes()
        return output_dir + '/truth.txt'

    def _generate_jar_output(self, run_params, input_dir, output_dir):
        jar_execution_command = 'java -jar ' + self.microhaplotype_jar_location \
                                + ' -bam ' + input_dir + '/' + run_params['Bam File'] \
                                + ' -bed ' + input_dir + '/' + run_params['Bed File'] \
                                + ' -info ' + input_dir + '/' + run_params['Info File'] \
                                + ' -minCov ' + run_params['Minimum Coverage'] \
                                + ' -o ' + output_dir + '/jar_output.txt'
        jar_execution_command_list = shlex.split(jar_execution_command)
        subprocess.check_output(jar_execution_command_list)
        return output_dir + '/jar_output.txt'

    def _create_current_scenario_dir(self, sub_dirname):
        current_scenario_dir = self.scenarios_output_dir + '/' + sub_dirname
        self._create_dir(current_scenario_dir)
        return current_scenario_dir

    def _get_scenarios_dir(self, parent_output_dir):
        scenarios_output_dir = parent_output_dir + "/mh_nightly_output-" + datetime.now().strftime('%d_%m_%Y_%H_%M_%S')
        self._create_dir(scenarios_output_dir)
        return scenarios_output_dir

    @staticmethod
    def _create_dir(dir_path):
        if not os.path.isdir(dir_path):
            os.mkdir(dir_path)

    @staticmethod
    def _get_run_params(readme_file_location):
        with open(readme_file_location) as readme_file:
            return {readme_line.split(':')[0]: readme_line.split(':')[1] for readme_line in
                    map(lambda x: x.rstrip('\n'), readme_file.readlines())}


if __name__ == '__main__':
    logging.basicConfig(format='%(asctime)s %(levelname)s : %(message)s', level=logging.INFO)
    args_obj = ArgumentParser()
    nightly_obj = ExecuteNightlyTest(scenarios_dir=args_obj.args.scenariosdir,
                                     microhaplotype_jar_location=args_obj.args.mhjarlocation,
                                     output_dir=args_obj.args.outdir)
    nightly_obj.test_scenarios()
    nightly_obj.output_results()
