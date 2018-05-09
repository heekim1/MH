import argparse

import utils


class ArgumentParser(object):
    def __init__(self):
        self.args = self.get_commandline_arguments()

    @staticmethod
    def get_commandline_arguments():
        description = "Script to compare different microhaplotype outputs"
        epilog = "EXAMPLE: python main.py --file1 file1.txt --file2 file2.json"
        parser = argparse.ArgumentParser(description=description, epilog=epilog)

        required_args_group = parser.add_argument_group('required arguments')
        required_args_group.add_argument('-f', '--file1', dest='file1_path', required=True,
                                         type=lambda x: utils.is_valid_file(parser, x))
        required_args_group.add_argument('-i', '--file2', dest='file2_path', required=True,
                                         type=lambda x: utils.is_valid_file(parser, x))
        return parser.parse_args()
