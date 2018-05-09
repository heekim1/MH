import argparse
import utils


class ArgumentParser(object):
    def __init__(self):
        self.args = self.get_commandline_arguments()

    @staticmethod
    def get_commandline_arguments():
        description = "Script to get microhaplotype allele counts"
        epilog = "EXAMPLE: python main.py --bam in.bam --bed in.bed --info info.txt --out out.txt --mincov 0.03"
        parser = argparse.ArgumentParser(description=description, epilog=epilog)

        required_args_group = parser.add_argument_group('required arguments')
        required_args_group.add_argument('-b', '--bam', dest='bam_file_path', required=True,
                                         type=lambda x: utils.is_valid_file(parser, x))
        required_args_group.add_argument('-e', '--bed', dest='bed_file_path', required=True,
                                         type=lambda x: utils.is_valid_file(parser, x))

        parser.add_argument('-i', '--info', dest='info_file_path', type=lambda x: utils.is_valid_file(parser, x))
        parser.add_argument('-o', '--out', dest='out_file_path')
        parser.add_argument('-m', '--mincov', help='Allowed values 0-1', dest='min_coverage', default=0.02,
                            type=lambda x: utils.is_valid_min_cov_value(parser, x))

        return parser.parse_args()
