import csv
import re

from allele import Allele
from base_file import File
from microhaplotype import Microhaplotype


class TxtFile(File):
    UNAVAILABLE_DATA = []

    def __init__(self, txt_file_path):
        super(TxtFile, self).__init__(file_path=txt_file_path, file_type='txt',
                                      unavailable_data=self.UNAVAILABLE_DATA)

    def parse(self):
        """Parse the txt file."""
        with open(self.file_path, 'rU') as txt_file_handle:
            txt_file_handle.readline()  # Get out the filters line
            number_of_contributors_line = txt_file_handle.readline()
            self.number_of_contributors = self._get_number_of_contributors(number_of_contributors_line)
            txt_lines = csv.DictReader(txt_file_handle, delimiter='\t')
            for single_txt_line_dict in txt_lines:
                microhap_name = single_txt_line_dict['MarkerId']
                self.mh_dict[microhap_name] = self._get_mh_obj(single_txt_line_dict)

    @staticmethod
    def _get_number_of_contributors(number_of_contributors_line):

        """ Parse the number of contributors line in the txt file.
        ::param number_of_contributors_line: Example #The minimum number of contributor : 7
        ::return Integer of number of contributors Example: 7"""

        return int(number_of_contributors_line.split(' : ')[-1])

    def _get_mh_obj(self, mh_single_line_dict):

        """ Get an object of class Microhaplotype populated with its attributes.
        :param mh_single_line_dict: A csv.DictReader dict of a single mh txt file line
        :return: An object of class Microhaplotype
        """
        ret_mh_obj = Microhaplotype()
        ret_mh_obj.name = mh_single_line_dict['MarkerId']
        ret_mh_obj.ae = float(mh_single_line_dict['AvgAe'])
        ret_mh_obj.prob_det_mix = float(mh_single_line_dict['ProbOfDetMix'])
        ret_mh_obj.number_of_alleles = int(mh_single_line_dict['NumOfAlleles'])
        ret_mh_obj.number_of_contributors = int(mh_single_line_dict['NumOfContributors'])
        ret_mh_obj.analytical_threshold = int(mh_single_line_dict['MinimumCoverage'])
        ret_mh_obj.snps_rsid_list = mh_single_line_dict['SNPs']
        ret_mh_obj.alleles = self._get_alleles(mh_single_line_dict['AlleleCoverageDetail(-strand:+strand:avgMQ)'])
        return ret_mh_obj

    @staticmethod
    def _get_alleles(allele_details_string):

        """ Get an populated allele dict.
        :param allele_details_string: Example [AA]296(143:153:69),[GT]806(410:396:70)
        :return: A dict of type key:allele seq value: populated object of Allele class
        """

        ret_allele_dict = {}
        for single_allele in filter(lambda x: x != "", allele_details_string.split(',')):
            allele_match = re.search(r'\[(.*?)\](\d+)\((\d+):(\d+):(\d+)\)', single_allele)
            allele_name = allele_match.group(1)
            ret_allele_dict.setdefault(allele_name, Allele())
            ret_allele_dict[allele_name].sequence = allele_match.group(1)
            ret_allele_dict[allele_name].coverage_total = int(allele_match.group(2))
            ret_allele_dict[allele_name].coverage_reverse = int(allele_match.group(3))
            ret_allele_dict[allele_name].coverage_forward = int(allele_match.group(4))
            ret_allele_dict[allele_name].avg_mapping_quality = int(allele_match.group(5))
        return ret_allele_dict


# =============Test================
if __name__ == '__main__':
    test_txt = TxtFile('./test_compare_microhaplotypes/test_output.txt')
    test_txt.parse()
    print test_txt.number_of_contributors
    print test_txt.mh_dict
