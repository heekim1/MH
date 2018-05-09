import yaml  # Using YAML instead of JSON to parse JSON and yaml will load it as strings instead of unicode

from allele import Allele
from base_file import File
from microhaplotype import Microhaplotype


class JsonFile(File):
    UNAVAILABLE_DATA = ['number_of_contributors_total', 'avg_mapping_quality']

    def __init__(self, json_file_path):
        super(JsonFile, self).__init__(file_path=json_file_path, file_type='json',
                                       unavailable_data=self.UNAVAILABLE_DATA)

    def parse(self):

        """Parse the json file"""

        with open(self.file_path, 'rU') as json_file_handle:
            json_results_list = yaml.safe_load(json_file_handle)['result']
            for single_microhap_dict in json_results_list:
                microhap_name = single_microhap_dict['locus']
                self.mh_dict[microhap_name] = self._get_mh_obj(single_microhap_dict)

    def _get_mh_obj(self, single_microhap_dict):
        ret_mh_obj = Microhaplotype()
        ret_mh_obj.name = single_microhap_dict['locus']
        ret_mh_obj.ae = float(single_microhap_dict['avgAe'])
        ret_mh_obj.prob_det_mix = float(single_microhap_dict['probOfDetMix'])
        ret_mh_obj.number_of_alleles = int(single_microhap_dict['numOfAlleles'])
        ret_mh_obj.number_of_contributors = int(single_microhap_dict['numOfCont'])
        ret_mh_obj.analytical_threshold = int(single_microhap_dict['minCoverage'])
        ret_mh_obj.snps_rsid_list = '/'.join(single_microhap_dict['snps'])
        ret_mh_obj.alleles = self._get_alleles(single_microhap_dict['profile'])
        return ret_mh_obj

    @staticmethod
    def _get_alleles(alleles_list):
        ret_allele_dict = {}
        for single_allele_details in alleles_list:
            allele_name = single_allele_details['displayValue']
            ret_allele_dict.setdefault(allele_name, Allele())
            ret_allele_dict[allele_name].sequence = allele_name
            ret_allele_dict[allele_name].coverage_total = int(single_allele_details['coverage']['total'])
            ret_allele_dict[allele_name].coverage_reverse = int(single_allele_details['coverage']['minus'])
            ret_allele_dict[allele_name].coverage_forward = int(single_allele_details['coverage']['plus'])
        return ret_allele_dict
