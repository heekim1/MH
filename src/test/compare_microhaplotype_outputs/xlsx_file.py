import pandas as pd

from allele import Allele
from base_file import File
from microhaplotype import Microhaplotype


class XlsxFile(File):
    UNAVAILABLE_DATA = ['number_of_contributors_total', 'ae', 'prob_det_mix', 'snps_rsid_list', 'analytical_threshold',
                        'coverage_forward', 'coverage_reverse', 'avg_mapping_quality']

    EXCEL_COLUMN_NAMES = ['SNo', 'Locus', 'Genotype', 'Allele Count', 'Min Contributors', 'Minus Coverage']
    ALLELE_COLUMN_MAPPING = {'Allele Name': 'Locus', 'Coverage Total': 'Allele Count'}

    def __init__(self, xlsx_file_path):
        super(XlsxFile, self).__init__(file_path=xlsx_file_path, file_type='xlsx',
                                       unavailable_data=self.UNAVAILABLE_DATA)

    def parse(self):
        pandas_df = self._get_df_from_excel()
        pandas_df.columns = self.EXCEL_COLUMN_NAMES
        active_microhaplotype_name = ''
        for index, single_row in pandas_df.iterrows():
            if single_row['Locus'] == 'Allele':
                continue
            elif not pd.isnull(single_row['SNo']):
                active_microhaplotype_name = single_row['Locus']
                self.mh_dict[single_row['Locus']] = self._get_mh_obj(single_row)
            else:
                allele_name = single_row[self.ALLELE_COLUMN_MAPPING['Allele Name']]
                self.mh_dict[active_microhaplotype_name].alleles[allele_name] = self._get_allele_obj(single_row)

    def _get_df_from_excel(self):
        return pd.read_excel(self.file_path)

    @staticmethod
    def _get_mh_obj(genotype_row):
        ret_mh_obj = Microhaplotype()
        ret_mh_obj.name = str(genotype_row['Locus'])
        ret_mh_obj.number_of_alleles = int(genotype_row['Allele Count'])
        ret_mh_obj.number_of_contributors = int(genotype_row['Min Contributors'])
        return ret_mh_obj

    def _get_allele_obj(self, allele_row):
        ret_allele_obj = Allele()
        ret_allele_obj.sequence = str(allele_row[self.ALLELE_COLUMN_MAPPING['Allele Name']])
        ret_allele_obj.coverage_total = int(allele_row[self.ALLELE_COLUMN_MAPPING['Coverage Total']])
        return ret_allele_obj


# =========== Test =========
if __name__ == '__main__':
    xlsx_file = XlsxFile('./test_compare_microhaplotypes/inputs_test_compare_microhaplotypes/ui.xlsx')
    xlsx_file.parse()
