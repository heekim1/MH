import csv


class Info(object):
    def __init__(self, info_file_path):
        self.info_file_path = info_file_path
        self.info_file_headers = ('MHLongName', 'MHShortName', 'ProbOfDetMix', 'Ae', 'Unknown1', 'Unknown2', 'Chr',
                                  'Pos', 'RSIDList', 'LocusName', 'NumberOfSNPs')
        self.info_dict = {}

    def get_info_dict(self):
        with open(self.info_file_path, 'rU') as info_file_handle:
            info_file_lines = csv.DictReader(info_file_handle, fieldnames=self.info_file_headers, delimiter='\t')
            self.info_dict = {info_file_single_line['MHShortName']: self._transform_data_types(info_file_single_line)
                              for info_file_single_line in info_file_lines}

    def add_info(self, mh_obj):
        if mh_obj.name in self.info_dict.keys():
            mh_obj.ae = self.info_dict[mh_obj.name]['Ae']
            mh_obj.prob_det_mix = self.info_dict[mh_obj.name]['ProbOfDetMix']
            mh_obj.snps_rsid_list = self.info_dict[mh_obj.name]['RSIDList']

    @staticmethod
    def _transform_data_types(mh_info_dict):
        mh_info_dict['Ae'] = float(mh_info_dict['Ae'])
        mh_info_dict['ProbOfDetMix'] = float(mh_info_dict['ProbOfDetMix'])
        return mh_info_dict
