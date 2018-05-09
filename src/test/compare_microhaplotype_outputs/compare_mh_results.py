import logging
import sys
from deepdiff import DeepDiff

from argument_parser import ArgumentParser
from comparison_details import ComparisonDetails
from json_file import JsonFile
from txt_file import TxtFile
from xlsx_file import XlsxFile


class CompareMHResults(object):
    def __init__(self, file1_path, file2_path):
        self.file1_path = file1_path
        self.file2_path = file2_path
        self.file1_obj = self._get_file_obj(self.file1_path)
        self.file1_obj.parse()
        self.file2_obj = self._get_file_obj(self.file2_path)
        self.file2_obj.parse()
        self.excluded_comparison_types = self.file1_obj.unavailable_data.union(self.file2_obj.unavailable_data)
        self.no_of_contributors_comparison = ComparisonDetails()
        self.mh_results_comparison = ComparisonDetails()

    @staticmethod
    def _get_file_obj(file_path):
        if file_path.lower().endswith('.json'):
            return JsonFile(file_path)
        elif file_path.lower().endswith('.txt'):
            return TxtFile(file_path)
        elif file_path.lower().endswith('.xlsx'):
            return XlsxFile(file_path)
        else:
            logging.error("Unknown file type for file %s" % file_path)
            sys.exit(1)

    def compare(self):
        if 'number_of_contributors_total' not in self.excluded_comparison_types:
            self._compare_number_of_contributors()
        self._compare_mh_results()

    def _compare_number_of_contributors(self):
        if self.file1_obj.number_of_contributors != self.file2_obj.number_of_contributors:
            self.no_of_contributors_comparison.flag = "FAIL"
            self.no_of_contributors_comparison.differences = '%d|%d' % (self.file1_obj.number_of_contributors,
                                                                        self.file2_obj.number_of_contributors)
        else:
            self.no_of_contributors_comparison.flag = "PASS"

    def _compare_mh_results(self):
        mh_dict_differences = DeepDiff(self.file1_obj.mh_dict, self.file2_obj.mh_dict)
        mh_dict_differences = self._filter_type_differences(mh_dict_differences)
        if mh_dict_differences:
            self.mh_results_comparison.flag = "FAIL"
            self.mh_results_comparison.differences = str(mh_dict_differences)
        else:
            self.mh_results_comparison.flag = "PASS"

    def _filter_type_differences(self, mh_dict_differences):
        if 'type_changes' in mh_dict_differences:
            post_filtering_changes = {}
            for change_key in mh_dict_differences['type_changes']:
                change_type = change_key.split('.')[-1]
                if change_type not in self.excluded_comparison_types:
                    post_filtering_changes[change_key] = mh_dict_differences['type_changes'][change_key]
            if post_filtering_changes:
                mh_dict_differences['type_changes'] = post_filtering_changes
            else:
                del mh_dict_differences['type_changes']

        return mh_dict_differences

    def __str__(self):
        return ','.join([self.no_of_contributors_comparison.flag, self.no_of_contributors_comparison.differences,
                         self.mh_results_comparison.flag, self.mh_results_comparison.differences])


if __name__ == '__main__':
    args_obj = ArgumentParser()
    compare_obj = CompareMHResults(file1_path=args_obj.args.file1_path, file2_path=args_obj.args.file2_path)
    compare_obj.compare()
    print compare_obj
