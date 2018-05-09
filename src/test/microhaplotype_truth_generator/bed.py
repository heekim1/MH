import csv
import logging
import pybedtools
import sys


class Bed(object):
    def __init__(self, bed_file_path):
        self.bed_file_path = bed_file_path
        self.bed_dict = {}

    def get_bed_dict(self):
        """
        Convert the mh bed file into a dict of the format
        {mh_name:[(tuple of chr, pos)]}
        """
        with open(self.bed_file_path, 'rU') as bed_file_handle:
            for bed_file_line in bed_file_handle:
                if self._bed_header_line_check(bed_file_line):
                    mh_name, chromosome, start_position = self._parse_bed_line(bed_file_line)
                    self.bed_dict.setdefault(mh_name, [])
                    if (chromosome, start_position) not in self.bed_dict[mh_name]:
                        self.bed_dict[mh_name].append((chromosome, start_position))

    @staticmethod
    def _parse_bed_line(bed_line):
        bed_line_list = bed_line.rstrip('\n').split('\t')
        try:
            chromosome = bed_line_list[0]
            start_position = int(bed_line_list[1])
            mh_name = bed_line_list[7]
        except (IndexError, ValueError):
            logging.error('Bed file incorrectly formatted at line %s' % bed_line)
            sys.exit(1)
        return mh_name, chromosome, start_position

    @staticmethod
    def _bed_header_line_check(bed_line):
        """ Checks if given line is a bed header """

        return not (bed_line.startswith('#') or bed_line.startswith('track'))


# ================ Test ==================
if __name__ == '__main__':
    bed_obj = Bed('/Users/daklep/Documents/Projects/8_HID_NGS/2_Data_Simulation/'
                  '13_microhaplotypeTestSuite/inputs/hotspotBedFiles/base/mh38_hotspot_sorted.bed')
    bed_obj.get_bed_dict()
    print bed_obj.bed_dict
