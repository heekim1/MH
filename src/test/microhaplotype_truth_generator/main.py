import pysam

from argument_parser import ArgumentParser
from bed import Bed
from filters import Filters
from info import Info
from microhaplotype import Microhaplotype
from output_microhaplotypes import OutputMicrohaplotypes


class Microhaplotyper(object):
    ANALYTICAL_THRESHOLD = 0.02
    ABSOLUTE_ANALYTICAL_THRESHOLD = 10
    MIN_MAPPING_QUALITY = 60
    MAX_MAPPING_QUALITY = 254
    REMOVE_DELETIONS = True

    def __init__(self, bam, bed, info, out, analytical_threshold=ANALYTICAL_THRESHOLD):

        self.bam = bam
        self.bed = bed
        self.info = info
        self.out = out
        self.analytical_threshold = analytical_threshold

        # Init
        self.bed_obj = Bed(self.bed)
        self.info_obj = Info(self.info)
        self.mh_dict = {}
        self.filters_obj = Filters(analytical_threshold=self.analytical_threshold,
                                   abs_analytical_threshold=self.ABSOLUTE_ANALYTICAL_THRESHOLD,
                                   min_mapq=self.MIN_MAPPING_QUALITY,
                                   max_mapq=self.MAX_MAPPING_QUALITY,
                                   remove_deletions=self.REMOVE_DELETIONS)
        self.short_reads_count = 0
        self.total_reads_count = 0
        self.number_of_contributors_overall = 0

    def get_microhaplotypes(self):
        self._process_bed_file()
        if self.info:
            self._process_info_file()
        self._process_microhaplotypes()

    def _process_bed_file(self):
        self.bed_obj.get_bed_dict()

    def _process_info_file(self):
        self.info_obj.get_info_dict()

    def _process_microhaplotypes(self):
        self._get_mh_count()
        self._get_number_of_contributors_overall()

    def _get_mh_count(self):
        bam_file_obj = pysam.AlignmentFile(self.bam, 'rb')
        for mh_name in self.bed_obj.bed_dict.keys():

            self.mh_dict.setdefault(mh_name, Microhaplotype(name=mh_name))
            mh_number_of_bases_targeted = len(self.bed_obj.bed_dict[mh_name])

            for chromosome, start_position in self.bed_obj.bed_dict[mh_name]:
                self.mh_dict[mh_name].process_overlapping_reads(chromosome, start_position, bam_file_obj)

            reads_processed, short_reads = self.mh_dict[mh_name].get_alleles(mh_number_of_bases_targeted)
            self.total_reads_count += reads_processed
            self.short_reads_count += short_reads
            self.mh_dict[mh_name].calculate_avg_mapping_quality()
            self.mh_dict[mh_name].calculate_total_coverage()
            self.mh_dict[mh_name].filter_alleles_and_set_analytical_threshold(self.filters_obj)
            self.mh_dict[mh_name].calculate_number_of_alleles()
            self.mh_dict[mh_name].calculate_number_of_contributors()
            if self.info:
                self.info_obj.add_info(self.mh_dict[mh_name])

    def _get_number_of_contributors_overall(self):
        self.number_of_contributors_overall = max([self.mh_dict[mh_name].number_of_contributors for mh_name in
                                                   self.mh_dict.keys()])


if __name__ == '__main__':
    args_obj = ArgumentParser()
    mh_obj = Microhaplotyper(bam=args_obj.args.bam_file_path,
                             bed=args_obj.args.bed_file_path,
                             info=args_obj.args.info_file_path,
                             out=args_obj.args.out_file_path,
                             analytical_threshold=args_obj.args.min_coverage)
    mh_obj.get_microhaplotypes()
    output_obj = OutputMicrohaplotypes(mh_obj)
    output_obj.output_microhaplotypes()
