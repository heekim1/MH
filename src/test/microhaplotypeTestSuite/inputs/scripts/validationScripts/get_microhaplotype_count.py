import pysam
from microhaplotype import Microhaplotype


class GetMicrohaplotypeCount:

    def __init__(self):
        pass

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
            self.mh_dict[mh_name].filter_alleles(self.filters_obj)
            if self.info:
                self.info_obj.add_info(self.mh_dict[mh_name])





