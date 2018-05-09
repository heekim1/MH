from __future__ import division

import math

from allele import Allele
from read import Read


class Microhaplotype(object):
    def __init__(self, name=''):
        self.name = name
        self.reads = {}
        self.alleles = {}
        self.alleles_filtered = {}
        self.ae = None
        self.prob_det_mix = None
        self.snps_rsid_list = None
        self.total_coverage = 0
        self.number_of_alleles = 0
        self.number_of_contributors = 0
        self.analytical_threshold = 0

    def __str__(self):
        return "%s\t%f\t%f\t%d\t%d\t%s\t%s\t%s\t%d\n" % (self.name, self.ae, self.prob_det_mix,
                                                         self.number_of_alleles, self.number_of_contributors,
                                                         self.get_allele_coverage_str(),
                                                         self.get_allele_coverage_str_detailed(), self.snps_rsid_list,
                                                         self.analytical_threshold)

    def process_overlapping_reads(self, chromosome, start_position, bam_file_obj):
        end_position = start_position + 1
        for pileup_column in bam_file_obj.pileup(chromosome, start_position, end_position, truncate=True,
                                                 max_depth=800000):
            self._process_pileup_column(pileup_column)

    def _process_pileup_column(self, pileup_column):
        for pileup_read in pileup_column.pileups:
            read_name = pileup_read.alignment.query_name
            nucleotide = '_' if pileup_read.query_position is None else pileup_read.alignment.query_sequence[
                pileup_read.query_position]
            self.reads.setdefault(read_name, Read(name=read_name))
            self.reads[read_name].mh_sequence += nucleotide
            self.reads[read_name].is_reverse = pileup_read.alignment.is_reverse
            self.reads[read_name].mapping_quality = pileup_read.alignment.mapping_quality

    def get_alleles(self, mh_number_of_bases_targeted):
        ret_short_reads_processed, ret_total_reads_processed = 0, 0

        for read_name in self.reads.keys():
            if len(self.reads[read_name].mh_sequence) == mh_number_of_bases_targeted:
                allele_name = self.reads[read_name].mh_sequence
                self.alleles.setdefault(allele_name, Allele(sequence=allele_name))
                self.alleles[allele_name].coverage_total += 1
                if self.reads[read_name].is_reverse:
                    self.alleles[allele_name].coverage_reverse += 1
                else:
                    self.alleles[allele_name].coverage_forward += 1
                self.alleles[allele_name].total_mapping_quality += self.reads[read_name].mapping_quality
            else:
                ret_short_reads_processed += 1
            ret_total_reads_processed += 1

        return ret_total_reads_processed, ret_short_reads_processed

    def calculate_avg_mapping_quality(self):
        [self.alleles[allele_name].compute_avg_mapping_quality() for allele_name in self.alleles.keys()]

    def calculate_total_coverage(self):
        self.total_coverage = sum([self.alleles[allele_name].coverage_total for allele_name in self.alleles.keys()])

    def filter_alleles_and_set_analytical_threshold(self, filter_obj):
        for allele_name in self.alleles.keys():
            if self._check_allele(self.alleles[allele_name], filter_obj):
                self.alleles_filtered[allele_name] = self.alleles[allele_name]

        self.analytical_threshold = filter_obj.get_analytical_threshold(self.total_coverage)

    def _check_allele(self, allele_obj, filter_obj):
        return filter_obj.check_analytical_threshold(self.total_coverage, allele_obj.coverage_total) and \
               filter_obj.check_absolute_analytical_threshold(allele_obj.coverage_total) and \
               filter_obj.check_min_mapq(allele_obj.avg_mapping_quality) and \
               filter_obj.check_max_mapq(allele_obj.avg_mapping_quality) and \
               filter_obj.check_deletions(allele_obj.sequence)

    def calculate_number_of_alleles(self):
        self.number_of_alleles = len(self.alleles_filtered.keys())

    def calculate_number_of_contributors(self):
        self.number_of_contributors = int(math.ceil(len(self.alleles_filtered.keys()) / 2))

    def get_allele_coverage_str(self):
        return ",".join([self.alleles_filtered[allele_name].get_simple_representation()
                         for allele_name in self.alleles_filtered])

    def get_allele_coverage_str_detailed(self):
        return ",".join([self.alleles_filtered[allele_name].get_detailed_representation()
                         for allele_name in self.alleles_filtered])
