from __future__ import division


class Filters(object):
    def __init__(self, analytical_threshold, abs_analytical_threshold, min_mapq, max_mapq, remove_deletions):
        self.analytical_threshold = analytical_threshold
        self.abs_analytical_threshold = abs_analytical_threshold
        self.min_mapq = min_mapq
        self.max_mapq = max_mapq
        self.remove_deletions = remove_deletions

    def get_analytical_threshold(self, total_coverage):
        return int(self.analytical_threshold * total_coverage)

    def check_analytical_threshold(self, total_coverage, allele_coverage):
        return allele_coverage >= int(self.analytical_threshold * total_coverage)

    def check_absolute_analytical_threshold(self, allele_coverage):
        return allele_coverage >= self.abs_analytical_threshold

    def check_min_mapq(self, allele_mapping_quality):
        return allele_mapping_quality >= self.min_mapq

    def check_max_mapq(self, allele_mapping_quality):
        return allele_mapping_quality <= self.max_mapq

    def check_deletions(self, allele_sequence):
        return False if self.remove_deletions is True and '_' in allele_sequence else True

    def __str__(self):
        return "##Filters : (Analytical Threshold=%f, Absolute Analytical Threshold=%d, " \
               "Mapping Quality Range=[%d,%d], Remove Deletions=%s)\n" % \
               (self.analytical_threshold, self.abs_analytical_threshold, self.min_mapq, self.max_mapq,
                self.remove_deletions)
