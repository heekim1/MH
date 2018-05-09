from __future__ import division
import logging
import sys

class Allele(object):
    def __init__(self, sequence='', coverage_total=0, coverage_forward=0, coverage_reverse=0, total_mapping_quality=0):
        self.sequence = sequence
        self.coverage_total = coverage_total
        self.coverage_forward = coverage_forward
        self.coverage_reverse = coverage_reverse
        self.total_mapping_quality = total_mapping_quality
        self.avg_mapping_quality = 0.0

    def compute_avg_mapping_quality(self):
        try:
            self.avg_mapping_quality = int(self.total_mapping_quality/self.coverage_total)
        except ZeroDivisionError:
            logging.error('Total coverage is zero for allele call %s' % self.sequence)
            sys.exit(1)

    def get_detailed_representation(self):
        return "[" + self.sequence + "]" + str(self.coverage_total) + "(" + ":".join(map(str, [self.coverage_reverse,
                                                                                               self.coverage_forward,
                                                                                               self.avg_mapping_quality]
                                                                                         )) + ")"

    def get_simple_representation(self):
        return "[" + self.sequence + "]" + str(self.coverage_total)



