class Allele(object):
    def __init__(self, sequence=''):
        self.sequence = sequence
        self.coverage_total = None
        self.coverage_forward = None
        self.coverage_reverse = None
        self.avg_mapping_quality = None
