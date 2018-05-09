class Microhaplotype(object):
    def __init__(self, name=''):
        self.name = name
        self.alleles = {}
        self.ae = None
        self.prob_det_mix = None
        self.snps_rsid_list = None
        self.number_of_alleles = None
        self.number_of_contributors = None
        self.analytical_threshold = None
