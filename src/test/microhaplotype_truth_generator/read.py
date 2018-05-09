class Read(object):
    def __init__(self, name=None, mh_sequence='', is_reverse=None, mapping_quality=None):
        self.name = name
        self.mh_sequence = mh_sequence
        self.is_reverse = is_reverse
        self.mapping_quality = mapping_quality
