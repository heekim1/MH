import abc


class File(object):
    __metaclass__ = abc.ABCMeta

    def __init__(self, file_path, file_type, unavailable_data):
        self.file_path = file_path
        self.file_type = file_type
        self.unavailable_data = set(unavailable_data)
        self.number_of_contributors = None
        self.mh_dict = {}

    @abc.abstractmethod
    def parse(self):
        pass
