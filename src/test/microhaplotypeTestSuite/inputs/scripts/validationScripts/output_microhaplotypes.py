class OutputMicrohaplotypes:
    def __init__(self, mh_obj):
        self.mh_obj = mh_obj

    def output_microhaplotypes(self):
        # self._output_to_stdout()
        if self.mh_obj.out:
            self._output_to_file()

    def _output_to_file(self):
        with open(self.mh_obj.out, 'w') as mh_out:
            mh_out.write(str(self.mh_obj.filters_obj))
            mh_out.write("#The minimum number of contributors : %d\n" % self.mh_obj.number_of_contributors_overall)
            mh_out.write('MarkerId\tAvgAe\tProbOfDetMix\tNumOfAlleles\tNumOfContributors\tAlleleCoverage\t'
                         'AlleleCoverageDetail(-strand:+strand:avgMQ)\tSNPs\tMinimumCoverage\n')
            [mh_out.write(str(self.mh_obj.mh_dict[mh_name])) for mh_name in self.mh_obj.mh_dict.keys() if
             self.mh_obj.mh_dict[mh_name].ae]

