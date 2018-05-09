#### Execute Nightly Test
This is a set of scripts that will depend on output from:-
- microhaplotypeTestSuite - To get scenario bam files for execution
- microhaplotype_truth_generator - To get the expected results for a given simulated or real bam file
- compare_microhaplotype_outputs - To compare actual jar output with the expected one from microhaplotype_truth_generator.
This module can compare any of the two output formats for microhaplotype - txt, json, xlsx(ui) with each other.

### Usage
python execute_nightly_test.py --scenariosdir /path/to/scenarios/dir/from/microhaplotypeTestSuite
                 --mhjarlocation /path/to/microhaplotyper.jar --outdir /path/to/output/dir/