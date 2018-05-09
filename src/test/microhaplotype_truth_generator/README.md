## Microhaplotype truth generator
This script will generate the expected truth for given bam files. This can be used to generate the expected output file.

## Usage.
```
from microhaplotype_truth_generator.main import Microhaplotyper
from microhaplotype_truth_generator.output_microhaplotypes import OutputMicrohaplotypes

        truth_obj = Microhaplotyper(bam=bamfile,
                                    bed=bedfile,
                                    info=infofile,
                                    analytical_threshold=min_coverage),
                                    out=output_truth_file_path)
        truth_obj.get_microhaplotypes()
        output_obj = OutputMicrohaplotypes(truth_obj)
        output_obj.output_microhaplotypes()
```
