#### Compare Microhaplotype Output
Module to compare microhaplotype outputs i.e. txt, json and xlsx with each other and give the difference if any between them.

### Usage
Import the class and given the paths to any two microhaplotype output formats
```
from compare_microhaplotype_outputs.compare_mh_results import CompareMHResults
CompareMHResults('/path/to/truth.txt', '/path/to/converge_results.json')
```