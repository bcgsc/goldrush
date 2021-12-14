# Detecting hairpins in Nanopore reads

Uses minimizers to identify and flag putative hairpin artifact reads in Nanopore sequencing data

## Dependencies

* [btllib](https://github.com/bcgsc/btllib)  python wrappers
  * Compile wrappers as described in README, then add path to `btllib/python` directory to `PYTHONPATH`
* [statsmodels](https://www.statsmodels.org/stable/index.html)
* scipy
* numpy

## Usage

```commandline
usage: detect_hairpins.py [-h] -i INDEX [--perc PERC] [-e E] [--upper_slope UPPER_SLOPE] [--lower_slope LOWER_SLOPE] [-c C] [--corr CORR] [-b BINS]
                          [-m MAPPED_BIN_THRESHOLD] [-o O] [-r R] [-v] -k K -w W
                          FA

Detect hairpin artifacts in nanopore reads

positional arguments:
  FA                    Input fasta file, or '-' if piping to standard in

optional arguments:
  -h, --help            show this help message and exit
  -i INDEX, --index INDEX
                        samtools faidx index for input reads
  --perc PERC           Percentage error allowed for yintercept
  -e E                  Length of ends to consider
  --upper_slope UPPER_SLOPE
                        Upper threshold for slope
  --lower_slope LOWER_SLOPE
                        Lower threshold for slope
  -c C                  Threshold for correlation
  --corr CORR           Correlation coefficient to use. Valid values are pearson or spearman
  -b BINS, --bins BINS  Number of bins for minimizer distribution check
  -m MAPPED_BIN_THRESHOLD, --mapped-bin-threshold MAPPED_BIN_THRESHOLD
                        Threshold number of bins with mapped minimizers
  -o O                  Output file for hairpin classifications [stdout]
  -r R                  Path to random forest models
  -v                    Verbose logging of filtered minimizers
  -k K                  Kmer size
  -w W                  Window size
```

## Credits

Written by Lauren Coombe