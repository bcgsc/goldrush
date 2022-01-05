# Detecting hairpins in Nanopore reads

Uses minimizers to identify and flag putative hairpin artifact reads in Nanopore sequencing data

## Dependencies

* [btllib](https://github.com/bcgsc/btllib)  python wrappers
  * Compile wrappers as described in README, then add path to `btllib/python` directory to `PYTHONPATH`
* [statsmodels](https://www.statsmodels.org/stable/index.html)
* scipy
* numpy
* pandas

## Usage

```commandline
usage: detect_hairpins.py [-h] -i INDEX -k K -w W [--perc PERC] [-e E] [--upper_slope UPPER_SLOPE] [--lower_slope LOWER_SLOPE] [-c C] [--corr CORR] [-b BINS]
                          [-m MAPPED_BIN_THRESHOLD] [-o O] [-r R] [-v]
                          FA

Detect hairpin artifacts in nanopore reads

positional arguments:
  FA                    Input fasta file, or '-' if piping to standard in

optional arguments:
  -h, --help            show this help message and exit
  -i INDEX, --index INDEX
                        samtools faidx index for input reads
  -k K                  Kmer size
  -w W                  Window size
  --perc PERC           Percentage error allowed for yintercept [10]
  -e E                  Length of ends to consider (bp) [5000]
  --upper_slope UPPER_SLOPE
                        Upper threshold for slope [-0.75]
  --lower_slope LOWER_SLOPE
                        Lower threshold for slope [-1.25]
  -c C                  Threshold for correlation [-0.75]
  --corr CORR           Correlation coefficient to use. Valid values are pearson or spearman [spearman]
  -b BINS, --bins BINS  Number of bins for minimizer distribution check [10]
  -m MAPPED_BIN_THRESHOLD, --mapped-bin-threshold MAPPED_BIN_THRESHOLD
                        Threshold number of bins with mapped minimizers [5]
  -o O                  Output file for hairpin classifications [stdout]
  -r R                  Path to random forest models
  -v                    Verbose logging of filtered minimizers
```

## Credits

Written by Lauren Coombe