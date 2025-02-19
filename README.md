![GitHub release (latest by date)](https://img.shields.io/github/v/release/bcgsc/goldrush)
[![Downloads](https://img.shields.io/github/downloads/bcgsc/goldrush/total?logo=github)](https://github.com/bcgsc/goldrush/releases)
![Conda](https://img.shields.io/conda/dn/bioconda/goldrush?label=Conda)

![Logo](https://github.com/bcgsc/GoldRush/blob/main/img/GoldRush-logo.png)

# GoldRush: memory-efficient _de novo_ assembly of long reads

## Description of the algorithm
GoldRush iterates through the input long reads to produce a "golden path" of reads comprising ~1-fold coverage of the target genome. These "golden path" reads, or "goldtigs", are then polished, corrected for misassemblies, and finally scaffolded to generate the final genome assembly. GoldRush is memory-efficient, and has a linear time complexity in the number of reads.


### General steps in the algorithm:
1. **GoldPath** (aka GoldRush-Path): selecting the golden path reads
2. **[GoldPolish](https://github.com/bcgsc/goldpolish)** (aka GoldRush-Edit): polishing the genome
3. **[Tigmint-long](https://github.com/bcgsc/tigmint)**: correcting the genome
4. **[GoldChain](https://github.com/bcgsc/ntlink)** (aka GoldRush-Link): scaffolding the genome
5. **[GoldPolish-Target](https://github.com/bcgsc/goldpolish)**: targeted polishing of the genome


## Credits
Concept: Inanc Birol and Rene L. Warren

Design and implementation: Johnathan Wong, Vladimir Nikolic, and Lauren Coombe

Logo Design: Rene L. Warren

## Usage
```
GoldRush
v1.2.2

Usage: goldrush [COMMAND] [OPTION=VALUE]…

For example, to run the default pipeline on reads 'reads.fq' with a genome size of gsize:
goldrush run reads=reads G=gsize

        Commands:

        run                             run default GoldRush pipeline: GoldRush-Path + Polisher (GoldPolish by default) + Tigmint-long + ntLink (default 5 rounds) + GoldPolish-Target
        goldrush-path                   run GoldRush-Path
        path-polish                     run GoldRush-Path, then GoldPolish
        path-tigmint                    run GoldRush-Path, then GoldPolish, then Tigmint-long
        path-tigmint-ntLink             run GoldRush-Path, then GoldPolish, then Tigmint-long, then ntLink (default 5 rounds)
        path-tigmint-ntLink-target      run GoldRush-Path, then GoldPolish, then Tigmint-long, then ntLink (default 5 rounds), then GoldPolish-Target

        General options (required):
        reads                           read name [reads]. File must have .fq or .fastq extension, but do not include the suffix in the supplied read name
        G                               haploid genome size (bp) (e.g. '3e9' for human genome)

        General options (optional):
        t                               number of threads [48]
        z                               minimum size of contig (bp) to scaffold [1000]
        track_time                      If 1 then track the run time and memory usage, if 0 then don't [0]

        GoldRush-Path options:
        k                               base k value to generate hash [22]
        w                               weight of spaced seed (number of 1's) [16]
        tile                            tile size [1000]
        b                               during insertion, number of consecutive tiles to be inserted with the same ID [10]
        u                               minimum number of unassigned tiles for the read to be considered unassigned [5]
        a                               maximum number of tiles that can be assigned, minimum number of overlapping tiles kept after trimming [1]
        o                               occupancy of the miBF [0.1]
        x                               threshold for number of hits in miBF for a given frame to be considered assigned [10]
        h                               number of seed patterns to use [3]
        m                               minimum read length [20000]
        M                               maximum number of silver paths to generate [5]
        r                               ratio of full genome in golden path [0.9]
        P                               minimum average phred score for each read, or 0 to calculate threshold automatically [0]
        d                               remove reads with greater or equal than d difference between average phred quality of first half and second half of the read [5]
        p                               prefix to use for the output paths [goldrush_asm]

        Tigmint-long options:
        span                            min number of spanning molecules [2]
        dist                            maximum distance between alignments to be considered the same molecule [500]

        ntLink options:
        k_ntLink                        k-mer size for minimizers [40]
        w_ntLink                        window size for minimizers [250]
        rounds                          number of rounds of ntLink [5]

        GoldPolish options:
        shared_mem		        Shared memory path where polishing occurs [/dev/shm]

Notes:
        - GoldRush-Path generates silver paths before generating the golden path
        - Ensure that all input files are in the current working directory, making soft-links if needed
        - The input reads must be in random order. If they are sorted by chromosome position, shuffle them prior to GoldRush assembly.
```

Running `goldrush help` prints the help documentation.

* Input reads files must be in uncompressed fastq format


### Example
Input files:
* long read file `long_reads.fq`

Required Parameters:
* genome size `3e9`

GoldRush command:
```
goldrush run reads=long_reads G=3e9
```

**For more information about the GoldRush algorithm and tips for running GoldRush see our [wiki](https://github.com/bcgsc/goldrush/wiki)** 

# System Requirements

## Hardware Requirements

GoldRush does not require any specialized hardware. For the assembly of a human genome with ~60-fold coverage, GoldRush requires a computer with at least 64 GB of RAM. We recommend running GoldRush with at least 48 threads.

The runtime and RAM usage benchmarks below were generated using a server-class system (144 Intel(R) Xeon(R) Gold 6254 CPU @ 3.1 GHz with 2.9 TB RAM) with 48 threads specified on three different *H. sapiens* Oxford Nanopore Technology genomic long read datasets.

| Coverage       | Time (h)    | RAM (GB)    |
| ----------- | ----------- | ----------- |
| 67X  | 16.6       | 51.9        |
| 63X  | 20.8       | 53.9        |
| 71X  | 20.8       | 54.5        |

## Software Requirements

### OS Requirements

GoldRush has been tested on *Linux* operating systems (centOS7, ubuntu-20.04)

### Dependencies
 * [GCC 7+](https://gcc.gnu.org/) with [OpenMP](https://www.openmp.org/)
 * [python 3.9+](https://www.python.org/)
 * [zlib](https://zlib.net/)
 * [meson](https://mesonbuild.com/Getting-meson.html)
 * [ninja](https://github.com/ninja-build/ninja/)
 * [tcmalloc](https://google.github.io/tcmalloc/quickstart.html)
 * [sdsl-lite](https://github.com/simongog/sdsl-lite)
 * [boost](https://www.boost.org/doc/libs/1_61_0/more/getting_started/unix-variants.html)
 * [libdivsufsort](https://github.com/y-256/libdivsufsort)
 * [sparsehash](https://github.com/sparsehash/sparsehash)
 * [btllib 1.6.2+](https://github.com/bcgsc/btllib)
 * [Tigmint 1.2.6+](https://github.com/bcgsc/tigmint)
 * [ntLink 1.3.3+](https://github.com/bcgsc/ntlink)
 * [minimap2](https://github.com/lh3/minimap2)
 * [snakemake](https://github.com/snakemake/snakemake)
 * [intervaltree](https://github.com/chaimleib/intervaltree) 

## Installation
### Installing using conda:
```
conda install -c bioconda -c conda-forge goldrush
```

### Installing from source code:

#### Github repository main branch
 ```
  git clone https://github.com/bcgsc/goldrush.git
  cd goldrush
  git submodule init
  git submodule update
  meson --prefix /path/to/install build
  cd build
  ninja install
 ```
#### Downloading release tarball

 ```
  tar -xf goldrush-x.y.z.tar.xz
  cd goldrush-x.y.z
  meson --prefix /path/to/install build
  cd build
  ninja install
 ```
Compiling GoldRush from the source code takes ~2.5min on a typical machine.

### Testing Installation

 ```
 goldrush help
 cd tests
 ./goldrush_test_demo.sh
 ```
Running the above `goldrush_test_demo.sh` script will automatically download a small set of long reads from a ~1Mbp segment of *C. elegans* chromosome 3, and run GoldRush. It will also check that the final assembly file has an L50 of 1, as expected (ie. at least half of the assembly is in a single piece). The test should run in <2min on a typical machine.

The final assembly file for the test demo can be found in: `goldrush_test_golden_path.goldpolish-polished.span2.dist500.tigmint.fa.k40.w250.ntLink-5rounds.fa`

## Citation
If you use GoldRush in your research, please cite:

Wong J, Coombe L, Nikolić V, Zhang E, Nip KM, Sidhu P, Warren RL and Birol I (2023). Linear time complexity de novo long read genome assembly with GoldRush. Nature Communications, 14(1), 2906. https://doi.org/10.1038/s41467-023-38716-x

[![link](https://img.shields.io/badge/GoldRush-manuscript-brightgreen)](https://doi.org/10.1038/s41467-023-38716-x)

## Presentations
Wong, J., Nikolic, V., Coombe, L., Zhang, E., Warren, R., & Birol, I. (2022, July 10–14). GoldRush-Path: A de novo assembler for long reads with linear time complexity [Conference presentation]. Intelligent Systems for Molecular Biology 2022, Madison, WI, United States.  

Nikolic, V., Coombe, L., Wong, J., Birol, I., & Warren, R. (2022, July 10–14). GoldRush-Edit : A targeted, alignment-free polishing & finishing pipeline for long read assembly, using long read k-mers [Conference presentation]. Intelligent Systems for Molecular Biology 2022, Madison, WI, United States. 

Coombe, L., Warren, R., Nikolic, V., Wong, J., & Birol, I. (2022, July 10–14). GoldRush-Link: Integrating minimizer-based overlap detection and gap-filling to the ntLink long read scaffolder [Conference presentation]. Intelligent Systems for Molecular Biology 2022, Madison, WI, United States. 

## License
GoldRush Copyright (c) 2022 British Columbia Cancer Agency Branch. All rights reserved.

GoldRush is released under the GNU General Public License v3

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, version 3.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.

For commercial licensing options, please contact Patrick Rebstein (prebstein@bccancer.bc.ca).
