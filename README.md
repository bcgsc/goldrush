![GitHub release (latest by date)](https://img.shields.io/github/v/release/bcgsc/goldrush)
![Conda](https://img.shields.io/conda/dn/bioconda/goldrush?label=Conda)

![Logo](https://github.com/bcgsc/GoldRush/blob/main/img/GoldRush-logo.png)

# GoldRush: memory-efficient _de novo_ assembly of long reads

## Description of the algorithm
GoldRush iterates through the input long reads to produce a "golden path" of reads comprising ~1-fold coverage of the target genome. These "golden path" reads, or "goldtigs", are then polished, corrected for misassemblies, and finally scaffolded to generate the final genome assembly. GoldRush is memory-efficient, and has a linear time complexity in the number of reads.


### General steps in the algorithm:
1. GoldRush-Path: selecting the golden path reads
2. GoldRush-Edit: polishing the genome
3. Tigmint-long: correcting the genome
4. GoldRush-Link: scaffolding the genome



## Credits
Concept: Inanc Birol and Rene L. Warren

Design and implementation: Johnathan Wong, Vladimir Nikolic, and Lauren Coombe

Logo Design: Rene L. Warren

## Presentations
Wong, J., Nikolic, V., Coombe, L., Zhang, E., Warren, R., & Birol, I. (2022, July 10–14). GoldRush-Path: A de novo assembler for long reads with linear time complexity [Conference presentation]. Intelligent Systems for Molecular Biology 2022, Madison, WI, United States.  

Nikolic, V., Coombe, L., Wong, J., Birol, I., & Warren, R. (2022, July 10–14). GoldRush-Edit : A targeted, alignment-free polishing & finishing pipeline for long read assembly, using long read k-mers [Conference presentation]. Intelligent Systems for Molecular Biology 2022, Madison, WI, United States. 

Coombe, L., Warren, R., Nikolic, V., Wong, J., & Birol, I. (2022, July 10–14). GoldRush-Link: Integrating minimizer-based overlap detection and gap-filling to the ntLink long read scaffolder [Conference presentation]. Intelligent Systems for Molecular Biology 2022, Madison, WI, United States. 

## News
#### October 31, 2022
We're looking for C++ bioinformatics developers to work on GoldRush and other exciting bioinformatics tools with the Bioinformatics Technology Lab in Vancouver, BC! For more information, check out the [BCGSC website](https://bcgsc.ca/careers/research-programmer-birol-labbioinformatics-technology-lab) and our [Birol Lab website](http://www.birollab.ca/jobs/2022/08/29/ProgrammerJob.html).   

## Usage
```
GoldRush
v1.0.2

Usage: goldrush [COMMAND] [OPTION=VALUE]…

For example, to run the default pipeline on reads 'reads.fq' with a genome size of gsize:
goldrush run reads=reads G=gsize

        Commands:

        run                     run default GoldRush pipeline: GoldRush-Path + Polisher (GoldRush-Edit by default) + Tigmint-long + ntLink (
default 5 rounds)
        goldrush-path           run GoldRush-Path
        path-polish             run GoldRush-Path, then GoldRush-Edit
        path-tigmint            run GoldRush-Path, then GoldRush-Edit, then Tigmint-long
        path-tigmint-ntLink     run GoldRush-Path, then GoldRush-Edit, then Tigmint-long, then ntLink (default 5 rounds)

        General options (required):
        reads                   read name [reads]. File must have .fq or .fastq extension, but do not include the suffix in the supplied read name
        G                       haploid genome size (bp) (e.g. '3e9' for human genome)

        General options (optional):
        t                       number of threads [48]
        z                       minimum size of contig (bp) to scaffold [1000]
        track_time              If 1 then track the run time and memory usage, if 0 then don't [0]

        GoldRush-Path options:
        k                       base k value to generate hash [22]
        w                       weight of spaced seed (number of 1's) [16]
        tile                    tile size [1000]
        b                       during insertion, number of consecutive tiles to be inserted with the same ID [10]
        u                       minimum number of unassigned tiles for the read to be considered unassigned [5]
        a                       maximum number of tiles that can be assigned, minimum number of overlapping tiles kept after trimming [1]
        o                       occupancy of the miBF [0.1]
        x                       threshold for number of hits in miBF for a given frame to be considered assigned [10]
        h                       number of seed patterns to use [3]
        s                       spaced seed design [1011011110110111101101]
        m                       minimum read length [20000]
        M                       maximum number of silver paths to generate [5]
        r                       ratio of full genome in golden path [0.9]
        P                       minimum average phred score for each read [15]
        d                       remove reads with greater or equal than d difference between average phred quality of first half and second half of the read [5]
        p                       prefix to use for the output paths [w16_x10]

        Tigmint-long options:
        span                    min number of spanning molecules [2]
        dist                    maximum distance between alignments to be considered the same molecule [500]

        ntLink options:
        k_ntLink                k-mer size for minimizers [40]
        w_ntLink                window size for minimizers [250]
        rounds                  number of rounds of ntLink [5]

        GoldRush-Edit options:
        polisher_mapper         Whether to use ntlink or minimap2 for mappings [minimap2]

Notes:
        - GoldRush-Path generates silver paths before generating the golden path
        - Ensure that all input files are in the current working directory, making soft-links if needed

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

## Installation
### Installing using conda:
```
conda install -c bioconda goldrush
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
### Testing Installation

 ```
 goldrush help
 cd tests
 ./goldrush_test_demo.sh
 ```
 
## Dependencies
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
 * [btllib](https://github.com/bcgsc/btllib)
 * [Tigmint](https://github.com/bcgsc/tigmint)
 * [ntLink 1.3.3+](https://github.com/bcgsc/ntlink)
 * [minimap2](https://github.com/lh3/minimap2)

## Citation
If you use GoldRush in your research, please cite:

Johnathan Wong, Lauren Coombe, Vladimir Nikolić, Emily Zhang, Ka Ming Nip, Puneet Sidhu, René L Warren, and Inanç Birol. 2022. ‘GoldRush: A *de novo* Long Read Genome Assembler with Linear Time Complexity’. BioRxiv. https://doi.org/10.1101/2022.10.25.513734.


[![link](https://img.shields.io/badge/GoldRush-manuscript-brightgreen)](https://doi.org/10.1101/2022.10.25.513734)

## License
GoldRush Copyright (c) 2022 British Columbia Cancer Agency Branch. All rights reserved.

GoldRush is released under the GNU General Public License v3

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, version 3.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.

For commercial licensing options, please contact Patrick Rebstein (prebstein@bccancer.bc.ca).
