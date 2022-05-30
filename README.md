![GitHub release (latest by date)]


![Logo](https://github.com/bcgsc/GoldRush/blob/readme/img/GoldRush-logo.png)

# GoldRush: memory-efficient _de novo_ assembly of long reads

## Description of the algorithm
GoldRush iterates through the input long reads to produce a "golden path" of reads comprising ~1-fold coverage of the target genome. These "golden path" reads, or "goldtigs" are then polished, corrected for misassemblies, and finally scaffolded to generate the final genome assembly. GoldRush is memory-efficient, and has a linear time complexity in the number of reads.

### General steps in the algorithm:
1. GoldRush-Path: selecting the golden path reads
2. GoldRush-Edit: polishing the genome
3. Tigmint-long: correcting the genome
4. GoldRush-Link: scafolding the genome


## Credits
Concept: Inanc Birol and Rene L. Warren

Design and implementation: Johnathan Wong, Vladimir Nikolic, and Lauren Coombe

Logo Design: Rene L. Warren

## Presentations
Wong, J., Nikolic, V., Coombe, L., Warren, R., & Birol, I. (2022, July 10–14). GoldRush-Path: A de novo assembler for long reads with linear time complexity [Conference presentation]. Intelligent Systems for Molecular Biology 2022, Madison, WI, United States.  

Nikolic, V., Coombe, L., Wong, J., Birol, I., & Warren, R. (2022, July 10–14). GoldRush-Edit : A targeted, alignment-free polishing & finishing pipeline for long read assembly, using long read k-mers [Conference presentation]. Intelligent Systems for Molecular Biology 2022, Madison, WI, United States. 

Coombe, L., Warren, R., Nikolic, V., Wong, J., & Birol, I. (2022, July 10–14). GoldRush-Link: Integrating minimizer-based overlap detection and gap-filling to the ntLink long read scaffolder [Conference presentation]. Intelligent Systems for Molecular Biology 2022, Madison, WI, United States. 

## Usage
```
GoldRush

Usage: ./goldrush [COMMAND] [OPTION=VALUE]…

For example, to run the default pipeline on reads reads.fa.gz and a genome size of gsize:
goldrush run reads=reads G=gsize

        Commands:

        run                   run default GoldRush pipeline: GoldRush-Path + Polisher (GoldRush-Edit by default) + Tigmint + ntLink

        goldrush-path             run GoldRush-Path
        path-polish           run GoldRush-Path, then goldrush-edit
        path-ntLink           run GoldRush-Path, then goldrush-edit, then ntLink
        path-tigmint          run GoldRush-Path, then goldrush-edit, then tigmint
        path-tigmint-ntLink   run GoldRush-Path, then goldrush-edit, then tigmint, then ntLink (default 5 rounds)
        path-tigmint-ntJoin   run GoldRush-Path, then goldrush-edit, then tigmint, then ntJoin

        General options (required):
        reads                   read name [reads]. File must have .fq.gz or .fa.gz extension
        G                       haploid genome size (bp) (e.g. '3e9' for human genome)

        General options (optional):
        t                       number of threads for scaffolding and correction tools [8]
        z                       minimum size of contig (bp) to scaffold [1000]
        track_time              If 1 then track the run time and memory usage, if 0 then don't [0]

        GoldRush-Path options:
        k                       base k value to generated hash [22]
        w                       weight of spaced seed (number of 1's) [16]
        tile                    tile size to use in GoldRush-Path [1000]
        u                       minimum number of unassigned tiles for the read to be considered unassigned [5]
        a                       maximum number of tiles that can be assigned, minimum number of overlapping tiles kept after trimming [5]
        l                       number fo golden paths to produce [1]
        o                       occupancy of the miBF [0.1]
        x                       threshold for number of hits in miBF for a given frame to be considered assigned [10]
        h                       the number of seed patterns to use [3]
        j                       number of threads to use [48]
        s                       spaced seed design [1011011110110111101101]
        M                       maximum number of silver paths to generate [5]
        r                       ratio of full genome in golden path [0.9]
        p1                      prefix to use for the silver paths [w16_x10]

        Tigmint options:
        span                    min number of spanning molecules to be considered correctly assembled [2]
        dist                    maximum distance between alignments to be considered the same molecule [1000]

        ntLink options:
        k_ntLink                k-mer size for minimizers [64]
        w_ntLink                window size for minimizers [1000]
        rounds                  number of rounds of ntLink [5]

        ntJoin options:
        k_ntJoin                k-mer size for minimizers [24]
        w_ntJoin                window size for minimizers [50]
        no_cut                  If True, will not cut contigs at putative misassemblies [False]
        ref_weight              Weight of the reference assemblies [2]

Notes:
        - GoldRush-Path runs the new protocol which generates silver paths before generating the golden path
        - ntJoin runs with the two silver paths instead of the full five
        - Ensure that all input files are in the current working directory, making soft-links if needed
```

Running `goldrush help` prints the help documentation.

* Input reads files can be gzipped(not yet implemented) or not, and in fastq 

### Example
Input files:
* long read file `long_reads.fq`
Required Parameters:
* genome size `3e9`

GoldRush command:
```
goldrush run reads=long_reads G=3e9
```

See the wiki page for more details.

**For more information about the GoldRush algorithm and tips for running GoldRush see our [wiki](https://github.com/bcgsc/goldrush/wiki)** 

 ## Installation
 GoldRush is available from github.
  
 Installing from source code:
 ```
  git clone https://github.com/bcgsc/goldrush.git
  cd goldrush
  meson build
  cd build
  ninja
 ```

## License
GoldRush Copyright (c) 2022 British Columbia Cancer Agency Branch. All rights reserved.

GoldRush is released under the GNU General Public License v3

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, version 3.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.

For commercial licensing options, please contact Patrick Rebstein (prebstein@bccancer.bc.ca).
