![GitHub release (latest by date)](PH)
![Conda](PH)


![Logo](https://github.com/bcgsc/GoldRush/blob/readme/img/GoldRush-logo.png)

# Long Read Genome Assembler

## Memory-efficient genome assembly with a linear time complexity in number of reads 

## Description of the algorithm
GoldRush iterates through the reads and collect reads to constitute ~1X of the genome. These reads are then polished, miassembly corrected, and scaffolded. 

### General steps in the algorithm:
1. GR-Path: assembling the genome
    1. Iterates through the reads to generate a Multi-index Bloom Filter with no IDs inserted
    2. Iterates through the reads and insert reads that are not found in the MiBF
![alt text](https://github.com/bcgsc/GoldRush/blob/readme/img/GR-Path.png)
2. GR-Edit: polishing the genome
    1. Identify regions of the assembly that requires polishing using ntEdit
    2. Map reads to these erroneous regions
    3. Correct these regions individually with their mapped reads using Sealer
![alt text](https://github.com/bcgsc/GoldRush/blob/readme/img/GR-Edit.png)
3. Tigmint-long: correcting the genome
    1. Break regions in the genome unsupported by other reads
4. GR-Link: scafolding the genome
    1. Scaffold the genome using minimizers
    2. Overlapping regions are trimmed by mapping minizers of the overlaped region together
    3. Fill in gaps by identifying reads that support the scaffold joins
![alt text](https://github.com/bcgsc/GoldRush/blob/readme/img/GR-Link.png)


## Credits
Concept: Inanc Birol

Design and implementation: Johnathan Wong, Vladimir Nikolic, and Lauren Coombe

Logo Design: Rene L. Warren

## Publications
Wong, J., Nikolic, V., Coombe, L., Warren, R., & Birol, I. (2022, July 10–14). GoldRush-Path: A de novo assembler for long reads with linear time complexity [Conference presentation]. Intelligent Systems for Molecular Biology 2022, Madison, WI, United States. PH
Nikolic, V., Coombe, L., Wong, J., Birol, I., & Warren, R. (2022, July 10–14). GoldRush-Edit : A targeted, alignment-free polishing & finishing pipeline for long read assembly, using long read k-mers [Conference presentation]. Intelligent Systems for Molecular Biology 2022, Madison, WO, United States. PH
Coombe, L., Warren, R., Nikolic, V., Wong, J., & Birol, I. (2022, July 10–14). GoldRush-Link: Integrating minimizer-based overlap detection and gap-filling to the ntLink long read scaffolder [Conference presentation]. Intelligent Systems for Molecular Biology 2022, Madison, WO, United States. PH

## Citing GoldRush
If you use GoldRush in your research, please cite:

PH

## Usage
```
GoldRush

Usage: ./goldrush [COMMAND] [OPTION=VALUE]…

For example, to run the default pipeline on reads reads.fa.gz and a genome size of gsize:
goldrush run reads=reads G=gsize

        Commands:

        run                     run default GoldRush pipeline: GoldRush-Path + Racon + Tigmint + ntLink

        goldrush-path           run GoldRush-Path
        path-racon              run GoldRush-Path, then racon
        path-ntLink             run GoldRush-Path, then racon, then ntLink
        path-tigmint            run GoldRush-Path, then racon, then tigmint
        path-tigmint-ntLink     run GoldRush-Path, then racon, then tigmint, then ntLink (default 5 rounds)
        path-tigmint-ntJoin     run GoldRush-Path, then racon, then tigmint, then ntJoin

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
* long read file `long_reads.fq.gz`
* genome size `3e9`

GoldRush command:
```
GoldRush run reads=long_reads G=3e9
```

The final assembly will have the suffix `*goldrush.assembly.fa` (will implement this)


See the wiki page for more details.

**For more information about the GoldRush algorithm and tips for running GoldRush see our [wiki](https://github.com/bcgsc/goldrush/wiki)** (will make this])

 ## Installation
 GoldRush is available from conda and homebrew package managers.
 
 Installing using conda:
 ```
 conda install -c bioconda goldrush #PH
 ```
 
 Installing using brew:
 ```
 brew install brewsci/bio/goldrush #PH
 ```
 
 Installing from source code:
 ```
  git clone https://github.com/bcgsc/goldrush.git
  cd goldrush
  meson build
  cd build
  ninja
 ```

#### Testing your installation
To test your GoldRush installation: #add tests
```
cd tests
./test_installation.sh
```
The expected output files can be found in: `tests/expected_outputs`

## Dependencies
will compile a list of dependencies

Python dependencies can be installed with:
```
pip3 install -r requirements.txt
```

## License
GoldRush Copyright (c) 2020 British Columbia Cancer Agency Branch. All rights reserved.

GoldRush is released under the GNU General Public License v3

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, version 3.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.

For commercial licensing options, please contact Patrick Rebstein (prebstein@bccancer.bc.ca).
