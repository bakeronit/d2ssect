# d2ssect

![conda install test badge](https://github.com/iracooke/d2ssect/actions/workflows/conda.yml/badge.svg)
![linux install test badge](https://github.com/iracooke/d2ssect/actions/workflows/linux.yml/badge.svg)
![macos install test badge](https://github.com/iracooke/d2ssect/actions/workflows/macos.yml/badge.svg)

d2ssect (pronounced dissect) calculates an alignment-free distance between samples based on frequencies of shared kmers. Specifically, it provides a fast implementation of the [D2S statistic](https://www.liebertpub.com/doi/10.1089/cmb.2009.0198) which can be used as a standalone distance measure, or as input to a range of methods (eg see [these tools](https://github.com/chanlab-genomics/alignment-free-tools)) for phylogenetic and network analysis.


## Installation

`d2ssect` is available via [pypi](https://pypi.org/project/d2ssect/).  Installation requires python 3.7 or greater as well as the [jellyfish](https://github.com/gmarcais/Jellyfish) program and libraries.  We recommend installation into a conda environment as follows
```bash
conda create -n d2ssect python=3.7 kmer-jellyfish
conda activate d2ssect
pip install d2ssect
d2ssect -h
```

Alternatively, you may use an existing Jellyfish installation, or install Jellyfish without using conda. If using this method please note that;
	- Jellyfish version 2 is required (Jellyfish 1 will not work)
	- Installation of Jellyfish via linux package managers will not work as this installs the jellyfish binary but not libraries and headers needed by `d2ssect`

Once Jellyfish is installed you should then be able to install `d2ssect` using pip or pip3 as follows
```bash
pip install d2ssect
```


## Usage

Lets say we have a collection of fastq files corresponding to sequencing reads from different samples. We want to compare these with `d2ssect`.  First count kmers in these files using `jellyfish`

```bash
for f in *.fastq;do jellyfish count -m 21 -s 10000000 $f -o ${f%.fastq}.jf ;done
```

Note that the command above will create a corresponding `.jf` file for every `.fastq` file in the current directory. By keeping the base names of the `jf` and `fastq` files identical we can then run `d2ssect` as follows;

```bash
d2ssect -l *.jf -f *.fastq
```

## Outputs

`d2ssect` provides information on progress (sent to stderr) and will eventually produce a matrix of pairwise D2S values (one for each pair of samples) sent to stdout. 



