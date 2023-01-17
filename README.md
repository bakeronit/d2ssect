# d2ssect

![conda install test badge](https://github.com/iracooke/d2ssect/actions/workflows/conda.yml/badge.svg)
![linux install test badge](https://github.com/iracooke/d2ssect/actions/workflows/linux.yml/badge.svg)
![macos install test badge](https://github.com/iracooke/d2ssect/actions/workflows/macos.yml/badge.svg)

d2ssect (pronounced dissect) calculates an alignment-free distance between samples based on frequencies of shared kmers. It provides a fast implementation d2s statistics used by [the chan lab](https://github.com/chanlab-genomics/alignment-free-tools) for [alignment-free phylogenetic analysis](https://pubmed.ncbi.nlm.nih.gov/33961218/). 



## Installation

`d2ssect` relies heavily on [jellyfish](https://github.com/gmarcais/Jellyfish).  You need the jellyfish program and also the jellyfish libraries.  To check that jellyfish is installed you can do;
```bash
jellyfish --version
```
Which should return a version > 2. In addition, you need the jellyfish libraries and headers. If you installed jellyfish via `conda` or by compiling from source these will be present in the right locations.  If you installed it your linux package manager they probably won't be present. 

If you do not want to use `conda` we recommend installing Jellyfish from source.  Once done you should then be able to install `d2ssect` using pip

```bash
pip3 install d2ssect
```



## Usage

Lets say we have a collection of fasta files corresponding to sequencing reads from samples that we want to compare with `d2ssect`.  First count kmers in these files using `jellyfish`

```bash
for f in *.fasta;do jellyfish count -m 21 -s 10000000 $f -o ${f%.fasta}.jf ;done
```

Note that the command above will create a corresponding `.jf` file for every `.fasta` file in the current directory. By keeping the base names of the `jf` and `fasta` files identical we can then run `d2ssect` as follows;

```bash
python3 ../d2ssect/d2ssect/main.py -l *.jf -f *.fasta
```


## Building from source

```
CC=g++ pip install .
```