# d2ssect
A tool to calculate d2s scores using short fastq reads
This repo will test and benchmark the existing [alignment-free tools](https://github.com/chanlab-genomics/alignment-free-tools) and the improving versions.

The originally version of this pipeline including three big steps:
1. get jellyfish count results
2. calculate d2s using jellyfish dump results of every pair of samples
3. generate a matrix

Our goal is to integrate these three steps and try to increase the speed of d2s calculation.
