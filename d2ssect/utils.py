from operator import mul
from functools import reduce
from typing import Tuple
from collections import Counter
import math, logging

import multiprocessing, subprocess
from itertools import combinations
import numpy as np
import math, logging, os, argparse, sys
from os.path import exists


def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line[1:], []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))


def get_seqinfo(seqfile: str):
    """
    parse sequence files and generate info of seq length, number of seqs, and frequencies of nucleotides
    """
    logging.info(f'parse seqinfo for {seqfile}')

    #if seqfile.endswith(('fq','fastq')):
    #    f_type = 'fastq'
    #elif seqfile.endswith(('fa','fasta')):
    #    f_type = 'fasta'
    #else:
    #    f_type = 'fasta'
    #    logging.warning('Could not determine the seqfile type, will parse as fasta file')
    
    n_seq = 0
    char_counter = Counter()
    with open(seqfile, 'r') as fh:
        for name, seq in read_fasta(fh):
            n_seq += 1
            char_counter += Counter(seq)
    
    char_counter.pop('N', None)
    total_len = sum(char_counter.values())
    char_freq = {k: v/total_len for k, v in char_counter.items()}

    logging.debug(f"charfreq of {seqfile}: A:{char_freq['A']}; C:{char_freq['C']}; G:{char_freq['G']}; T:{char_freq['T']}; n_seq: {n_seq}; total_len: {total_len}")

    return n_seq, total_len, char_freq

def generate_matrix(d2s_combinations_list, n_sample):
    d2s_matrix = np.zeros((n_sample,n_sample),dtype='float')
    triu = np.triu_indices(n_sample,k=1)
    tril = np.tril_indices(n_sample, -1)
    d2s_matrix[triu] = d2s_combinations_list
    d2s_matrix[tril] = d2s_matrix.T[tril]
    
    return d2s_matrix