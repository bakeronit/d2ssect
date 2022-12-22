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

def get_num_possible_kmer(kmer: str, seqinfo: Tuple[dict, int]) -> float:
    n_seq, total_len, char_freq = seqinfo
    
    k_len = len(kmer)   
    n_kmer = total_len - (n_seq * (k_len -1))

    p_kmer = reduce(mul, [char_freq[i] for i in kmer])

    return n_kmer * p_kmer

def d2s_calculation(kmc1, kmc2, np1, np2) -> float:

    norm_kmc1 = kmc1 - np1
    norm_kmc2 = kmc2 - np2

    return norm_kmc1 * norm_kmc2 / math.sqrt((norm_kmc1 ** 2 + norm_kmc2 ** 2))


bin_path = os.path.realpath(__file__)
count_share_in_file = os.path.abspath(f'{bin_path}/../../bin/count_share_in_file')

def cal_d2s_self(jf, seqinfo):
    """
    this step prepare seq frequencies and self-mathcing d2s values for each sample
    """
    logging.info(f'Calculating self-matching d2s score for {jf}')
    d2s_self = 0.0
    cmd = f'jellyfish dump -ct {jf}'
    dump_ct = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1)
    dump_table = dump_ct.stdout

    for line in dump_table:
        line = line.decode().strip()
        kmer = line.split()[0]
        kmc = int(line.split()[1])

        np = get_num_possible_kmer(kmer, seqinfo)
        d2s_self += d2s_calculation(kmc, kmc, np, np)
    
    logging.debug(f'self-matching d2s score for {jf}: {d2s_self}')

    return d2s_self


def cal_d2s_and_transform(jf1, jf2, seqinfo1, seqinfo2,d2s_self1, d2s_self2):
    """
    calculate the d2s of two samples based on shared kmer table,
    do the logorithm transformation to get final distance score
    |ln(Sab/sqrt(Saa*Sbb))|
    """
    logging.info(f'Calculating d2s score for {jf1} vs {jf2}')
    d2s = 0.0
    cmd = f'{count_share_in_file} {jf1} {jf2}'
    kmer_table = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1).stdout
    
    for line in kmer_table:
        line = line.decode().strip()
        kmer = line.split()[0]
        kmc1, kmc2 = [ int(i) for i in line.split()[1:3] ]

        np1 = get_num_possible_kmer(kmer, seqinfo1)
        np2 = get_num_possible_kmer(kmer, seqinfo2)

        d2s += d2s_calculation(kmc1, kmc2, np1, np2)
    
    logging.debug(f'D2S score for {jf1} - {jf2}: {d2s}')
    logging.info(f'Distance for {jf1} - {jf2} is {abs(math.log(d2s/math.sqrt(d2s_self1*d2s_self2)))}')

    return abs(math.log(d2s/math.sqrt(d2s_self1*d2s_self2)))


def generate_matrix(d2s_combinations_list, n_sample):
    d2s_matrix = np.zeros((n_sample,n_sample),dtype='float')
    triu = np.triu_indices(n_sample,k=1)
    tril = np.tril_indices(n_sample, -1)
    d2s_matrix[triu] = d2s_combinations_list
    d2s_matrix[tril] = d2s_matrix.T[tril]
    
    return d2s_matrix