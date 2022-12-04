from operator import mul
from functools import reduce
from typing import Tuple
from collections import Counter
from Bio.SeqIO.FastaIO import SimpleFastaParser
import math, logging

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
        for name, seq in SimpleFastaParser(fh):
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