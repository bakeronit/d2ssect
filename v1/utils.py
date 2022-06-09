from typing import Tuple, TextIO
from Bio import SeqIO
from collections import Counter
from functools import reduce
from operator import mul
import math

## TODO: add fastq support
def get_seqinfo(seqfile: str) -> Tuple[dict, int]:
    n_seq = 0
    char_counter = Counter()
    with open(seqfile,'r') as fh:
        for record in SeqIO.parse(fh, "fasta"):
            n_seq += 1
            char_counter += Counter(record.seq)
    
    return n_seq, char_counter

def get_num_possible_kmer(kmer: str, seqinfo: Tuple[int, dict]) -> float:
    n_seq, char_counter = seqinfo
    k_len = len(kmer)
    char_counter.pop('N', None) # to get the same char freq as v0, uncomment this line
    total_len = sum(char_counter.values())
    total_num_kmer = total_len - (n_seq * (k_len - 1))
    freq = {k: round(v/total_len,12) for k, v in char_counter.items()}
    p_kmer = reduce(mul, [freq[i] for i in kmer]) # get product of a list, faster

    return total_num_kmer * p_kmer

def d2s_per_kmer(kmer: str, num_x: int, num_y: int, seqinfo_x: Tuple[int, dict], seqinfo_y: Tuple[int, dict]) -> float:
    np_x = get_num_possible_kmer(kmer, seqinfo_x)
    np_y = get_num_possible_kmer(kmer, seqinfo_y)

    norm_x = num_x - np_x
    norm_y = num_y - np_y

    return norm_x * norm_y / math.sqrt((norm_x ** 2 + norm_y ** 2))

#print(d2s_per_kmer("ACGATCACAT", 4, 9, "../fasta/DI-1-1_S6.fasta","../fasta/FI-2-21_S28.fasta"))
