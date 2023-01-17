import math, logging
from os.path import exists
from collections import Counter
import re

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

def read_fastq(fp):
    dnare = re.compile("^[NGATCgatcn]")
    name, seq = None, None
    for line in fp:
        line = line.rstrip()
        if line.startswith("@"):
            if name: yield (name, seq)
            name, seq = line[1:], ""
        if dnare.match(line):
            seq=line
    if name: yield (name, seq)  

def get_seqinfo(seqfile: str):
    """
    parse sequence files and generate info of seq length, number of seqs, and frequencies of nucleotides
    """
    logging.info(f'parse seqinfo for {seqfile}')

    seq_reader=None

    if seqfile.endswith(('fq','fastq')):
       f_type = 'fastq'
       seq_reader=read_fastq
       logging.debug("Reading seqfiles as fastq")
    elif seqfile.endswith(('fa','fasta')):
       f_type = 'fasta'
       seq_reader=read_fasta
       logging.debug("Reading seqfiles as fasta")
    else:
       f_type = 'fasta'
       logging.warning('Could not determine the seqfile type, will parse as fasta file')
       seq_reader=read_fasta
    


    n_seq = 0
    char_counter = Counter()
    with open(seqfile, 'r') as fh:
        for name, seq in seq_reader(fh):
            n_seq += 1
            char_counter += Counter(seq)
    
    char_counter.pop('N', None)
    total_len = sum(char_counter.values())
    char_freq = {k: v/total_len for k, v in char_counter.items()}

    logging.debug(f"charfreq of {seqfile}: A:{char_freq['A']}; C:{char_freq['C']}; G:{char_freq['G']}; T:{char_freq['T']}; n_seq: {n_seq}; total_len: {total_len}")

    return n_seq, total_len, char_freq
