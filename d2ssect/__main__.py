from d2ssect.utils import *
import d2ssect.jellyfish as jellyfish
import argparse, sys
from operator import mul
from functools import reduce
from typing import Tuple
import math, logging
from os.path import basename,splitext

import multiprocessing, subprocess
from itertools import combinations



def cal_d2s_cpp(jf1, jf2, seqinfo1, seqinfo2):
    n_seq1, total_len1, char_freq1 = seqinfo1
    n_seq2, total_len2, char_freq2 = seqinfo2

    return jellyfish.d2s(jf1,jf2,n_seq1,n_seq2,total_len1,total_len2, char_freq1,char_freq2)

def main():
    parser = argparse.ArgumentParser(
        description='Calculate d2s distance from jellyfish kmer count results',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('-l','--jf_files', nargs='+', help='<Required> List of files from jellyfish count', required = True)
    parser.add_argument('-f','--seq_files',nargs='+', help='<Required> List of seq files(fasta/fastq), one for each sample, must aligned with jf files', required=True)
    parser.add_argument('-t','--threads',type=int, default=1,help='Number of threads to use, when sample size are big (>100), using more threads is suggested')
    parser.add_argument('-p','--precision',type=int, default=6,help='Number of significant figures for numerical output')
    parser.add_argument('-o','--output',nargs='?',type=argparse.FileType('w'),default=sys.stdout)
    parser.add_argument('--debug',default=False, action='store_true',help="output intermedia results")
    args = parser.parse_args()

    if args.debug:
        logging.basicConfig(level=logging.DEBUG, format='#%(levelname)s :: %(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')
        logging.info(f'Logging intermediate files for debugging')
    else:
        logging.basicConfig(level=logging.INFO, format='#%(levelname)s :: %(asctime)s - %(message)s', datefmt='%d-%b-%y %H:%M:%S')


    if args.threads >= multiprocessing.cpu_count():
        threads = multiprocessing.cpu_count()
        logging.warning(f'Using all or more cpus than the limit: {threads}, try use less..')
    else:
        threads = args.threads
    
    logging.info(f'jf files list: {" ".join(args.jf_files)}')
    logging.info(f'seq files list: {" ".join(args.seq_files)}')

    if len(args.jf_files) != len(args.seq_files):
        logging.error(f'Yikes, different number of jf files and seq files. Exit..')
        sys.exit()
    
    n_sample = len(args.jf_files)
    index_combinations = list(combinations(range(n_sample),2))
    logging.info(f'Number of pairwise d2s calculation: {len(index_combinations)}')
    d2s_combinations_list = []

    sample_names = [splitext(basename(seqfile))[0] for seqfile in args.seq_files]

    if args.threads == 1:
        seqinfo_list = [ get_seqinfo(seqfile) for seqfile in args.seq_files ]

        for i,j in index_combinations:
            transformed_d2s = cal_d2s_cpp(args.jf_files[i], args.jf_files[j], 
                                        seqinfo_list[i], seqinfo_list[j])
            d2s_combinations_list.append(transformed_d2s)
    
    else:
        threads1 = n_sample if threads>n_sample else threads
        pool1 = multiprocessing.Pool(threads1)
        logging.info(f'Using {threads1} of cpus for parallisation')

        seqinfo_list = pool1.map(get_seqinfo, args.seq_files)
        
        pool1.close()

        threads2 = len(index_combinations) if threads>len(index_combinations) else threads
        pool2 = multiprocessing.Pool(threads2)
        logging.info(f'Using {threads2} of cpus for cross-sample comparison parallisation')

        args_s2 = [(args.jf_files[i],args.jf_files[j],seqinfo_list[i],seqinfo_list[j]) for i,j in index_combinations]
        d2s_combinations_list = pool2.starmap(cal_d2s_cpp, args_s2)        
        
        pool2.close()

    logging.info("All comparison done! Generating a matrix...")

    d2s_matrix = [([0]*n_sample) for i in range(n_sample)]
    for p in range(len(index_combinations)):
        i,j = index_combinations[p]
        v = d2s_combinations_list[p]
        d2s_matrix[i][j] = v
        d2s_matrix[j][i] = v

    for r in range(n_sample):
        row = d2s_matrix[r]
        sample = sample_names[r]
        args.output.write(sample + '\t' + '\t'.join(f'{score:.{args.precision}f}' for score in row) + '\n')

    # for row in d2s_matrix:
    #     args.output.write('\t'.join(f'{score:.8f}' for score in row) + '\n')

if __name__ == "__main__":
    main()

