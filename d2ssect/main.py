import multiprocessing, subprocess
from itertools import combinations
import numpy as np
from utils import *
import math, logging, os, argparse, sys
from os.path import exists


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

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Calculate d2s distance from jellyfish kmer count results',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('-l','--jf_files', nargs='+', help='<Required> List of files from jellyfish count', required = True)
    parser.add_argument('-f','--seq_files',nargs='+', help='<Required> List of seq files(fasta/fastq), one for each sample, must aligned with jf files', required=True)
    parser.add_argument('-t','--threads',type=int, default=1,help='Number of threads to use, when sample size are big (>100), using more threads is suggested')
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

    if args.threads == 1:
        seqinfo_list = [ get_seqinfo(seqfile) for seqfile in args.seq_files ]
        d2s_self_list = [cal_d2s_self(jf, seqinfo) for jf, seqinfo in zip(args.jf_files, seqinfo_list)]

        for i,j in index_combinations:
            transformed_d2s = cal_d2s_and_transform(args.jf_files[i], args.jf_files[j], 
                                        seqinfo_list[i], seqinfo_list[j], 
                                        d2s_self_list[i], d2s_self_list[j])
            d2s_combinations_list.append(transformed_d2s)

        d2s_matrix = generate_matrix(d2s_combinations_list, n_sample)
    
    else:
        threads1 = n_sample if threads>n_sample else threads
        pool1 = multiprocessing.Pool(threads1)
        logging.info(f'Using {threads1} of cpus for parallisation')

        seqinfo_list = pool1.map(get_seqinfo, args.seq_files)
        args_s1 = zip(args.jf_files, seqinfo_list)
        d2s_self_list = pool1.starmap(cal_d2s_self, args_s1)
        
        pool1.close()

        threads2 = len(index_combinations) if threads>len(index_combinations) else threads
        pool2 = multiprocessing.Pool(threads2)
        logging.info(f'Using {threads2} of cpus for cross-sample comparison parallisation')

        args_s2 = [(args.jf_files[i],args.jf_files[j],seqinfo_list[i],seqinfo_list[j],d2s_self_list[i],d2s_self_list[j]) for i,j in index_combinations]
        d2s_combinations_list = pool2.starmap(cal_d2s_and_transform, args_s2)
        d2s_matrix = generate_matrix(d2s_combinations_list, n_sample)
        
        pool2.close()

    logging.info("All comparison done! Generating a matrix...")
    for row in d2s_matrix:
        args.output.write('\t'.join(f'{score:.8f}' for score in row) + '\n')
