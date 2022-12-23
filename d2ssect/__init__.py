from d2ssect.utils import *
from d2ssect import jellyfish
import argparse, sys

def cal_d2s_cpp(jf1, jf2, seqinfo1, seqinfo2,d2s_self1, d2s_self2):
    n_seq1, total_len1, char_freq1 = seqinfo1
    n_seq2, total_len2, char_freq2 = seqinfo2

    d2s = jellyfish.d2s(jf1,jf2,n_seq1,n_seq2,total_len1,total_len2, char_freq1,char_freq2)
    return abs(math.log(d2s/math.sqrt(d2s_self1*d2s_self2)))


def cal_d2s_self_cpp(jf, seqinfo):
    n_seq, total_len, char_freq = seqinfo

    return jellyfish.d2s(jf,jf,n_seq,n_seq,total_len,total_len, char_freq,char_freq)

def d2ssectmain():
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

        d2s_self_list = [cal_d2s_self_cpp(jf, seqinfo) for jf, seqinfo in zip(args.jf_files, seqinfo_list)]

        for i,j in index_combinations:
            transformed_d2s = cal_d2s_cpp(args.jf_files[i], args.jf_files[j], 
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
        d2s_self_list = pool1.starmap(cal_d2s_self_cpp, args_s1)
        
        pool1.close()

        threads2 = len(index_combinations) if threads>len(index_combinations) else threads
        pool2 = multiprocessing.Pool(threads2)
        logging.info(f'Using {threads2} of cpus for cross-sample comparison parallisation')

        args_s2 = [(args.jf_files[i],args.jf_files[j],seqinfo_list[i],seqinfo_list[j],d2s_self_list[i],d2s_self_list[j]) for i,j in index_combinations]
#        d2s_combinations_list = pool2.starmap(cal_d2s_and_transform, args_s2)
        d2s_combinations_list = pool2.starmap(cal_d2s_cpp, args_s2)        
        d2s_matrix = generate_matrix(d2s_combinations_list, n_sample)
        
        pool2.close()

    logging.info("All comparison done! Generating a matrix...")
    for row in d2s_matrix:
        args.output.write('\t'.join(f'{score:.8f}' for score in row) + '\n')
