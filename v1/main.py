from itertools import combinations
from collections import defaultdict
import numpy as np
from utils import *
import sys, glob
import time

#./count_in_file *.jf | python main.py 
# these codes needs to be changed
jf_list = glob.glob("*.jf")
fasta_list = [f"../fasta/{i.rstrip('.jf')}.fasta" for i in jf_list]
n_sample = len(jf_list)
index_combinations = list(combinations(range(n_sample),2))

seqinfo_list = [get_seqinfo(i) for i in fasta_list]

def parse_kmer_table(kmer_table):
    d2s_matrix = np.zeros((n_sample,n_sample),dtype='float')

    for line in kmer_table:
        fields = line.strip().split()
        kmer = line.split()[0]

        for i, j in index_combinations:
            if fields[i+1] == "0" or fields[j+1] == "0":
                continue
            d2s_matrix[i,j] += d2s_per_kmer(kmer, int(fields[i+1]), int(fields[j+1]), seqinfo_list[i], seqinfo_list[j])
        
    return d2s_matrix
        

#start_time = time.time()
#print(parse_kmer_table(sys.stdin.readlines()))
#print("follow formula in paper--- %s seconds ---" % (time.time() - start_time))


def parse_kmer_table_v0(kmer_table):
    d2s_matrix = np.zeros((n_sample,n_sample),dtype='float')
    diagonal = [0] * n_sample

    for line in kmer_table:
        fields = line.strip().split()
        kmer = line.split()[0]
        not_calculate = [True] * n_sample # initialise as True, to avoid repeat calculation of kmerset1_vs_kmerset1..

        for i, j in index_combinations:
            if fields[i+1] == "0" and fields[j+1] == "0":
                continue
            elif fields[i+1] == "0" and not_calculate[j]:
                diagonal[j] += d2s_per_kmer(kmer, int(fields[j+1]), int(fields[j+1]), seqinfo_list[j], seqinfo_list[j])
                not_calculate[j] = False
            elif fields[j+1] == "0" and not_calculate[i]:
                diagonal[i] += d2s_per_kmer(kmer, int(fields[i+1]), i nt(fields[i+1]), seqinfo_list[i], seqinfo_list[i])
                not_calculate[i] = False
            else:
                if not_calculate[i]:
                    diagonal[i] += d2s_per_kmer(kmer, int(fields[i+1]), int(fields[i+1]), seqinfo_list[i], seqinfo_list[i])
                    not_calculate[i] = False
                if not_calculate[j]:
                    diagonal[j] += d2s_per_kmer(kmer, int(fields[j+1]), int(fields[j+1]), seqinfo_list[j], seqinfo_list[j])
                    not_calculate[j] = False
                d2s_matrix[i,j] += d2s_per_kmer(kmer, int(fields[i+1]), int(fields[j+1]), seqinfo_list[i], seqinfo_list[j])

    print(d2s_matrix)
    print(diagonal)
    for i, j in index_combinations:
        d2s_matrix[i,j] = abs(math.log(d2s_matrix[i,j]/math.sqrt((diagonal[i]*diagonal[j]))))
    

    ## to fill the lower triangle in matrix
    low = np.tril_indices(n_sample, -1)
    d2s_matrix[low] = d2s_matrix.T[low]
    np.savetxt("bbb",d2s_matrix,fmt='%.8f')

    return d2s_matrix

start_time = time.time()
print(parse_kmer_table_v0(sys.stdin.readlines()))
print("v0 method --- %s seconds ---" % (time.time() - start_time))