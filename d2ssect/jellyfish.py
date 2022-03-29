# jellyfish count -m 19 -t 5 -s 1G -o data/jf/DI-1-1_S6.jf data/fasta/DI-1-1_S6.fasta

import subprocess
from sys import stdout

def count(input,output):
	print(input)
	jfresult = subprocess.run(["jellyfish","count","-m 19", "-t 5", "-s 1G","-o",output,input])
	jfresult


def dump(input, output):
	print(f"dumping results from {input}")
	dump_result = subprocess.Popen(["jellyfish","dump","-ct", input], stdout=subprocess.PIPE)
	dump_sort = subprocess.Popen(["sort","-k1,1"],stdin=dump_result.stdout, stdout=subprocess.PIPE)
	dump_result.stdout.close()
	nkz_result = subprocess.run(["python","v0/Kmers_2_NumbericRepresentation_py3.py","-o", output], stdin=dump_sort.stdout)
	nkz_result