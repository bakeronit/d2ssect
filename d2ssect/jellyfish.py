# jellyfish count -m 19 -t 5 -s 1G -o data/jf/DI-1-1_S6.jf data/fasta/DI-1-1_S6.fasta

import subprocess

def count(input,output):
	print(input)
	jfresult = subprocess.run(["jellyfish","count","-m 19", "-t 5", "-s 1G","-o",output,input])
	jfresult