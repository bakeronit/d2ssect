#!/usr/bin/env python3

import argparse
from posixpath import basename
import tempfile
import shutil
import os
import d2ssect.jellyfish as jellyfish

parser = argparse.ArgumentParser(description='Calculate d2s statistics')

parser.add_argument('inputfiles',type=str,nargs='+',help='input files in fastq or fastq format')


args  = parser.parse_args()


# Setup temporary directory where we will store intermediate outputs
#
temp_path = tempfile.mkdtemp(prefix="d2ssect_",dir=".")
print("Storing intermediate files in ",temp_path)


# Do jellyfish counts and jellyfish dump for each file
#
for f in args.inputfiles:
	of=tempfile.mkstemp(dir=temp_path,prefix=os.path.basename(f))
	print(of)
	result = jellyfish.count(f,of[1])
	print(result)




# If cleanup
#shutil.rmtree(temp_path)