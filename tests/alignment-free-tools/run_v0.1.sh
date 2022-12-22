#!/usr/bin/bash

# Generates benchmark data

#jellyfish count -m 19 -t 5 -s 1G -o data/jf/DI-1-1_S6.jf data/fasta/DI-1-1_S6.fasta
#jellyfish count -m 19 -t 5 -s 1G -o data/jf/FI-2-21_S28.jf data/fasta/FI-2-21_S28.fasta
#jellyfish count -m 19 -t 5 -s 1G -o data/jf/MI-1-19_S9.jf data/fasta/MI-1-19_S9.fasta
#jellyfish count -m 19 -t 5 -s 1G -o data/jf/TAY_9_S28.jf data/fasta/TAY_9_S28.fasta

jellyfish dump -ct data/jf/DI-1-1_S6.jf | sort -k1,1 | python2 v0/Kmers_2_NumbericRepresentation.py -o data/nkz/DI-1-1_S6.nkz
python2 v0/Composition_of_InputSeqs.py --fasta data/fasta/DI-1-1_S6.fasta --freq data/charfreq/DI-1-1_S6.CharFreq
jellyfish dump -ct data/jf/FI-2-21_S28.jf | sort -k1,1 | python2 v0/Kmers_2_NumbericRepresentation.py -o data/nkz/FI-2-21_S28.nkz
python2 v0/Composition_of_InputSeqs.py --fasta data/fasta/FI-2-21_S28.fasta --freq data/charfreq/FI-2-21_S28.CharFreq
jellyfish dump -ct data/jf/MI-1-19_S9.jf | sort -k1,1 | python2 v0/Kmers_2_NumbericRepresentation.py -o data/nkz/MI-1-19_S9.nkz
python2 v0/Composition_of_InputSeqs.py --fasta data/fasta/MI-1-19_S9.fasta --freq data/charfreq/MI-1-19_S9.CharFreq
jellyfish dump -ct data/jf/TAY_9_S28.jf | sort -k1,1 | python2 v0/Kmers_2_NumbericRepresentation.py -o data/nkz/TAY_9_S28.nkz
python2 v0/Composition_of_InputSeqs.py --fasta data/fasta/TAY_9_S28.fasta --freq data/charfreq/TAY_9_S28.CharFreq

python2 v0/Calculate_D2S.py --kmerset1 data/nkz/DI-1-1_S6.nkz --kmerset1_freq data/charfreq/DI-1-1_S6.CharFreq \
        --kmerset2 data/nkz/FI-2-21_S28.nkz --kmerset2_freq data/charfreq/FI-2-21_S28.CharFreq --D2S_out data/d2s/DI-1-1_S6-FI-2-21_S28.txt
python2 v0/Calculate_D2S.py --kmerset1 data/nkz/DI-1-1_S6.nkz --kmerset1_freq data/charfreq/DI-1-1_S6.CharFreq \
	--kmerset2 data/nkz/MI-1-19_S9.nkz --kmerset2_freq data/charfreq/MI-1-19_S9.CharFreq --D2S_out data/d2s/DI-1-1_S6-MI-1-19_S9.txt
python2 v0/Calculate_D2S.py --kmerset1 data/nkz/DI-1-1_S6.nkz --kmerset1_freq data/charfreq/DI-1-1_S6.CharFreq \
	--kmerset2 data/nkz/TAY_9_S28.nkz --kmerset2_freq data/charfreq/TAY_9_S28.CharFreq --D2S_out data/d2s/DI-1-1_S6-TAY_9_S28.txt

python2 v0/Calculate_D2S.py --kmerset1 data/nkz/FI-2-21_S28.nkz --kmerset1_freq data/charfreq/FI-2-21_S28.CharFreq \
	--kmerset2 data/nkz/MI-1-19_S9.nkz --kmerset2_freq data/charfreq/MI-1-19_S9.CharFreq --D2S_out data/d2s/FI-2-21_S28-MI-1-19_S9.txt
python2 v0/Calculate_D2S.py --kmerset1 data/nkz/FI-2-21_S28.nkz --kmerset1_freq data/charfreq/FI-2-21_S28.CharFreq \
	--kmerset2 data/nkz/TAY_9_S28.nkz --kmerset2_freq data/charfreq/TAY_9_S28.CharFreq --D2S_out data/d2s/FI-2-21_S28-TAY_9_S28.txt

python2 v0/Calculate_D2S.py --kmerset1 data/nkz/MI-1-19_S9.nkz --kmerset1_freq data/charfreq/MI-1-19_S9.CharFreq \
	--kmerset2 data/nkz/TAY_9_S28.nkz --kmerset2_freq data/charfreq/TAY_9_S28.CharFreq --D2S_out data/d2s/MI-1-19_S9-TAY_9_S28.txt

python3 v0/phylip_amalg.py --data data/d2s --matrix data/matrix.txt