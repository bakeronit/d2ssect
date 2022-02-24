#!/usr/bin/bash

jellyfish count -m 19 -t 5 -s 1G -o jf/DI-1-1_S6.jf fasta/DI-1-1_S6.fasta
jellyfish count -m 19 -t 5 -s 1G -o jf/FI-2-21_S28.jf fasta/FI-2-21_S28.fasta
jellyfish count -m 19 -t 5 -s 1G -o jf/MI-1-19_S9.jf fasta/MI-1-19_S9.fasta
jellyfish count -m 19 -t 5 -s 1G -o jf/TAY_9_S28.jf fasta/TAY_9_S28.fasta

jellyfish dump -ct jf/DI-1-1_S6.jf | sort -k1,1 | python2 jackknifing/jf_scripts/Kmers_2_NumbericRepresentation.py -o nkz/DI-1-1_S6.nkz
python2 jackknifing/jf_scripts/Composition_of_InputSeqs.py --fasta fasta/DI-1-1_S6.fasta --freq charfreq/DI-1-1_S6.CharFreq
jellyfish dump -ct jf/FI-2-21_S28.jf | sort -k1,1 | python2 jackknifing/jf_scripts/Kmers_2_NumbericRepresentation.py -o nkz/FI-2-21_S28.nkz
python2 jackknifing/jf_scripts/Composition_of_InputSeqs.py --fasta fasta/FI-2-21_S28.fasta --freq charfreq/FI-2-21_S28.CharFreq
jellyfish dump -ct jf/MI-1-19_S9.jf | sort -k1,1 | python2 jackknifing/jf_scripts/Kmers_2_NumbericRepresentation.py -o nkz/MI-1-19_S9.nkz
python2 jackknifing/jf_scripts/Composition_of_InputSeqs.py --fasta fasta/MI-1-19_S9.fasta --freq charfreq/MI-1-19_S9.CharFreq
jellyfish dump -ct jf/TAY_9_S28.jf | sort -k1,1 | python2 jackknifing/jf_scripts/Kmers_2_NumbericRepresentation.py -o nkz/TAY_9_S28.nkz
python2 jackknifing/jf_scripts/Composition_of_InputSeqs.py --fasta fasta/TAY_9_S28.fasta --freq charfreq/TAY_9_S28.CharFreq

python2 jackknifing/calc_d2s/Calculate_D2S.py --kmerset1 nkz/DI-1-1_S6.nkz --kmerset1_freq charfreq/DI-1-1_S6.CharFreq \
        --kmerset2 nkz/FI-2-21_S28.nkz --kmerset2_freq charfreq/FI-2-21_S28.CharFreq --D2S_out d2s/DI-1-1_S6-FI-2-21_S28.txt
python2 jackknifing/calc_d2s/Calculate_D2S.py --kmerset1 nkz/DI-1-1_S6.nkz --kmerset1_freq charfreq/DI-1-1_S6.CharFreq \
	--kmerset2 nkz/MI-1-19_S9.nkz --kmerset2_freq charfreq/MI-1-19_S9.CharFreq --D2S_out d2s/DI-1-1_S6-MI-1-19_S9.txt
python2 jackknifing/calc_d2s/Calculate_D2S.py --kmerset1 nkz/DI-1-1_S6.nkz --kmerset1_freq charfreq/DI-1-1_S6.CharFreq \
	--kmerset2 nkz/TAY_9_S28.nkz --kmerset2_freq charfreq/TAY_9_S28.CharFreq --D2S_out d2s/DI-1-1_S6-TAY_9_S28.txt

python2 jackknifing/calc_d2s/Calculate_D2S.py --kmerset1 nkz/FI-2-21_S28.nkz --kmerset1_freq charfreq/FI-2-21_S28.CharFreq \
	--kmerset2 nkz/MI-1-19_S9.nkz --kmerset2_freq charfreq/MI-1-19_S9.CharFreq --D2S_out d2s/FI-2-21_S28-MI-1-19_S9.txt
python2 jackknifing/calc_d2s/Calculate_D2S.py --kmerset1 nkz/FI-2-21_S28.nkz --kmerset1_freq charfreq/FI-2-21_S28.CharFreq \
	--kmerset2 nkz/TAY_9_S28.nkz --kmerset2_freq charfreq/TAY_9_S28.CharFreq --D2S_out d2s/FI-2-21_S28-TAY_9_S28.txt

python2 jackknifing/calc_d2s/Calculate_D2S.py --kmerset1 nkz/MI-1-19_S9.nkz --kmerset1_freq charfreq/MI-1-19_S9.CharFreq \
	--kmerset2 nkz/TAY_9_S28.nkz --kmerset2_freq charfreq/TAY_9_S28.CharFreq --D2S_out d2s/MI-1-19_S9-TAY_9_S28.txt

python jackknifing/PHYLIP/phylip_amalg.py --data d2s --matrix matrix.txt