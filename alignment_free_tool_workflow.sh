#!/usr/bin/env bash

jellyfish_count() {
	OUTDIR=$1
	jellyfish count -m 19 -t 5 -s 1G -o ${OUTDIR}/DI-1-1_S6.jf data/fasta/DI-1-1_S6.fasta
	jellyfish count -m 19 -t 5 -s 1G -o ${OUTDIR}/FI-2-21_S28.jf data/fasta/FI-2-21_S28.fasta
	jellyfish count -m 19 -t 5 -s 1G -o ${OUTDIR}/MI-1-19_S9.jf data/fasta/MI-1-19_S9.fasta
	jellyfish count -m 19 -t 5 -s 1G -o ${OUTDIR}/TAY_9_S28.jf data/fasta/TAY_9_S28.fasta
}

jellyfish_dump() {
	OUTDIR=$1
	jellyfish dump -ct ${OUTDIR}/DI-1-1_S6.jf | sort -k1,1 | python2 v0/Kmers_2_NumbericRepresentation.py -o ${OUTDIR}/DI-1-1_S6.nkz
	python2 v0/Composition_of_InputSeqs.py --fasta data/fasta/DI-1-1_S6.fasta --freq ${OUTDIR}/DI-1-1_S6.CharFreq
	jellyfish dump -ct ${OUTDIR}/FI-2-21_S28.jf | sort -k1,1 | python2 v0/Kmers_2_NumbericRepresentation.py -o ${OUTDIR}/FI-2-21_S28.nkz
	python2 v0/Composition_of_InputSeqs.py --fasta data/fasta/FI-2-21_S28.fasta --freq ${OUTDIR}/FI-2-21_S28.CharFreq
	jellyfish dump -ct ${OUTDIR}/MI-1-19_S9.jf | sort -k1,1 | python2 v0/Kmers_2_NumbericRepresentation.py -o ${OUTDIR}/MI-1-19_S9.nkz
	python2 v0/Composition_of_InputSeqs.py --fasta data/fasta/MI-1-19_S9.fasta --freq ${OUTDIR}/MI-1-19_S9.CharFreq
	jellyfish dump -ct ${OUTDIR}/TAY_9_S28.jf | sort -k1,1 | python2 v0/Kmers_2_NumbericRepresentation.py -o ${OUTDIR}/TAY_9_S28.nkz
	python2 v0/Composition_of_InputSeqs.py --fasta data/fasta/TAY_9_S28.fasta --freq ${OUTDIR}/TAY_9_S28.CharFreq
}

cal_d2s() {
	OUTDIR=$1
	mkdir -p ${OUTDIR}/d2s
	python2 v0/Calculate_D2S.py --kmerset1 ${OUTDIR}/DI-1-1_S6.nkz --kmerset1_freq ${OUTDIR}/DI-1-1_S6.CharFreq \
    --kmerset2 ${OUTDIR}/FI-2-21_S28.nkz --kmerset2_freq ${OUTDIR}/FI-2-21_S28.CharFreq --D2S_out ${OUTDIR}/d2s/DI-1-1_S6-FI-2-21_S28.txt
	
	python2 v0/Calculate_D2S.py --kmerset1 ${OUTDIR}/DI-1-1_S6.nkz --kmerset1_freq ${OUTDIR}/DI-1-1_S6.CharFreq \
	--kmerset2 ${OUTDIR}/MI-1-19_S9.nkz --kmerset2_freq ${OUTDIR}/MI-1-19_S9.CharFreq --D2S_out ${OUTDIR}/d2s/DI-1-1_S6-MI-1-19_S9.txt
	
	python2 v0/Calculate_D2S.py --kmerset1 ${OUTDIR}/DI-1-1_S6.nkz --kmerset1_freq ${OUTDIR}/DI-1-1_S6.CharFreq \
	--kmerset2 ${OUTDIR}/TAY_9_S28.nkz --kmerset2_freq ${OUTDIR}/TAY_9_S28.CharFreq --D2S_out ${OUTDIR}/d2s/DI-1-1_S6-TAY_9_S28.txt

	python2 v0/Calculate_D2S.py --kmerset1 ${OUTDIR}/FI-2-21_S28.nkz --kmerset1_freq ${OUTDIR}/FI-2-21_S28.CharFreq \
	--kmerset2 ${OUTDIR}/MI-1-19_S9.nkz --kmerset2_freq ${OUTDIR}/MI-1-19_S9.CharFreq --D2S_out ${OUTDIR}/d2s/FI-2-21_S28-MI-1-19_S9.txt

	python2 v0/Calculate_D2S.py --kmerset1 ${OUTDIR}/FI-2-21_S28.nkz --kmerset1_freq ${OUTDIR}/FI-2-21_S28.CharFreq \
	--kmerset2 ${OUTDIR}/TAY_9_S28.nkz --kmerset2_freq ${OUTDIR}/TAY_9_S28.CharFreq --D2S_out ${OUTDIR}/d2s/FI-2-21_S28-TAY_9_S28.txt

	python2 v0/Calculate_D2S.py --kmerset1 ${OUTDIR}/MI-1-19_S9.nkz --kmerset1_freq ${OUTDIR}/MI-1-19_S9.CharFreq \
	--kmerset2 ${OUTDIR}/TAY_9_S28.nkz --kmerset2_freq ${OUTDIR}/TAY_9_S28.CharFreq --D2S_out ${OUTDIR}/d2s/MI-1-19_S9-TAY_9_S28.txt
}

generate_matrix() {
	OUTDIR=$1
	python3 v0/phylip_amalg.py --data ${OUTDIR}/d2s --matrix ${OUTDIR}/matrix.txt
}