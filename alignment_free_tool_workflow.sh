#!/usr/bin/env bash

jellyfish_count() {
	OUTDIR=$1
	SAMPLE=$2
	jellyfish count -m 19 -t 5 -s 1G -o ${OUTDIR}/${SAMPLE}.jf data/fasta/${SAMPLE}.fasta
}

jellyfish_dump() {
	OUTDIR=$1
	SAMPLE=$2
	jellyfish dump -ct ${OUTDIR}/${SAMPLE}.jf | sort -k1,1 | python2 v0/Kmers_2_NumbericRepresentation.py -o ${OUTDIR}/${SAMPLE}.nkz
	python2 v0/Composition_of_InputSeqs.py --fasta data/fasta/${SAMPLE}.fasta --freq ${OUTDIR}/${SAMPLE}.CharFreq
}

cal_d2s() {
	OUTDIR=$1
	SAMPLE1=$2
	SAMPLE2=$3
	mkdir -p ${OUTDIR}/d2s
	python2 v0/Calculate_D2S.py --kmerset1 ${OUTDIR}/${SAMPLE1}.nkz --kmerset1_freq ${OUTDIR}/${SAMPLE1}.CharFreq \
    --kmerset2 ${OUTDIR}/${SAMPLE2}.nkz --kmerset2_freq ${OUTDIR}/${SAMPLE2}.CharFreq --D2S_out ${OUTDIR}/d2s/${SAMPLE1}-${SAMPLE2}.txt
}

generate_matrix() {
	OUTDIR=$1
	python3 v0/phylip_amalg.py --data ${OUTDIR}/d2s --matrix ${OUTDIR}/matrix.txt
}