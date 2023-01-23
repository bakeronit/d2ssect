for f in fasta/*.fasta;do
	bn=$(basename $f)
	sample=${bn%.fasta}
	cat $f | bioawk -c fastx '{printf("@%s\n%s\n+\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n",$name,$seq)}' > fastq/$sample.fastq
done
