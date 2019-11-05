#!/bin/bash

#get fastq file name for read 1 and read two and store into variables
for file in ./*.bam
do
	
	echo "==========================================================================================================================================================="
	echo "GlennTools version 0.1 Bash Script" 
	echo "==========================================================================================================================================================="
	echo "Processing: " $file

	
	
	htseq-count -m union --stranded=no -f bam $file /gpfs/group/databases/Homo_sapiens/UCSC/hg19/Annotation/Archives/archive-2013-03-06-11-23-03/Genes/genes.gtf > $file.counts
done