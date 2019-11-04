#!/bin/bash

#get fastq file name for read 1 and read two and store into variables
for file in ./*rmdup.sam
do
	echo "Filtering: " $file
	sampleID=${file%%.rmdup.sam}
	samtools view -bSh $file > $sampleID.q40.rmdup.bam
done
