#get fastq file name for read 1 and read two and store into variables
for bamfile in ./*srt.bam
do
	sampleID=${bamfile%%.srt.bam}
	echo "=========================================================================================================================================================================================================="
	let pair++
	echo "remove Dups for $file:"
	echo "Sample ID: $sampleID:"
	echo $r2file
	echo "GlennTools version 0.1: GlennDUP!!!!"
	echo "=========================================================================================================================================================================================================="
			
	zcat ../${sampleID}_R1.fastq.gz | paste - - - - > ../pastedfq/${sampleID}_R1.fastq.pst
	zcat ../${sampleID}_R2.fastq.gz | paste - - - - > ../pastedfq/${sampleID}_R2.fastq.pst

	samtools view -h $bamfile | python UMIremovedup.py ../pastedfq/${sampleID}_R1.fastq.pst ../pastedfq/${sampleID}_R2.fastq.pst | samtools view -bhS - | samtools sort - > ${sampleID}.srt.rmdup.bam

done