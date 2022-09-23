# Step 1: index the reference genome
bwa index sacCer3.fa 

for SAMPLE in *.fastq
do
	# Step 2: align the reads with the reference genome
	bwa mem -t 4 -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}" -o ${SAMPLE}.sam sacCer3.fa ${SAMPLE}.fastq
	
	# Step 3a: create sorted bam files
	samtools sort -O bam -o ${SAMPLE}.bam ${SAMPLE}.sam
	
	# Step 3b: create an index for each sorted bam file
	samtools index ${SAMPLE}.bam