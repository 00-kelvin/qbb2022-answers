echo "Question 2.1:" 
grep -c '>' contigs.fasta

samtools faidx contigs.fasta

echo "Question 2.2:" 
awk '{sum+=$2} END {print sum}' "contigs.fasta.fai"

echo "Question 2.3:"

cut -f2 contigs.fasta.fai | sort -n -r | head -1
