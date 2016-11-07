### Trimmomatic 
```
#$ -N trimmomatic
#$ -cwd
#$ -S /bin/bash
#$ -l mem=7G
#$ -q medium*


module load sra-tools trimmomatic bwa fastqc

for i in $(cat fastq.txt)

do
	trimmomatic SE $1 $1.trimmed.fastq ILLUMINACLIP:/data/apps/trimmomatic/0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

done
```
