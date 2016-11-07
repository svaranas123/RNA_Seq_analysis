### Fastq files
```
#$ -N fastq_download
#$ -cwd
#$ -S /bin/bash
#$ -l mem=7G
#$ -q medium*


module load sra-tools/2.7.0



for i in $(cat srr.txt)
do 
```
