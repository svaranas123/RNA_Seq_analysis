### HTseq counts

```
#$ -N HTseq_count   
#$ -q medium*
#$ -cwd
#$ -pe threads 4
#$ -S /bin/bash


mkdir Counts

for sorted_bam_path in $(find -name *.bam)
do
    counts_file=$(echo $sorted_bam_path | grep -o "SRR3663[0-9]*")_$1_ct
    echo "The target bam file is: "$sorted_bam_path
    echo "==================================================="
        htseq-count -f bam \
                    -t gene \
                    -i Name  \
                    -r pos \
                    $sorted_bam_path \
                    ./GCF_000001635.25_GRCm38.p5_genomic.gff  \
                    | \
                    grep -v '^__' > ./Counts/$counts_file
    echo "The count data has been written into: $counts_file"
    echo "==================================================="
done


```
