###STAR Index
```
#$ -N STAR_index
#$ -q short*
#$ -cwd
#$ -pe threads 2
#$ -S /bin/bash

mkdir genomeDir

STAR --runMode genomeGenerate     \
    --genomeDir ./genomeDir       \
    --genomeFastaFiles ../../GCF_000001635.25_GRCm38.p5_genomic.fna  \
    --runThreadN 4      \
    --sjdbGTFfile GCF_000001635.25_GRCm38.p5_genomic.gff   \
    --sjdbOverhang 101
```

