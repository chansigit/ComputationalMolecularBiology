#!/bin/bash
#$ -S /bin/bash
#$ -q week.q
#$ -j yes
#$ -N sjBlast1
#$ -o /data/home/schen/scripts/Nov2016/
/data/software/blast/blast-2.2.27/bin/blastn  -db /data/home/CMB/homework4/blastdb/hg18 -query /data/home/CMB/homework4/sequence2.fa -evalue 1e-9 -outfmt 6 
 
