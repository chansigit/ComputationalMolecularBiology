#!/bin/bash
#$ -S /bin/bash
#$ -q week.q
#$ -j yes
#$ -N sjBowtie
#$ -o /data/home/schen/scripts/Nov2016/
/data/software/bowtie/bowtie-0.12.7/bin/bowtie -f /data/home/CMB/homework4/bowtie_index_hg18/hg18 /data/home/CMB/homework4/sequence2.fa
