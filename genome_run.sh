#!/bin/csh

#$ -M alaguna@nd.edu	 # Email address for job notification
#$ -m ae		 # Send mail when job begins, ends and aborts
#$ -pe mpi-24 24	 # Specify parallel environment and legal core size
#$ -q long		 # Specify queue
#$ -N genome_exact	 # Specify job name

module load python3      # Required modules

python3 genome.py GCF_000005845.2_ASM584v2_genomic.fna bp100/bp100 
