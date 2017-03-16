#!/bin/bash -login
#PBS -l walltime=2:00:00,nodes=1:ppn=1,mem=12Gb
#mdiag -A ged
#PBS -m abe			#send email to myself
#PBS -N tblastn		#give name to the job


module load BLAST+/2.2.30
#export BLASTDB=/mnt/research/common-data/Bio/blastdb:$BLASTDB

cd $PBS_O_WORKDIR

tblastn -query $input \
        -db $DB \
        -out $input.tbn \
        -outfmt "6 std qlen slen" \
        -num_threads 1 \
        -evalue 1e-5 \
        -max_target_seqs 50 \
        -seg no

qstat -f ${PBS_JOBID}

