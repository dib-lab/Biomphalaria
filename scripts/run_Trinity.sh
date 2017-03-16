#!/bin/bash
#PBS -l walltime=18:00:00,nodes=1:ppn=6,mem=128Gb
#mdiag -A ged
#PBS -m abe			#send email to myself
#PBS -N trinity		#give name to the job


module load trinity/2.2.0
#module load GNU/4.9
#module load trinity/2.3.2
#module load bowtie2/2.2.6
#module load Java/1.8.0_31 ## load automatically
#module load SAMTools/0.1.19 ## load automatically


cd $PBS_O_WORKDIR

Trinity --seqType fq --max_memory 100G --CPU 6 \
--left $(echo ${lf[*]} | tr ' ' ',') \
--right $(echo ${rt[*]} | tr ' ' ',')


qstat -f ${PBS_JOBID}

