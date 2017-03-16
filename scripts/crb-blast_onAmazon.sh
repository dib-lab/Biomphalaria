## Functional annotation of whole transcriptome assembly
## start m3.medium instance with 50GB extra-storage (EBS)
## log on to your machine and prepare the storage drive
lsblk #Check for the directory structure
sudo file -s /dev/xvdb #output should be "data", which means no file system on the device
sudo mkfs -t ext4 /dev/xvdb
sudo mkdir workDir
sudo mount /dev/xvdb workDir
lsblk
## install required software
sudo su
apt-get update
\curl -sSL https://get.rvm.io | bash -s stable --ruby
source /usr/local/rvm/scripts/rvm
gem install crb-blast
#wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.29/ncbi-blast-2.2.29+-x64-linux.tar.gz
#tar zxvpf ncbi-blast-2.2.29+-x64-linux.tar.gz 
#export PATH=$PATH:$HOME/ncbi-blast-2.2.29+/bin

## download required sequences 
cd workDir
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz
#wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz
scp mansourt@hpc.msu.edu:/mnt/home/mansourt/temp/Trinity.clean.201.exp2.UniVec.fasta biomph.fasta

## run crb-blast
crb-blast --query biomph.fasta --target uniprot_sprot.fasta --threads 4 --output biomph_crbBlast.tsv

