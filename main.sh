#!/bin/sh

## construction of the basic diretory structure
git clone https://github.com/drtamermansour/Biomphalaria.git
cd Biomphalaria
Biomphalaria=$(pwd)

## create working directory and define paths for raw data and scripts
mkdir -p $Biomphalaria/{prepdata,QC_raw,adap_remove}
script_path=$Biomphalaria/scripts
prepData=$Biomphalaria/prepdata
QC_raw=$Biomphalaria/QC_raw
trimmed_data=$Biomphalaria/adap_remove

## Download the data from GSAF website
cd $prepData
bash gsaf_download.sh
chmod u-w $prepData/*.fastq.gz

## QC check for the original data files 
cd $QC_raw
for f in ${prepData}/*.gz; do qsub -v INPUT_FILE="$f" ${script_path}/fastqc.sh; done
mv ${prepData}/{*.zip,*.html} $QC_raw/.
#============================================================================================
## adaptor removal only 
cd ${trimmed_data}
bash ${script_path}/run_adapter_trimmer.sh ${prepData} ${trimmed_data} ${script_path}/adapter_removal.sh

## combine single reads 
for f in ${trimmed_data}/*; do if [ -d "${f}" ]; then echo "${f}"; cat "${f}"/s1_se.fq "${f}"/s2_se.fq > "${f}"/$(basename "$f").s_se.fq; fi; done;

## QC check
for dir in ${trimmed_data}/*; do cd $dir; for f in *_pe.fq *.s_se.fq; do qsub -v INPUT_FILE="$f" ${script_path}/fastqc.sh; done; done;
#===========================================================================================
## prepare the files for subsequent analysis
for dir in ${trimmed_data}/*; do if [ -d "${dir}" ]; then echo $dir; cd $dir; qsub -v output=$(basename "$dir").s_pe.fq ${script_path}/interleave.sh; fi; done;

## filter abundance
mkdir $Biomphalaria/abundFilter
cd $Biomphalaria/abundFilter
input_files=()
for f in ${trimmed_data}/*/*.s_[ps]e.fq; do input_files+=($f); done
qsub -v graph=$"countgraph_k20.kt",input_files="${input_files[*]}" ${script_path}/load-into-counting.sh
qsub -v input=$"countgraph_k20.kt",files="${input_files[*]}" ${script_path}/filter_abund.sh
for f in filter_abund.e*;do grep -B2 '^output in' "$f";done > filter_abund.summary
## remove short reads
module load FASTX/0.0.14
for f in *.abundfilt;do cat $f | fastx_clipper -l 25 > $f.long;done
## break out the orphaned and still-paired reads and rename files (this step ends with .s_pe.fq and .s_se.fq for each sample)
#for i in *.s_pe.*.abundfilt.long; do extract-paired-reads.py $i; done
for i in *.s_pe.*.abundfilt.long; do qsub -v input=$i $script_path/extract-paired-reads.sh; done
##  combine the orphaned reads into a single file & rename pe files
for i in *.s_se.fq.abundfilt.long; do
   pe_orphans=$(basename $i .s_se.fq.abundfilt.long).s_pe.fq.abundfilt.long.se
   cat $i $pe_orphans > $(basename $i .abundfilt.long)
done
#rm *.abundfilt.se
for i in *.abundfilt.long.pe; do mv $i $(basename $i .abundfilt.long.pe); done

## read count after filtration
echo PE $(($(cat *.s_pe.fq | wc -l)/8)) > read_count_after_filtration
echo SE $(($(cat *.s_se.fq | wc -l)/4)) >> read_count_after_filtration

## QC check
for f in *.s_pe.fq *.s_se.fq; do qsub -v INPUT_FILE="$f" ${script_path}/fastqc.sh; done;

## split interleaved files & ## merge the single reads to the end of left reads (this step reform the data into .s_pe1_se.fq & .s_pe2.fq for each sample )
for i in *.s_pe.fq; do qsub -v input=$i $script_path/split-paired-reads.sh; done
for f in *.s_pe.fq; do
  base=$(basename $f);
  file_end=$(tail -n 8 $f | head -n 4);
  file_endA=$(tail -n 4 $base.1);
  if [ "$file_end" != "$file_endA" ];then echo $f.1;fi
  file_end=$(tail -n 4 $f);
  file_endB=$(tail -n 4 $base.2);
  if [ "$file_end" != "$file_endB" ];then echo $f.2;fi
done

for f in *.s_pe.fq.1; do echo $f; base=$(basename $f .s_pe.fq.1);
 cat $f $base.s_se.fq | awk '{if (NR % 4 == 1) {print $1"/1"} else {print $0} }' > $base.s_pe1_se.fq; done

for f in *.s_pe.fq.2; do echo $f; base=$(basename $f .s_pe.fq.2); 
 cat $f | awk '{if (NR % 4 == 1) {print $1"/2"} else {print $0} }' > $base.s_pe2.fq; done

## check fastq files
module load SAMStat/20130521
for f in *.s_pe1_se.fq *.s_pe2.fq;do samstat $f;done

## Trinity
lf_files=()
for f in *.s_pe1_se.fq; do lf_files+=($f); done;
rt_files=()
for f in *.s_pe2.fq; do rt_files+=($f); done;
qsub -v lf="${lf_files[*]}",rt="${rt_files[*]}" ${script_path}/run_Trinity.sh
echo $(grep "^>" trinity_out_dir/Trinity.fasta | wc -l) ## 149545

## SeqClean (trim polyA tails and remove dust "low complexity seq")
module load SeqClean/20130718
seqclean_dir=$Biomphalaria/abundFilter/trinity_out_dir/seqclean
mkdir $seqclean_dir
cd $seqclean_dir
seqclean $Biomphalaria/abundFilter/trinity_out_dir/Trinity.fasta
## no of transcripts
echo $(($(grep "^>" $Biomphalaria/abundFilter/trinity_out_dir/Trinity.fasta | wc -l) - $(grep "^>" Trinity.fasta.clean | wc -l))) ## 38

## remove small transcripts
perl ${script_path}/removesmalls.pl 201 Trinity.fasta.clean > Trinity.fasta.clean.201
echo $(($(grep "^>" Trinity.fasta.clean | wc -l) - $(grep "^>" Trinity.fasta.clean.201 | wc -l))) ## 745
trinity_transcriptome=$seqclean_dir/Trinity.fasta.clean.201
##########
## get the longest isoform from my asselbly (edit the Ctr+v,tab at coomand line)
cat Trinity.fasta.clean.201 | awk '/^>/ {if(N>0) printf("\n"); printf("%s\t",$0);N++;next;} {printf("%s",$0);} END {if(N>0) printf("\n");}' | sed 's/_i/_i\t/;s/ len=/\t/;s/ path=/\tpath=/' | sort -t $'\t' -k1,1 -k3,3nr | sort -t $'\t' -k1,1 -u -s | cut -f 1,2,5 | sed 's/Ctr+v,tab/./' | tr "\t" "\n"  | fold -w 60 > Trinity.fasta.clean.201.longest
###########
## Abundance estimation
cd $seqclean_dir
qsub -v index="salmon_index",transcriptome="$trinity_transcriptome" ${script_path}/salmonIndex.sh

cd $Biomphalaria/abundFilter
for f in $Biomphalaria/abundFilter/*.s_pe.fq.1; do if [ -f $f ]; then
 identifier=$(basename ${f%.s_pe.fq.1}); echo $identifier;
fi;done | uniq > identifiers.txt
identifiers=$Biomphalaria/abundFilter/identifiers.txt

while read identifier;do
 ls ${identifier}.s_pe.fq.1 ${identifier}.s_pe2.fq ${identifier}.s_se.fq
 qsub -v index="$seqclean_dir/salmon_index",identifier=$identifier ${script_path}/salmonQuant_PE.sh
 qsub -v index="$seqclean_dir/salmon_index",identifier=$identifier ${script_path}/salmonQuant_SE.sh
done < $identifiers
find ./*.quant -name *.sf -exec grep -H "mapping rate" {} \; | sort > salmonQuant_summary.txt
python $script_path/gather-counts2.py -i "$(pwd)"
echo "transcript"$'\t'"length" > transcripts.lengthes
sf=$(find ./*.quant -name \*.sf | head -n1)
cat $sf | grep -v "^#" | awk -F "\t" -v OFS='\t' '{print $1,$2}' >> transcripts.lengthes

# generation of gene map
perl ${script_path}/get_Trinity_gene_to_trans_map.pl $trinity_transcriptome > gene_transcript_map

#module load R/3.0.1
> targets_list
while read identifier;do
  echo $(pwd) $identifier
  bash $script_path/run_calcTPM.sh "$(pwd)" "$identifier" "transcripts.lengthes" "gene_transcript_map" ${script_path}/calcTPM2.R
  #Rscript ${script_path}/calcTPM_tis.R "$(pwd)" "$identifier" "transcripts.lengthes" "gene_transcript_map" >> targets_list
done < $identifiers
bash $script_path/abund_est.sh

## exclude the unexpressed transcripts from the transcriptome
cd $seqclean_dir
module load QIIME/1.8.0
filter_fasta.py --input_fasta_fp $trinity_transcriptome --output_fasta_fp Trinity.clean.201.exp.fasta --seq_id_fp $Biomphalaria/abundFilter/unexp_isoformTPM --negate
exp_transcriptome=$seqclean_dir/Trinity.clean.201.exp.fasta  
echo $(($(grep "^>" Trinity.fasta.clean.201 | wc -l) - $(grep "^>" Trinity.clean.201.exp.fasta | wc -l))) ## 4445
###################
## repeat the Salmon filteration
cd $Biomphalaria/abundFilter
mkdir SalmonRun1
mv *.quant salmonQuant_summary.txt *.quant.counts transcripts.lengthes gene_transcript_map *.dataSummary_comp targets_list *TPM SalmonRun1/.

cd $seqclean_dir
qsub -v index="salmon_index2",transcriptome="$exp_transcriptome" ${script_path}/salmonIndex.sh

cd $Biomphalaria/abundFilter
while read identifier;do
 ls ${identifier}.s_pe.fq.1 ${identifier}.s_pe2.fq ${identifier}.s_se.fq
 qsub -v index="$seqclean_dir/salmon_index2",identifier=$identifier ${script_path}/salmonQuant_PE.sh
 qsub -v index="$seqclean_dir/salmon_index2",identifier=$identifier ${script_path}/salmonQuant_SE.sh
done < $identifiers
find ./*.quant -name *.sf -exec grep -H "mapping rate" {} \; | sort > salmonQuant_summary.txt
python $script_path/gather-counts2.py -i "$(pwd)"
echo "transcript"$'\t'"length" > transcripts.lengthes
sf=$(find ./*.quant -name \*.sf | head -n1)
cat $sf | grep -v "^#" | awk -F "\t" -v OFS='\t' '{print $1,$2}' >> transcripts.lengthes

# generation of gene map
perl ${script_path}/get_Trinity_gene_to_trans_map.pl $exp_transcriptome > gene_transcript_map

#module load R/3.0.1
> targets_list
while read identifier;do
  echo $(pwd) $identifier
  bash $script_path/run_calcTPM.sh "$(pwd)" "$identifier" "transcripts.lengthes" "gene_transcript_map" ${script_path}/calcTPM2.R
  #Rscript ${script_path}/calcTPM_tis.R "$(pwd)" "$identifier" "transcripts.lengthes" "gene_transcript_map" >> targets_list
done < $identifiers
bash $script_path/abund_est.sh

## exclude the unexpressed transcripts from the transcriptome
cd $seqclean_dir
module load QIIME/1.8.0
filter_fasta.py --input_fasta_fp $exp_transcriptome --output_fasta_fp Trinity.clean.201.exp2.fasta --seq_id_fp $Biomphalaria/abundFilter/unexp_isoformTPM --negate
exp_transcriptome=$seqclean_dir/Trinity.clean.201.exp2.fasta  
echo $(($(grep "^>" Trinity.fasta.clean.201 | wc -l) - $(grep "^>" Trinity.clean.201.exp2.fasta | wc -l))) ## 4526

################
## remove univec contaminants
## The UniVec Database: http://www.ncbi.nlm.nih.gov/tools/vecscreen/univec/
## Contamination in Sequence Databases: http://www.ncbi.nlm.nih.gov/tools/vecscreen/contam/
## About VecScreen (parameters and categories): http://www.ncbi.nlm.nih.gov/tools/vecscreen/about/
## Interpretation of VecScreen Results: http://www.ncbi.nlm.nih.gov/tools/vecscreen/interpretation/
mkdir -p $Biomphalaria/UniVec/db
cd $Biomphalaria/UniVec/db
wget ftp://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec
mv UniVec UniVec.fasta
#wget ftp://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec_Core
#mv UniVec_Core UniVec_Core.fasta
module load BLAST+/2.2.30
makeblastdb -in UniVec.fasta -dbtype nucl
#makeblastdb -in UniVec_Core.fasta -dbtype nucl
cd $Biomphalaria/UniVec
input=$exp_transcriptome
DB=$Biomphalaria/UniVec/db/UniVec.fasta
blastn -task blastn -reward 1 -penalty -5 -gapopen 3 -gapextend 3 -dust yes -soft_masking true -evalue 700 -searchsp 1750000000000 -query $input -db $DB -num_threads 1 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send evalue bitscore score" -out univec_blastn.out
cat univec_blastn.out | awk '{print $1}' | sort | uniq | wc -l ## 30791

## choose all hit with score more than the weak hit for terminal match
cat univec_blastn.out | awk 'BEGIN{OFS="\t";} {if($14>=16) print $1,$7,$8,$9,$14;}' | sort -k1,1 -k2,2n -k3,3nr -k5,5nr > univec_blastn.out.sig
cat univec_blastn.out.sig | awk '{print $1}' | sort | uniq | wc -l ## 30791

## remove contained matches and merge overlapping or neighbering matches(assign highest score)
b_qseqid=""
while read qseqid qstart qend qlen score;do
if [ "$qseqid" == "$b_qseqid" ];then
  max=$(( b_score > score ? b_score : score ))
  if [ "$qend" -le "$b_qend" ]; then b_score=$max; ## keep first with highest score
  elif [ "$qstart" -le "$((b_qend+24))" ];then
#    echo -e "$b_qseqid\t$b_qstart\t$b_qend\t$b_qlen\t$max"; ## echo first with highest score
    b_qend=$qend;b_score=$max;## keep second with  highest score
  else
    echo -e "$b_qseqid\t$b_qstart\t$b_qend\t$b_qlen\t$b_score"; ## echo first
    b_qstart=$qstart;b_qend=$qend;b_score=$score;## keep second
  fi
elif [ "$b_qseqid" != "" ];then
  echo -e "$b_qseqid\t$b_qstart\t$b_qend\t$b_qlen\t$b_score"; ## echo first
  b_qseqid=$qseqid;b_qstart=$qstart;b_qend=$qend;b_qlen=$qlen;b_score=$score; ## keep second
else b_qseqid=$qseqid;b_qstart=$qstart;b_qend=$qend;b_qlen=$qlen;b_score=$score; ## keep first
fi; done < univec_blastn.out.sig > univec_blastn.out.sig2 # 147863
if [ "$(tail -n1 univec_blastn.out.sig2 | awk '{print $1}')" != "$b_qseqid" ];then 
 echo -e "$b_qseqid\t$b_qstart\t$b_qend\t$b_qlen\t$b_score" >> univec_blastn.out.sig2; fi ## echo last record 

## select moderate and strong hits
cat univec_blastn.out.sig2 | awk '(((($2<=25)||(($4-$3)<=25))&&$5>=19)||(($2>25)&&(($4-$3)>25)&&($5>=25)))' | sort -k1,1 -k2,2n -k3,3nr -k5,5nr > univec_blastn.out.sig3 ## Biom: 122
#confirmation#cat univec_blastn.out.sig3 | sort -u -k1,1 -k2,2 -k3,3 --merge | wc -l ## Biom: 122
awk '{print $1,$2,$3,$4;}' univec_blastn.out.sig3 > cont.report.len

cat cont.report.len | awk '($4-$3) < 201 && ($2-1) < 201 {print $1}' > excludeIDs  ## Biom: 23
## remove excludeIDs (this would remove their duplicates as well if exist)
cat excludeIDs | awk '{print $1}' | grep -v -F -w -f - cont.report.len > cont.report.len.keep ## 883 (i.e. we removed 192 transcript, none of them has dup matches)
#sort cont.report.len.keep | awk '{print $1}' | uniq -c | awk '($1>1) {print $2}' | grep -F -w -f - cont.report.len | sort > cont.report.len.dup ## 2
cat cont.report.len.keep | awk '($4-$3) <= ($2-1) {print $1,$2}' | sort -k1,1 -k2,2n | sort -u -k1,1 --merge > suffix_cut ## 462 ## cut from $2-len
cat cont.report.len.keep | awk '($4-$3) > ($2-1) {print $1,$3}' | sort -k1,1 -k2,2nr | sort -u -k1,1 --merge > prefix_cut ## 421 ## cut from 1-$3

module load QIIME/1.8.0
filter_fasta.py --input_fasta_fp $exp_transcriptome --output_fasta_fp suffix_cut.fa --seq_id_fp suffix_cut
while read line; do
  if [[ $line == \>* ]]; then
    suffix=$(echo $line | awk -F '[> ]' '{print $2}' | grep -w -f - suffix_cut | awk '{print $2}')
    echo $line;
  else echo ${line:0:$suffix-1}; fi
done < suffix_cut.fa > suffix_cut.fa.fixed ## Biom: 54

while read contig;do mark=$(echo $contig | cut -d" " -f 1); grep $mark suffix_cut.fa.fixed | sed 's/>//';done < prefix_cut > suffix_cut_pass
filter_fasta.py --input_fasta_fp $exp_transcriptome --output_fasta_fp temp_prefix_cut.fa --seq_id_fp suffix_cut --negate
filter_fasta.py --input_fasta_fp temp_prefix_cut.fa --output_fasta_fp prefix_cut.1.fa --seq_id_fp prefix_cut
filter_fasta.py --input_fasta_fp suffix_cut.fa.fixed --output_fasta_fp prefix_cut.2.fa --seq_id_fp suffix_cut_pass
cat prefix_cut.1.fa prefix_cut.2.fa > prefix_cut.fa
while read line; do
  if [[ $line == \>* ]]; then
    prefix=$(echo $line | awk -F '[> ]' '{print $2}' | grep -w -f - prefix_cut | awk '{print $2}')
    echo $line;
  else echo ${line:$prefix}; fi
done < prefix_cut.fa > prefix_cut.fa.fixed ## Biom: 44
filter_fasta.py --input_fasta_fp $exp_transcriptome --output_fasta_fp ${exp_transcriptome%.fasta}.UniVec.fasta --seq_id_fp cont.report.len --negate
univec_exp_tran=${exp_transcriptome%.fasta}.UniVec.fasta

filter_fasta.py --input_fasta_fp suffix_cut.fa.fixed --output_fasta_fp suffix_cut.fa.fixed_final --seq_id_fp suffix_cut_pass --negate
cat suffix_cut.fa.fixed_final >> $univec_exp_tran
cat prefix_cut.fa.fixed >> $univec_exp_tran  ## Biom: 144213

echo $(($(grep "^>" $exp_transcriptome | wc -l) - $(grep "^>" $univec_exp_tran | wc -l))) ## 23 (and trimmed 54+44)
################
## http://www.bioinformatics.org/cd-hit/
#cd $seqclean_dir
#module load CDHIT/4.6.1c
#awk '{print $1}' $univec_exp_tran | sed 's/TRINITY_//' > univec_exp_tranHD.fasta
#cd-hit-est -i univec_exp_tranHD.fasta -o univec_exp_tran95 -c 0.95 -n 8 -T 6 &> univec_exp_tran95.log  ## 118489  clusters
################
## repeat salmon analysis to generate backmapping stats and generate final isoform expression table
cd $Biomphalaria/abundFilter
mkdir SalmonRun2
mv *.quant salmonQuant_summary.txt *.quant.counts transcripts.lengthes gene_transcript_map *.dataSummary_comp targets_list *TPM SalmonRun2/.

cd $seqclean_dir
qsub -v index="salmon_index3",transcriptome="$univec_exp_tran" ${script_path}/salmonIndex.sh

cd $Biomphalaria/abundFilter
while read identifier;do
 ls ${identifier}.s_pe.fq.1 ${identifier}.s_pe2.fq ${identifier}.s_se.fq
 qsub -v index="$seqclean_dir/salmon_index3",identifier=$identifier ${script_path}/salmonQuant_PE.sh
 qsub -v index="$seqclean_dir/salmon_index3",identifier=$identifier ${script_path}/salmonQuant_SE.sh
done < $identifiers
## finish the abundance exp module to generate isoform expression table for the later annotation step
find ./*.quant -name *.sf -exec grep -H "mapping rate" {} \; | sort > salmonQuant_summary.txt
python $script_path/gather-counts2.py -i "$(pwd)"
echo "transcript"$'\t'"length" > transcripts.lengthes
sf=$(find ./*.quant -name \*.sf | head -n1)
cat $sf | grep -v "^#" | awk -F "\t" -v OFS='\t' '{print $1,$2}' >> transcripts.lengthes
# generation of gene map
perl ${script_path}/get_Trinity_gene_to_trans_map.pl $univec_exp_tran > gene_transcript_map  ## no of gene 128739 / isoforms 144213
#module load R/3.0.1
> targets_list
while read identifier;do
  echo $(pwd) $identifier
  bash $script_path/run_calcTPM.sh "$(pwd)" "$identifier" "transcripts.lengthes" "gene_transcript_map" ${script_path}/calcTPM2.R
  #Rscript ${script_path}/calcTPM_tis.R "$(pwd)" "$identifier" "transcripts.lengthes" "gene_transcript_map" >> targets_list
done < $identifiers
bash $script_path/abund_est.sh
## print a new version of isoform expression table with rounded values
cat allTissues_isoformTPM | awk '{x=sprintf("%.2f", $2);print $1,x;}' >  allTissues_isoformTPM.rounded
module load R/3.0.1
Rscript ${script_path}/exp.R "allTissues_isoformTPM.rounded" > exp_assessment
## generate assembly statistics
cd $seqclean_dir
module load Bioperl/1.6.923
perl ${script_path}/seq_stats.pl $univec_exp_tran >  $univec_exp_tran.MatzStat

module load trinity/2.2.0
TrinityStats.pl $univec_exp_tran > $univec_exp_tran.TrinityStat
#################
## check of assembly completness 
cd ~
git clone https://gitlab.com/ezlab/busco.git
cd busco/
wget http://busco.ezlab.org/datasets/metazoa_odb9.tar.gz
tar -zxvf metazoa_odb9.tar.gz
module load BLAST+/2.2.31
module load HMMER/3.1b2
module load augustus/3.2.2
mkdir augustus
cp -r /opt/software/augustus/3.2.2--GCC-4.4.5/config augustus/.
chmod -R 755 augustus/config
export AUGUSTUS_CONFIG_PATH="/busco/augustus/config/"
cd $seqclean_dir
python ~/busco/BUSCO.py -i Trinity.clean.201.exp2.UniVec.fasta -l ~/busco/metazoa_odb9 -o busco_repot -m tran --cpu 4 

#################
## Functional annotation of whole transcriptome assembly
#################
## install crb-blast
wget https://pypi.python.org/packages/d4/0c/9840c08189e030873387a73b90ada981885010dd9aea134d6de30cd24cb8/virtualenv-15.1.0.tar.gz
tar xzf virtualenv*.gz
rm virtualenv*.gz
cd virtualenv-15*
python2.7 virtualenv.py ../crbBlastenv
source $HOME/crbBlastenv/bin/activate                 ## may be the virtualenv is not needed
cd ../
command curl -sSL https://rvm.io/mpapis.asc | gpg2 --import -
\curl -sSL https://get.rvm.io | bash -s stable --ruby
source /mnt/home/mansourt/.rvm/scripts/rvm
gem install crb-blast
# WARNING: You have '~/.profile' file, you might want to load it,
#    to do that add the following line to '/mnt/home/mansourt/.bash_profile':
#     source ~/.profile
#################
## a) generation of gene map
#cd $seqclean_dir
#perl ${script_path}/get_Trinity_gene_to_trans_map.pl $univec_exp_tran > gene_trans_map
#gene_transcript_map=$seqclean_dir/gene_trans_map

## b) Generation of possible ORFs
module load TransDecoder/2.0.1
TransDecoder.LongOrfs -t $univec_exp_tran
LongOrfs=$seqclean_dir/$(basename $univec_exp_tran).transdecoder_dir/longest_orfs.pep ## supplementary file 1
while read line; do
  if [[ $line == \>* ]]; then
    echo $line | awk -F'[>|: ]' '{print $2,$7","$9","$13}';fi
done < $LongOrfs > $LongOrfs.summary ## 48327
cat $LongOrfs.summary | awk '{if(a[$1]==""){a[$1] = $2}else{a[$1] = a[$1]";"$2;}};END{for(i in a){print i,a[i]}}' > $LongOrfs.summary2 ## 40408
grep "complete" $LongOrfs.summary2 | wc -l ## 8678

## c) Run blast to search Trinity transcripts & Transdecoder-predicted proteins
module load BLAST+/2.2.30
## prepare reference database
mkdir ${Biomphalaria}/uniprot
cd ${Biomphalaria}/uniprot
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz
makeblastdb -in uniprot_sprot.fasta -dbtype prot
target="uniprot_sprot"
#wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz
#gunzip uniref90.fasta.gz
#makeblastdb -in uniref90.fasta -dbtype prot

## praper denovo transcriptome database
cd $seqclean_dir
makeblastdb -in $univec_exp_tran -dbtype nucl
query=$(basename ${univec_exp_tran%.fasta})

## blast denovo transcriptome aganist target database
mkdir $seqclean_dir/blastx_dir
cd $seqclean_dir/blastx_dir
cp $univec_exp_tran trans.fa
perl ${script_path}/splitFasta.pl trans.fa 300  ## 500 for uniprot_uniref90

#DB=/mnt/ls15/scratch/users/mansourt/Tamer/Biomphalaria/uniprot/uniprot_sprot.fasta
for f in subset*_trans.fa; do
 ## header of blast output follows the crb-blast format= qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen
 qsub -v input=$f,DB=${Biomphalaria}/uniprot/uniprot_sprot.fasta ${script_path}/blastx_targetDB.sh;
 # blastx -query $f -db $DB -out $f.bx -outfmt "6 std qlen slen" -num_threads 4 -evalue 1e-5 -max_target_seqs 50 -seg no
done
cat subset*_trans.fa.bx > ../${query}_into_${target}.1.blast

## blast reference proteome aganist denovo transcriptome database
mkdir ${Biomphalaria}/uniprot/tblastn_dir
cd ${Biomphalaria}/uniprot/tblastn_dir
cp ../uniprot_sprot.fasta uniprot.fa
perl ${script_path}/splitFasta.pl uniprot.fa 150  ## 500 for uniprot_uniref90

#DB=$univec_exp_tran
for f in subset*_uniprot.fa; do
 ## header of blast output follows the crb-blast format= qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen
 qsub -v input=$f,DB=$univec_exp_tran ${script_path}/tblastn_targetDB.sh;
 #tblastn -query $f -db $DB -out $f.tbn -outfmt "6 std qlen slen" -num_threads 1 -evalue 1e-5 -max_target_seqs 50 -seg no
done
cat subset*_uniprot.fa.tbn > ../${target}_into_${query}.2.blast

## run crb-blast
mkdir $seqclean_dir/crb_blast
cd $seqclean_dir/crb_blast
head $univec_exp_tran > query.fasta
head ${Biomphalaria}/uniprot/uniprot_sprot.fasta > target.fasta
touch query.nin query.nhr query.nsq
touch target.fasta target.psq target.pin target.phr
cp $seqclean_dir/${query}_into_${target}.1.blast query_into_target.1.blast
cp ${Biomphalaria}/uniprot/${target}_into_${query}.2.blast target_into_query.2.blast
#source /mnt/home/mansourt/.rvm/scripts/rvm
crb-blast --query query.fasta --target target.fasta --threads 4 --output output.tsv   ## supplementary file 2

## generate a liftover homolog list
sort -k1,1 -k12,12nr -k11,11n  query_into_target.1.blast | sort -u -k1,1 --merge > query_into_target.1.blast.best
awk '{print $1}' output.tsv | grep -vFwf - query_into_target.1.blast.best > homologs.blast  ## supplementary file 3

## generate a liftover ORF list
awk '{print $1}' query_into_target.1.blast.best | grep -vFwf - $LongOrfs.summary2 > ORF.noBlast

## stats
awk '{print $1}' query_into_target.1.blast | sort | uniq | wc -l ## no of transcripts with significant blast hit  ## 32056
wc -l output.tsv      ## no of transcripts with orthologs                  ## 15246 
wc -l homologs.blast  ## no of transcripts with homologs                   ## 16810
wc -l ORF.noBlast     ## no of transcripts with long ORF but no blast hits ## 15161
grep "complete" ORF.noBlast | wc -l ## 2627

## transfer the annotation to the fasta files
cd $seqclean_dir
grep "^>" ${Biomphalaria}/uniprot/uniprot_sprot.fasta | sed 's/>//' > ${Biomphalaria}/uniprot/uniprot_sprot.info
ref_info=${Biomphalaria}/uniprot/uniprot_sprot.info
exp_report=$Biomphalaria/abundFilter/allTissues_isoformTPM.rounded
crb_report=crb_blast/output.tsv
homolog_report=crb_blast/homologs.blast
ORF_report=crb_blast/ORF.noBlast
while read line; do
  if [[ $line == \>* ]]; then
    id=$(echo $line | awk -F '[> ]' '{print $2}')
    len=$(echo $line | awk -F '[> ]' '{print $3}')
    exp=$(grep -w $id $exp_report | awk '{print $2}')
    ortholog=$(grep -w $id $crb_report | awk -F '\t' '{print $2}' | grep -wf - $ref_info)
    if [ "$ortholog" != "" ];then
      echo ">"$id $len "exp="$exp "ortholog=["$ortholog"]";
    else 
      homolog=$(grep -w $id $homolog_report | awk -F '\t' '{print $2}' | grep -wf - $ref_info)
      if [ "$homolog" != "" ];then
        echo ">"$id $len "exp="$exp "homolog=["$homolog"]";
      else
        ORF=$(grep -w $id $ORF_report | awk '{print $2}')
        if [ "$ORF" != "" ];then
         echo ">"$id $len "exp="$exp "ORF=["$ORF"]";
        else
         echo ">"$id $len "exp="$exp; fi; fi; fi; 
  else echo $line; fi
done < $univec_exp_tran > ${univec_exp_tran%.fasta}.ann.fasta ## supplementary file 4
ann_exp_transcriptome=${univec_exp_tran%.fasta}.ann.fasta

#############################################################################
## annotation of high expression transcripts
mkdir $seqclean_dir/hiExp
cdmcd $seqclean_dir/hiExp
awk '$2>100' $exp_report > hiExp
awk '{print $1}' hiExp | grep -Fwf - $LongOrfs.summary2 > LongORFs.hiExp  
awk '{print $1}' hiExp | grep -Fwf - ../${query}_into_${target}.1.blast > blast.hiExp
cat hiExp | awk '{print $1}' | grep -Fwf - ../$crb_report > ortholog.hiExp
cat hiExp | awk '{print $1}' | grep -Fwf - ../$homolog_report > homolog.hiExp
cat hiExp | awk '{print $1}' | grep -Fwf - ../$ORF_report > ORF.hiExp
wc -l *hiExp ## 848 hiExp /  487 LongORFs.hiExp / 336 ortholog.hiExp / 92 homolog.hiExp / 145 ORF.hiExp
grep "complete" LongORFs.hiExp | wc -l ## 258
grep "complete" ORF.hiExp | wc -l ## 88
awk '{print $1}' blast.hiExp | sort | uniq | wc -l ## no of transcripts with significant blast hit  ## 428

module load QIIME/1.8.0
filter_fasta.py --input_fasta_fp $ann_exp_transcriptome --output_fasta_fp hiExp.fasta --seq_id_fp hiExp

makeblastdb -in hiExp.fasta -dbtype nucl

## blast denovo transcriptome aganist target database
mkdir $seqclean_dir/hiExp/blastx_dir && cd $seqclean_dir/hiExp/blastx_dir
cp ../hiExp.fasta trans.fa
perl ${script_path}/splitFasta.pl trans.fa 100  ## 500 for uniprot_uniref90

#module load BLAST+/2.2.30
for f in subset*_trans.fa; do
 qsub -v input=$f,DB="/mnt/research/common-data/Bio/blastdb/nr" ${script_path}/blastx_targetDB.sh;
 #blastx -query $f -db /mnt/research/common-data/Bio/blastdb/nr -out $f.bx -outfmt "6 std qlen slen" -num_threads 8 -evalue 1e-5 -max_target_seqs 50 -seg no
done
cat subset*_trans.fa.bx > ../query_into_target.1.blast

## blast reference proteome aganist denovo transcriptome database
mkdir $seqclean_dir/hiExp/tblastn_dir && cd $seqclean_dir/hiExp/tblastn
#ln -s /mnt/research/common-data/Bio/blastdb/FASTA/nr.fa .
awk '{print $1}' /mnt/research/common-data/Bio/blastdb/FASTA/nr.fa > nr.fa
perl ${script_path}/splitFasta.pl nr.fa 100 

DB="$seqclean_dir/hiExp/hiExp.fasta"
for f in subset1[6-9]_nr.fa; do
 qsub -v input=$f,DB=$DB ${script_path}/tblastn_targetDB.sh;
 #tblastn -query $f -db $DB -out $f.tbn -outfmt "6 std qlen slen" -num_threads 4 -evalue 1e-5 -max_target_seqs 50 -seg no
done
cat subset*_nr.fa.tbn > ../target_into_query.2.blast

## run crb-blast
mkdir $seqclean_dir/hiExp/crb_blast
cd $seqclean_dir/hiExp/crb_blast
head ../hiExp.fasta > query.fasta
head /mnt/research/common-data/Bio/blastdb/FASTA/nr.fa | awk '{print $1}'  > target.fasta
touch query.nin query.nhr query.nsq
touch target.fasta target.psq target.pin target.phr
ln -s ../query_into_target.1.blast .
ln -s ../target_into_query.2.blast .
#source /mnt/home/mansourt/.rvm/scripts/rvm
crb-blast --query query.fasta --target target.fasta --threads 4 --output hiExp_output.tsv  ## possible supplementary file

## generate a liftover homolog list
sort -k1,1 -k12,12nr -k11,11n  query_into_target.1.blast | sort -u -k1,1 --merge > query_into_target.1.blast.best
awk '{print $1}' hiExp_output.tsv | grep -vFwf - query_into_target.1.blast.best > hiExp_homologs.blast ## possible supplementary file

## generate a liftover ORF list
awk '{print $1}' query_into_target.1.blast.best | grep -vFwf - ../LongOrfs.hiExp > hiExp_ORF.noBlast

## stats
awk '{print $1}' query_into_target.1.blast | sort | uniq | wc -l ## no of transcripts with significant blast hit  ## 542
wc -l hiExp_output.tsv      ## no of transcripts with orthologs                  ## 469  
wc -l hiExp_homologs.blast  ## no of transcripts with homologs                   ## 73
wc -l hiExp_ORF.noBlast     ## no of transcripts with long ORF but no blast hits ## 61
grep "complete" hiExp_ORF.noBlast | wc -l ## 41

## transfer the annotation to the fasta files
cd $seqclean_dir/hiExp
grep "^>" /mnt/research/common-data/Bio/blastdb/FASTA/nr.fa | sed 's/^>//' > nr.info
nr_ref_info=nr.info
#exp_report=$Biomphalaria/abundFilter/allTissues_isoformTPM.rounded
nr_crb_report=crb_blast/hiExp_output.tsv
nr_homolog_report=crb_blast/hiExp_homologs.blast
nr_ORF_report=crb_blast/hiExp_ORF.noBlast
while read line; do
  if [[ $line == \>* ]]; then
    id=$(echo $line | awk -F '[> ]' '{print $2}')
    len=$(echo $line | awk -F '[> ]' '{print $3}')
    exp=$(echo $line | awk -F '[> ]' '{print $4}')
    #exp=$(grep -w $id $exp_report | awk '{print $2}')
    ortholog=$(grep -w $id $nr_crb_report | awk -F '\t' '{print $2}' | grep -wf - $nr_ref_info)
    if [ "$ortholog" != "" ];then
      echo ">"$id $len "exp="$exp "ortholog=["$ortholog"]";
    else
      homolog=$(grep -w $id $nr_homolog_report | awk -F '\t' '{print $2}' | grep -wf - $nr_ref_info)
      if [ "$homolog" != "" ];then
        echo ">"$id $len "exp="$exp "homolog=["$homolog"]";
      else
        ORF=$(grep -w $id $nr_ORF_report | awk '{print $2}')
        if [ "$ORF" != "" ];then
         echo ">"$id $len "exp="$exp "ORF=["$ORF"]";
        else
         echo ">"$id $len "exp="$exp; fi; fi; fi;
  else echo $line; fi
done < hiExp.fasta > hiExp.nr_ann.fasta  ## supplementary file 5

## upload files to iPlant
mkdir -p ~/temp/Biomphalaria
cd ~/temp/Biomphalaria
cp $LongOrfs supp1_LongOrfs.fa
cp $seqclean_dir/crb_blast/output.tsv supp2_allTrans_orthologs.tab
cp $seqclean_dir/crb_blast/homologs.blast supp3_allTrans_homologs.tab
cp $ann_exp_transcriptome supp4_allTrans_ann.fasta
cp $seqclean_dir/hiExp/hiExp.nr_ann.fasta supp5_hiTrans_ann.fasta

iinit
imkdir /iplant/home/drtamermansour/Biomphalaria
icd /iplant/home/drtamermansour/Biomphalaria
iput * .

