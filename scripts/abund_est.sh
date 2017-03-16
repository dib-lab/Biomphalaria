#script_path=$(dirname "${BASH_SOURCE[0]}")

## merge the TPM for tissue specific files
targets=()
i=1
rm -f temp.*
for f in *.dataSummary_comp;do
  cat $f | tail -n+2 | awk '{print $2,$7}' | uniq > $f.gene
  cat $f | tail -n+2 | awk '{print $3,$6}' > $f.isoform
  target=${f%.dataSummary_comp}
  targets+=($target)
  if [ $i -eq 1 ];then cat $f.gene > temp.$i;else join -t" " --nocheck-order temp.$((i-1)) $f.gene > temp.$i;fi
  if [ $i -eq 1 ];then cat $f.isoform > isotemp.$i;else join -t" " --nocheck-order isotemp.$((i-1)) $f.isoform > isotemp.$i;fi
  ((i+=1))
done
#echo "geneName" "${targets[@]}" > allTissues_geneTPM
mv temp.$((i-1)) allTissues_geneTPM                    
#echo "isoformName" "${targets[@]}" > allTissues_isoformTPM
mv isotemp.$((i-1)) allTissues_isoformTPM                
rm -f *.gene *.isoform *temp.*

echo "geneName" "${targets[@]}" > unexp_geneTPM
cat allTissues_geneTPM | awk '$2==0' >> unexp_geneTPM ## 7777
echo "geneName" "${targets[@]}" > unexp_isoformTPM
cat allTissues_isoformTPM | awk '$2==0' >> unexp_isoformTPM ## 11380

## merge the rawCounts for tissue specific files
#targets=()
#i=1
#rm temp.*
#for f in *.dataSummary_comp;do
#  cat $f | tail -n+2 | awk '{print $3,$5}' > $f.isoform
#  target=${f%.dataSummary_comp}
#  targets+=($target)
#  if [ $i -eq 1 ];then cat $f.isoform > isotemp.$i;else join -t" " --nocheck-order isotemp.$((i-1)) $f.isoform > isotemp.$i;fi
#  ((i+=1))
#done
#mv isotemp.$((i-1)) allTissues_isoformRaw                   ## 880285 (11379 has no exp)
#rm *.isoform *temp.*


