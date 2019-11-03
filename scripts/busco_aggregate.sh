ls *_*/ | sort | uniq | grep ".fna" > uniq_dna_BUSCO
ls *_*/ | sort | uniq | grep ".faa" > uniq_prot_BUSCO
cat uniq_dna_BUSCO uniq_prot_BUSCO > uniq_BUSCO
all_lines=`cat uniq_BUSCO`
mkdir ../busco_clusters
for bu in $all_lines
do
echo $bu
cat *_*/$bu > ../busco_clusters/$bu
done
