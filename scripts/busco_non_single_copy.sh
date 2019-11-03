#! /bin/bash
#Separate single copy BUSCO genes fron non-single
mkdir non_single_copy_busco_sequences
#Check total number of loci
ls single_copy_busco_sequences/*.fna | wc -l
#Check total number of sequences
grep ">" single_copy_busco_sequences/*.fna | wc -l
#Move DNA and PROT fastas with more than one entry to non_single_copy_busco_sequences 
for i in single_copy_busco_sequences/*
do
l=$(grep ">" $i | wc -l)
if (( $l > 1 ))
then
mv $i non_single_copy_busco_sequences
fi
done
