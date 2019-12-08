# Drosophila Phylogeny Pipeline 

:red_circle: **BUSCO processing**
1) Separate single copy BUSCOs from non single copy BUSCOs    
```bash
busco_non_single_copy.sh
```
2) Create directories with "clean" IDs     
```bash
busco_copy_rename.sh
```     
3) Aggregate BUSCOs into clusters of 1:1 orthologs 
```bash
busco_aggregate.sh
```
:red_circle: **Generation of DNA MSAs using [MAFFT](https://mafft.cbrc.jp/alignment/software/) v7.427 L-INS-i strategy (most accurate) with minimal SLURM requirements**

#SBATCH --time=05:00:00 --ntasks=1 --nodes=1 --mem-per-cpu=3G
```bash
mafft --localpair --maxiterate 1000 [in_fasta] > out_fasta.aln
```
Note: 
1) for large MSAs minimal SLUM requirements 

#SBATCH --time=100:00:00 --ntasks=15 --nodes=1 --mem-per-cpu=5G
```bash
mafft --localpair --maxiterate 1000 --thread 15 [in_fasta] > out_fasta.aln
```
2) for ultra-large MSAs minimal SLUM requirements (-- auto option will most likely choose the [E-INS-i](https://mafft.cbrc.jp/alignment/software/manual/manual.html) model)

#SBATCH --time=100:00:00 --ntasks=30 --nodes=1 --mem-per-cpu=10G
```bash
mafft --auto --thread 30 [in_fasta] > out_fasta.aln
```
:red_circle: **Trimming locus MSAs by removing sites that have less than 3 non-gap characters using [fasta_site_trim.py](https://github.com/SchriderLab/drosophila_phylogeny/blob/master/scripts/fasta_site_trim.py)**
```bash
fasta_site_trim.py --Nbase 3 --input [in_fasta]
```
:red_circle: **Gene (locus) tree inference using [IQTREE](http://www.iqtree.org/) v1.6.5**
```bash
iqtree -s [in_fasta.trimmed] -bb 1000 -nt 2 -m GTR+I+G -blmin 1e-300 -safe
```


:red_circle: **Concatenation of BUSCO MSAs into a Supermatrix using geneSticher.py from [Utensils](https://github.com/ballesterus/Utensils)**  
```bash
python2.7 geneSticher.py -in *.aln
```

:red_circle: **ML phylogenetic inference using [IQTREE](http://www.iqtree.org/) v1.6.5**

#SBATCH --time=100:00:00 --ntasks=15 --nodes=1 --mem-per-cpu=10G
```bash
iqtree -s SuperMatrix.trim.al -nt 15 -m GTR+I+G -bb 1000 -pre ML -safe -bb 1000 -alrt 1000 -abayes
```
:red_circle: **Generation of CODON MSAs using [PRANK](http://wasabiapp.org/software/prank/) v.170427 with minimal SLURM requirements**

1) Ectract CDS from BUSCO-extracted DNA sequences in fasta that correspond to BUSCO-translated Protein sequences with [dna2cds_busco.py](https://github.com/SchriderLab/drosophila_phylogeny/blob/master/scripts/dna2cds_busco.py)
```bash
dna2cds_busco.py -p [in_fasta] 
```
2) Align 
```bash
prank -codon -F -d=[in_fasta] -o=[out_fasta]
```





