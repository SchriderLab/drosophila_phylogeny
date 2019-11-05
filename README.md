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
:red_circle: **Generation of MSAs using [MAFFT](https://mafft.cbrc.jp/alignment/software/) v7.427 L-INS-i strategy (most accurate) with minimal SLURM requirements**

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
:red_circle: **Concatenation of BUSCO MSAs into a Supermatrix using geneSticher.py from [Utensils](https://github.com/ballesterus/Utensils)**  
```bash
python2.7 geneSticher.py -in *.aln
```

:red_circle: **Trimming Supermatrix by removing sites that have only one non-gap character using [trimAl](http://trimal.cgenomics.org/introduction)**
```bash
trimal -in SuperMatrix.al -out SuperMatrix.trim.al -gt 0.0164
```
:red_circle: **ML phylogenetic inference using [IQTREE](http://www.iqtree.org/) v1.6.5**

#SBATCH --time=100:00:00 --ntasks=15 --nodes=1 --mem-per-cpu=10G
```bash
iqtree -s SuperMatrix.trim.al -nt 15 -m GTR+I+G -bb 1000 -pre ML -safe -bb 1000 -alrt 1000 -abayes
```
