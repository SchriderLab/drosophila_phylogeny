# Drosophila Phylogeny Pipeline 

**BUSCO processing**
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
**Generation of MSAs using [MAFFT](https://mafft.cbrc.jp/alignment/software/) v7.427 L-INS-i strategy (most accurate)**
```bash
mafft --localpair --maxiterate 1000 [in_fasta]
```
Note: 
1) for large MSAs minimal SLUM requirements 

#SBATCH --time=100:00:00 --ntasks=15 --nodes=1 --mem-per-cpu=5G
```bash
mafft --localpair --maxiterate 1000 --thread 15 [in_fasta]
```
2) for ultra-large MSAs minimal SLUM requirements
#SBATCH --time=100:00:00 --ntasks=30 --nodes=1 --mem-per-cpu=10G
```bash
mafft --auto --thread 30 [in_fasta]
```
