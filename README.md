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
