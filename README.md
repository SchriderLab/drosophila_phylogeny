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
iqtree -s [supermatrix.aln] -nt 15 -m GTR+I+G -pre ML -safe -bb 1000 -alrt 1000 -abayes
```
:red_circle: **ASTRAL phylogenetic inference using [ASTRAL](https://github.com/smirarab/ASTRAL) v5.6.3**
```bash
java -jar astral.5.6.3.jar -i [gene_trees] -t 3
```
# Drosophila Introgression Analysis 
:red_circle: **DCT/BLT** using [blt_dct_test.r](https://github.com/SchriderLab/drosophila_phylogeny/blob/master/scripts/blt_dct_test.r) 
```bash
Rscript blt_dct_test.r -t [gene_trees] -s [species tree] -n [node number] -c [cores] -p [prefix] -o [outgroup species]
```
```
Options:
	-t CHARACTER, --trees=CHARACTER
		input gene trees in newick

	-s CHARACTER, --species_tree=CHARACTER
		rooted species tree in newick

	-n NUMERIC, --node=NUMERIC
		node number of a tested clade in a species tree

	-c NUMERIC, --cores=NUMERIC
		number of cores for parallel computing

	-p CHARACTER, --prefix=CHARACTER
		clade prefix

	-o CHARACTER, --outgroup=CHARACTER
		outgroup species (only one allowed)

	-h, --help
		Show this help message and exit
```    
Example for clade 1 (use ```nodelabels()``` function in [ape](https://cran.r-project.org/web/packages/ape/index.html) package to obtain node number):
```bash
Rscript blt_dct_test.r -t gene_trees_noboot_dna_mafft -s iqtree_rooted.tre -n 307 -c 10 -p C1 -o Anopheles_gambiae
```
Interpretaton of the output for a given triplet with taxa A,B,C  
```clade```: clade ID   
```P1out```: outgroup species 1 e.g. if P1out = A, then topology is (A,(B,C))   
```P2out```: outgroup species 2 e.g. if P1out = B, then topology is (B,(A,C))  
```P3out```: outgroup species 3 e.g. if P1out = C, then topology is (C,(A,B))  
```CountP1```:  gene tree counts for triplet with the outgroup species A, i.e. (A,(B,C))  
```CountP2```:  gene tree counts for triplet with the outgroup species B, i.e. (B,(A,C))  
```CountP3```:  gene tree counts for triplet with the outgroup species C, i.e. (C,(A,B))   
```PvalueChi```: Chi-square test P value form comparison of two discordant gene tree counts (e.g. if CountP1 = 140, CountP2 = 23 and CountP3 = 1569, then P value comes from the comparison of CountP1 and CountP2)    
```meanT_concord```: mean value of the genetic distances between sister species obtined from all concordant gene trees (e.g. if CountP1 = 140, CountP2 = 23 and CountP3 = 1569, then meanT_concord comes from the gene trees with the P3out outgroup)    
```meanT_discord1```: mean value of the genetic distances between sister species obtined from all discordant gene trees with the smallest count (e.g. if CountP1 = 140, CountP2 = 23 and CountP3 = 1569, then meanT_discord1 comes from the gene trees with the P2out outgroup)    
```meanT_discord2```: mean value of the genetic distances between sister species obtined from all discordant gene trees with the intermediate count (e.g. if CountP1 = 140, CountP2 = 23 and CountP3 = 1569, then meanT_discord2 comes from the gene trees with the P1out outgroup)  
```PvalueWCOMC1```: Mann–Whitney U test P value from comparison of meanT_concord and meanT_discord1  
```PvalueWCOMC2```: Mann–Whitney U test P value from comparison of meanT_concord and meanT_discord2  
```PvalueWC1C2```: Mann–Whitney U test P value from comparison of meanT_discord1 and meanT_discord2  


