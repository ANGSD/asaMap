# asaMap

Article
http://biorxiv.org/content/early/2015/01/22/014001

The implementation is still a work in progress. More features will be added at some point

# Install
git clone https://github.com/angsd/asaMap.git;

cd asaMap

make
#Options
```
   -p '(null)'		plink prefix filename
   -o '(null)'		output filename
   -c '(null)'		covariance matrix filename
   -y '(null)'		phenotypes
   -a '(null)'		admixproportions (for source pop1)
   -f '(null)'		allele frequencies

 optional arguments:
   -m '0'		model 0=add 1=rec
   -b '(null)'		file containing the start
   -i '10'		max number of iterations
   -0 '0'		full 1:M1,2:M2,3:M3
   -1 '0'		null 1:M1,2:M2,3:M3, 4:M4 5:M5
   -r '100'		random seed
   -t '1.000000e-08'	float for breaking EM update
   -P '1'		number of threads
```
All files must be specified: -p -c -y -a -f -o

#Run example using ADMIXTURE 

```
PH=pheno.files
COV=cov.file
PF=plinkFile #without postfix e.g. no .bim / .fam / .bed
#run admixture
admixture $PF.bed 2

#get admixture proportions for population 1
cut -f1 -d" " $PF.2.Q > Q

#run asaMap
../asaMap -p $PF -o out -c $COV -y $PH -a Q -f $PF.2.P
```

### p-values
Easy to optain in R

```
res <- read.table("out",head=T)
pval <- function(a,b,df)
  1- pchisq(-2*(a-b),df=df)
#pvalues for M1 vs. M5
pval(res$llh.M1.,res$llh.M5.,df=2)

#pvalues for M1 vs. M4
pval(res$llh.M1.,res$llh.M4.,df=1)
```

###models

| model | parameters | notes |
| --- | --- | --- |
| M1 |(\beta_1,\beta_2)\in R^2  | population specific effects |
| M2 | \beta_1=0,\beta_2\in R | no effect in population  |
| M3 | \beta_1\in R, \beta_2=0 | no effect in population 2 |
| M4 |\beta_1=\beta_2\in R | same effect in both populations |
| M5 | \beta_1=\beta_2=0 | no effect in any population |


#Input files
###Genotypes
plink binary (.bim .bam .fam)

###Phenotypes (response)
A file with each individuals phenotypes on each line. e.g. 
```
>head pheno 
-0.712027291121767
-0.158413122435864
-1.77167888612947
-0.800940619551485
0.3016297021294
0.596892506547882
-0.661786423692485
-0.405728519330873
-1.04224674183241
0.0881848860116932
```
### extra covariates (in addition to the intersept and genotypes)
A file where each column is a covariate and each row is an individual
```
>head cov
0.0127096117618385 -0.0181281029917176 -0.0616739439849275 -0.0304606694443973
0.0109944672768584 -0.0205785925514037 -0.0547523583405743 -0.0208813157640705
0.0128395346453956 -0.0142116856067135 -0.0471689997039534 -0.0266186436009881
0.00816783754598649 -0.0189271733933446 -0.0302259313905976 -0.0222247658768436
0.00695928218989132 -0.0089960963981644 -0.0384886176827146 -0.0126490197701687
0.00908359304129912 -0.019562503526549 -0.0276058491506046 -0.0202388414332682
0.0193657006952317 -0.0219605099975189 -0.050537627417191 -0.0236411635865132
0.015862252334236 -0.0134969241244036 -0.0336244748700029 -0.0222294652006281
0.0194100156955457 -0.0371103372950621 0.00813012568415838 -0.015311879434991
0.0190516629849255 -0.0194012185542486 -0.0413589828106922 -0.0292318169458017
```
###admixture proportions
Right now the implementation is just for two populations and the admixture porportions is a file with a column consisting or the admixture proportions for population 1. 
```
>head Q
0.398834
0.491314
0.399352
0.634947
0.000010
0.467470
0.497885
0.800272
0.999990
0.032108
```

### ancestral allele frequencies
A files where each column is a population and each row is a SNP. 
```
>head P
0.907268 0.903723
0.498913 0.543486
0.865254 0.819332
0.737799 0.738876
0.812931 0.823395
0.783896 0.785741
0.711677 0.865518
0.592533 0.587921
0.834116 0.850966
0.837011 0.852711
```
