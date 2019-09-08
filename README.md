# eQTLMAPT 
**eQTLMAPT** - eQTL Mediaiton Analysis with accelerated Permutation Testing approaches

## Introduction
`eQTLMAPT` uses R language to achieve fast and accurate eQTL mediation analysis with efficient permutation procedures to control for multiple testing. eQTLMAPT provides adaptive permutation scheme which prunes the permutation process opportunely, and it models the null distribution using generalized Pareto distribution (GPD) trained from a few permutations. The adjusted P-values can be estimated at any significance level based on the null distribution with much less computational cost compared with traditional fixed number permutation strategies. In addition, eQTLMAPT provides flexible options for users to choose proper permutation schemes combining with various confounding factors adjustment methods based on their practical needs. Experiments on real eQTL dataset demonstrate that eQTLMAPT provides higher resolution of estimated significance and is order of magnitude faster than compared methods with similar accuracy.

## Installation
    install.packages("devtools", dependencies = T)  
    devtools::install_github("hitbc/eQTLMAPT")

## Input
`snp.dat`  A eQTL genotype matrix.  Each row is an eQTL, each column is a sample.  
`fea.dat`  A feature profile matrix. Each row is for one feature profile, each column is a sample.  
`conf/known.conf`  A known confounders matrix which is adjusted in all mediation tests. Each row is a confounder, each column is a sample.   
`cov.pool`  The pool of candidate confounding variables from which potential confounders are adaptively selected to adjust for each mediation test. Each row is a covariate, each column is a sample.  
`trios.idx`  The matrix of selected trios indexes (row numbers) for mediation tests. Each row consists of the index (i.e., row number) of the eQTL in eQTL genotype matrix, the index of cis-gene transcript in feature profile matrix, and the index of trans-gene in feature profile matrix.  
`cl`  If parallel computing is required, cluster information needs to be provided.  
For other parameter information, refer specifically to [eQTLMAPT](https://github.com/hitbc/eQTLMAPT).  

## Output
`nperm`  If adaptive permutation scheme is adopted, the actual permutation number is output.  
`nominal.p`  The nominal P-values.  
`empirical.p`  The empirical P-Value.  
`empirical.p.gpd`  The empirical P-value obtained by GPD fitting.  
`std.error`  The std.error of cis-gene in liner regression.  
`t_stat`  The t_stat of cis-gene in liner regression.  
`beta`  The beta of SNP in liner regression.  
`beta.total`  The beta of cis-gene in liner regression.  
`beta.change`  The proportions mediated.  
`pc.matrix`  PCs will be returned if the PCs based on expression data are used as the pool of potential confounders. Each column is a PC.  
`sel.conf.ind`  An indicator matrix with dimension of the number of trios by the number of covariates in `cov.pool` or `pc.matrix` if the principal components (PCs) based on expression data are used as the pool of potential confounders.  


## Citation
Paper are published at [paper]().

## Contact
If you need help, please contact **ydwang@hit.edu.cn** or **jiajiepeng@nwpu.edu.cn**。
