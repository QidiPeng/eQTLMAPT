# eQTLMAPT 
**eQTLMAPT** - eQTL Mediaiton Analysis with accelerated Permutation Testing approaches

## Introduction
eQTLMAPT uses R language to achieve fast and accurate eQTL mediation analysis with efficient permutation procedures to control for multiple testing. eQTLMAPT provides adaptive permutation scheme which prunes the permutation process opportunely, and it models the null distribution using generalized Pareto distribution (GPD) trained from a few permutations. The adjusted P-values can be estimated at any significance level based on the null distribution with much less computational cost compared with traditional fixed number permutation strategies. In addition, eQTLMAPT provides flexible options for users to choose proper permutation schemes combining with various confounding factors adjustment methods based on their practical needs. Experiments on real eQTL dataset demonstrate that eQTLMAPT provides higher resolution of estimated significance and is order of magnitude faster than compared methods with similar accuracy.

## Installation
    install.packages("devtools", dependencies = T)  
    devtools::install_github("hitbc/eQTLMAPT")

## Input data
* The eQTL genotype matrix consists of one eQTL per row and one sample per column.  
* The feature profile matrix consists of one feature profile per row and one sample per column.  
* The confounders matrix which consists of one confounder per row and ne sample per column. There are two types of confounders matrix: known confounders matrix which is adjusted in all mediation tests, and a pool of candidate confounding variables from which potential confounders are adaptively selected to adjust for each mediation test.  
* The matrix of selected trios indexes (row numbers) for mediation tests. Each row consists of the index (i.e., row number) of the eQTL in eQTL genotype matrix, the index of cis-gene transcript in feature profile matrix, and the index of trans-gene in feature profile matrix.  
* If parallel computing is required, cluster information needs to be provided.  
For other parameter information, refer specifically to [eQTLMAPT](https://github.com/hitbc/eQTLMAPT)

## Running


## Results/Output


## Citation
Our paper was published in 论文地址。

## Contact
If you need help, please contact 邮箱地址。
