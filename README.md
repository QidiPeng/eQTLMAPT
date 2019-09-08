# eQTLMAPT 
**eQTLMAPT** - eQTL Mediaiton Analysis with accelerated Permutation Testing approaches

## Introduction
eQTLMAPT uses R language to achieve fast and accurate eQTL mediation analysis with efficient permutation procedures to control for multiple testing. eQTLMAPT provides adaptive permutation scheme which prunes the permutation process opportunely, and it models the null distribution using generalized Pareto distribution (GPD) trained from a few permutations. The adjusted P-values can be estimated at any significance level based on the null distribution with much less computational cost compared with traditional fixed number permutation strategies. In addition, eQTLMAPT provides flexible options for users to choose proper permutation schemes combining with various confounding factors adjustment methods based on their practical needs. Experiments on real eQTL dataset demonstrate that eQTLMAPT provides higher resolution of estimated significance and is order of magnitude faster than compared methods with similar accuracy.

## Dependences
1. parallel  
2. readr

R version 3.2.3

## Building xQTLImp
'<install.packages("devtools", dependencies = T)>'  
'devtools::install_github("QidiPeng/eQTLMAPT")'

## 输入


## 运行



## 输出


## Citation
Our paper was published in 论文地址。

## Contact
If you need help, please contact 邮箱地址。
