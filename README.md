# eQTLMAPT 
**eQTLMAPT**:fast and accurate eQTL mediation analysis with efficient permutation testing approaches.

## Introduction
`eQTLMAPT` implements fast and accurate eQTL mediation analysis with efficient permutation procedures to control for multiple testing. eQTLMAPT provides adaptive permutation scheme which prunes the permutation process opportunely, and it models the null distribution using generalized Pareto distribution (GPD) trained from a few permutations. The adjusted P-values can be estimated at any significance level based on the null distribution with much less computational cost compared with traditional fixed number permutation strategies. In addition, eQTLMAPT provides flexible options for users to choose proper permutation schemes combining with various confounding factors adjustment methods based on their practical needs. Experiments on real eQTL dataset demonstrate that eQTLMAPT provides higher resolution of estimated significance and is order of magnitude faster than compared methods with similar accuracy.

## Installation
    install.packages("devtools", dependencies = T)  
    devtools::install_github("QidiPeng/eQTLMAPT")

## Input
```bash
`snp.dat`          The genotype matrix.  Each row is an eQTL, each column is a sample.  
`fea.dat`          The gene expression profile matrix. Each row is a gene\'s expression profile, each column is a sample.  
`known.conf/conf`  The known confounder matrix which is adjusted in all mediation tests. Each row is a confounder, each column is a sample. 
`cov.pool`         The pool of candidate confounding variables from which potential confounders are adaptively selected to adjust for each mediation test. Each row is a covariate, each column is a sample.  Only need to sepcify in adaptive confounder selection mode, if not given in this mode, principal components will be calculated instead.
`trios.idx`        The trios matrix of 3 columns. Each row represents a trio (eQTL, cis-gene, trans-gene). The first element represents the index of the eQTL in `snp.dat`; The second element represents the index of cis-gene in `fea.dat`,  and the third element represents the index of the trans-gene in `fea.dat`.  
`cl`               If parallel computing is required, cluster information needs to be provided.  
```
For other parameter information, refer specifically to help function.  

## Output
```bash
`nperm`            The executed permutation times, for adaptive permutation scheme only.  
`nominal.p`        The nominal P-value by testing the significance of `beta2` in the regression formula `trans_gene ~ beta1 * SNP + beta2 * cis_gene + err`, using t-test.  
`empirical.p`      The permutation P-value by testing the significance of `beta2` using permutation test.  
`empirical.p.gpd`  The permutation P-value estimated by GPD approximation.  
`std.error`        The standard error of `beta2`.  
`t_stat`           The `t-statistics` in testing the significance of `beta2`.  
`beta`             The `beta1` in the regression formula `trans_gene ~ beta1 * SNP + beta2 * cis_gene + err`.  
`beta.total`       The `beta.total` in the regression formula `trans_gene ~ beta.total * SNP + err`.  
`beta.change`      The proportions mediated, calculated by `(beta.total-beta)/beta.total`.  
`pc.matrix`        The principal components (PCs) matrix of expression profiles. This will be returned if the PCs are used as the pool of potential confounders. Each column is a PC. Returned only in adaptive confounder selection mode.
`sel.conf.ind`     The indicator matrix indicating which confounders are selected during mediation analysis. Returned only in adaptive confounder selection mode. Dimension of `sel.conf.ind` is the number of trios by the number of covariates in `cov.pool` or `pc.matrix`.  
```

## Demo data
    ## generate a cluster with 4 nodes for parallel computing  
    cl = makeCluster(4)    
    ## Using the first `pc.num` = 10 PCs as candidate confunding variable pools.  
    ## The maximum number of permutation is `Maxperm` = 10000 in the adaptive permutation scheme. And when permutation number better than original statistics upon `Minperm` = 100 stop.  
    ## When the empirical P-value is less than `gpd.perm` = 0.01, a more accurate empirical P-value is estimated using the GPD fit.  
    output <- gmap.ac.gpd(snp.dat = dat$snp.dat, fea.dat = dat$fea.dat, known.conf = dat$known.conf, trios.idx = dat$trios.idx[1:10,], cl = cl, cov.pool = NULL, pc.num = 10, Minperm = 100, Maxperm = 10000, gpd.perm = 0.01)  
 The `dat` object used here is packaged in [example.rda](https://github.com/QidiPeng/eQTLMAPT/blob/master/data/example.rda).


## Contact
If you need help, please contact **ydwang@hit.edu.cn** or **jiajiepeng@nwpu.edu.cn**.
