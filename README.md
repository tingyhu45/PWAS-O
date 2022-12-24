# PWAS-O: Omnibus proteome-wide association study

### An omnibus proteome-wide association study (PWAS) integrating multiple statistical methods

![image](https://github.com/tingyhu45/PWAS-O/blob/main/PWAS-O_Framework.png)

For this Omnibus PWAS, we：
* use [FUSION software](http://gusevlab.org/projects/fusion/) to train protein imputation model and estimate pQTL weights with the most predictive model out of penalized linear regression model with Elastic-Net penalty and LASSO penalty, regular linear regression model with best unbiased linear predictor (BLUP), single variant model with Top pQTL (Top 1).
* use [TIGAR tool](https://github.com/yanglab-emory/TIGAR) to train protein imputation model and estimate pQTL weights with a nonparamtric Bayesian latent Dirchlet Process Regression model (DPR) as well as penalized linear regression model with Elastic-Net penalty (as implemented by PrediXcan).
* use TIGAR to conduct the summary-level association test integrating pQTL weights and GWAS summary statistics.
* apply R pacakage [ACAT](https://github.com/yaowuliu/ACAT) to combine the PWAS p-values from different tools based on the Cauchy Association test.


# Getting started
 ### 1. Download [TIGAR](https://github.com/yanglab-emory/TIGAR) and complete its software setup
 
* [BGZIP](http://www.htslib.org/doc/bgzip.html)
* [TABIX](http://www.htslib.org/doc/tabix.html) 
* Python 3.5 modules/libraries: pandas, numpy, scipy, sklearn, statsmodels
 
 
