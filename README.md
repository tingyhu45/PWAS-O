# PWAS-O: Omnibus proteome-wide association study


For this Omnibus PWAS, we：
* use [FUSION software](http://gusevlab.org/projects/fusion/) to train protein imputation model with the most predictive model out of penalized linear regression model with Elastic-Net penalty and LASSO penalty, regular linear regression model with best unbiased linear predictor (BLUP), single variant model with Top pQTL (Top 1).
* use [TIGAR tool](https://github.com/yanglab-emory/TIGAR) to train protein imputation model with a nonparamtric Bayesian latent Dirchlet Process Regression model (DPR) as well as penalized linear regression model with Elastic-Net penalty (as implemented by PrediXcan).
* apply R pacakage [ACAT](https://github.com/yaowuliu/ACAT) to combine the PWAS p-values from different tools based on the Cauchy Association test.
