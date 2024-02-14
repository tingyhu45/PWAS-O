# PWAS-O: Omnibus proteome-wide association study

### An omnibus proteome-wide association study (PWAS) integrating multiple statistical methods

![image](https://github.com/tingyhu45/PWAS-O/blob/main/PWAS-O_Framework.png)

For this Omnibus PWAS, we：
* use [FUSION software](http://gusevlab.org/projects/fusion/) to train protein imputation model and estimate pQTL weights with the most predictive model out of penalized linear regression model with Elastic-Net penalty and LASSO penalty, regular linear regression model with best unbiased linear predictor (BLUP), single variant model with Top pQTL (Top 1).
* use [TIGAR tool](https://github.com/yanglab-emory/TIGAR) to train protein imputation model and estimate pQTL weights with a nonparamtric Bayesian latent Dirchlet Process Regression model (DPR) as well as penalized linear regression model with Elastic-Net penalty (as implemented by PrediXcan).
* use TIGAR to conduct the summary-level association test integrating pQTL weights and GWAS summary statistics.
* apply R package [ACAT](https://github.com/yaowuliu/ACAT) to combine the PWAS p-values from different tools based on the Cauchy Association test.

### Reference
[Hu T et al. Omnibus Proteome-Wide Association Study (PWAS-O) Identified 43 Risk Genes for Alzheimer’s Disease Dementia. medrxiv; 2022.](https://www.medrxiv.org/content/10.1101/2022.12.25.22283936v2)


# Getting started
 ### 1. Download [TIGAR](https://github.com/yanglab-emory/TIGAR) and complete its software setup
 
* [BGZIP](http://www.htslib.org/doc/bgzip.html)
* [TABIX](http://www.htslib.org/doc/tabix.html) 
* Python 3.5 modules/libraries: pandas, numpy, scipy, sklearn, statsmodels
	* using conda
	```
	# create the environment tigarenv
	conda create --name tigarenv python=3.5 pandas numpy scipy scikit-learn statsmodels
	# deactivate the conda environment
	conda deactivate
	# activate the environment
	conda activate tigarenv
	# set the PYTHONPATH
	export PYTHONPATH=${CONDA_PREFIX}/lib/python3.5/site-packages/:$PYTHONPATH
	```
	* using pip
	```
	# install pip
	# install virtualenv
	python3 -m pip install --user virtualenv
	# cd to preferred install_directory
	cd ${install_dir}
	# create the virtual environment tigarenv in the current directory
	python3 -m virtualenv tigarenv --python=python3.5
	# activate the environment
	source ${install_dir}/tigarenv/bin/activate
	# install the packages
	python3 -m pip install numpy==1.15.2 pandas==0.23.4 scikit-learn==0.20.0 scipy==1.1.0 statsmodels==0.9.0
	# deactivate the environment
	deactivate
 	```
 
### 2. Download [FUSION tool](http://gusevlab.org/projects/fusion/) and complete its software setup

* FUSION:
```
wget https://github.com/gusevlab/fusion_twas/archive/master.zip
unzip master.zip
```
* add the bundled GCTA binary `gcta_nr_robust` to path
* download and install [GEMMA software](https://github.com/genetics-statistics/GEMMA/releases) 

* Launch R and install required libraries
```
install.packages(c('optparse','RColorBrewer'))
install.packages('plink2R-master/plink2R/',repos=NULL)
install.packages(c('glmnet','methods'))
```

# Example usage

## Step 1. Train protein imputation models 

### TIGAR tool for DPR and Elastic-Net model

* (1) Genotype data: vcf or dosage file
* (2) Training Sample ID file: headerless, single-column file containing sampleIDs to use
* (3) protein abundance file:
*  First 5 columns are _Chromosome number, Gene start position, Gene end position, Target gene ID, Gene name (optional, could be the same as Target gene ID)_

#### train nonparametric Bayesian DPR model
```
# Setup input file paths
Protein_Exp_train_file="${TIGAR_dir}/ExampleData/protein_exp.txt"
train_sample_ID_file="${TIGAR_dir}/ExampleData/sampleID.txt"
genofile="${TIGAR_dir}/ExampleData/example.vcf.gz"
out_dir="${TIGAR_dir}/ExampleData/output"

# Call TIGAR model training shell script
${TIGAR_dir}/TIGAR_Model_Train.sh \
--model DPR \
--gene_exp ${Gene_Exp_train_file} \
--train_sampleID ${train_sample_ID_file} \
--chr 1 \
--genofile ${genofile} \
--genofile_type vcf \
--format GT \
--maf 0.01 \
--hwe 0.0001 \
--cvR2 1 \
--dpr 1 \
--ES fixed \
--thread 2 \
--out_dir ${out_dir} \
--TIGAR_dir ${TIGAR_dir}
```

#### train Elastic-Net model
```
Protein_Exp_train_file="${TIGAR_dir}/ExampleData/protein_exp.txt"
train_sample_ID_file="${TIGAR_dir}/ExampleData/sampleID.txt"
genofile="${TIGAR_dir}/ExampleData/example.vcf.gz"
out_dir="${TIGAR_dir}/ExampleData/output"

${TIGAR_dir}/TIGAR_Model_Train.sh \
--model elastic_net \
--gene_exp ${Gene_Exp_train_file} \
--train_sampleID ${train_sample_ID_file} \
--chr 1 \
--genofile ${genofile} \
--genofile_type vcf \
--format GT \
--maf 0.01 \
--hwe 0.0001 \
--cvR2 1 \
--alpha 0.5 \
--thread 2 \
--out_dir ${out_dir} \
--TIGAR_dir ${TIGAR_dir}
```


### FUSION tool
FUSION computes pQTL weights based on binary file per gene
```
Rscript ${FUSION_dir}/FUSION.compute_weights.R \
	--bfile $b_dir \
	--tmp $tmp_dir \
	--out $out_dir \
	--PATH_plink ${pato_to_plink}/plink \
	--PATH_gcta ${pato_to_gcta}gcta_nr_robust \
	--PATH_gemma ${pato_to_gemma}/gemma-linux \
	--models top1,lasso,enet,blup
 ```



## Step 2. Conduct summary-level PWAS using TIGAR 

#### We integrate pQTL weights from protein imputation models with GWAS summary statistics to conduct the gene-based association test

* `--asso`: `2` summary-level TWAS using GWAS summary Z-score statistics and reference LD
* `--weight`: Path to SNP weight (eQTL effect size) file
* `--Zscore`: Path to GWAS summary Zscore statistics
* `--LD`: Path to reference LD (SNP genotype covariance matrix) that should be bgzipped and tabixed (generated by using `${TIGAR_dir}/TWAS/Get_LD.py` script with individual-level genotype data from a reference sample)
* `--window`: Window size (in base pairs) around gene region from which to include SNPs (default: `1000000` [+- 1MB region around gene region])
* `--test_stat`: burden Z test statistic to calculate: `FUSION`, `SPrediXcan`, or `both` (default both)

```
gene_anno_file="${TIGAR_dir}/ExampleData/gene_anno.txt"
pQTL_ES_file="${TIGAR_dir}/ExampleData/eQTLweights.txt.gz"
Zscore_file="${TIGAR_dir}/ExampleData/CHR1_GWAS_Zscore.txt.gz"
LD_file="${TIGAR_dir}/ExampleData/CHR1_reference_cov.txt.gz"

${TIGAR_dir}/TIGAR_TWAS.sh \
--asso 2 \
--gene_anno ${gene_anno_file} \
--chr 1 \
--weight ${pQTL} \
--Zscore ${Zscore_file} \
--LD ${LD_file} \
--thread 2 \
--out_dir ${out_dir} \
--TIGAR_dir ${TIGAR_dir}
```


Here reference LD genotype covariance file can be generated using reference genotype data, which should be generated per chromosome with the corresponding genome segmentation block annotation file.
LD files for PWAS should be made with the same reference data used for model training.

LD files for TWAS should be made with the same reference data used for model training.
* `--genom_block`: Path to genome segmentation block annotation based on LD
* `--sampleID`: Path to a file with sampleIDs to be used for generating reference LD files.
* `--genofile`: Path to the reference genotype file (bgzipped and tabixed). It is recommended that genotype files with data for multiple chromosomes be split into per-chromosome files.
```
block_annotation="${TIGAR_dir}/ExampleData/example_genome_block_CHR1.txt"
sample_id="${TIGAR_dir}/sampleID.txt"

${TIGAR_dir}/TIGAR_LD.sh \
--genome_block ${block_annotation} \
--sampleID ${sample_id}
--chr 1 \
--genofile ${genofile} \
--genofile_type vcf \
--format GT \
--maf 0.01 \
--thread 2 \
--out_dir ${out_dir} \
--TIGAR_dir ${TIGAR_dir}
```

## Step 3. Combine p-values from three methods using ACAT-O test
Install the ACAT package in R

```
library(devtools)
devtools::install_github("yaowuliu/ACAT")
```

Launch R and load the package
```
library(ACAT)

ACAT_withNA = function(p_vec){
  p_vec_noNA = p_vec[is.na(p_vec) == F]
  ACAT(p_vec_noNA)
  
}
```

## Step 4. Analyze the PWAS-O results
* Generate the Manhattan plot
* Q-Q plot
* pQTL weight plots by different methods
* Use [GIFT](https://yuanzhongshang.github.io/GIFT/) tool to further conduct PWAS risk gene fine mapping

#### The example scripts and sample data can be found at https://github.com/tingyhu45/PWAS-O/tree/main/Example


# Data availability

refer to this [page](https://doi.org/10.7303/syn53498818)

(1) Reference LD covariance files 

(2) pQTL weight files from protein imputation models, and the PWAS results

