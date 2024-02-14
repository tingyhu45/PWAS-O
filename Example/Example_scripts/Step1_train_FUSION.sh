#!/usr/bin/bash


# Load modules 
module load R/4.2.2


# genotype directory: binary files
# protein abundance should be written to geno.fam file
cd ${PWAS-O_dir}/Example


#### train FUSION and format output file

# 1. run FUSION

# required tools
plink_dir=$path_to_plink
gcta_dir=$path_to_gcta
gemma_dir=$path_to_gemma

# input file
b_dir=Exp.geno
tmp_dir=Exp.tmp
out_file=Exp_fusion_out


echo "start fusion training"

Rscript ${path-to-fusion}/fusion_twas-master/FUSION.compute_weights.R \
--bfile $b_dir \
--tmp $tmp_dir \
--out $out_file \
--hsq_p 1.5 \
--hsq_set 0.01\
--PATH_plink ${plink_dir} \
--PATH_gcta ${gcta_dir} \
--PATH_gemma ${gemma_dir} \
--models top1,lasso,enet,blup,bslmm


# 2. format FUSION pQTL weights for Stage II association test
out_dir=${PWAS-O_dir}/Example

Rscript Example_scripts/Step1_train_FUSION_aggregate_wt.R ${out_dir}






