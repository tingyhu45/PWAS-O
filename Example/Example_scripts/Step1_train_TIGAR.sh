#!/usr/bin/bash


# activate the environment
conda activate tigarenv  
export PYTHONPATH=${CONDA_PREFIX}/lib/python3.5/site-packages/:$PYTHONPATH

# set up number of threads
nthread=1

### arguments for all scripts
TIGAR_dir=${path-to-TIGAR}

# output directory
out_dir=${path-for-training-output}


### DPR Training

## input:
cd ${PWAS-O_dir}/Example


# file of sampleID in protein and genotype file to use
train_sampleID=Exp_train_sample_ID.txt

# genotype file
genofile=Exp_geno.vcf.gz

# protein abundance file
protein=Exp_Protein_Abundance.txt

# output prefix
out_prefix_dpr=Exp_DPR_train
out_prefix_elas=Exp_Elastic_train

# chr
chr=$chromsome

############################################
#### TRAINING PROTEIN ABUNDANCE PREDICTION MODELS

# DPR model
${TIGAR_dir}/TIGAR_Model_Train.sh \
--model DPR \
--gene_exp ${protein} \
--train_sampleID ${train_sampleID} \
--genofile ${genofile} \
--genofile_type vcf \
--format GT \
--maf 0.01 \
--chr ${chr} \
--hwe 0.00001 \
--cvR2 0 \
--cvR2_threshold 0.005 \
--dpr 1 \
--ES fixed \
--thread ${nthread} \
--out_prefix ${out_prefix_dpr} \
--TIGAR_dir ${TIGAR_dir} \
--sub_dir 0 \
--out_dir ${out_dir}


# Elastic-net model
${TIGAR_dir}/TIGAR_Model_Train.sh \
--model elastic_net \
--gene_exp ${protein} \
--train_sampleID ${train_sampleID} \
--genofile ${genofile} \
--chr ${chr} \
--genofile_type vcf \
--format GT \
--maf 0.01 \
--hwe 0.0001 \
--cvR2 1 \
--cvR2_threshold 0.005 \
--alpha 0.5 \
--thread ${nthread} \
--out_prefix ${out_prefix_elas} \
--job_suf ${job_suf} \
--TIGAR_dir ${TIGAR_dir} \
--sub_dir 0 \
--out_dir ${out_dir}




