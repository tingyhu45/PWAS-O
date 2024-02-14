#!/usr/bin/bash



############
### ACTIVATE TIGAR PYTHON ENVIRONMENT
module load Anaconda3/4.2.0
source activate tigarenv  
export PYTHONPATH=${CONDA_PREFIX}/lib/python3.5/site-packages/:$PYTHONPATH

############
### LOAD TABIX
## load modules
module load tabix/0.2.6

# TIGAR directory
TIGAR_dir=${path-to-TIGAR}

############ PWAS using pQTL weights from DPR, EN model and FUSION
## input files
cd ${PWAS-O_dir}/Example

# pQTL weight file
dpr_weight=Exp_DPR_train_pQTLweights.txt.gz
dpr_gene_info=Exp_DPR_train_GeneInfo.txt


EN_weight=Exp_Elastic_train_pQTLweights.txt.gz
EN_gene_info=Exp_Elastic_train_GeneInfo.txt


FUSION_weight=Exp_fusion_weight.txt.gz
FUSION_gene_info=Exp_fusion_train_GeneInfo.txt


# GWAS summary level statistics (from publicly available GWAS studies)
Zscore=Exp_GWAS_Zscore.txt

############
## Require input LD covariance file
LD=$LD_covariance_file


########
## OUTPUT
dpr_out_dir=${PWAS-O_dir}/Example/DPR_PWAS_result
EN_out_dir=${PWAS-O_dir}/Example/EN_PWAS_result
FUSION_out_dir=${PWAS-O_dir}/Example/FUSION_PWAS_result

# chr
chr=$chromsome


############
### RUN PWAS for DPR
${TIGAR_dir}/TIGAR_TWAS.sh \
--asso 2 \
--gene_anno ${dpr_gene_info} \
--Zscore ${Zscore} \
--weight ${dpr_weight} \
--LD ${LD} \
--chr ${chr} \
--thread ${NSLOTS:-1} \
--TIGAR_dir ${TIGAR_dir} \
--sub_dir 0 \
--out_dir ${dpr_out_dir}



############
### RUN PWAS for EN
${TIGAR_dir}/TIGAR_TWAS.sh \
--asso 2 \
--gene_anno ${EN_gene_info} \
--Zscore ${Zscore} \
--weight ${EN_weight} \
--LD ${LD} \
--chr ${chr} \
--thread ${NSLOTS:-1} \
--TIGAR_dir ${TIGAR_dir} \
--sub_dir 0 \
--out_dir ${EN_out_dir}


############
### RUN PWAS for FUSION
${TIGAR_dir}/TIGAR_TWAS.sh \
--asso 2 \
--gene_anno ${FUSION_gene_info} \
--Zscore ${Zscore} \
--weight ${FUSION_weight} \
--LD ${LD} \
--chr ${chr} \
--thread ${NSLOTS:-1} \
--TIGAR_dir ${TIGAR_dir} \
--sub_dir 0 \
--out_dir ${FUSION_out_dir}


############
## Unload modules
module unload tabix/0.2.6
conda deactivate






