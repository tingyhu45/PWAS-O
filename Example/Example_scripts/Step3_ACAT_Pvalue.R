
# load packages
library(data.table)
library(readr)
library(dplyr)
library(ACAT)

# read in PWAS output files: reaplce with output directories
setwd("Example/")
PWAS_DPR <- fread("Exp_DPR_PWAS_result/Exp_PWAS_result.txt")
PWAS_EN <- fread("Exp_EN_PWAS_result/Exp_PWAS_result.txt")
PWAS_fusion <- fread("Exp_FUSION_PWAS_result/Exp_PWAS_result.txt") 

# combine 3 sets of PWAS outputs:  DPR, EN and FUSION
PWAS_acat <- PWAS_DPR %>% full_join(PWAS_EN,by="GeneName") %>% 
  full_join(PWAS_fusion,by="GeneName") 

# ACAT function
ACAT_withNA = function(p_vec){
  p_vec_noNA = p_vec[is.na(p_vec) == F]
  p_vec_noZero = p_vec_noNA[which(p_vec_noNA != 0)]
  ACAT(p_vec_noZero)
}

# apply ACAT function to 3 pvalues
pvalue_index = which(names(PWAS_acat) %in% c("SPred_PVAL.x","SPred_PVAL.y","SPred_PVAL"))
PWAS_acat$P_ACAT <- apply(PWAS_acat[,..pvalue_index],1,ACAT_withNA)

# adjust the ACAT p-value for FDR
PWAS_acat$P_ACAT_fdr <- p.adjust(PWAS_acat$P_ACAT, method = 'fdr', n=length(PWAS_acat$P_ACAT))

# annotate gene info: chromosome, position, and Ensemble ID
chrom_index = which(names(PWAS_acat) %in% c("CHROM.x","CHROM.y","CHROM"))
PWAS_acat$CHROM <- apply(PWAS_acat[,..chrom_index], 1, function(i) ifelse(all(is.na(i)), NA, i[!is.na(i)]))

start_index = which(names(PWAS_acat) %in% c("GeneStart.x","GeneStart.y","GeneStart"))
PWAS_acat$GeneStart <- apply(PWAS_acat[,..start_index], 1, function(i) ifelse(all(is.na(i)), NA, i[!is.na(i)]))

end_index = which(names(PWAS_acat) %in% c("GeneEnd.x","GeneEnd.y","GeneEnd"))
PWAS_acat$GeneEnd <- apply(PWAS_acat[,..end_index], 1, function(i) ifelse(all(is.na(i)), NA, i[!is.na(i)]))

TargetID_index = which(names(PWAS_acat) %in% c("TargetID.x","TargetID.y","TargetID"))
PWAS_acat$TargetID <- apply(PWAS_acat[,..TargetID_index], 1, function(i) ifelse(all(is.na(i)), NA, i[!is.na(i)]))

# select and rename columns 
PWAS_acat = PWAS_acat %>% select(CHROM,GeneName,TargetID,GeneStart,GeneEnd,
                                 SPred_Z.x,SPred_PVAL.x,SPred_Z.y,SPred_PVAL.y,SPred_Z,SPred_PVAL,
                                 P_ACAT,P_ACAT_fdr)

names(PWAS_acat) = c("Chrom","GeneName","EnsembleID","GeneStart","GeneEnd",
                     "Zscore_DPR","Pvalue_DPR","Zscore_EN","Pvalue_EN","Zscore_FUSION","Pvalue_FUSION",
                     "Pvalue_ACAT","Pvalue_ACAT_FDR")
