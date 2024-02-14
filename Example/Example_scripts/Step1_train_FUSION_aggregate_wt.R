#!/usr/bin/env Rscript

Sys.setlocale("LC_ALL", "C")
options(stringsAsFactors=F)

library(dplyr)
library(data.table)

### load in arguments
args=(commandArgs(TRUE))
print(args)
if(length(args)==0) {
  stop("Error: No arguments supplied!")
} else {
  out_dir = as.character(args[[1]])

}

# directories
fusion.out = file.path(out_dir,"Exp.wgt.RDat")


# select the best model: most significant among model with non-all-zero
if (file.exists(fusion.out)){
  print("start formatting fusion output")
  
  load(fusion.out)
  
  # Remove NAs (these should not be here) and models with all zero weights
  wgt.matrix[is.na(wgt.matrix)] = 0
  
  mod.nonzero <- which(!apply(wgt.matrix==0,2,all))
  wgt.matrix <- wgt.matrix[ ,mod.nonzero]
  cv.performance <- cv.performance[ ,mod.nonzero]
  
  # which rows have rsq
  row.rsq = grep( "rsq" , rownames(cv.performance) )
  # which rows have p-values
  row.pval = grep( "pval" , rownames(cv.performance) )
  
  # Identify the best model
  # get the most significant model
  mod.best = which.min(apply(cv.performance[row.pval,,drop=F],2,min,na.rm=T))
  
  print(cv.performance[row.pval,])
  
  if ( length(mod.best) == 0 ) {
    cat( "WARNING : " , unlist(wgtlist[w,]) , " did not have a predictive model ... skipping entirely\n" )
    FAIL.ctr = FAIL.ctr + 1
    next
  }
  
    
  # Format snps file 
  colnames(snps) = c("CHROM", "SNP", "BP", "POS", "ALT", "REF")
  
  # Create a data frame including SNPs and estimated pQTL weights (non-zero)
  wgt.df = data.frame(POS = snps$POS, ES = wgt.matrix[, mod.best],
                      TargetID = geneID)
  
  
  # Merge snps with wgt.df to get estimated weights for each SNP
  df_out = merge(snps, wgt.df,by="POS")
  df_out = df_out %>% arrange(as.numeric(POS))
  
  # create output files
  out_cols = c("CHROM", "POS","snpID", "REF","ALT", "TargetID", "ES")
  # MyDf = data.frame(matrix(nrow = 0, ncol = length(out_cols))) 
  # colnames(MyDf) = out_cols
  
  df_out$snpID = paste(df_out$CHROM,
                       df_out$POS,
                       df_out$REF,
                       df_out$ALT,sep=":")
  
  df_out = df_out[, out_cols]
  
  out_file = file.path(out_dir,"Exp_fusion_weight.txt")
  
  # Save the weights
  write.table(df_out, 
              file = out_file, 
              row.names = F,
              col.names = T,
              append = F,
              quote = F,
              sep = "\t")
  
}


