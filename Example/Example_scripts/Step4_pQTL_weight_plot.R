
# load package
library(ggplot2)
library(ggpubr)


# Plot cis-weights vs. position: TIGAR and FUSION 
weight_pos_plot <- function(data_TIGAR,
                            data_fusion,
                            point_alpha = 0.5,
                            xlab="Position on chr11 (Mb)",
                            ylab="pQTL Weights for MRPL16",
                            GWAS_sig_level=1e-05
){
  p_TIGAR <- ggplot(
    data_TIGAR,
    aes(x=position,
        y=ES)) +
    geom_point(size=2.5,
               alpha=point_alpha) +
    
    # plot significant genes based on GWAS
    geom_point(
      data=subset(data_TIGAR, GWAS_p < GWAS_sig_level), # highlight significant GWAS signals
      aes(
        x=position, 
        y=ES,
        color=-log10(GWAS_p)), 
      size=2.5, 
      alpha=point_alpha)+ 
    scale_color_gradient(low="#F0E442",high="#E74C3C") +
    
    theme(axis.title.x=element_blank(),axis.title.y=element_blank())+  # hide axis labels
    ggtitle("TIGAR")+
    geom_hline(
      yintercept=0, 
      linetype='dashed', 
      size=0.3, 
      color="black") +
    labs(color="-log10(GWAS p-value)") +
    theme(plot.title=element_text(size=14,face="bold"),
          legend.title = element_text(color="black",size=10,face="bold")) 
  #+ guides(color=guide_legend(title=expression((paste("-log10","(GWAS p-value)")))))
  
  
  
  #   p_EN <- ggplot(
  #     data_EN,
  #     aes(x=position,
  #         y=ES)) +
  #       geom_point(size=2.5,
  #         alpha=point_alpha) +
  #      #plot significant genes based on GWAS
  #     geom_point(
  # 			data=subset(data_EN, GWAS_p < GWAS_sig_level),
  # 			aes(
  # 				x=position, 
  # 				y=ES,
  # 				color=GWAS_p), 
  # 	      size=2.5, 
  # 			alpha=point_alpha) +
  #     
  #     theme(axis.title.x=element_blank(),axis.title.y=element_blank()) +
  #     scale_y_continuous(n.breaks = 7)+
  #     ggtitle("PrediXcan")+
  #     geom_hline(
  # 			yintercept=0, 
  # 			linetype='dashed', 
  # 			size=0.3, 
  # 			color="black") +
  #     theme(plot.title=element_text(size=14,face="bold"),legend.position = "bottom")
  
  p_fusion <-  ggplot(
    data_fusion,
    aes(x=position,
        y=ES)) +
    geom_point(size=2.5,
               alpha=point_alpha) +
    geom_point(
      data=subset(data_fusion, GWAS_p < GWAS_sig_level),
      aes(
        x=position, 
        y=ES,
        color=GWAS_p), 
      size=2.5, 
      alpha=point_alpha) +
    scale_color_gradient(low="#F0E442",high="#E74C3C") +
    
    theme(axis.title.x=element_blank(),axis.title.y=element_blank()) +
    scale_y_continuous(n.breaks = 7)+
    ggtitle("FUSION")+
    geom_hline(
      yintercept=0, 
      linetype='dashed', 
      size=0.3, 
      color="black") +
    theme(plot.title=element_text(size=14,face="bold"),legend.position = "bottom")
  
  
  p_all <- ggarrange(p_TIGAR,p_fusion,ncol =2,common.legend = T, legend="bottom")
  
  
  p_all_annotate <- annotate_figure(p_all, 
                                    bottom = text_grob(xlab,face="bold",size=14),
                                    left= text_grob(ylab,face="bold",size=14,rot=90))
  
  return(p_all_annotate)
}

