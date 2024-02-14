
# load package
library(ggplot2)
library(ggrepel)
library(gg.gap)

# Manhattan plot
manhattan_plot1 <- function(
  ## data: 
  data, # dataframe (columns: 'CHR','POS','Pvalue','label_text')
  ## genes:
  point_colors=c('#44B5AD', '#3A948E', '#36807A', '#2f615d'), # color for non-sig genes
  sig_color1='#FD923F', # color for sig. genes without labels
  sig_color2='#D92B26', # color for sig. genes with labels
  point_alpha=0.9, # transparency for genes
  ## significance level:
  sig_level=0.05, # significance level value
  sig_level_line_col='black', # line color
  sig_linetype='dashed', # linetype
  sig_line_size=1, # line size
  ## plot theme, other plot variables:
  chr_vec=1:22, # chromosomes to plot
  chr_gap=80, # gap between chromosomes in x-axis
  theme=theme_grey(), # ggplot2 theme (can use custom theme)
  plot_bg_col=NULL, # background color if different from theme
  panel_border=element_blank(), # ggplot panel border (default: blank)				
  text_size=13, # text size
  ## point labelling:
  geom_label_size=2.5, # label text size
  label_fill='white', # label background color
  label_col='black', # label border color
  label_seg_col='black', # color of line from label to point
  min_segment_length=0.01, # minimum length of line from label to point
  segment_size=0.2, # line from label to point
  label_force=2, # force of repulsion between overlapping text labels
  point_padding=1e-06, # padding around genes,
  label.padding = 0.25,
  box.padding=0.3, #padding around boxes
  seed=NA, # optional seed for generating label positions
  max_iter=20000, # number of iterations to use to generate label,positions
  max.overlaps=5,
  title = "Manhattan Plot of the AD PWAS",
  plot_broken_axis = TRUE, # If p-value with large magnitude: break y-axis
  broken_segments = list(c(20,30),c(40,60)),
  broken_tick_width = c(5,10,10),
  broken_rel_heights = c(0.9,0,0.1,0,0.1),
  broken_y_lim = c(0,70)
){
  
  
  # setup dataframe for plotting
  plot_dat <- NULL # dataframe
  endPos <- 0  # place on x-axis where last chromosome ended
  x_axis_chr_breaks <- NULL # chromosome label positions
  x_axis_chr_labels <- NULL # list of chromosomes to label
  
  # get plot positions from chromosome positions
  for (chr in chr_vec) {
    # get data for chr
    temp <- data[data$CHR==chr, ]
    
    if (nrow(temp) > 0) {
      # append chromosome to list of chromosomes to label
      x_axis_chr_labels <- c(x_axis_chr_labels, chr)
      
      # get unique positions for this chr
      uniq_pos <- sort(unique(temp$POS))
      uniq_pos <- setNames(order(uniq_pos), uniq_pos)
      
      # set POS to order value of the unique positions
      temp$POS <- uniq_pos[as.character(temp$POS)]
      
      # get plot positions for genes on this chr
      temp$plotPos <- (temp$POS - min(temp$POS, na.rm=TRUE) ) + endPos + 1
      
      # get new end position based on max position
      endPos <- max(temp$plotPos, na.rm=TRUE) + chr_gap
      
      # append label position
      x_axis_chr_breaks <- c(x_axis_chr_breaks, mean(temp$plotPos, na.rm=TRUE) )
      
      # add rows to plot_dat
      plot_dat <- rbind(plot_dat, temp)
    }
  }
  
  # set min, max values for axes
  min_x <- min(plot_dat$plotPos)
  max_x <- max(plot_dat$plotPos)
  max_y <- max(-log10(plot_dat$Pvalue))
  
  # plot
  p <- ggplot(
    data=plot_dat, 
    aes(
      x=plotPos, 
      y=-log10(Pvalue), 
      label=label_text)) + 
    # plot non-sig. genes
    geom_point(
      data=subset(plot_dat, Pvalue >= sig_level),
      aes(
        x=plotPos, 
        y=-log10(Pvalue), 
        color=factor(CHR)),
      size=2, 
      alpha=point_alpha) + 
    # color non-sig. genes
    scale_color_manual(
      values=rep(point_colors, 22)) +
    # plot sig. genes
    geom_point(
      data=subset(plot_dat, Pvalue < sig_level),
      aes(
        x=plotPos, 
        y=-log10(Pvalue), 
        fill=factor(CHR)),
      size=ifelse(subset(plot_dat, Pvalue < sig_level)$label_text=='', 2.25, 2.5),
      color=ifelse(subset(plot_dat, Pvalue < sig_level)$label_text=='', sig_color1, sig_color2),
      alpha=point_alpha) +
    
    # add labels
    geom_label_repel(
      data=subset(plot_dat, Pvalue < sig_level ), 
      min.segment.length=min_segment_length,
      segment.size=segment_size,
      segment.color=label_seg_col,
      box.padding=box.padding,
      size=geom_label_size,
      alpha=1,
      ylim=c(-log10(sig_level), max_y+20),
      xlim = c(min_x, max_x),
      force =label_force,
      point.padding=point_padding,
      label.padding = label.padding,
      max.iter=max_iter,
      colour=label_col,
      fill = label_fill,
      force_pull=10,
      seed=seed,
      direction = c("both", "y", "x"),
      #direction = c("both"),
      #max.overlaps = getOption("ggrepel.max.overlaps", default = 10)) 
      max.overlaps = max.overlaps) +
    
    
    # significance level line
    geom_hline(
      yintercept=-log10(sig_level), 
      linetype=sig_linetype, 
      size=sig_line_size, 
      color=sig_level_line_col) +
    # remove legend
    guides(
      color='none', 
      fill='none') + 
    # set axes titles
    labs(
      x='Chromosome', 
      y=expression(bold(-"log"[10]("q-value")))) + 
    # x-axis labels, breaks
    scale_x_continuous(
      breaks=x_axis_chr_breaks, 
      labels=x_axis_chr_labels, 
      expand=c(0.01, 0)) + 
    # don't clip drawing to extent of plot panel
    coord_cartesian(
      clip='off') +
    # pad y-axis
    #scale_y_continuous(
     # expand=c(0.05, 0)) +
    # theme variable; (ie, theme_grey() or theme_bw())
    theme +
    # convenient theme options for a manhattan plot, 
    # 	sets text size
    # 	removes minor grids
    # 	sets panel border from variable
    # 	sets plot background from variable
    # this will overwrite /contradicting/ user-submitted values from 'theme'
    theme(
      text=element_text(size=text_size,face='bold'), 
      axis.title.y = element_text(face='bold',size = text_size),
      axis.text.x = element_text(
        face='bold', 
        size=text_size-1, 
        angle=-90, 
        vjust=0.5, 
        hjust=0),
      panel.grid.major.x=element_blank(),
      panel.grid.minor.x=element_blank(),
      panel.border=panel_border,
      plot.background=element_rect(fill=plot_bg_col))  
  #+ggtitle(title)
  
  if (plot_broken_axis==F){
    return(p)
  } else{
    p2 <- gg.gap(plot=p,segments = broken_segments,
                 tick_width = broken_tick_width,
                 rel_heights = broken_rel_heights,
                 ylim = broken_y_lim,
                 margin = c(top = 0, right = 0, bottom = 0, left = 0.8)) +
      ggtitle(title)
    return(p2)
  }
}






