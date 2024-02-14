# plot Q-Q plot


library(ggplot2)


############

myqq <- function(pvector, title="Quantile-Quantile Plot of -log10(P-values)", size = 24,point_size=3) {
  ci <- 0.95
  pvector = pvector[!is.na(pvector)]
  pvector = pvector[!is.nan(pvector)]
  n <- length(pvector)
  plotdata <- data.frame(
    observed = -log10(sort(pvector)),
    expected = -log10(1:n/n),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = seq(n), shape2 = rev(seq(n)))),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = seq(n), shape2 = rev(seq(n))))
  )

  plotdata <- plotdata[!is.infinite(plotdata$observed),]
  
  ggplot(plotdata, aes(x = expected, y = observed)) +
    geom_point(size=point_size) +
    geom_ribbon(aes(ymax = cupper, ymin = clower), fill = "grey30", alpha = 0.5) +
    xlim(c(0, max(plotdata$expected))) +
    #  ylim(c(0, max_p)) +
    labs(x = "Expected", y = "Observed", title = title) +
    geom_abline(intercept = 0, slope = 1, colour = "red") +
    theme(text = element_text(size = size, face = "bold"))
}

# apply Q-Q plot function
qq_ACAT <- myqq(P_ACAT,size=18,point_size=2.6) # P_ACAT: a vector of ACAT p-values