
#' ---
#' title: "Example of statistical assessment of empirical results"
#' author: "Christian Blum <christian.blum@iiia.csic.es>
#' date: "Today"
#' output: pdf_document
#' ---
#' 

#' This document contains the code used for the statistical assessment of empirical results
#' of algorithm comparision
#' 
#' ## Data loading and variable initializations
#' 
#' Load the package and create directories for the results of the analysis

library(scmamp)
#library(reshape2)
#library(ggplot2)

data.file <- "results_table.dat"
sep <- "\t"
plot.dir<-"./"

#' First of all, we will load the results.
#' 

data <- read.csv(data.file, sep=sep)

#' Now we create a summarization table with the average results per instance
#' 

average.function <- mean
instance.descriptors <- c("length", "sigma", "t")
to.ignore <- c("inst")

summary.data <- summarizeData(data=data, fun=average.function, 
                              group.by=instance.descriptors, ignore=to.ignore
)

#' ## Statistical assessment
#' 
#' First, verify that there is at least one algorithm that is different from the rest.
#' For that we can use Friedman's test
#' 
alg.columns <- c("BRKGA", "BRKGALlama3InstructFree","BRKGAGemini2.5Flash","BRKGAGPT4.1Mini","BRKGALlama4Maverick")
friedmanTest(data[, alg.columns])


#' Now, two-fold post-hoc comparison. For the table, the best in each instance vs. the rest
#' 

all.vs.best.results <- postHocTest(data=data, algorithms=alg.columns, 
                                   group.by=instance.descriptors, test="friedman",
                                   control="max", use.rank=FALSE, sum.fun=average.function,
                                   correct="finner", alpha=0.05)

#' Now, all vs. all for the plots.
#' 

# Set decreasing=TRUE when the best results are the smallest.
all.vs.all <- postHocTest(data=data, algorithms=alg.columns, test="friedman", 
                          control=NULL, use.rank=TRUE, sum.fun=average.function,
                          correct="finner", alpha=0.05, decreasing=TRUE)

#' Same but grouped 
#' 

#all.vs.all.grouped <- postHocTest(data, algorithms=alg.columns, test="friedman",
#                                  group.by=c("nstr"), control=NULL, use.rank=TRUE,
#                                  sum.fun=average.function, correct="finner", alpha=0.05, decreasing=TRUE)

#' ## Plots
#' 
#' Create a plot similar to Demsar's critical difference plot.
#' 
#' Non-grouped plot
#' 

pdf(file=paste0(plot.dir, "CD_plot_all2.pdf"), width=8, height=2.2)
plotRanking(all.vs.all$corrected.pval, summary=all.vs.all$summary, alpha=0.05)
dev.off()

#' Grouped plots for one column in group.by
#' 

#for (i in 1:dim(all.vs.all.grouped$corrected.pval)[3]){
#  pdf(file=paste0(plot.dir, "CD_plot_nstr_",
#                  gsub("/","-",all.vs.all.grouped$summary[i,1]),
#                  ".pdf"), width=6, height=2.2)
#  plotRanking(all.vs.all.grouped$corrected.pval[, , i], summary=all.vs.all.grouped$summary[i,-1], alpha=0.05)
#  dev.off()
#}


#' ## References
#' Demšar, J. (2006) Statistical Comparisons of Classifiers over Multiple Data Sets. _Journal of Machine Learning Research_, 7, 1-30. 

