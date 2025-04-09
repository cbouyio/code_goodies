# Source file with useful functions for preprocessing, differential, enrichment and cluster analysis of RNA-seq data.
# cbouyio, since 2019, UMR 7216 & UFR 8512, ERL 11933

# Libraries
library(mclust)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(diptest)


# Global functions
my_palette <- function(n, nc = "RdYlGn"){
  return(colorRampPalette(brewer.pal(11, name = nc))(n))
}

jaccard_similarity <- function(a, b) {
  # Return the Jaccard similarity between two sets A and B. Useful to compare clusterings
  inters <- length(intersect(a, b))
  uni <- length(union(a, b))
  return(inters/uni)
}

# Processing functions ----------------
discretise <- function(x, t) {
  # Low level function to return a 3 level discretisation
  if (x <= -t) return(-1)
  else if (x >= t) return(1)
  else return(0)
}

discretiseList <- function(v, thres, fc = 1, ...) {
  # Low level function which discretises a list, either by a 'hard' threshold, or by MAD
  med <- median(v)
  md <- mad(v)
  #med = mean(v)
  #md = sd(v)
  if (missing(thres)) {
    thres2 = med + fc*md
  } else {
    thres2 = thres
  }
  return(sapply(v, FUN = discretise, t = thres2, ...))
}

discretiseMatrix <- function(df, ...){
  # Function to discretise a data.frame
  ll <- apply(df, 2, FUN = discretiseList, ...)
  return(as.data.frame(ll))
}

export_cytoscape <- function(m, filename, ...) {
  # Function to export a correlation matrix 'm' to a cytoscape network file
  # m: correlation matrix
  fn <- filename
  write("from\tto\tedgeWeight\tsign", file = fn)
  for (i in 1:nrow(m)) {
    for (j in i:ncol(m)) {
      if (m[i,j] != 0.0) {
        write(sprintf("%s\t%s\t%f\t%s", colnames(m)[i], rownames(m)[j], abs(m[i,j]), sign(m[i,j])), file = fn, append = TRUE)
      }
    }
  }
}


# Filtering functions -----------------
remove_all_zeros <- function(df, ...){
  # Function to remove rows that contain nothing but zeros in a data frame
  return(df[rowSums(df[]) > 0,])
}

filter_low_counts <- function(gem, exps, g = 1, t = 5, ...){
  # Function to filter out lowly expressed (i.e. counts) genes
  # gem  : Gene Expression Matrix. A data frame with the conditions as columns and the gene expression profiles as rows.
  # exps : Experiments. A factor specifying the grouping of experiments. MUST have length the number of columns of gem and MUST specify the different levels (i.e. treatments) of the data.
  # g    : Groups. The number of groups which we ask to have more that the threshold reads.
  # t    : Counts threshold. A threshold of counts that all the replicates in the specified number of groups must be above.
  rows <- vector()
  # THIS is the vector of number of replicates per experiment.
  nR <- tabulate(exps)
  for (i in 1:nrow(gem)) {
    row <- as.numeric(gem[i,])
    # Calculate how many times you get more counts than the threshold in each experiment.
    agT <- aggregate(row~exps, FUN = function(v){return(sum(v >= t))})$row
    # This condition is checking for the counts to be higher than the threshold t, in at least g experiments.
    if ( sum(agT == nR) >= g ) {
      rows <- append(rows, rownames(gem)[i])
    }
  }
  return(gem[rows, ])
}

filter_noisy_counts <- function(gem, exps, c = 0.5, ...){
  # Function to filter out noisy expressed genes
  # gem  : Gene Expression Matrix. A data frame with the conditions as columns and the gene expression profiles as rows.
  # exps : Experiments. A factor specifying the grouping of experiments. MUST have equal length with the columns of gem and MUST specify the different levels (i.e. treatments) of the data.
  # c    : Coefficient of variation threshold.
  rows <- vector()
  for (i in 1:nrow(gem)) {
    row <- as.numeric(gem[i,])
    means <- aggregate(row~exps, FUN = mean)$row
    sds <- aggregate(row~exps, FUN = sd)$row
    cvs <- sds/means
    cvs[is.na(cvs)] <- Inf
    lg <- length(levels(exps))
    # Condition to filter all genes whose coefficient of variation is more than c in at least one experiment.
    if (sum(cvs <= c) == lg) {
      rows <- append(rows, rownames(gem)[i])
    }
  }
  return(gem[rows, ])
}

filter_low_expression <- function(e, groups, thres = 3, samples = 1, coefVar = 0.5, ...){
  # Function to filter lowly expressed genes
  # e       : raw counts data.frame (or cpm, or tpm matrix).
  # groups  : factor designating the grouping(s) (conditions, treatments etc.) it MUST be of equal length to the columns of e and its levels must represent the different groups (conditions).
  # thres   : the threshold of the *mean* expression whithin a group.
  # samples : the minimum number of groups (i.e. conditions) that we want the mean to be higher than the threshold 'thres'.
  # coefVar : The coefficient of variation threshold, within replicates, used to remove "noisily" measured genes.
  rows <- vector()
  for (i in 1:nrow(e)) {
    row <- as.numeric(e[i,])
    means <- aggregate(row~groups, FUN = mean)$row
    sds <- aggregate(row~groups, FUN = sd)$row
    cvs <- sds/means
    lg <- length(levels(groups))
    # This condition filters for "noisy" genes AND for lowly expressed genes.
    if ( (sum(cvs <= coefVar) == lg) & (sum(means >= thres) >= samples) ) { # The first part of this condition checks for the coefficient of variability in ALL groups by grouping sums in the number of samples and the second part is checking for the experimental threshold we specify.
      rows <- append(rows, rownames(e)[i])
    }
  }
  return(e[rows,])
}

filter_informative_genes <- function(e, grouping, test = "t", thres = 0.01, ...){
  # Function to filter informative genes between two conditions
  # e      : gene expression data frame (or cpm, or tpm matrix).
  # groups : factor designating the grouping (conditions, treatment etc.) it MUST be of equal length to the columns of e and it MUST have only two levels.
  # test   : One of "t", "d" or "w" for a parametric t-Student test, a Hartigan dip bimodality or a non-parametric Wilcoxon test (default t-test)
  # thres  : the threshold of the t-test p-value.
  rows <- vector()
  for (i in 1:nrow(e)) {
    dft <- data.frame(gexpr = as.numeric(e[i,]), cond = grouping)
    if (test == "t") {
      tt <- t.test(gexpr~cond, data = dft)
    } else if (test == "w") {
      tt <- wilcox.test(gexpr~cond, data = dft)
    } else if (test == "d") {
      tt <- dip.test(dft$gexpr)
    } else {
      print("Wrong test option, select only one of 'w', 'd' or 't' for Wilcoxon, DIP or t-Student test.")
      q()
    }
    # This condition filters for "informative" genes.
    if ( tt$p.value <= thres ) {  # Looks in the UNADJUSTED p-value TODO perhaps some adjustment for multiple testing.
      rows <- append(rows, rownames(e)[i])
    }
  }
  return(e[rows,])
}


# Microarray functions ------------------------------------
select_max_probset <- function(df, ...){
  # Select the one probeset for each gene. The probeset with the highest average value among all experiments
  # df MUST contain two ID columns one for probeset names and a second one with gene/transcripts or any other entity names
  # First we calculate the average intensities for all probesets
  avgAll <- rowSums(df[, 3:ncol(df)])
  df["AvgInt"] <- avgAll
  #Then use the awesome AVE function SUPERB solution!!!
  dfA <- df[as.logical(ave(df$AvgInt, df[,2], FUN = function(x) x == max(x))),]
  dfA["AvgInt"] <- NULL
  return(dfA)
}


## Ploting functions ------------------
plotMatrixBoxplot <- function(df, ...) {
  # Function that plots a violin - jitter plot of a numeric matrix
  dff <- df %>% rownames_to_column(var = "GeneID") %>% gather(Experiment, LogFC, -GeneID, factor_key = TRUE)
  p <- ggplot(dff, aes(x = Experiment, y = LogFC, color = Experiment)) + geom_violin(trim = FALSE) + geom_jitter(aes(alpha = 0.5), position = position_jitter(0.25)) + stat_summary(fun.data= mean_sdl, geom = "pointrange", color = "dimgrey", size = 1) + stat_summary(fun.y = median, geom = "point", shape = 95, size = 10, color = "black") + scale_color_brewer(palette = "Dark2") + theme_linedraw()
  return(p)
}

geneBoxplotCond <- function(matrix, name, experiments, treatments, jit_width = 0.1, point_size = 2, ...){
  # Function to plot expression of individual gene
  # Experiment : are the different fractions.
  # Treatment  : is High or Low glucose.
  ge <- data.frame(t(matrix[name,]));
  ge$exp <- experiments;
  ge$treat <- treatments;
  colnames(ge)[1] <- "TPM";
  p <- ggplot(ge, aes(exp, TPM));
  p + geom_jitter(aes(color = treat), width = jit_width, size = point_size) + ggtitle(name);
}

plot_semiSupervised_clust <- function(data, k, method, scale = FALSE, title = "", ...){
  # Nicely plots a k-means and other semi-supervised clusterings
  # Scaling.
  if (scale == TRUE) {
    data <- as.data.frame(scale(data))
  }
  # Calculate the clustering.
  clustMeth <- match.fun(method)
  clustRes <- clustMeth(data, k, ...)
  # Append id and cluster
  dfcall <- cbind(data, id = seq(nrow(data)), cluster = clustRes$cluster)
  # Add idsort, the id number ordered by cluster
  dfcall$idsort <- dfcall$id[order(dfcall$cluster)]
  dfcall$idsort <- order(dfcall$idsort)
  # Generate cluster colours.
  clusterCols <- as.character(sort(clustRes$cluster))
  # Title
  if (title == "") {
    ti = paste(method, " clustering of ", deparse(substitute(data)), ", scale:", as.character(scale))
  } else {
    ti = title
  }
  # Plotting
  heatmap(as.matrix(data)[order(clustRes$cluster),], Rowv = NA, Colv = NA, scale = "none", labRow = NA, cexCol = 0.75, col = my_palette(32), RowSideColors = clusterCols, ylab = "Genes", main = ti)
  invisible(list(res = clustRes, df = dfcall))
}


plot_unSupervised_clust <- function(data, method, scale = FALSE, title = TRUE, ...){
  # Nicely plots a k-means clustering or other unsupervised clustering
  # Scaling.
  if (scale == TRUE) {
    data <- as.data.frame(scale(data))
  }
  # Calculate the clustering.
  clustMeth <- match.fun(method)
  clustRes <- clustMeth(data, ...)
  # Append id and cluster
  dfcall <- cbind(data, id = seq(nrow(data)), cluster = clustRes$classification)
  # Add idsort, the id number ordered by cluster
  dfcall$idsort <- dfcall$id[order(dfcall$cluster)]
  dfcall$idsort <- order(dfcall$idsort)
  # Generate cluster colours.
  noClust <- max(clustRes$classification)
  #clusterCols <- as.character(sort(clustRes$classification))
  clusterCols <- brewer.pal(n = noClust, name = "Dark2")[as.factor(as.character(sort(clustRes$classification)))]
  # Title
  if (title == TRUE) {
    ti <- paste(method, " clustering of, ", deparse(substitute(data)), " scale: ", as.character(scale))
  } else {
    ti <- NULL
  }
  # Plotting
  heatmap(as.matrix(data)[order(clustRes$classification),], scale = "none", cexCol = 1, col = my_palette(32), RowSideColors = clusterCols, main = ti, cexRow = 1.5)
  invisible(list(res = clustRes, df = dfcall))
}


# Enrichment functions ----------------
comp_allEGO <- function(degsUP, degsDOWN, org, key, uni, pv = 0.01){
  degsALL <- c(degsUP, degsDOWN)
  ego1 <- enrichGO(gene = degsALL, OrgDb = org, keyType = key, ont = "MF", pvalueCutoff = pv, universe = uni, readable = FALSE, pool = TRUE)
  ego2 <- enrichGO(gene = degsUP, OrgDb = org, keyType = key, ont = "MF", pvalueCutoff = pv, universe = uni, readable = FALSE, pool = TRUE)
  ego3 <- enrichGO(gene = degsDOWN, OrgDb = org, keyType = key, ont = "MF", pvalueCutoff = pv, universe = uni, readable = FALSE, pool = TRUE)
  ego4 <- enrichGO(gene = degsALL, OrgDb = org, keyType = key, ont = "BP", pvalueCutoff = pv, universe = uni, readable = FALSE, pool = TRUE)
  ego5 <- enrichGO(gene = degsUP, OrgDb = org, keyType = key, ont = "BP", pvalueCutoff = pv, universe = uni, readable = FALSE, pool = TRUE)
  ego6 <- enrichGO(gene = degsDOWN, OrgDb = org, keyType = key, ont = "BP", pvalueCutoff = pv, universe = uni, readable = FALSE, pool = TRUE)
  namesEgos <- c("All_MF", "Up_MF", "Down_MF", "All_BP", "Up_BP", "Down_BP")
  allEgos <- list(ego1, ego2, ego3, ego4, ego5, ego6)
  names(allEgos) <- namesEgos
  return(allEgos)
}

plot_allEGO <- function(allEgo, cats = 20, ti = ""){
  return()
  #for eg in allEgo...
  #dotplot(eg, cats, )
}


# Low level function to add counts on a boxplot.
n_fun <- function(x){
  return(data.frame(y = max(x), label = length(x)))
}
