# Source file with useful functions for preprocessing, differential, enrichment and cluster analysis of RNA-seq data.
# cbouyio, since 2019, UMR 7216 & UFR 8512, ERL 11933

# Libraries
library(mclust)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(diptest)
library(enrichplot)
library(clusterProfiler)


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

compute_allEGOs <- function(degsUP, degsDOWN, org, key, uni, pv = 0.05, ...) {
  # Computes all GO enrichments for list of genes of interest. UP and DOWN and union.
  # Returns a list of 6 enrichResult instances.
  degsALL <- c(degsUP, degsDOWN)
  ego1 <- enrichGO(gene = degsALL, OrgDb = org, keyType = key, ont = "MF", pvalueCutoff = pv, universe = uni, readable = FALSE, pool = TRUE, ...)
  ego2 <- enrichGO(gene = degsUP, OrgDb = org, keyType = key, ont = "MF", pvalueCutoff = pv, universe = uni, readable = FALSE, pool = TRUE, ...)
  ego3 <- enrichGO(gene = degsDOWN, OrgDb = org, keyType = key, ont = "MF", pvalueCutoff = pv, universe = uni, readable = FALSE, pool = TRUE, ...)
  ego4 <- enrichGO(gene = degsALL, OrgDb = org, keyType = key, ont = "BP", pvalueCutoff = pv, universe = uni, readable = FALSE, pool = TRUE, ...)
  ego5 <- enrichGO(gene = degsUP, OrgDb = org, keyType = key, ont = "BP", pvalueCutoff = pv, universe = uni, readable = FALSE, pool = TRUE, ...)
  ego6 <- enrichGO(gene = degsDOWN, OrgDb = org, keyType = key, ont = "BP", pvalueCutoff = pv, universe = uni, readable = FALSE, pool = TRUE, ...)
  namesEgos <- c("All_MF", "Up_MF", "Down_MF", "All_BP", "Up_BP", "Down_BP")
  allEgos <- list(ego1, ego2, ego3, ego4, ego5, ego6)
  names(allEgos) <- namesEgos
  return(allEgos)
}


compute_allGSEGOs <- function(degsLFCList, org, key, pv = 0.5, ...) {
  # Computes BP and MF GO GSEA enrichments for a LogFC sorted list of genes of interest.
  # Returns a list of 2 enrichResult instances.
  gse1 <- gseGO(geneList = degsLFCList, OrgDb = org, keyType = key, ont = "MF", pvalueCutoff = pv, ...)
  gse2 <- gseGO(gene = degsLFCList, OrgDb = org, keyType = key, ont = "BP", pvalueCutoff = pv, ...)
  return(list("MF_GSEA" =  gse1, "BP_GSEA" = gse2))
}


compute_allKEGGs <- function(degsUP, degsDOWN, org, uni, pv = 0.05, ...){
  # Computes all KEGG enrichments for list of genes of interest. UP and DOWN and union. Org must be a KEGGdb name.
  # Returns a list of 3 enrichResult instances.
  degsALL <- c(degsUP, degsDOWN)
  ekeg1 <- enrichKEGG(gene = degsALL, organism = org, pvalueCutoff = pv, universe = uni, ...)
  ekeg2 <- enrichKEGG(gene = degsUP, organism = org, pvalueCutoff = pv, universe = uni, ...)
  ekeg3 <- enrichKEGG(gene = degsDOWN, organism = org, pvalueCutoff = pv, universe = uni, ...)
  return(list("KEGG_ALL" = ekeg1, "KEGG_UP" = ekeg2, "KEGG_DOWN" = ekeg3))
}


compute_allMKEGGs <- function(degsUP, degsDOWN, org, uni, pv = 0.05, ...){
  # Computes all KESS enrichments for list of genes of interest. UP and DOWN and union. Org must be a KEGGdb name.
  # Returns a list of 3 enrichResult instances.
  degsALL <- c(degsUP, degsDOWN)
  ekeg1 <- enrichMKEGG(gene = degsALL, organism = org, pvalueCutoff = pv, universe = uni, ...)
  ekeg2 <- enrichMKEGG(gene = degsUP, organism = org, pvalueCutoff = pv, universe = uni, ...)
  ekeg3 <- enrichMKEGG(gene = degsDOWN, organism = org, pvalueCutoff = pv, universe = uni, ...)
  return(list("MKEGG_ALL" = ekeg1, "MKEGG_UP" = ekeg2, "KEGG_MDOWN" = ekeg3))
}


compute_wikiPaths <- function(degsLFCList, org, orgDB, key, pv = 0.05, ...){
  # Compute the WikiPathways enrichments
  # Transform the key types to ENTREZ...
  nn <- bitr(names(degsLFCList), fromType = key, toType = "ENTREZID", OrgDb = orgDB)
  # Select the ones that we have a key...
  degsLFCList2 <- degsLFCList[nn[[key]]]
  names(degsLFCList2) <- nn$ENTREZID
  return(list("WikiPaths" = gseWP(degsLFCList2, organism = org, pvalueCutoff = pv, ...)))
}


compute_Reactome <- function(degs, org, orgDB, key, uni, pv = 0.05, ...){
  # Compute Reactome Paths enrichments
  # Transform the key types to ENTREZ...
  nn <- bitr(degs, fromType = key, toType = "ENTREZID", OrgDb = orgDB)
  return(list("REACTOME" = enrichPathway(nn$ENTREZID, organism = org, pvalueCutoff = pv, universe = uni, ...)))
}


plot_allEnrich <- function(allEgos, cats = 20, ti = "", ...) {
  # Manage the plotting of an allEGOs result.
  nms <- names(allEgos)
  for (nm in nms) {
    ego <- allEgos[[nm]]
    if (dim(data.frame(ego))[1] == 0) {
      plot.new()
      text(x = 0.5, y = 0.5, paste("No enrichment found for", nm, "in", ti))
    } else {
      nEnriched <- dim(data.frame(ego))[1]
      plot(dotplot(ego, showCategory = cats, title = paste(ti, nEnriched, "Enrichments of", nm, sep = " "), ...))
      # R error... FIXME plot(goplot(ego), ...)
    }
  }
}

plotAdvanced_Enrich <- function(ego, degsDF, ct1 = 25, cl = 5, ct2 = 4, ct3 = 25, ti = "...", ...){
  # Tree of GOs
  egoPT <- pairwise_termsim(ego)
  p1 <- tryCatch({treeplot(egoPT, showCategory = ct1, cluster.params = list(n = cl), ...)}, error = function(e) {plot.new(); text("Too few enrichments to construct a tree")}) # Enriched  GO clusters
  # Cnetplot
  # prepare the logFC object
  lfc <- sort(setNames(degsDF$logFC, rownames(degsDF)), decreasing = TRUE)
  p2 <- cnetplot(ego, layout = "randomly", foldChange = lfc, color_edge = "category", size_edge = 0.1, shadowtext = "category", showCategory = ct2,  ...)  # Coloured with the fold change
  # GO plot DAG
  p3 <- goplot(ego, showCategory = ct3, max.overlaps = Inf, ...)  # A way to look at the lower level of enrichments
  ggarrange(p1, p2, p3, ncol = 1, nrow = 3)
}



# Low level function to add counts on a boxplot.
n_fun <- function(x){
  return(data.frame(y = max(x), label = length(x)))
}
