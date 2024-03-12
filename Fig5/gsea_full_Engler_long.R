library(tidyverse)
library(ggplot2)
library(glue)
library(DESeq2)
library(biomaRt)
library(fgsea)
library(ggpubr)

# Read in a .csv with mouse/human genes name conversions
mouse_to_human <- read.csv("")

# Read in a csv for IAM gene signature (mouse gene names)
gsets <- read.csv("")

# Lets convert mouse to human genes
for (i in colnames(gsets)){
  genes <- gsets[[i]]
  new_genes <- sapply(genes, function(x){
    if (x != ""){
    genes <- mouse_to_human %>%
      dplyr::filter(MGI.symbol == x) %>%
      .$HGNC.symbol %>%
      .[1] # takes the first one...might not be the best choice
  } else {
    return("")
  }
  })
  new_genes[is.na(new_genes)] <- ""
  gsets[[i]] <- new_genes
}

# write out a file with converted mouse->human gene names 
write.csv(gsets, "")

# Read in the pseudobulk stats...these are the stats for our "OVA" analysis, where we compare each cluster to all other clusters
# This was generated previously for PMID: 37146132
pseudo_stats = read.csv("/Users/tk857/Dropbox (Partners HealthCare)/Lab/Scripts/Asthma Brush Project/myeloid_harmonized_pseudobulk_stats.csv")

# Get the gene set to a named list
mac_gene_sets <- lapply(colnames(gsets), function(n){
  genes <- gsets[[n]]
  genes <- genes[genes != ""]
  genes <- genes[!is.na(genes)]
  return(genes)
})
names(mac_gene_sets) <- colnames(gsets)

# The data has old cluster assignments, here are the corresponding clusters
annotations <- c("11" = "MC1 (CXCL10)",
                  "5" = "MC2 (SPP1)",
                  "6" = "MC3 (AREG)",
                  "7" = "Mac1 (FABP4)",
                  "13" = "quiesMac",
                  "2" = "quiesMC",
                  "14" = "Cycling (PCLAF)",
                  "1" = "MC4 (CCR2)",
                  "8" = "Mac2 (A2M)",
                  "4" = "pDC (TCF4)",
                  "9" = "migDC (CCR7)",
                  "10" = "DC1 (CLEC9A)",
                  "3" = "DC2 (CD1C)",
                  "12" = "AS DC (AXL)"
)

# A dataframe for which clusters to look at with which gene sets
pair_df = data.frame("cluster" = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14"),
                     "gene_set" = c("Engler", "Engler", "Engler", "Engler", "Engler", "Engler", "Engler", "Engler", "Engler", "Engler", "Engler", "Engler","Engler", "Engler"))

# initialize empty vector to store list of p values
pval_list = c()
name_list = c()

# duplicate to generate p-values

plot_list <- list()
for (i in 1:nrow(pair_df)){
  # Get your cluster, gene set
  cluster <- pair_df$cluster[i]
  gs <- pair_df$gene_set[i]

  # Get the relevant stats for that cluster
  pseudo_clust <- pseudo_stats[pseudo_stats$cluster == glue("cluster_{cluster}"),]

  # Create ranking based on log-fc, limit to genes with >5% expression
  ranks <- pseudo_clust %>%
    dplyr::filter(percent_cells > 5) %>%
    dplyr::arrange(desc(log_fc)) %>%
    dplyr::select(feature, log_fc) %>%
    deframe(.)

  # Run GSEA
  fgsea_res <- fgsea(pathways = mac_gene_sets, stats = ranks, nperm = 10000)

  annot <- annotations[[cluster]]

  # Get GSEA stats
  nes <- round(fgsea_res$NES[fgsea_res$pathway == gs], 3)
  pval <- round(fgsea_res$pval[fgsea_res$pathway == gs], 4)
  n_genes <- fgsea_res$size[fgsea_res$pathway == gs]
  pval_list <- append(pval_list, pval)
  name_list <- append(name_list, annotations[[cluster]])
}

# the saved file still needs FDR multiple test correction which can be used instad of p-values
# the fdr remains index matched to the cluster dataframes
fdr <- p.adjust(pval_list, method = "fdr")

plot_list <- list()
for (i in 1:nrow(pair_df)){
  # Get your cluster, gene set
  cluster <- pair_df$cluster[i]
  gs <- pair_df$gene_set[i]

  # Get the relevant stats for that cluster
  pseudo_clust <- pseudo_stats[pseudo_stats$cluster == glue("cluster_{cluster}"),]

  # Create ranking based on log-fc, limit to genes with >5% expression
  ranks <- pseudo_clust %>%
    dplyr::filter(percent_cells > 5) %>%
    dplyr::arrange(desc(log_fc)) %>%
    dplyr::select(feature, log_fc) %>%
    deframe(.)

  # Run GSEA
  fgsea_res <- fgsea(pathways = mac_gene_sets, stats = ranks, nperm = 10000)

  annot <- annotations[[cluster]]

  # Get GSEA stats
  nes <- round(fgsea_res$NES[fgsea_res$pathway == gs], 3)
  #pval <- round(fgsea_res$pval[fgsea_res$pathway == gs], 4)
  # replace with FDR from above
  fdr_val <- round(fdr[i], 4)
  n_genes <- fgsea_res$size[fgsea_res$pathway == gs]
  #pval_list <- append(pval_list, pval)

  # Make the GSEA plot
  rnk <- rank(-ranks)
  ord <- order(rnk)

  statsAdj <- ranks[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj) ^ 1)
  statsAdj <- statsAdj / max(abs(statsAdj))

  pathway <- unname(as.vector(na.omit(match(mac_gene_sets[[gs]], names(statsAdj)))))
  pathway <- sort(pathway)

  gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway,
                            returnAllExtremes = TRUE)

  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops

  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x=c(0, xs, n + 1), y=c(0, ys, 0))

  diff <- (max(tops) - min(bottoms)) / 8

  x=y=NULL

  p <- ggplot(toPlot, aes(x = x, y = y)) +
    # geom_point(color="blue", size=0.1) +
    geom_line(color="blue") +
    geom_hline(yintercept=0, colour="black") +
    geom_segment(data=data.frame(x=pathway),
                     mapping=aes(x=x, y=-0.15,
                                 xend=x, yend=-0.25),
                     size=0.4) +
    scale_y_continuous(expand = c(0.05,0.05)) +
    xlab("Rank") + ylab("Enrichment score") +
    geom_text(aes(label = "")) +
    annotate("text", label = glue("NES : {nes}"), x = length(ranks) - 1000, y  =0.9) +
    annotate("text", label = glue("FDR : {fdr_val}"), x = length(ranks) - 1000, y = 0.8) +
    annotate("text", label = glue("# genes : {n_genes}"), x = length(ranks) - 1000, y = 0.7) +
    ggtitle(glue("IAM signature : {annot}")) +
    theme_classic(base_size = 12)

  # p <- plotEnrichment(sig_gene_list[[gs]], ranks) +
  #   ggtitle(glue("{gs} signature : {annot}")) +
  #   geom_text(aes(label = "")) +
  #   annotate("text", label = glue("NES : {nes}"), x = length(ranks) - 1000, y  =1) +
  #   annotate("text", label = glue("p-value : {pval}"), x = length(ranks) - 1000, y = 0.93) +
  #   theme_classic()
  plot_list <- c(plot_list, list(p))
}

plots <- ggarrange(plotlist = plot_list, ncol = 3, nrow = 5)
# specify output directory and name below
ggsave("GSEA_plots_full.pdf", plots, useDingbats=FALSE, width = 16, height = 20, units="in", dpi = 300)

