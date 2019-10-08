
rm(list=ls())

######################################################################
# helper functions
######################################################################


####################################################################
# binquire() and inquire()
####################################################################

source("~/inquire.R")
# not yet tested

biocUpgrade <- function(){
  source("https://bioconductor.org/biocLite.R")
  biocLite("BiocUpgrade")
}

detach_package <- function(pkg, character.only = FALSE) {
  if(!character.only) {
    pkg <- deparse(substitute(pkg))
  }
  search_item <- paste("package", pkg, sep = ":")
  while(search_item %in% search()) {
    detach(search_item, unload = TRUE, character.only = TRUE)
  }
}

####################################################################
# required packages
####################################################################

binquire(DESeq2)     # DE analysis
# inquire(xlsx)        # xlsx printing
binquire(edgeR)      # DE analysis

binquire(biomaRt)     # for Annotations
# for this necessary:
# sudo apt install libssl-dev
binquire(org.Mm.eg.db)
binquire(org.Hs.eg.db)

binquire(GO.db)       # for GO analysis
binquire(GOstats)     # for GO analysis

binquire(gage)        # for KEGG analysis
binquire(gageData)
binquire(KEGG.db)
binquire(annotate)
binquire(genefilter)
binquire(vsn)

# inquire(Cairo)       # pdf output
# # necessary:
# # sudo apt install libcairo2-dev libxt-dev
# CairoFonts(
#   regular="Arial:style=Regular",
#   bold="Arial:style=Bold",
#   italic="Arial:style=Italic",
#   bolditalic="Arial:style=Bold Italic,BoldItalic",
#   symbol="Symbol")

inquire(ggplot2)     # for plots
binquire(Glimma)      # for interactive plots

#####################################################################
inquire(SummarizedExperiment)

inquire(ggplot2)               # for graphics
inquire(reshape)               # for graphics

inquire(gplots)                # for heatmap
inquire(RColorBrewer)          # for heatmap

binquire(pheatmap)             # for heatmap
require(gskb)                  # for gskb 
require(clusterProfiler)
require(GenomicFeatures)
require(openxlsx)
require(clusterProfiler)
require(org.Mm.eg.db)
require(plotly)
require(magrittr)

options(java.parameters = "-Xmx3000m")
# https://stackoverflow.com/questions/21937640/handling-java-lang-outofmemoryerror-when-writing-to-excel-from-r

####################################################################
# pseudoannotation
####################################################################

pseudodf <- function(vec) {
  df <- data.frame(GeneID = as.character(vec), symbol = as.character(vec), stringsAsFactors=FALSE)
  rownames(df) <- as.character(vec)
  colnames(df) <- "symbol"
  df
}

####################################################################
# non-overlapping gene lengths from gtf
# (required for DGEList obj and normalizations)
# thanks to Irsan # https://www.biostars.org/p/83901/
####################################################################

gtf2gene.length <- function(gtf.path) {
  gc()
  require(GenomicFeatures)
  txdb <- makeTxDbFromGFF(gtf.path, format = "gtf")
  exons.list.per.gene <- exonsBy(txdb, by="gene")
  exonic.gene.sizes <- as.data.frame(sum(width(reduce(exons.list.per.gene))))
  colnames(exonic.gene.sizes) <- "length"
  exonic.gene.sizes
}



#####################################################################
# Load openxlsx
#####################################################################

require(openxlsx)

#####################################################################
# list flattener
#####################################################################

once.flatten <- function(obj) unlist(obj, recursive = FALSE)
k.flatten <- function(obj, k) repeated(obj, k, once.flatten)
flatten <- function(obj) unlist(obj)

#####################################################################
# utility function look list elements by head
#####################################################################

lhead <- function(l) lapply(l, head)

#####################################################################
# read xlsx files to dfs list
#####################################################################

xlsx2df.list <- function(xlsx.path, rowNames = TRUE, colNames = TRUE, ...) {
  wb <- loadWorkbook(xlsx.path)
  sheetNames <- names(wb)
  res <- lapply(sheetNames, function(sheetName) {
    read.xlsx(wb, sheet = sheetName, rowNames = rowNames, colNames = colNames, ...)
  })
  names(res) <- sheetNames
  res
}


#####################################################################
# printing dfs to xlsx files
#####################################################################

withNames <- function(...) {
  # Returns a list constructed by using 
  # alterning name, obj, name, obj arguments
  p.l <- list(...)
  len <- length(p.l)
  if (len %% 2 == 1) {
    stop("withNames call with odd numbers of arguments")
    print()
  }
  seconds <- p.l[seq(2, len, 2)]
  firsts <- p.l[seq(1, len, 2)]
  names(seconds) <- unlist(firsts)
  seconds
}

write.dfs <- function(df.list, fpath) {
  wb <- createWorkbook()
  Map(function(data, name) {
    addWorksheet(wb, name)
    writeData(wb, name, data, rowNames = TRUE, colNames = TRUE)
  }, df.list, names(df.list))
  saveWorkbook(wb, file = fpath, overwrite = TRUE)
}


#####################################################################
# select a df by list of rownames -> df.list
#####################################################################

select.by.names.vec.list <- function(data.df, list.of.names.vecs) {
  lapply(list.of.names.vecs, function(names.vec) data.df[names.vec, ])
}

#####################################################################
# Return for a df, its sortings by l2FC, p-val and rownames
#####################################################################

df.sortings <- function(DE.res.df) {
  list(l2FC.srt = DE.res.df[order(DE.res.df$log2FoldChange, decreasing = TRUE), ],
       p.srt = DE.res.df[order(DE.res.df$pvalue, decreasing = FALSE), ],
       names.srt = DE.res.df[order(rownames(DE.res.df), decreasing = FALSE), ])
}

#####################################################################
# print DE sortings
#####################################################################

print.DE.sortings <- function(DE.list, fname, dir = ".") {
  DE.list.with.sortings <- lapply(once.flatten(lapply(DE.list, df.sortings)), as.data.frame)
  write.dfs(DE.list.with.sortings, file.path(dir, fname))
}

#####################################################################
# print cnts selections of DE 
#####################################################################

print.cnts.DE.sortings <- function(cnts, DE.list, fname, dir = ".") {
  DE.list.with.sortings <- lapply(once.flatten(lapply(DE.list, df.sortings)),
                                  as.data.frame)
  DE.list.with.sortings.names <- lapply(DE.list.with.sortings, rownames)
  DE.cnts.list.with.sortings <- select.by.names.vec.list(as.data.frame(cnts),
                                                         DE.list.with.sortings.names)
  write.dfs(DE.cnts.list.with.sortings, file.path(dir, fname))
}



####################################################################
# helper functions for timestamp
####################################################################

time.now <- function() format(Sys.time(), "%y%m%d%H%M%S")

#####################################################################
# helper functions for averaging tables
#####################################################################

counts.avg.build <- function(cnts.df, cols.list, groups){
  cnts.df <- as.data.frame(cnts.df)
  ncol_old <- length(colnames(cnts.df))
  ncol_limit <- length(groups) + ncol_old
  new_col_names <- c(colnames(cnts.df), groups)
  cnts.df <- cbind(cnts.df,
                   lapply(cols.list,
                          function(idxs) rowMeans(cnts.df[, idxs, drop = FALSE])))
  colnames(cnts.df) <- new_col_names
  cnts.df[, (ncol_old + 1):ncol_limit]
}

counts.std.build <- function(cnts.df, cols.list, groups){
  cnts.df <- as.data.frame(cnts.df)
  rowSds <- function(df) apply(df, 1, sd)
  ncol_old <- length(colnames(cnts.df))
  ncol_limit <- length(groups) + ncol_old
  new_col_names <- c(colnames(cnts.df), groups)
  cnts.df <- cbind(cnts.df,
                   lapply(cols.list,
                          function(idxs) rowSds(cnts.df[, idxs, drop = FALSE])))
  colnames(cnts.df) <- new_col_names
  cnts.df[, (ncol_old + 1):ncol_limit]
}

# ####################################################################
# Read-in meta.df and indirpath to count-matrix
# ####################################################################

read.tab <- function(fpath) {
  read.delim(fpath, sep = "\t", head = F, row.names = 1, stringsAsFactors = F)
}

read.dfs2table <- function(meta.df, indirpath) {
  files <- file.path(indirpath, as.character(meta.df$fileName))
  df.list <- lapply(files, read.tab)
  res.df <- Reduce(cbind, df.list)
  names(res.df) <- meta.df$sampleName
  res.df
}


# ####################################################################
# Translate vector values
# ####################################################################

translate.vec <- function(values.vec, from.vec, to.vec) {
  tr.vec <- to.vec
  names(tr.vec) <- from.vec
  tr.vec[values.vec]
}

# ####################################################################
# cluster df sorting helper functions
# ####################################################################

# df to list and back ################################################

select_df <- function(df, val, col.selector) {
  df[ df[, col.selector] == val, ]
}

df2dflist <- function(df, col.selector) { # actually it is split()
  col.vals <- unique(df[, col.selector])
  dfl <- lapply(seq(col.vals),
                function(i) select_df(df,
                                      val = col.vals[i],
                                      col.selector))
  names(dfl) <- col.vals
  dfl
}

dflist2df <- function(dfl) {
  Reduce(rbind, dfl)
}

# ordering df lists  #################################################

order.list.by.col.mean <- function(dfl, col.selector) {
  dfl[order(unlist(lapply(dfl, function(df) mean(df[, col.selector]))))]
}

order.list.by.df.function <- function(dfl, df.function, ...) {
  dfl[order(unlist(lapply(dfl, function(df) df.function(df, ...))))]
}

list.order.by.col.mean <- function(dfl, col.selector) {
  order(unlist(lapply(dfl, function(df) mean(df[, col.selector]))))
}

list.order.manually.by.vec <- function(dfl, vec) {
  dfl[vec]
}

repeated.values.into.df.list <- function(dfl, col.selector, vals) {
  Map(function(df, val) {df[, col.selector] <- val; df}, dfl, vals)
}

list.order.by.col.mean.diffs <-
  function(dfl, col.selector.one, col.selector.two, decreasing = FALSE) {
    dfl[order(unlist(lapply(dfl, function(df) mean(df[, col.selector.one]))) -
                unlist(lapply(dfl, function(df) mean(df[, col.selector.two]))),
              decreasing = decreasing)]
  }

#######################################################################


#######################################################################
# create a  DESeq2 obj out of
# - path to counts
# - file names in meta.txt
# - (outputpath)
#######################################################################

meta2DESeq2.obj <- function(meta.df, indirpath, normalized=FALSE) {
  # require(DESeq2)
  count.matrix <- read.dfs2table(meta.df, indirpath)
  DESeq2.obj <- DESeqDataSetFromMatrix(
    countData = count.matrix,
    colData = meta.df,
    design = ~ 0 + condition)
  if (normalized) {
    DESeq2.obj <- estimateSizeFactors(DESeq2.obj)
    DESeq2.obj <- estimateDispersions(DESeq2.obj)
  }
  DESeq2.obj
}



#######################################################################
# create a raw count table out of
# - path to counts
# - file names in meta.txt
# - (outputpath)
# create a DESeq2-normalized table out of
# - path to counts
# - file names in meta.txt
# - (outputpath)
# create an averaged DESeq2-normalized table out of
# - path to counts
# - file names in meta.txt
# - (outputpath)
#######################################################################

meta2cnts <- function(meta.df, DESeq2.obj, outdirpath=".",
                      dataname = dataname,
                      printp=FALSE, normalized=FALSE, averaged=FALSE,
                      sheetName) {
  cnts <- counts(DESeq2.obj, normalized=normalized)
  if (averaged) {
    cond <- unique(as.character(meta.df$condition))
    cond.cols <- lapply(cond,
                        function(cd) which(as.character(meta.df$condition) == cd))
    names(cond.cols) <- cond
    cnts <- counts.avg.build(cnts, cond.cols, cond)
  }
  if (printp) {
    filename <- paste0(ifelse(normalized, "nrm-counts-", "raw-counts-"),
                       ifelse(averaged, "avg-", ""),
                       dataname, "-", core.name(meta.df), ".xlsx")
    write.dfs(withNames(ifelse(normalized, "nrm-counts", "raw-counts"), cnts),
              file.path(outdirpath, filename))
  }
  cnts
}




#######################################################################
# create list of differentially expressed genes
# create list of differentially upregulated genes
# create list of differentially downregulated genes
# print them out
# out of
# path to counts
# file name in meta.txt - metapath
# also for single-rep RNAseq analysis! ("DESeq2")
#######################################################################


DEanalysis <- function(meta.df, DESeq2.obj.disp, outdirpath=".", dataname="",
                       printp=FALSE, prior.count=0,
                       alpha=0.05, lFC=1, filterp=FALSE) {
  dds <- DESeq(DESeq2.obj)
  res <- results(dds, contrast=c("condition", num(meta.df), denom(meta.df)),
                 cooksCutoff = Inf,
                 independentFiltering = FALSE) # those two avoid NA!
  if (filterp) {
    res <- subset(subset(res, padj < alpha), abs(log2FoldChange) > lFC)
  }
  
  up <- subset(res, res$log2FoldChange > 0)
  down <- subset(res, res$log2FoldChange < 0)
  res <- list(all=res, up=up, down=down)
  
  if (printp) {
    filecorename <- paste0("DE_", ifelse(filterp, "sig_", ""),
                           dataname, "_", core.name(meta.df), "_", time.now())
    if (filterp) { filecorename <- paste0(filecorename, "_", alpha, "_", lFC) }
    filename <- paste0(filecorename, ".xlsx")
    print.DE.sortings(res, fname = filename, dir = DE.outdir)
  }
  res
}




#######################################################################
# create an averaged heatmap and group-pictures out of
# - k
# - (DESeq2 obj OR DEGList obj OR raw count table OR cpm-normalized table OR
#   DESeq2-normalized table)
#   metapath and indirpath
# - (outputpath)
#######################################################################

meta2heatmap <- function(meta.df, cnts.avg.nrm, resSig, outdirpath=".",
                         dataname = dataname, selected.genes = NULL, name.add = "", # for showing in name
                         k = k, printp=FALSE,
                         alpha = alpha, lFC= lFC, filterp=FALSE,
                         xlsxp=TRUE, csvp=FALSE, tsvp=FALSE) {
  res.names <- rownames(resSig$all)
  
  ## if sth given in 'selected.genes': if "up" or "down" use them from resSig, else the vector given:
  if (length(selected.genes) > 0) {
    if (selected.genes == "up") {
      res.names <- rownames(resSig$up)
    } else if (selected.genes == "down") {
      res.names <- rownames(resSig$down)
    } else {
      res.names <- selected.genes
    }
  }
  
  cnts.res <- cnts.avg.nrm[res.names, ]
  
  # gene-wise normalization
  scaledata <- t(scale(t(cnts.res)))
  scaledata <- scaledata[complete.cases(scaledata), ]
  
  # k means clustering
  kClust <- kmeans(scaledata, centers = k, nstart = 1000, iter.max = 30)
  kClusters <- kClust$cluster
  
  # function to find centroid (cluster core) in cluster i
  clust.centroid <- function(i, dat, clusters) {
    ind = (clusters == i)
    colMeans(dat[ind, ])
  }
  kClustcentroids <- sapply(levels(factor(kClusters)),
                            clust.centroid, scaledata, kClusters) ## is a matrix
  
  # plot centroids
  Kmolten <- melt(kClustcentroids)
  colnames(Kmolten) <- c("sample", "cluster", "value")
  # ensure correct factorizing
  Kmolten$sample <- factor(Kmolten$sample,
                           levels = unique(meta.df$condition))
  {
    p1 <- ggplot(Kmolten, aes(x=factor(sample, levels = unique(meta.df$condition)),
                              y=value, group = cluster,
                              colour=as.factor(cluster))) +
      geom_point() +
      geom_line() +
      xlab("Time") +
      ylab("Expression") +
      labs(title = "Cluster Expression by Group", color = "Cluster")
    png(paste0(outdirpath, "/k", k, "_", core.name(meta.df), paste0("-ClusterAll_", alpha, "_", lFC, "_", time.now(), ".png")))
    print(p1)
    dev.off()
  }
  
  # check similarity of centroids
  print(cor(kClustcentroids))
  
  for (i in 1:k) {
    # subset cores molten df to plot core
    assign(paste0("core", i), Kmolten[Kmolten$cluster == i, ])
    eval(parse(text = paste0("core", i, "$sample <- factor(core", i,
                             "$sample, levels = unique(meta.df$condition))")))
    # get clusters
    assign(paste0("K", i), scaledata[kClusters == i, ])
    
    # calculate correlation with core
    assign(paste0("corescore", i),
           eval(parse(text = paste0("function(x) {cor(x, core", i, "$value)}"))))
    assign(paste0("score", i),
           eval(parse(text = paste0("apply(K", i, ", 1, corescore", i, ")"))))
    
    # get data frame into long format for plotting
    assign(paste0("K", i, "molten"),
           eval(parse(text = paste0("melt(K", i, ")"))))
    eval(parse(text = paste0("colnames(K", i, "molten) <- c('gene', 'sample', 'value')")))
    
    # add the score
    eval(parse(text = paste0("K", i, "molten <- merge(K", i, "molten, score", i, ", by.x = 'gene', by.y = 'row.names', all.x = T)")))
    eval(parse(text = paste0("colnames(K", i, "molten) <- c('gene', 'sample', 'value', 'score')")))
    
    # order
    eval(parse(text = paste0("K", i, "molten <- K", i, "molten[order(K", i, "molten$score), ]")))
  }
  
  # plot cluster groups
  for (i in 1:k) {
    text = paste0("sp", i, " <- ggplot(K", i, "molten, aes(x=factor(sample,",
                  " levels=unique(meta.df$condition)), y=value)) + ",
                  "geom_line(aes(colour=score, group=gene)) + ",
                  "scale_color_gradientn(colours=c('blue1', 'red2')) + ",
                  # this adds the core
                  "geom_line(data=core", i, ", aes(sample,value,group=cluster), ",
                  "color='black', inherit.aes=FALSE) + ",
                  "xlab('Time') + ",
                  "ylab('Expression') + ",
                  "labs(title='Cluster ", i, " Expression by Group', color = 'Score'); ",
                  "png('", outdirpath, "/k", k, "_", core.name(meta.df), "-Cluster", i, "_", dataname, "_", alpha, "_", lFC, name.add, ".png'); ",
                  "print(sp", i, "); dev.off()"
    )
    eval(parse(text = text))
  }
  
  # prepare heatmap
  colors.kr <- colorRampPalette(c("black", "red"))(100)
#  colors.kgr <- colorRampPalette(c("black", "grey88", "red"))(100)
  
  eval(parse(text = paste0("scores <- c(", paste(paste("score", 1:k, sep = ''),
                                                 collapse = ", "), ")")))
  # add cluster number and score for sorting the data
  scaledata.k <-cbind(scaledata,
                      cluster = kClust$cluster,
                      score = scores[rownames(scaledata)])
  scaledata.k <- scaledata.k[order(scaledata.k[, "cluster"], scaledata.k[, "score"]), ]
  
  # outpath
  outname <- file.path(outdirpath, paste0("k", k, "_khmap_", dataname, "_", time.now(), "_UO_", alpha, "_", lFC, name.add))
  
  # for gaps
  gaps.idxs <- cumsum(table(scaledata.k[, "cluster"])) # scaledata.k is one matrix! only column 'cluster'
  # dim(scaledata.k)
  # unordered clusters
  
  ##################
  # black to red
  ##################
  
  {
    svg(paste0(outname, ".svg"))
    pheatmap(scaledata.k[, 1:length(unique(meta.df$condition))],
             cluster_rows = F,
             cluster_cols = F,
             cellwidth = 40,
             col = colors.kr,
             fontsize_row = 0.5,
             border_color = NA,
             gaps_row = gaps.idxs          # gap after each block
    )
    dev.off()
  }
  
  {
    setEPS()
    postscript(paste0(outname, ".eps"))
    pheatmap(scaledata.k[, 1:length(unique(meta.df$condition))],
             cluster_rows = F,
             cluster_cols = F,
             cellwidth = 40,
             col = colors.kr,
             fontsize_row = 0.5,
             border_color = NA,
             gaps_row = gaps.idxs          # gap after each block
             
    )
    dev.off()
  }
  

#   ##################
#   # black grey red
#   ##################
#   
#   {
#     svg(paste0(outname, "kgr.svg"))
#     pheatmap(scaledata.k[, 1:length(unique(meta.df$condition))],
#              cluster_rows = F,
#              cluster_cols = F,
#              cellwidth = 40,
#              col = colors.kgr,
#              fontsize_row = 0.5,
#              border_color = NA,
#              gaps_row = gaps.idxs          # gap after each block
#     )
#     dev.off()
#   }
#   
#   {
#     setEPS()
#     postscript(paste0(outname, "kgr.eps"))
#     pheatmap(scaledata.k[, 1:length(unique(meta.df$condition))],
#              cluster_rows = F,
#              cluster_cols = F,
#              cellwidth = 40,
#              col = colors.kgr,
#              fontsize_row = 0.5,
#              border_color = NA,
#              gaps_row = gaps.idxs          # gap after each block
#              
#     )
#     dev.off()
#   }
# 


  # print scaledata
  scaledata_list <- df2dflist(scaledata.k, "cluster")
  outfname <- paste0(outname, '.xlsx')
  names(scaledata_list) <- paste("cluster", 1:k, sep = "")
  write.dfs(scaledata_list, outfname)
  
  eval(parse(text = paste0("Kmolten.list <- list(", paste(paste("K", 1:k, "molten", sep = ''), collapse = ", "), ")")))
  names(Kmolten.list) <- paste("K", 1:k, "molten", sep = '')
  
  eval(parse(text = paste0("core.list <- list(", paste(paste("core", 1:k, sep = ''), collapse = ", "), ")")))
  names(Kmolten.list) <- paste("core", 1:k, sep = '')
  
  eval(parse(text = paste0("score.list <- list(", paste(paste("score", 1:k, sep = ''), collapse = ", "), ")")))
  names(Kmolten.list) <- paste("score", 1:k, sep = '')
  
  data <- list(scaledata_list, Kmolten.list, core.list, score.list)
  names(data) <- c("scaledata_list", "Kmolten.list", "core.list", "score.list")
  saveRDS(data, file = paste0(outdirpath, "/scaledata.Kmolten.core.score.list.", dataname, ".", time.now(), ".", name.add, ".rds"))
  data
}


#######################################################################
# create an averaged heatmap and group-pictures out of
# unordered list (xlsx) and an ordervector
# - path to unordered .xlsx file
# - ordervector
# - outpath
#######################################################################

unordered2orderedheatmap <- function(UO.heatmap.path, order.vec, outdirpath=".",
                                     sheet="all", gaps.after.blocks=c(),
                                     alpha=alpha, lFC=lFC, dataname=dataname) {
  
  dfl <- xlsx2df.list(UO.heatmap.path)
  # remove 'all' sheet
  all.idx <- which(names(dfl) == "all")
  if (length(all.idx > 0)) {
    dfl <- dfl[-all.idx]
  }
  # within cluster sort by score
  dflsc <- lapply(dfl, function(df) df[order(df[, "score"]), ])
  dflo <- list.order.manually.by.vec(dflsc, order.vec)
  k <- length(dfl)
  names(dflo) <- 1:k
  dflo.corr <- repeated.values.into.df.list(dflo, "cluster", 1:k)
  scaledata.k.ordered <- dflist2df(dflo.corr)
  if (length(gaps.after.blocks) > 0) {
    gaps_idxs <- cumsum(unlist(lapply(df2dflist(scaledata.k.ordered, "cluster"),
                                      function(df) dim(df)[1])))
    gaps_idxs <- gaps_idxs[gaps.after.blocks]
  }
  dir.create(outdirpath, recursive = T, showWarnings = F)
  outname <- file.path(outdirpath, paste0("k", k, "_khmap_", dataname, "_", time.now(), "_O_", alpha, "_", lFC))
  outfname <- paste0(outname, ".xlsx")
  colfunc <- colorRampPalette(c("black", "red"))
  
  {
    setEPS()
    postscript(paste0(outname, ".eps"))
    # svg(paste0(outname, ".svg"))
    pheatmap(scaledata.k.ordered[, 1:(length(colnames(dfl[[1]]))-2)], # without cluster scores
             cluster_rows = F,
             cluster_cols = F,
             cellwidth = 40,
             gaps_row = ifelse(length(gaps.after.blocks) > 0, list(gaps_idxs), list(dim(df)[1]))[[1]],
             border_color = NA,
             col = colfunc(100),
             fontsize_row = 0.5
    )
    dev.off()
  }
  
  {
    svg(paste0(outname, ".svg"))
    pheatmap(scaledata.k.ordered[, 1:(length(colnames(dfl[[1]]))-2)], # without cluster scores
             cluster_rows = F,
             cluster_cols = F,
             cellwidth = 40,
             gaps_row = ifelse(length(gaps.after.blocks) > 0, list(gaps_idxs), list(dim(df)[1]))[[1]],
             border_color = NA,
             col = colfunc(100),
             fontsize_row = 0.3
    )
    dev.off()
  }
  
  scaledata_list <- df2dflist(scaledata.k.ordered, "cluster")
  names(scaledata_list) <- paste("cluster", 1:k, sep = "")
  write.dfs(scaledata_list, outfname)
}



UO2O.hm <- function(filepath, s.K.c.s.list, order.vec, outdirpath=".",
                    sheet="all", gaps.after.blocks=c(),
                    alpha=alpha, lFC=lFC, dataname=dataname) {
  
  colfunc <- colorRampPalette(c("black", "red"))
  
  if (filepath == "") {
    scaledata.Kmolten.core.score.list <- s.K.c.s.list
  } else if (endsWith(filepath, suffix = ".rds")) {
    scaledata.Kmolten.core.score.list <- readRDS(filepath)
  } else {
    load_obj <- function(f) {
      env <- new.env()
      nm <- load(f, env)[1]
      env[[nm]]
    }
    scaledata.Kmolten.core.score.list <- load_obj(filepath)
  }
  scaledata_list <- scaledata.Kmolten.core.score.list$scaledata_list
  new.scaledata <- list.order.manually.by.vec(scaledata_list, order.vec)
  k <- length(scaledata_list)
  names(new.scaledata) <- 1:k
  new.scaledata.corr <- repeated.values.into.df.list(new.scaledata, "cluster", 1:k)
  
  scaledata.k.ordered <- dflist2df(new.scaledata.corr)
  if (length(gaps.after.blocks) > 0) {
    gaps_idxs <- cumsum(unlist(lapply(df2dflist(scaledata.k.ordered, "cluster"),
                                      function(df) dim(df)[1])))
    gaps_idxs <- gaps_idxs[gaps.after.blocks]
  }
  outname <- file.path(outdirpath, paste0("k", k, "_khmap_", dataname, "_", time.now(), "_O_", alpha, "_", lFC))
  outfname <- paste0(outname, ".xlsx")
  df <- scaledata_list[[1]] # exemplarisches df
  {
    # setEPS()
    # postscript(paste0(outname, ".eps"))
    svg(paste0(outname, ".svg"))
    pheatmap(scaledata.k.ordered[, 1:(length(colnames(df))-2)], # without cluster scores
             cluster_rows = F,
             cluster_cols = F,
             cellwidth = 40,
             gaps_row = ifelse(length(gaps.after.blocks) > 0, list(gaps_idxs), list(dim(df)[1]))[[1]],
             border_color = NA,
             col = colfunc(100),
             fontsize_row = 0.5
    )
    dev.off()
  }
  
  {
    svg(paste0(outname, ".svg"))
    pheatmap(scaledata.k.ordered[, 1:(length(colnames(df))-2)], # without cluster scores
             cluster_rows = F,
             cluster_cols = F,
             cellwidth = 40,
             gaps_row = ifelse(length(gaps.after.blocks) > 0, list(gaps_idxs), list(dim(df)[1]))[[1]],
             border_color = NA,
             col = colfunc(100),
             fontsize_row = 0.3
    )
    dev.off()
  }
  scaledata_list <- df2dflist(scaledata.k.ordered, "cluster")
  names(scaledata_list) <- paste("cluster", 1:k, sep = "")
  write.dfs(scaledata_list, outfname)
  
  Kmolten.list <- scaledata.Kmolten.core.score.list$Kmolten.list
  new.Kmolten.list <- Kmolten.list[order.vec]
  names(new.Kmolten.list) <- paste("K", 1:k, "molten", sep = '')
  for (i in 1:k) {
    eval(parse(text = paste0("K", i, "molten <- new.Kmolten.list[[", i, "]]")))
  }
  
  core.list <- scaledata.Kmolten.core.score.list$core.list
  new.core.list <- core.list[order.vec]
  names(new.core.list) <- paste("core", 1:k, sep = '')
  for (i in 1:k) {
    eval(parse(text = paste0("core", i, " <- new.core.list[[", i, "]]")))
  }
  
  if (endsWith(filepath, suffix = ".rds")) {
    score.list <- scaledata.Kmolten.core.score.list$score.list
    new.score.list <- score.list[order.vec]
    names(new.score.list) <- paste("score", 1:k, sep = '')
    for (i in 1:k) {
      eval(parse(text = paste0("score", i, " <- new.score.list[[", i, "]]")))
    }
  }
  
  # plot cluster groups
  for (i in 1:k) {
    text = paste0("sp", i, " <- ggplot(K", i, "molten, aes(x=factor(sample,",
                  " levels=unique(meta.df$condition)), y=value)) + ",
                  "geom_line(aes(colour=score, group=gene)) + ",
                  "scale_color_gradientn(colours=c('blue1', 'red2')) + ",
                  # this adds the core
                  "geom_line(data=core", i, ", aes(sample,value,group=cluster), ",
                  "color='black', inherit.aes=FALSE) + ",
                  "xlab('Time') + ",
                  "ylab('Expression') + ",
                  "labs(title='Cluster ", i, " Expression by Group', color = 'Score'); ",
                  "png('", outdirpath, "/k", k, "_", core.name(meta.df), "-Cluster", i, "_", dataname, "_", alpha, "_", lFC, ".png'); ",
                  "print(sp", i, "); dev.off()"
    )
    eval(parse(text = text))
  }
  
}

#######################################################################
# MDS plot, iMDS plot
# - metapath
# - indirpath
# - lFC, alpha
# - method
# - outdirpath
#######################################################################

meta2iMDSplot <- function(meta.df, DESeq2.obj.disp, outdirpath=".",
                          dataname = dataname, top=500,
                          launch = TRUE) {
  title <- paste0("MDS Plot ", dataname, " ", time.now())
  filename <- paste0("iMDS_", dataname, "_", core.name(meta.df), "_", time.now(), "_DESeq2", collapse = "_")
  glMDSPlot(DESeq2.obj,
            top = top,
            path = outdirpath,
            main = title,
            html = filename,
            launch = launch
  )
}


#######################################################################
# volcano plot, iVolcano Plot, MD plot, iMD plot
# - metapath
# - indirpath
# - lFC, alpha
# - method
# - outdirpath
#######################################################################

meta2iVolcano <- function(meta.df, DESeq2.obj, DESeq2.obj.disp, outdirpath=".",
                          dataname= dataname,
                          lFC=lFC, alpha=alpha, launch = TRUE) {
  wald.test <- nbinomWaldTest(DESeq2.obj.disp)
  res.DESeq2 <- results(wald.test, alpha=alpha, pAdjustMethod = "BH",
                        contrast = c("condition", num(meta.df), denom(meta.df)),
                        cooksCutoff = Inf,
                        independentFiltering = FALSE) # avoid NA in padj
  title <- paste0("Volcano Plot ", dataname, " ", time.now())
  filename <- paste0("iVolcano_", dataname, "_", core.name(meta.df), "_", time.now(), "_", alpha, "_", lFC, "_DESeq2")
  {
    glXYPlot(x=res.DESeq2$log2FoldChange, y=-log10(res.DESeq2$pvalue),
             counts = counts(DESeq2.obj)[rownames(res.DESeq2), ],
             anno = pseudodf(rownames(counts(DESeq2.obj))),
             groups = meta.df$condition,
             samples = meta.df$sampleName,
             xlab = "log2FC",
             ylab = "log10padj",
             main = title,
             status = ifelse(abs(res.DESeq2$log2FoldChange) > lFC &
                               res.DESeq2$padj < alpha, 1, 0),
             side.main = "symbol",
             side.xlab = "Group",
             side.ylab = "Counts",
             path = outdirpath,
             html = filename,
             launch = launch)
  }
  title <- paste0("iMD Plot ", dataname, " ", time.now())
  filename <- paste0("iMD_", dataname, "_", core.name(meta.df), "_", time.now(), "_DESeq2")
  {
    glMDPlot(x = res.DESeq2,
             counts = counts(DESeq2.obj)[rownames(res.DESeq2), ],
             anno = pseudodf(rownames(counts(DESeq2.obj))), # GeneID and symbol as col
             groups = factor(meta.df$condition, levels=unique(meta.df$condition)),
             samples = meta.df$sampleName,
             ylab = "log2FC",
             xlab = "Average log10 CPM",
             main = title,
             status = ifelse(abs(res.DESeq2$log2FoldChange) > lFC &
                               res.DESeq2$padj < alpha, 1, 0),
             side.xlab = "Group",
             side.ylab = "Counts",
             side.main = "symbol",
             path = outdirpath,
             html = filename,
             launch = launch)
  }
}



#####################################################################
# GOI
#####################################################################

bra.down <- c("Msgn1", "Osr1", "Rspo3", "Fgf8", "Wnt3a")
eo.down.repr <- c("Dmdx1", "Dpf3", "Foxa1", "Hey1", "Hhex", "Tcf7l2", "Tle2")
eo.down <- c("Mesp1", "Foxa2", "Sox17", "Lhx1", "Cer1")
pp.up <- c("Lefty1", "Lefty2", "Nodal", "Wnt8a", "Fgf5", "Otx2", "Cldn6")
episc.up <- c("Nanog", "Pou5f1", "Sox2", "L1td1", "Utf1")
ne.up <- c("Sox1", "Sox3", "Olig3", "Gbx2", "Pou3f1", "Msx3", "Slc7a3", "Zic1", "Zic2", "Nkx1-2", "Epha2", "Efnb1", "Nkx6-2")
me.down <- c("Pdgfra", "Foxc1", "Foxc2", "Isl1", "Kdr", "Mixl1", "Hand1", "Myh6")

gois <- c(me.down, ne.up, episc.up, pp.up, eo.down, eo.down.repr, bra.down)

other.gois <- c("T", "Wnt3", "Wnt3a", "Nodal", "Fgf8", "Mixl1")

#####################################################################
# List Genes
#####################################################################

down.genes <- c("Foxc2", "Gsc", "Mixl1", "Foxc1", "Pdgfra", "Kdr", "Tbx1", 
                "Amot", "Eya1", "Prrx1", "Lhx1", "Tnnt1", "Isl1", "Tbx20", 
                "Myh7", "Myh6", "Eya2", "Tbx18", "Snai1", "Snai2", "Hand1", 
                "Hand2", "Mesp1", "Pax2", "Mesp2", "Col2a1", "Myocd", "Pax3", 
                "Wt1", "Dll3", "Prickle1", "Msgn1", "Rspo3", "Osr1", "Twist1", 
                "Twist2", "Tbx6", "Meox1", "Gata6", "Gata4", "Foxa2", "Sox17", 
                "Lama1", "Tgfa", "Fn1", "Cer1", "Tgfb2", "Bmper", "Chrd", "Lgr5", 
                "Bmp2", "Fzd7", "Wnt3a", "Bmp7", "Cfc1", "Dkk1", "Fgf1", 
                "Tdgf1", "Wnt3", "Tgfb1", "Nodal", "Dact1", "Bmp4", "Dll1", 
                "Hey1", "Hhex", "Id1", "Dmbx1", "Dpf3", "Sall1", "Foxa1", 
                "Tcf7l2", "Hesx1", "Tcf7l1", "Hdac7", "Otx2", "Fbn2", "Otx1")

up.genes <- c("Six3", "Fabp7", "Ntrk2", "Zic1", "Mab21l2", "Nefl", "Olig3", 
              "Pou4f2", "Nkx1-2", "Sox3", "Sema3c", "Epha2", "Zic3", "Efnb1", 
              "Radil", "Syt11", "Slc7a3", "Nes", "Zic5", "Snph", "Nfasc", 
              "Pou4f1", "Elavl3", "Nrcam", "Grin1", "Gfra3", "Phox2a", "Nkx6-2",
              "Pax6", "Nsg2", "Msx3", "Ephb1", "Ncan", "Nova2", "Zic2", "Grik3",
              "Epha1", "Bcl11a", "Hoxa2", "Tubb3", "Sox1", "Neurod1", "Neurog1", 
              "Stmn2", "Atxn1", "Cntn2", "Neurl1a", "Sema3c", "Gap43", "Fgf5", 
              "Tbx3", "Cldn6", "Pou3f1", "Lefty2", "Gbx2", "Nanog", "L1td1", 
              "Wnt8a", "Pou5f1", "Lefty1", "Utf1", "Esrrb", "Sox2", "Lin28a", 
              "Dnmt3a", "Cbx7", "Kdm2b", "Atf7ip2")

meta2gois <- function(meta.df, cnts.avg.nrm, cnts.nrm, res, gois, outdirpath=".",
                      dataname = dataname, title.leader,
                      alpha = alpha, lFC= lFC) {
  cnts.gois <- cnts.avg.nrm[gois, ]
  
  # gene-wise normalization
  scaledata <- t(scale(t(cnts.gois)))
  scaledata <- scaledata[complete.cases(scaledata), ]
  
  # prepare heatmap
  colfunc <- colorRampPalette(c("black", "red"))
  
  # outpath
  outname <- file.path(outdirpath, paste0(title.leader, dataname, "_", time.now(), "_UO_", alpha, "_", lFC))
  
  # unordered clusters
  {
    setEPS()
    postscript(paste0(outname, ".eps"))
    # svg(paste0(outname, ".svg"))
    pheatmap(scaledata,
             cluster_rows = F,
             cluster_cols = F,
             cellwidth = 40,
             col = colfunc(100),
             fontsize_row = 1,
             border_color = NA
    )
    dev.off()
  }
  
  # print counts
  gois.present <- gois[gois %in% rownames(cnts.nrm)]
  res.all.gois <- gois[gois %in% rownames(res$all)]
  res.up.gois <- gois[gois %in% rownames(res$up)]
  res.down.gois <- gois[gois %in% rownames(res$down)]
  
  res.gois <- list(cnts.nrm[gois.present, ], cnts.gois, res.all.gois, res.up.gois, res.down.gois)
  names(res.gois) <- c("goi.counts", "goi.counts.avg", "goi.DE.all", "goi.DE.up", "goi.DE.down")
  write.dfs(res.gois, paste0(outname, ".xlsx"))
}


###########################
# Nice PCA plots
###########################

#######################################################################
# plotly SVG
#######################################################################
require(plotly)

ensure_package <- Vectorize(
  FUN =
    function(package,
             load = TRUE)
    {
      installed_packages <- rownames(installed.packages())
      
      if ( !(package %in% installed_packages) )
      {
        install.packages(package,
                         repos = "https://cloud.r-project.org/")
      }
      if ( load ) library(package, character.only = TRUE)
    }
)

export_plotly2SVG <- function(plotly_graph,
                              filename = NULL,
                              parent_path = paste0(getwd(), "/output"),
                              width = 800,
                              height = 600,
                              remove_title = FALSE,
                              font_family = "Arial",
                              incl_PDF_copy = FALSE,
                              incl_PNG_copy = FALSE,
                              png_scaling_factor = 1.8,
                              autocrop_png = TRUE)
{
  ensure_package(package = c("dplyr",
                             "plotly",
                             "readr",
                             "RSelenium",
                             "rsvg",
                             "stringr"),
                 load = FALSE)
  
  ensure_package("magrittr")
  
  # remove trailing slash in `parent_path`
  parent_path %<>% normalizePath()
  
  # ensure `parent_path` exists
  if ( !dir.exists(parent_path) ) dir.create(path = parent_path,
                                             recursive = TRUE)
  
  # generate sensible filename
  if ( is.null(filename) )
  {
    auto_name <- deparse(substitute(plotly_graph))
    
    filename <- dplyr::if_else(
      condition = auto_name == ".",
      true = "plotly_graph.svg",
      false = paste0(deparse(substitute(plotly_graph)), ".svg")
    )
  } else
  {
    filename %<>%
      stringr::str_replace(pattern = "([^\\.svg])$", 
                           replacement = "\\1.svg")
  }
  
  filepath <- paste0(parent_path, "/", filename)
  
  # delete old SVG file
  if ( file.exists(filepath) )
  {
    unlink(x = filepath)
  }
  
  if ( remove_title )
  {
    plotly_graph %<>%
      plotly::layout(title = "",
                     margin = list(t = 0))
  }
  
  if ( !is.null(font_family) )
  {
    plotly_graph %<>%
      plotly::layout(font = list(family = font_family))
  }
  
  # temporarily export plot to a HTML file
  tempfile <- tempfile(pattern = "plotly_temp_",
                       tmpdir = parent_path,
                       fileext = ".html")
  
  export_plotly2HTML(plotly_graph = plotly_graph,
                     filename = basename(tempfile),
                     parent_path = parent_path)
  
  on.exit(unlink(tempfile),
          add = TRUE)
  
  # get <div> ID of exported htmlwidget
  htmlwidget_id <-
    stringr::str_extract(string = readr::read_file(file = tempfile),
                         pattern = "(?<=<div id=\")htmlwidget-[^\"]+")
  
  # initialize Chrome as RSelenium driver
  selenium_driver <-
    RSelenium::rsDriver(browser = "chrome",
                        extraCapabilities = list(
                          chromeOptions = list(
                            prefs = list(
                              "profile.default_content_settings.popups" = 0L,
                              "download.prompt_for_download" = FALSE,
                              "download.default_directory" = parent_path
                            )
                          )
                        ),
                        verbose = FALSE)
  
  # navigate to temporary HTML file
  selenium_driver$client$navigate(url = paste0("file://", normalizePath(tempfile)))
  
  # download plot as SVG using the native
  # [`Plotly.downloadImage`](https://plot.ly/javascript/plotlyjs-function-reference/#plotlydownloadimage) function
  selenium_driver$client$executeScript(
    script = paste0("Plotly.downloadImage(document.getElementById('", htmlwidget_id, "'), ",
                    "{format: 'svg', width: ", width, ", height: ", height, ", filename: '",
                    tools::file_path_sans_ext(x = filename), "'});"),
    args = list(NULL)
  )
  
  # wait for SVG to be saved to disk
  Sys.sleep(time = 1)
  
  # convert to PDF
  if ( incl_PDF_copy )
  {
    rsvg::rsvg_pdf(svg = filepath,
                   file = paste0(tools::file_path_sans_ext(parent_path), ".pdf"))
  }
  
  # convert to PNG
  if ( incl_PNG_copy )
  {
    filepath_png <- paste0(tools::file_path_sans_ext(parent_path), ".png")
    
    rsvg::rsvg_png(svg = filepath,
                   file = filepath_png,
                   width = png_scaling_factor * width,
                   height = png_scaling_factor * height)
    
    if ( autocrop_png ) autocrop_png(path_to_png = filepath_png)
  }
}



export_plotly2HTML <- function(plotly_graph,
                               filename = NULL,
                               parent_path = paste0(getwd(), "/output"),
                               selfcontained = FALSE,
                               libdir = "plotly_files",
                               disable_legend_toggling = NULL,
                               add_web_font = "https://fonts.googleapis.com/css?family=Work+Sans:200,300,400,600,700")
{
  ensure_package(package = c("checkmate",
                             "dplyr",
                             "readr",
                             "stringr"),
                 load = FALSE)
  
  ensure_package("magrittr")
  
  # remove trailing slash in `parent_path`
  parent_path %<>% normalizePath()
  
  # ensure `parent_path` exists
  if ( !dir.exists(parent_path) ) dir.create(path = parent_path,
                                             recursive = TRUE)
  
  # generate sensible filename
  if ( is.null(filename) )
  {
    auto_name <- deparse(substitute(plotly_graph))
    
    filename <- dplyr::if_else(
      condition = auto_name == ".",
      true = "plotly_graph.html",
      false = paste0(deparse(substitute(plotly_graph)), ".html")
    )
  }
  
  filepath <- paste0(parent_path, "/", filename)
  
  htmlwidgets::saveWidget(
    widget = plotly_graph,
    file = filepath,
    selfcontained = selfcontained,
    libdir = libdir
  )
  
  if ( !is.null(disable_legend_toggling) )
  {
    test_char <- checkmate::check_choice(x = disable_legend_toggling,
                                         choices = "all")
    test_num <- checkmate::check_numeric(x = disable_legend_toggling,
                                         lower = 1,
                                         upper = length(plotly_graph$x$attrs),
                                         min.len = 1,
                                         max.len = length(plotly_graph$x$attrs),
                                         unique = TRUE,
                                         any.missing = FALSE,
                                         all.missing = FALSE)
    
    if ( !isTRUE(test_char) & !isTRUE(test_num) )
    {
      stop("Invalid argument provided: disable_legend_toggling\n",
           ifelse(!isTRUE(test_char) & is.character(disable_legend_toggling),
                  paste0(test_char, ". Or alternatively can also be a vector of integers >= 1 and <= number of traces."),
                  ""),
           ifelse(!isTRUE(test_num) & is.numeric(disable_legend_toggling),
                  paste0(test_num, ". Or alternatively can also be \"all\"."),
                  ""))
      
    } else if ( isTRUE(test_char) )
    {
      css_rules <-
        c("",
          "/* hides the svg dom element that has the click handler responsible for toggling */",
          ".legend .traces .legendtoggle {",
          "  display: none;",
          "}",
          "/* just for presentation: shows the default cursor instead of the text cursor */",
          ".legend .traces .legendtext {",
          "  cursor: default;",
          "}",
          "")
    } else
    {
      disable_legend_toggling %<>% as.integer()
      
      css_rules <-
        c("",
          "/* hides the svg dom element that has the click handler responsible for toggling */")
      
      for ( i in disable_legend_toggling )
      {
        css_rules %<>%
          c(paste0(".legend .groups:nth-of-type(", i, ") .legendtoggle",
                   dplyr::if_else(i == last(disable_legend_toggling),
                                  " {",
                                  ","), " "))
      }
      
      css_rules %<>%
        c("  display: none;",
          "}",
          "/* just for presentation: shows the default cursor instead of the text cursor */")
      
      for ( i in disable_legend_toggling )
      {
        css_rules %<>%
          c(paste0(".legend .groups:nth-of-type(", i, ") .legendtext",
                   dplyr::if_else(i == last(disable_legend_toggling),
                                  " {",
                                  ","), " "))
      }
      
      css_rules %<>%
        c("  cursor: default;",
          "}",
          "")
    }
    
    # write modified .css file
    plotly_dir <-
      list.dirs(path = paste0(parent_path, "/", libdir),
                full.names = TRUE,
                recursive = FALSE) %>%
      stringr::str_subset(pattern = "plotlyjs")
    
    readr::read_lines(file = paste0(plotly_dir, "/plotly-htmlwidgets.css")) %>%
      c(css_rules) %>%
      readr::write_lines(path = paste0(plotly_dir, "/plotly_htmlwidgets_custom.css"),
                         append = FALSE)
    
    # modify dependency path in HTML file
    readr::read_file(file = filepath) %>%
      stringr::str_replace(pattern = "plotly-htmlwidgets\\.css",
                           replacement = "plotly_htmlwidgets_custom.css") %>%
      readr::write_file(path = filepath,
                        append = FALSE)
  }
  
  if ( !is.null(add_web_font) )
  {
    webfont_tag <-
      "<link href=\"" %>%
      paste0(checkmate::assert_character(x = add_web_font,
                                         pattern = "^https?://\\w.*",
                                         ignore.case = TRUE,
                                         any.missing = FALSE,
                                         all.missing = FALSE,
                                         unique = TRUE)) %>%
      paste0("\" rel=\"stylesheet\" />")
    
    readr::read_file(file = filepath) %>%
      stringr::str_replace(pattern = "<link href=",
                           replacement = paste0(webfont_tag, "\n<link href=")) %>%
      readr::write_file(path = filepath,
                        append = FALSE)
  }
}

autocrop_png <- function(path_to_png)
{
  ensure_package("magrittr")
  
  imager::load.image(file = path_to_png) %>%
    imager::autocrop() %>%
    imager::pad(nPix = 4,
                axes = "xy",
                pos = 0) %>%
    imager::flatten.alpha() %>%
    imager::save.image(file = path_to_png)
}

## ## require(devtools)
## ## devtools::install_github("ropensci/RSelenium")

## # install.packages("XML")

# # usage:
# export_plotly2SVG(p, filename = "test.plotly.svg") # fantastic!!


#######################################################################
# MDS plot, iMDS plot for svg using plotly
# - metapath
# - indirpath
# - lFC, alpha
# - method
# - outdirpath
#######################################################################

meta2plyMDSplot <- function(meta.df, DESeq2.obj.disp, outdirpath=".",
                            dataname = dataname, top=500,
                            launch = TRUE,
                            DE = FALSE) {
  
  title <- paste0("PCA Plot ", dataname, if (DE) { "_DE" } else {""},  " ", time.now())
  filename <- paste0("PCA_", dataname, "_", core.name(meta.df), "_", time.now(), "_DESeq2", if (DE) { "_DE" } else {""}, collapse = "_")
  
  a1 <- glMDSPlot(if (DE) {as.data.frame(DESeq2.obj)} else {DESeq2.obj},
                  top = top,
                  path = outdirpath,
                  main = title,
                  html = filename,
                  launch = launch)
  
  df <- as.data.frame(a1$points)
  
  p <- plot_ly(df, x = df$V1, y = df$V2, text = rownames(df),
               mode = "markers", color = meta.df$condition, marker = list(size=11))
  p <- layout(p, title = title,
              xaxis = list(title = "PC 1"),
              yaxis = list(title = "PC 2"))
  
  export_plotly2SVG(p, 
                    filename = paste0(filename, "PC1_PC2.svg"), 
                    parent_path = outdirpath)
  
  p1 <- plot_ly(df, x = df$V2, y = df$V3, text = rownames(df),
                mode = "markers", color = meta.df$condition, marker = list(size=11))
  p1 <- layout(p1, title = title,
               xaxis = list(title = "PC2"),
               yaxis = list(title = "PC3"))
  
  export_plotly2SVG(p1, 
                    filename = paste0(filename, "PC2_PC3.svg"),
                    parent_path = outdirpath)
}


#######################################################################
# MDS plot, iMDS plot for svg using plotly only DE table
# - metapath
# - indirpath
# - lFC, alpha
# - method
# - outdirpath
#######################################################################

meta2plyMDSplotDE <- function(cntsDEfpath, outdirpath=".",
                              dataname = dataname, corename = corename, top=500,
                              launch = TRUE) {

  title <- paste0("PCA Plot ", dataname, "_DE ", time.now())
  filename <- paste0("PCA_DE_", dataname, "_", corename, "_", time.now(), "_DESeq2", collapse = "_")
  
  countsDE <- xlsx2df.list(cntsDEfpath)[[1]]

  a1 <- glMDSPlot(countsDE,
                  top = top,
                  path = outdirpath,
                  main = title,
                  groups = unique(colnames(countsDE)),
                  html = filename,
                  launch = launch)
  
  df <- as.data.frame(a1$points)
  
  p <- plot_ly(df, x = df$V1, y = df$V2, text = rownames(df),
               mode = "markers", color = meta.df$condition, marker = list(size=11))
  p <- layout(p, title = title,
              xaxis = list(title = "PC 1"),
              yaxis = list(title = "PC 2"))
  
  export_plotly2SVG(p, 
                    filename = paste0(filename, "PC1_PC2.svg"), 
                    parent_path = outdirpath)
  
  p1 <- plot_ly(df, x = df$V2, y = df$V3, text = rownames(df),
                mode = "markers", color = meta.df$condition, marker = list(size=11))
  p1 <- layout(p1, title = title,
               xaxis = list(title = "PC2"),
               yaxis = list(title = "PC3"))
  
  export_plotly2SVG(p1, 
                    filename = paste0(filename, "PC2_PC3.svg"),
                    parent_path = outdirpath)
}


#######################################################################
# PCA data frame to plotly svg
#######################################################################

df2plotlyPCA <- function(df, colorPalette, colorCol, title, outdir, outfname, pdfp = TRUE, pngp = TRUE) {
  dir.create(outdir, recursive=TRUE)
  
  p <- plot_ly(df, x = df[, 1], y = df[, 2], text = rownames(df),
               mode = "markers", color = colorPalette[df[, colorCol]],
               marker = list(size=11),
               incl_PDF_copy = pdfp,
               incl_PNG_copy = pngp)
  p <- layout(p, title = title,
              xaxis = list(title = "PC1"),
              yaxis = list(title = "PC2"))
  
  export_plotly2SVG(p, 
                    filename = paste0(outfname, "_PC1_PC2.svg"), 
                    parent_path = outdirpath,
                    incl_PDF_copy = pdfp,
                    incl_PNG_copy = pngp)
  
  #   p1 <- plot_ly(df, x = df$V2, y = df$V3, text = rownames(df),
  #                 mode = "markers", color = colorPalette[df[, colorCol]],
  #                 marker = list(size=11))
  #   p1 <- layout(p1, title = title,
  #                xaxis = list(title = "PC2"),
  #                yaxis = list(title = "PC3"))
  #   
  #   export_plotly2SVG(p1, 
  #                     filename = paste0(filename, "_PC2_PC3.svg"),
  #                     parent_path = file.path(outdir, "output"),
  #                     incl_PDF_copy = pdfp,
  #                     incl_PNG_copy = pngp)
  
}


# create color palette `length(unique(df$cluster))` should be simply the number of clusters
# require(RColorBrewer)
# col.palette <- colorRampPalette(brewer.pal(n=12,name="Set3"))(length(unique(df$cluster)))

# call the function:
# df2plotlyPCA(df, col.palette, "cluster", title="PCA", outdir="test/PCA", outfname = "test")



###########################
# GO analysis (stand 20180917)
###########################


#######################################
# load packages
#######################################

# source("~/Dropbox/R/central-scripts/io.xlsx.R")
require(clusterProfiler)
require(org.Mm.eg.db)
keytypes(org.Mm.eg.db)


# [1] "ACCNUM"       "ALIAS"        "ENSEMBL"     
# [4] "ENSEMBLPROT"  "ENSEMBLTRANS" "ENTREZID"    
# [7] "ENZYME"       "EVIDENCE"     "EVIDENCEALL" 
# [10] "GENENAME"     "GO"           "GOALL"       
# [13] "IPI"          "MGI"          "ONTOLOGY"    
# [16] "ONTOLOGYALL"  "PATH"         "PFAM"        
# [19] "PMID"         "PROSITE"      "REFSEQ"      
# [22] "SYMBOL"       "UNIGENE"      "UNIPROT"  



#######################################
# ALIAS to ENTREZID data frame rownames
#######################################

alias2entrezid.df <- function(DE.df) {
  # converts rownames from ALIAS to ENTREZID
  genes <- rownames(DE.df)
  genes2eg <- bitr(genes,
                   fromType = "ALIAS",
                   toType   = "ENTREZID",
                   OrgDb    = "org.Mm.eg.db")
  
  # remove duplicated rows of aliases
  genes2eg.unique <- genes2eg[!duplicated(genes2eg$ALIAS), ]
  
  # remove duplicated rows of entrezids
  genes2eg.unique <- genes2eg[!duplicated(genes2eg$ENTREZID), ]
  df.new <- DE.df[genes2eg.unique$ALIAS, ]
  
  # again remove rows with duplicated rownames
  df.new <- df.new[!duplicated(rownames(df.new)), ]
  rownames(df.new) <- genes2eg.unique$ENTREZID
  
  # save/hide not found genes in the object # sorted looks nicer!
  not.found.genes <- sort(setdiff(genes, genes2eg.unique$ALIAS))
  attr(df.new, "not.found.genes") <- not.found.genes
  
  # save found genes in the object
  found.genes <- genes2eg.unique$ALIAS # don't sort to match with df.new!
  attr(df.new, "found.genes") <- found.genes
  
  df.new
}


#######################################
# found and not founds
#######################################

dfListNotFounds <- function(eg.df.list) {
  lapply(eg.df.list, function(df) data.frame(alias = attr(df, "not.found.genes")))
}

dfListFounds    <- function(eg.df.list) {
  names.list <- lapply(eg.df.list, function(df) attr(df, "found.genes"))
  Map(f = function(df, names.vec) {df$alias <- names.vec; df}, eg.df.list, names.list)
}




#######################################
# GO overrepresentation
# simplify is discussed here:
#######################################
# https://github.com/GuangchuangYu/clusterProfiler/issues/28


dfEnrichGO <- function(eg.df, ont) {
  res <- enrichGO(gene          = rownames(eg.df),
                  OrgDb         = org.Mm.eg.db,
                  keyType       = "ENTREZID",
                  ont           = ont,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  res <- simplify(res,
                  cutoff = 0.7,
                  by = "p.adjust",
                  select_fun = min)
  res <- setReadable(res, org.Mm.eg.db, keytype = "ENTREZID")
  res
}

dfListEnrichGO <- function(eg.df.list, ont) {
  lapply(eg.df.list, function(df) dfEnrichGO(df, ont))
}



#######################################
# GO GSEA
#######################################
# prepare own geneList
# https://github.com/GuangchuangYu/DOSE/wiki/how-to-prepare-your-own-geneList


####################################################
# to make possible `simplify` on gseaResult object
####################################################

require(magrittr)
setMethod("simplify", signature(x="gseaResult"),
          function(x, cutoff=0.7, by="p.adjust", select_fun=min, measure="Wang", semData=NULL) {
            if (!x@setType %in% c("BP", "MF", "CC"))
              stop("simplify only applied to output from enrichGO...")
            x@result %<>% simplify_internal(., cutoff, by, select_fun,
                                            measure, x@setType, semData)
            return(x)
          }
)

# https://github.com/GuangchuangYu/clusterProfiler/blob/master/R/simplify.R
simplify_internal <- function(res, cutoff=0.7, by="p.adjust", select_fun=min, measure="Rel", ontology, semData) {
  if (missing(semData) || is.null(semData)) {
    if (measure == "Wang") {
      semData <- GOSemSim::godata(ont = ontology)
    } else {
      stop("godata should be provided for IC-based methods...")
    }
  } else {
    if (ontology != semData@ont) {
      msg <- paste("semData is for", semData@ont, "ontology, while enrichment result is for", ontology)
      stop(msg)
    }
  }
  
  sim <- GOSemSim::mgoSim(res$ID, res$ID,
                          semData = semData,
                          measure=measure,
                          combine=NULL)
  
  ## to satisfy codetools for calling gather
  go1 <- go2 <- similarity <- NULL
  
  sim.df <- as.data.frame(sim)
  sim.df$go1 <- row.names(sim.df)
  sim.df <- tidyr::gather(sim.df, go2, similarity, -go1)
  
  sim.df <- sim.df[!is.na(sim.df$similarity),]
  
  ## feature 'by' is attached to 'go1'
  sim.df <- merge(sim.df, res[, c("ID", by)], by.x="go1", by.y="ID")
  sim.df$go2 <- as.character(sim.df$go2)
  
  ID <- res$ID
  
  GO_to_remove <- character()
  for (i in seq_along(ID)) {
    ii <- which(sim.df$go2 == ID[i] & sim.df$similarity > cutoff)
    ## if length(ii) == 1, then go1 == go2
    if (length(ii) < 2)
      next
    
    sim_subset <- sim.df[ii,]
    
    jj <- which(sim_subset[, by] == select_fun(sim_subset[, by]))
    
    ## sim.df <- sim.df[-ii[-jj]]
    GO_to_remove <- c(GO_to_remove, sim_subset$go1[-jj]) %>% unique
  }
  
  res[!res$ID %in% GO_to_remove, ]
}

#############################################


#############################################
# GO GSEA
#############################################

dfGseGO <- function(eg.df, ont) {
  genes <- eg.df$log2FoldChange
  names(genes) <- rownames(eg.df)
  genes <- sort(genes, decreasing = TRUE) # ensure decreasing order
  res <- gseGO(gene          = genes,
               OrgDb         = org.Mm.eg.db,
               keyType       = "ENTREZID",
               ont           = ont,
               nPerm         = 1000,
               minGSSize     = 10,
               maxGSSize     = 1000,
               pAdjustMethod = "BH",
               pvalueCutoff  = 0.05)
  res <- simplify(res,              # maybe I should write that code! Did! (see above!) works!
                  cutoff = 0.7,
                  by = "p.adjust",
                  select_fun = min)
  res <- setReadable(res, org.Mm.eg.db, keytype = "ENTREZID")
  res
}

dfListGseGO <- function(eg.df.list, ont) {
  lapply(eg.df.list, function(df) dfGseGO(df, ont))
}






#######################################################
# GO BP,CC,MF entire analysis with or without selection
# CNTS removed, since not used in GO
#######################################################

do_GO_enrichment <- function(fpathDExlsx,  # fpathCNTSxlsx, 
                             outBase, 
                             ont,
                             fpathCLUSTERxlsx = "",
                             select_clusters = c(),
                             DE.list.mode.p=T, # take only c(1, 4, 7) of xlsx sheets
                             gsea.p=T) { # do gsea? log2FC needed!
  
  fname <- basename(fpathDExlsx)
  
  if (length(select_clusters > 0)) {
    clust.numbers <- paste0("_", paste(unique(select_clusters), collapse="-"), "_")
  } else {
    clust.numbers <- "_"
  }
  
  decorateFname <- function(textPiece) gsub(".xlsx", 
                                            paste0(clust.numbers, textPiece, ont, ".xlsx"), 
                                            fname)
  overrepFname_ont <- decorateFname("over_GO_")
  gseaFname_ont <- decorateFname("gsea_GO_")
  notFoundFnameDE <- decorateFname("not_found_")
  FoundFnameDE <- decorateFname("found_")
  
  if (!dir.exists(outBase)) {
    dir.create(outBase, recursive = TRUE)
  }
  
  # read-in data
  if (DE.list.mode.p) {
    xlsxDE_dfList   <- xlsx2df.list(fpathDExlsx)[c(1, 4, 7)]
  } else {
    xlsxDE_dfList   <- xlsx2df.list(fpathDExlsx)
  }
  
  if (fpathCLUSTERxlsx != "" && length(select_clusters) > 0) {
    # read-in cluster tables
    cluster.dfs <- xlsx2df.list(fpathCLUSTERxlsx)
    
    # extract gene names from clusters
    clusternames <- lapply(cluster.dfs, rownames)
    
    # extract genes from the clusters
    select_cluster_names <- unique(unlist(clusternames[select_clusters]))
    xlsxDE_dfList <- lapply(xlsxDE_dfList, function(df) df[intersect(rownames(df), select_cluster_names), ])
    
  }
  
  # translate to eg
  xlsxDE_dfListEg <- lapply(xlsxDE_dfList, alias2entrezid.df)
  
  # not founds
  not_founds_list_DE       <- dfListNotFounds(xlsxDE_dfListEg)
  
  # founds
  founds_list_DE           <- dfListFounds(xlsxDE_dfListEg)
  
  # do overrep ont
  xlsxDE_dfListEg_GO_over_ont <- dfListEnrichGO(xlsxDE_dfListEg, ont) 
  
  # do gsea ont
  if (gsea.p) {
    xlsxDE_dfListEg_GO_gsea_ont <- dfListGseGO(xlsxDE_dfListEg, ont) 
  }
  
  # write results
  write.dfs(not_founds_list_DE, file.path(outBase, notFoundFnameDE))
  write.dfs(founds_list_DE, file.path(outBase, FoundFnameDE))
  write.dfs(xlsxDE_dfListEg_GO_over_ont, file.path(outBase, overrepFname_ont))
  if (gsea.p) {
  write.dfs(xlsxDE_dfListEg_GO_gsea_ont, file.path(outBase, gseaFname_ont))
  }
}

do_all_GO_enrichment <- function(fpathDExlsx, 
                                 outBase,
                                 fpathCLUSTERxlsx = "",
                                 select_clusters = c(),
                                 DE.list.mode.p=T,
                                 gsea.p=T)  {
  do_GO_enrichment(fpathDExlsx, outBase, "BP", fpathCLUSTERxlsx, select_clusters, DE.list.mode.p, gsea.p)
  do_GO_enrichment(fpathDExlsx, outBase, "CC", fpathCLUSTERxlsx, select_clusters, DE.list.mode.p, gsea.p)
  do_GO_enrichment(fpathDExlsx, outBase, "MF", fpathCLUSTERxlsx, select_clusters, DE.list.mode.p, gsea.p)
}



#########################
# KEGG
#########################



#######################################
# KEGG overrepresentation 
#######################################


dfEnrichKEGG <- function(eg.df) {
  res <- enrichKEGG(gene          = rownames(eg.df),
                    keyType       = "kegg",
                    pAdjustMethod = 'BH',
                    organism      = 'mmu',
                    pvalueCutoff  =  0.05,
                    qvalueCutoff  =  0.05)
  res <- setReadable(res, org.Mm.eg.db, keytype = "ENTREZID")
  res
}

dfListEnrichKEGG <- function(eg.df.list) {
  lapply(eg.df.list, function(df) dfEnrichKEGG(df))
}

#######################################
# KEGG GSEA
#######################################


dfGseKEGG <- function(eg.df, pCutOff = 0.05, org = 'mmu') {
  genes <- eg.df$log2FoldChange
  names(genes) <- rownames(eg.df)
  genes <- sort(genes, decreasing = TRUE)
  res <- gseKEGG(gene          = genes,
                 keyType       = "kegg",
                 nPerm         = 1000,
                 minGSSize     = 10,
                 organism      = org,
                 pvalueCutoff  = pCutOff)
  res <- setReadable(res, org.Mm.eg.db, keytype = "ENTREZID")
  res
}

dfListGseKEGG <- function(eg.df.list) {
  lapply(eg.df.list, function(df) dfGseKEGG(df))
}


#############################
# MKEGG
#############################


#######################################
# KEGG Module overrepresentation test
#######################################


dfEnrichMKEGG <- function(eg.df, pCutOff = 0.05, qCutOff = 0.05, org = 'mmu') {
  res <- enrichMKEGG(gene          = rownames(eg.df),
                     keyType       = "kegg",
                     pAdjustMethod = 'BH',
                     organism      = org,
                     pvalueCutoff  = pCutOff,
                     qvalueCutoff  = qCutOff)
  res <- setReadable(res, org.Mm.eg.db, keytype = "ENTREZID")
  res
}

dfListEnrichMKEGG <- function(eg.df.list) {
  lapply(eg.df.list, function(df) dfEnrichMKEGG(df))
}

#######################################
# KEGG Module GSEA test
#######################################


dfGseMKEGG <- function(eg.df, pCutOff = 0.05, org = 'mmu') {
  genes <- eg.df$log2FoldChange
  names(genes) <- rownames(eg.df)
  genes <- sort(genes, decreasing = TRUE)
  res <- gseMKEGG(gene          = genes,
                  keyType       = "kegg",
                  pAdjustMethod = 'BH',
                  organism      = org,
                  pvalueCutoff  = pCutOff)
  res <- setReadable(res, org.Mm.eg.db, keytype = "ENTREZID")
  res
}

dfListGseMKEGG <- function(eg.df.list, pCutOff = 0.05, org = 'mmu') {
  lapply(eg.df.list, function(df) dfGseMKEGG(df, pCutOff = pCutOff, org = org))
}

do_KEGG_enrichment <- function(fpathDExlsx,
                               outBase, 
                               fpathCLUSTERxlsx = "",
                               select_clusters = c(),
                               DE.list.mode.p=T, # take only c(1, 4, 7) of xlsx sheets
                               gsea.p=T) { #gsea? log2FC column needed!
  
  fname <- basename(fpathDExlsx)
  
  if (length(select_clusters > 0)) {
    clustNumbers <- paste0("_", paste(unique(select_clusters), collapse="-"), "_")
  } else {
    clustNumbers <- "_"
  }
  
  if (!dir.exists(outBase)) {
    dir.create(outBase, recursive = TRUE)
  }
  
  decorateFname <- function(textPiece, clustNumbers_ = clustNumbers, fname_ = fname) {
    gsub(".xlsx", paste0(clustNumbers_, textPiece, ".xlsx"), fname_)
  }
  overKEGGFname <- decorateFname("over_KEGG")
  gseaKEGGFname <- decorateFname("gsea_KEGG")
  overMKEGGFname <- decorateFname("over_MKEGG")
  gseaMKEGGFname <- decorateFname("gsea_MKEGG")
  notFoundFnameDE <- decorateFname("not_found_")
  FoundFnameDE <- decorateFname("found_")
  
  # read-in data
  if (DE.list.mode.p) {
    xlsxDE_dfList   <- xlsx2df.list(fpathDExlsx)[c(1, 4, 7)]
  } else {
    xlsxDE_dfList   <- xlsx2df.list(fpathDExlsx)
  }
  
  if (fpathCLUSTERxlsx != "" && length(select_clusters) > 0) {
    # read-in cluster tables
    cluster.dfs <- xlsx2df.list(fpathCLUSTERxlsx)
    
    # extract gene names from clusters
    clusternames <- lapply(cluster.dfs, rownames)
    
    # extract genes from the clusters
    select_cluster_names <- unique(unlist(clusternames[select_clusters]))
    xlsxDE_dfList <- lapply(xlsxDE_dfList, function(df) df[intersect(rownames(df), select_cluster_names), ])
    
  }
  
  # translate to eg
  xlsxDE_dfListEg <- lapply(xlsxDE_dfList, alias2entrezid.df)
  
  # not founds
  not_founds_list_DE       <- dfListNotFounds(xlsxDE_dfListEg)
  
  # founds
  founds_list_DE           <- dfListFounds(xlsxDE_dfListEg)
  
  # do overrep KEGG
  xlsxDE_dfListEg_er_KEGG <- dfListEnrichKEGG(xlsxDE_dfListEg) 
  
  # do gsea KEGG
  if (gsea.p) {
    xlsxDE_dfListEg_gs_KEGG <- dfListGseKEGG(xlsxDE_dfListEg)
  }
  
  # do overrep MKEGG
  xlsxDE_dfListEg_er_MKEGG <- dfListEnrichMKEGG(xlsxDE_dfListEg) 
  
  # do gsea MKEGG
  if (gsea.p) {
    xlsxDE_dfListEg_gs_MKEGG <- dfListGseMKEGG(xlsxDE_dfListEg)
  }
  
  # write results
  write.dfs(not_founds_list_DE, file.path(outBase, notFoundFnameDE))
  write.dfs(founds_list_DE, file.path(outBase, FoundFnameDE))
  write.dfs(xlsxDE_dfListEg_er_KEGG, file.path(outBase, overKEGGFname))
  write.dfs(xlsxDE_dfListEg_er_MKEGG, file.path(outBase, overMKEGGFname))
  if (gsea.p) {
    write.dfs(xlsxDE_dfListEg_gs_KEGG, file.path(outBase, gseaKEGGFname))
    write.dfs(xlsxDE_dfListEg_gs_MKEGG, file.path(outBase, gseaMKEGGFname))
  }
}




################################
# gskb enrichment
################################

require(gskb)

# load gskb data which were prepared previously
all_ez_grouped <- readRDS(file = "/media/josephus/archive/gskb/all/mGSKB_Entrez_cl_grouped.list.gmt.RDS")

names(all_ez_grouped)
# [1] "MousePath_Co-expression_eg.gmt" "MousePath_GO_eg.gmt"            "MousePath_Location_eg.gmt"      "MousePath_Metabolic_eg.gmt"     "MousePath_miRNA_eg.gmt"        
# [6] "MousePath_Other_eg.gmt"         "MousePath_Pathway_eg.gmt"       "MousePath_TF_eg.gmt"  



#######################################
# MSigDB/GSKB overrepresentation test
#######################################

dfEnrichDB <- function(df, gmt, pCutOff = 0.05, pAdj = "BH") {
  if (typeof(gmt) == "character") {
    term2gene <- read.gmt(gmt)
  } else {
    term2gene <- gmt
  }
  
  res <- enricher(rownames(df), TERM2GENE=term2gene, pvalueCutoff = pCutOff, pAdjustMethod = pAdj)
  res <- setReadable(res, org.Mm.eg.db, keytype = "ENTREZID")
  res
}

dfListEnrichDB <- function(df.list, gmt, pCutOff = 0.05, pAdj = "BH") {
  if (typeof(gmt) == "character") {
    term2gene <- read.gmt(gmt)
  } else {
    term2gene <- gmt
  }
  lapply(df.list, function(df) dfEnrichDB(df, gmt, pCutOff = pCutOff, pAdj = pAdj))
}

#######################################
# MSigDB/GSKB GSEA
#######################################

dfGseDB <- function(df, gmt, pCutOff = 0.05, pAdj = "BH") {
  if (typeof(gmt) == "character") {
    term2gene <- read.gmt(gmt)
  } else {
    term2gene <- gmt
  }
  genes <- df$log2FoldChange
  names(genes) <- rownames(df)
  genes <- sort(genes, decreasing = TRUE)
  res <- GSEA(genes, 
              TERM2GENE=term2gene, 
              nPerm = 1000,
              minGSSize = 10,
              maxGSSize = 1000,
              pvalueCutoff = pCutOff, 
              pAdjustMethod = pAdj)
  res <- setReadable(res, org.Mm.eg.db, keytype = "ENTREZID")
  res
}

dfListGseDB <- function(df.list, gmt, pCutOff = 0.05, pAdj = "BH") {
  if (typeof(gmt) == "character") {
    term2gene <- read.gmt(gmt)
  } else {
    term2gene <- gmt
  }
  lapply(df.list, function(df) dfGseDB(df, gmt, pCutOff = pCutOff, pAdj = pAdj))
}



#######################
# gskb general analysis
#######################

do_gskb_enrichment <- function(fpathDExlsx,
                               gmt,
                               gmtname,
                               outBase, 
                               fpathCLUSTERxlsx = "",
                               select_clusters = c(),
                               DE.list.mode.p=T, # take only c(1, 4, 7) of xlsx sheets
                               gsea.p=T) { #gsea? log2FC needed!
  
  fname <- basename(fpathDExlsx)
  outBase <- file.path(outBase, gmtname)
  
  if (length(select_clusters > 0)) {
    clustNumbers <- paste0("_", paste(unique(select_clusters), collapse="-"), "_")
  } else {
    clustNumbers <- "_"
  }
  
  if (!dir.exists(outBase)) {
    dir.create(outBase, recursive = TRUE)
  }
  
  decorateFname <- function(textPiece, gmtname, clustNumbers_ = clustNumbers, fname_ = fname) {
    gsub(".xlsx", paste0(clustNumbers_, textPiece, ".xlsx"), fname_)
  }
  overDBFname <- decorateFname(paste0("over_gskb_", gmtname))
  gseaDBFname <- decorateFname(paste0("gsea_gskb_", gmtname))
  notFoundFnameDE <- decorateFname(paste0("not_found_", gmtname))
  FoundFnameDE <- decorateFname(paste0("found_", gmtname))
  
  # read-in data
  if (DE.list.mode.p) {
    xlsxDE_dfList   <- xlsx2df.list(fpathDExlsx)[c(1, 4, 7)]
  } else {
    xlsxDE_dfList   <- xlsx2df.list(fpathDExlsx)
  }
  
  if (fpathCLUSTERxlsx != "" && length(select_clusters) > 0) {
    # read-in cluster tables
    cluster.dfs <- xlsx2df.list(fpathCLUSTERxlsx)
    
    # extract gene names from clusters
    clusternames <- lapply(cluster.dfs, rownames)
    
    # extract genes from the clusters
    select_cluster_names <- unique(unlist(clusternames[select_clusters]))
    xlsxDE_dfList <- lapply(xlsxDE_dfList, function(df) df[intersect(rownames(df), select_cluster_names), ])
    
  }
  
  # translate to eg
  xlsxDE_dfListEg <- lapply(xlsxDE_dfList, alias2entrezid.df)
  
  # not founds
  not_founds_list_DE       <- dfListNotFounds(xlsxDE_dfListEg)
  
  # founds
  founds_list_DE           <- dfListFounds(xlsxDE_dfListEg)
  
  # do overrep gskb
  xlsxDE_dfListEg_er_DB <- dfListEnrichDB(xlsxDE_dfListEg, gmt = gmt) 
  
  # do gsea gskb
  if (gsea.p) {
    xlsxDE_dfListEg_gs_DB <- dfListGseDB(xlsxDE_dfListEg, gmt = gmt)
  }
  
  # write results
  write.dfs(not_founds_list_DE, file.path(outBase, notFoundFnameDE))
  write.dfs(founds_list_DE, file.path(outBase, FoundFnameDE))
  write.dfs(xlsxDE_dfListEg_er_DB, file.path(outBase, overDBFname))
  if (gsea.p) {
    write.dfs(xlsxDE_dfListEg_gs_DB, file.path(outBase, gseaDBFname))
  }
}



do_all_gskb_enrichments <- function(DEfpath, outBase, all_ez_grouped, xlsx.path = "", select_group = c(), DE.list.mode.p=T, gsea.p=T) {
  do_gskb_enrichment(DEfpath, 
                     all_ez_grouped$MousePath_Pathway_eg.gmt,
                     "mpathway",
                     outBase,
                     fpathCLUSTERxlsx = xlsx.path,
                     select_clusters = select_group,
                     DE.list.mode.p = DE.list.mode.p,
                     gsea.p = gsea.p)
  
  do_gskb_enrichment(DEfpath,
                     all_ez_grouped$MousePath_Metabolic_eg.gmt,
                     "mmetabolic",
                     outBase,
                     fpathCLUSTERxlsx = xlsx.path,
                     select_clusters = select_group,
                     DE.list.mode.p = DE.list.mode.p,
                     gsea.p = gsea.p)
  
  do_gskb_enrichment(DEfpath,
                     all_ez_grouped$MousePath_TF_eg.gmt,
                     "mTF",
                     outBase,
                     fpathCLUSTERxlsx = xlsx.path,
                     select_clusters = select_group,
                     DE.list.mode.p = DE.list.mode.p,
                     gsea.p = gsea.p)
  
  do_gskb_enrichment(DEfpath,
                     all_ez_grouped$MousePath_GO_eg.gmt,
                     "mGO",
                     outBase,
                     fpathCLUSTERxlsx = xlsx.path,
                     select_clusters = select_group,
                     DE.list.mode.p = DE.list.mode.p,
                     gsea.p = gsea.p)
  
  do_gskb_enrichment(DEfpath,
                     all_ez_grouped$MousePath_Other_eg.gmt,
                     "mOther",
                     outBase,
                     fpathCLUSTERxlsx = xlsx.path,
                     select_clusters = select_group,
                     DE.list.mode.p = DE.list.mode.p,
                     gsea.p = gsea.p)
  
  do_gskb_enrichment(DEfpath,
                     all_ez_grouped$MousePath_miRNA_eg.gmt,
                     "mmiRNA",
                     outBase,
                     fpathCLUSTERxlsx = xlsx.path,
                     select_clusters = select_group,
                     DE.list.mode.p = DE.list.mode.p,
                     gsea.p = gsea.p)
  
  do_gskb_enrichment(DEfpath,
                     all_ez_grouped$`MousePath_Co-expression_eg.gmt`,
                     "mCoexpr",
                     outBase,
                     fpathCLUSTERxlsx = xlsx.path,
                     select_clusters = select_group,
                     DE.list.mode.p = DE.list.mode.p,
                     gsea.p = gsea.p)
  
  do_gskb_enrichment(DEfpath,
                     all_ez_grouped$`MousePath_Location_eg.gmt`,
                     "mlocation",
                     outBase,
                     fpathCLUSTERxlsx = xlsx.path,
                     select_clusters = select_group,
                     DE.list.mode.p = DE.list.mode.p,
                     gsea.p = gsea.p)
}

