
rm(list=ls())




######################################################################
# helper functions
######################################################################


####################################################################
# binquire() and inquire()
####################################################################

binquire <- function(packagename) {
  pckgname <- toString(substitute(packagename))
  eval(substitute(
    if (!require(pckgname)) {
      source("https://bioconductor.org/biocLite.R")
      biocLite(pckgname)
      require(pckgname)
    }
  ))
}

inquire <- function(packagename) {
  pckgname <- toString(substitute(packagename))
  eval(substitute(
    if (!require(pckgname)) {
      install.packages(pckgname)
      require(pckgname)
    }
  ))
}

gitinquire <- function(gitrepo) {
  gitrep <- toString(substitute(gitrepo))
  eval(substitute(
    if (!require(gitrep)) {
      inquire(devtools)
      devtools::install_github(gitrep)
      require(gitrep)
    }
  ))
}
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
inquire(xlsx)        # xlsx printing
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

df.sortings.limma <- function(DE.res.df) {
  list(l2FC.srt = DE.res.df[order(DE.res.df$logFC, decreasing = TRUE), ],
       p.srt = DE.res.df[order(DE.res.df$adj.P.Val, decreasing = FALSE), ],
       names.srt = DE.res.df[order(rownames(DE.res.df), decreasing = FALSE), ],
       B.srt = DE.res.df[order(DE.res.df$B, decreasing = TRUE), ])
}

#####################################################################
# print DE sortings
#####################################################################

print.DE.sortings <- function(DE.list, fname, dir = ".") {
  DE.list.with.sortings <- lapply(once.flatten(lapply(DE.list, df.sortings)), as.data.frame)
  write.dfs(DE.list.with.sortings, file.path(dir, fname))
}

print.DE.sortings.limma <- function(DE.list, fname, dir = ".") {
  DE.list.with.sortings <- lapply(once.flatten(lapply(DE.list, df.sortings.limma)), as.data.frame)
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

print.cnts.DE.sortings.limma <- function(cnts, DE.list, fname, dir = ".") {
  DE.list.with.sortings <- lapply(once.flatten(lapply(DE.list, df.sortings.limma)),
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

df2dflist <- function(df, col.selector) {
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
# DE analysis with limma
#######################################################################
require(edgeR)

DEanalysis.limma <- function(meta.df, DESeq2.obj,
                             outdirpath=".",
                             dataname = "",
                             printp = FALSE, prior.count = 3,
                             alpha = 0.05, lFC = 1, filterp = FALSE) {
  cnts <- counts(DESeq2.obj, normalized=FALSE)
  dge <- DGEList(counts = cnts)
  design <- model.matrix(~ 0 + condition, data = meta.df)
  colnames(design) <- levels(meta.df$condition)
  
  # remove rows with zero or very low counts
  keep <- filterByExpr(dge, design)
  dge <- dge[keep, , keep.lib.sizes=FALSE]
  
  # apply scale normalization to RNA-seq read counts and TMM normalization method
  dge <- calcNormFactors(dge)
  
  # DE limma-trend
  logCPM <- cpm(dge, log = TRUE, prior.count = prior.count) # damp down variances of logs of low counts
  
  fit <- lmFit(logCPM, design)
  fit <- eBayes(fit, trend = TRUE)
  res <- topTable(fit, coef = ncol(design), number = Inf)
  
  if (filterp) {
    res <- subset(subset(res, adj.P.Val < alpha), abs(logFC) > lFC)
  }
  
  up <- subset(res, res$logFC > 0)
  down <- subset(res, res$logFC < 0)
  res <- list(all=res, up=up, down=down)
  
  if (printp) {
    filecorename <- paste0("DE_", ifelse(filterp, "sig_", ""),
                           paste0(dataname, "-limma"), "_", core.name(meta.df), "_", time.now())
    if (filterp) { filecorename <- paste0(filecorename, "_", alpha, "_", lFC) }
    filename <- paste0(filecorename, ".xlsx")
    print.DE.sortings.limma(res, fname = filename, dir = DE.outdir)
  }
  res
}


DEanalysis.limma.corrected <- function(meta.df, DESeq2.obj,
                             outdirpath=".",
                             dataname = "",
                             printp = FALSE, prior.count = 3,
                             alpha = 0.05, lFC = 1, filterp = FALSE) {
  cnts <- counts(DESeq2.obj, normalized=FALSE)
  dge <- DGEList(counts = cnts)
  design <- model.matrix(~ 0 + condition, data = meta.df)
  colnames(design) <- levels(meta.df$condition)
  
  contr.matrix <- makeContrasts(
    k2vswtn = k2 - wtn,
    levels = colnames(design))
  # remove heteroscedasticity from count data
  
  
  # remove rows with zero or very low counts
  keep <- filterByExpr(dge, design)
  dge <- dge[keep, , keep.lib.sizes=FALSE]
  
  # apply scale normalization to RNA-seq read counts and TMM normalization method
  dge <- calcNormFactors(dge, method = "TMM")
  
  v <- voom(dge, design)
  
  # fit linear models for comparisons of interest
  vfit <- lmFit(v, design)
  vfit <- contrasts.fit(vfit, contrasts = contr.matrix)
  efit <- eBayes(vfit, trend = TRUE)
  res <- topTable(efit, coef = ncol(1), number = Inf)
  
  if (filterp) {
    res <- subset(subset(res, adj.P.Val < alpha), abs(logFC) > lFC)
  }
  
  up <- subset(res, res$logFC > 0)
  down <- subset(res, res$logFC < 0)
  res <- list(all=res, up=up, down=down)
  
  if (printp) {
    filecorename <- paste0("DE_", ifelse(filterp, "sig_", ""),
                           paste0(dataname, "-limma"), "_", core.name(meta.df), "_", time.now())
    if (filterp) { filecorename <- paste0(filecorename, "_", alpha, "_", lFC) }
    filename <- paste0(filecorename, ".xlsx")
    print.DE.sortings.limma(res, fname = filename, dir = DE.outdir)
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
                            clust.centroid, scaledata, kClusters)
  
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
  colfunc <- colorRampPalette(c("black", "red"))
  eval(parse(text = paste0("scores <- c(", paste(paste("score", 1:k, sep = ''),
                                                 collapse = ", "), ")")))
  # add cluster number and score for sorting the data
  scaledata.k <-cbind(scaledata,
                      cluster = kClust$cluster,
                      score = scores[rownames(scaledata)])
  scaledata.k <- scaledata.k[order(scaledata.k[, "cluster"], scaledata.k[, "score"]), ]
  
  # outpath
  outname <- file.path(outdirpath, paste0("k", k, "_khmap_", dataname, "_", time.now(), "_UO_", alpha, "_", lFC, name.add))
  
  # unordered clusters
  {
    # setEPS()
    # postscript(paste0(outname, ".eps"))
    svg(paste0(outname, ".svg"))
    pheatmap(scaledata.k[, 1:length(unique(meta.df$condition))],
             cluster_rows = F,
             cluster_cols = F,
             cellwidth = 40,
             col = colfunc(100),
             fontsize_row = 0.5,
             border_color = NA
    )
    dev.off()
  }
  
  {
    setEPS()
    postscript(paste0(outname, ".eps"))
    # svg(paste0(outname, ".svg"))
    pheatmap(scaledata.k[, 1:length(unique(meta.df$condition))],
             cluster_rows = F,
             cluster_cols = F,
             cellwidth = 40,
             col = colfunc(100),
             fontsize_row = 0.5,
             border_color = NA
    )
    dev.off()
  }
  
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
  colfunc <- colorRampPalette(c("black", "red"))
  df <- read.xlsx(UO.heatmap.path, sheetName=sheet, row.names=1, stringsAsFactors=FALSE)
  dfl <- df2dflist(df, "cluster")
  # within cluster sort by score
  dflsc <- lapply(dfl, function(df) df[order(df[, "score"]), ])
  dflo <- list.order.manually.by.vec(dflsc, order.vec)
  k <- length(unique(df$cluster))
  names(dflo) <- 1:k
  dflo.corr <- repeated.values.into.df.list(dflo, "cluster", 1:k)
  scaledata.k.ordered <- dflist2df(dflo.corr)
  if (length(gaps.after.blocks) > 0) {
    gaps_idxs <- cumsum(unlist(lapply(df2dflist(scaledata.k.ordered, "cluster"),
                                      function(df) dim(df)[1])))
    gaps_idxs <- gaps_idxs[gaps.after.blocks]
  }
  outname <- file.path(outdirpath, paste0("k", k, "_khmap_", dataname, "_", time.now(), "_O_", alpha, "_", lFC))
  outfname <- paste0(outname, ".xlsx")
  {
    setEPS()
    postscript(paste0(outname, ".eps"))
    # svg(paste0(outname, ".svg"))
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






####################################################################
# DE settings
####################################################################
alpha <- 0.05       # padj BH limit
FDR   <- 0.05       # false discovery rate, q-value
lFC   <- 1       # log2FC limit
GO.alpha <- 0.01
GO.FDR   <- 0.05
GO.N     <- 30
GO.cl.N  <- 10



# ####################################################################
# # These are run at 0808
# ####################################################################

# go through the following and cycle where stands # [cycle-begin] [cycle-end]


####################################################################
# paths for JT-cells-56R-improved (big clustering)
####################################################################

indirpath  <- "/media/josephus/archive_big/kworkspace/count/JT-cells-hybrid-new"

gtf.path <- paste0("/media/josephus/Elements/genome/Mus_musculus/",
                   "UCSC/mm10/Annotation/Archives/",
                   "archive-2015-07-17-14-33-26/Genes/genes.gtf")



dataname <- "JT-cells-hybrid-new-0808"
outdirpath <- "/media/josephus/Elements/DEanalysis/JT-cells-hybrid-new-0808"
# limma.outdirpath <- "/media/josephus/Elements/DEanalysis/JT-cells-hybrid-new-0723/limma"
metapath <- "/media/josephus/archive_big/kworkspace/count/JT-cells-hybrid-new/rsubread/mm10/meta.txt"


print(dataname) # for checking

####################################################################
# get meta information
####################################################################

meta.df <- read.table(metapath, sep = '\t',
                      header = TRUE,
                      stringsAsFactors = FALSE)
meta.df$condition <- factor(meta.df$condition, # ensures preferred order
                            levels = unique(meta.df$condition))

denom <- function(meta.df) as.character(meta.df$condition[meta.df$testing == "denom"][1])
num   <- function(meta.df) as.character(meta.df$condition[meta.df$testing == "num"][1])
core.name <- function(meta.df) paste0(num(meta.df), "-vs-", denom(meta.df), collapse = '')

####################################################################
# inferred paths
####################################################################

outdirpath <- file.path(outdirpath, core.name(meta.df))
# limma.outdirpath <- file.path(outdirpath, "limma", core.name(meta.df))
meta.outdir <- file.path(outdirpath, "meta")
cnts.outdir <- file.path(outdirpath, "count-table")
DE.outdir <- file.path(outdirpath, "DE-table")
plot.outdir <- file.path(outdirpath, "glimma-plots")
hm.outdir <- file.path(outdirpath, "heatmap")
hm.final.outdir <- file.path(hm.outdir, "final")
hm.up.outdir <- file.path(hm.outdir, "up")
hm.down.outdir <- file.path(hm.outdir, "down")
hm.all.outdir <- file.path(hm.outdir, "all")
hm.up.list.outdir <- file.path(hm.outdir, "up-list")
hm.down.list.outdir <- file.path(hm.outdir, "down-list")
hm.all.list.outdir <- file.path(hm.outdir, "all-list")
GO.outdir <- file.path(outdirpath, "GO")
pathway.outdir <- file.path(outdirpath, "pathway")
goi.outdir <- file.path(outdirpath, "goi")
lgenes.outdir <- file.path(outdirpath, "listgenes")

#####################################################################
# ensure existence of output paths
#####################################################################

dir.create(path=outdirpath, recursive = TRUE, showWarnings = FALSE)
# dir.create(path=limma.outdirpath, recursive = TRUE, showWarnings = FALSE)
dir.create(path=meta.outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(path=cnts.outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(path=DE.outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(path=plot.outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(path=hm.final.outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(path=hm.up.outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(path=hm.down.outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(path=hm.all.outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(path=hm.up.list.outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(path=hm.down.list.outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(path=hm.all.list.outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(path=GO.outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(path=pathway.outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(path=goi.outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(path=lgenes.outdir, recursive = TRUE, showWarnings = FALSE)


# https://guangchuangyu.github.io/2016/01/go-analysis-using-clusterprofiler/
#######################################################################
# Create meta table
#######################################################################

write.table(meta.df,
            file = file.path(meta.outdir, basename(metapath)),
            sep = "\t", row.names = FALSE, quote = FALSE)

#######################################################################
# Create DESeq2.obj and DESeq2.obj.disp
#######################################################################
indirpath <- file.path(indirpath, "rsubread/mm10")
DESeq2.obj <- meta2DESeq2.obj(meta.df, indirpath = indirpath)
DESeq2.obj.disp <- meta2DESeq2.obj(meta.df, indirpath, normalized = TRUE)

#######################################################################
# Create the count tables
#######################################################################

cnts.raw <- meta2cnts(meta.df, DESeq2.obj, outdirpath = cnts.outdir,
                      dataname = dataname,
                      printp = TRUE, normalized = FALSE, averaged = FALSE,
                      sheetName = "raw.all")
cnts.nrm <- meta2cnts(meta.df, DESeq2.obj.disp, outdirpath = cnts.outdir,
                      dataname = dataname,
                      printp = TRUE, normalized = TRUE, averaged = FALSE,
                      sheetName = "normalized.all")

cnts.avg.nrm <- meta2cnts(meta.df, DESeq2.obj.disp, outdirpath = cnts.outdir,
                          dataname = dataname,
                          printp = TRUE, normalized = TRUE, averaged = TRUE,
                          sheetName = "avg.normalized.all")

#######################################################################
# Create DE table
#######################################################################

res <- DEanalysis(meta.df, DESeq2.obj.disp, outdirpath = DE.outdir,
                  dataname = dataname,
                  printp = TRUE)
resSig <- DEanalysis(meta.df, DESeq2.obj.disp, outdirpath = DE.outdir,
                     dataname = dataname,
                     printp = TRUE,
                     filterp = TRUE, alpha = alpha, lFC = lFC)
print.cnts.DE.sortings(cnts.nrm, 
                       resSig, 
                       file.path(paste0("DE-cnts-sig-", dataname, "_", 
                                        core.name(meta.df), "_",
                                        time.now(), "-", alpha, "-", lFC, ".xlsx")), 
                       dir = DE.outdir)

print.cnts.DE.sortings(cnts.avg.nrm, 
                       resSig, 
                       file.path(paste0("DE-cnts-avg-sig-", dataname, "_", 
                                        core.name(meta.df), "_",
                                        time.now(), "-", alpha, "-", lFC, ".xlsx")), 
                       dir = DE.outdir)

# collect from counts order according to 'padj' value


















################ Here I stop analysis, because I see that JT-cells-hybrid-new is already what I want to do
################ the newer anlysis were for limma!































#######################################################################
# Create limma DE table
#######################################################################

res.limma <- DEanalysis.limma(meta.df, DESeq2.obj, outdirpath = DE.outdir,
                              dataname = dataname,
                              printp = TRUE)
resSig.limma <- DEanalysis.limma(meta.df, DESeq2.obj, outdirpath = DE.outdir,
                                 dataname = dataname,
                                 printp = TRUE,
                                 filterp = TRUE, alpha = alpha, lFC = lFC)

dataname.limma <- paste0(dataname, "-limma")
print.cnts.DE.sortings.limma(cnts.nrm,
                             resSig.limma,
                             file.path(paste0("DE-cnts-sig-", dataname.limma, "_",
                                              core.name(meta.df), "_",
                                              time.now(), "-", alpha, "-", lFC, ".xlsx")),
                             dir = DE.outdir)

print.cnts.DE.sortings.limma(cnts.avg.nrm, 
                             resSig.limma, 
                             file.path(paste0("DE-cnts-avg-sig-", dataname.limma, "_", 
                                              core.name(meta.df), "_",
                                              time.now(), "-", alpha, "-", lFC, ".xlsx")), 
                             dir = DE.outdir)

#######################################################################
# Create limma.corr DE table
#######################################################################

res.limma.corr <- DEanalysis.limma.corrected(meta.df, DESeq2.obj, outdirpath = DE.outdir,
                              dataname = dataname,
                              printp = TRUE)
resSig.limma <- DEanalysis.limma.corrected(meta.df, DESeq2.obj, outdirpath = DE.outdir,
                                 dataname = dataname,
                                 printp = TRUE,
                                 filterp = TRUE, alpha = alpha, lFC = lFC)

dataname.limma <- paste0(dataname, "-limma")
print.cnts.DE.sortings.limma(cnts.nrm,
                             resSig.limma,
                             file.path(paste0("DE-cnts-sig-", dataname.limma, "_",
                                              core.name(meta.df), "_",
                                              time.now(), "-", alpha, "-", lFC, ".xlsx")),
                             dir = DE.outdir)

print.cnts.DE.sortings.limma(cnts.avg.nrm, 
                             resSig.limma, 
                             file.path(paste0("DE-cnts-avg-sig-", dataname.limma, "_", 
                                              core.name(meta.df), "_",
                                              time.now(), "-", alpha, "-", lFC, ".xlsx")), 
                             dir = DE.outdir)


#######################################################################
# Create plots
#######################################################################

meta2iMDSplot(meta.df, DESeq2.obj.disp, outdirpath,
              dataname, top=500, launch = FALSE)

meta2iVolcano(meta.df, DESeq2.obj, DESeq2.obj.disp, outdirpath,
              dataname, alpha = alpha, lFC = lFC,
              launch = FALSE)

#######################################################################
# Create heatmap
#######################################################################


scaledata.Kmolten.core.list <- meta2heatmap(meta.df, cnts.avg.nrm, resSig,
                                            outdirpath = hm.all.outdir, selected.genes = NULL, name.add = "all",
                                            dataname = dataname,
                                            k = 12, printp=TRUE,
                                            alpha=alpha, lFC=lFC, filterp=TRUE,
                                            xlsxp=TRUE, csvp=FALSE, tsvp=FALSE)

scaledata.Kmolten.core.list <- meta2heatmap(meta.df, cnts.avg.nrm, resSig,
                                            outdirpath = hm.all.outdir, selected.genes = NULL, name.add = "all",
                                            dataname = dataname,
                                            k = 8, printp=TRUE,
                                            alpha=alpha, lFC=lFC, filterp=TRUE,
                                            xlsxp=TRUE, csvp=FALSE, tsvp=FALSE)

scaledata.Kmolten.core.list <- meta2heatmap(meta.df, cnts.avg.nrm, resSig,
                                            outdirpath = hm.up.outdir, selected.genes = "up", name.add = "up",
                                            dataname = dataname,
                                            k = 12, printp=TRUE,
                                            alpha=alpha, lFC=lFC, filterp=TRUE,
                                            xlsxp=TRUE, csvp=FALSE, tsvp=FALSE)

scaledata.Kmolten.core.list <- meta2heatmap(meta.df, cnts.avg.nrm, resSig,
                                            outdirpath = hm.up.outdir, selected.genes = "up", name.add = "up",
                                            dataname = dataname,
                                            k = 8, printp=TRUE,
                                            alpha=alpha, lFC=lFC, filterp=TRUE,
                                            xlsxp=TRUE, csvp=FALSE, tsvp=FALSE)

scaledata.Kmolten.core.list <- meta2heatmap(meta.df, cnts.avg.nrm, resSig,
                                              outdirpath = hm.down.outdir, selected.genes = "down", name.add = "down",
                                              dataname = dataname,
                                              k = 12, printp=TRUE,
                                              alpha=alpha, lFC=lFC, filterp=TRUE,
                                              xlsxp=TRUE, csvp=FALSE, tsvp=FALSE)
scaledata.Kmolten.core.list <- meta2heatmap(meta.df, cnts.avg.nrm, resSig,
                                            outdirpath = hm.down.outdir, selected.genes = "down", name.add = "down",
                                            dataname = dataname,
                                            k = 8, printp=TRUE,
                                            alpha=alpha, lFC=lFC, filterp=TRUE,
                                            xlsxp=TRUE, csvp=FALSE, tsvp=FALSE)



scaledata.Kmolten.core.list <- meta2heatmap(meta.df, cnts.avg.nrm, resSig,
                                            outdirpath = hm.up.list.outdir, selected.genes = c(up.genes), name.add = "up.list.genes",
                                            dataname = dataname,
                                            k = 12, printp=TRUE,
                                            alpha=alpha, lFC=lFC, filterp=TRUE,
                                            xlsxp=TRUE, csvp=FALSE, tsvp=FALSE)

scaledata.Kmolten.core.list <- meta2heatmap(meta.df, cnts.avg.nrm, resSig,
                                            outdirpath = hm.up.list.outdir, selected.genes = c(up.genes), name.add = "up.list.genes",
                                            dataname = dataname,
                                            k = 8, printp=TRUE,
                                            alpha=alpha, lFC=lFC, filterp=TRUE,
                                            xlsxp=TRUE, csvp=FALSE, tsvp=FALSE)

scaledata.Kmolten.core.list <- meta2heatmap(meta.df, cnts.avg.nrm, resSig,
                                            outdirpath = hm.down.list.outdir, selected.genes = c(down.genes), name.add = "down.list.genes",
                                            dataname = dataname,
                                            k = 12, printp=TRUE,
                                            alpha=alpha, lFC=lFC, filterp=TRUE,
                                            xlsxp=TRUE, csvp=FALSE, tsvp=FALSE)

scaledata.Kmolten.core.list <- meta2heatmap(meta.df, cnts.avg.nrm, resSig,
                                          outdirpath = hm.down.list.outdir, selected.genes = c(down.genes), name.add = "down.list.genes",
                                          dataname = dataname,
                                          k = 8, printp=TRUE,
                                          alpha=alpha, lFC=lFC, filterp=TRUE,
                                          xlsxp=TRUE, csvp=FALSE, tsvp=FALSE)


scaledata.Kmolten.core.list <- meta2heatmap(meta.df, cnts.avg.nrm, resSig,
                                          outdirpath = hm.all.list.outdir, selected.genes = c(up.genes, down.genes), name.add = "list.genes",
                                          dataname = dataname,
                                          k = 12, printp=TRUE,
                                          alpha=alpha, lFC=lFC, filterp=TRUE,
                                          xlsxp=TRUE, csvp=FALSE, tsvp=FALSE)

scaledata.Kmolten.core.list <- meta2heatmap(meta.df, cnts.avg.nrm, resSig,
                                            outdirpath = hm.all.list.outdir, selected.genes = c(up.genes, down.genes), name.add = "list.genes",
                                            dataname = dataname,
                                            k = 8, printp=TRUE,
                                            alpha=alpha, lFC=lFC, filterp=TRUE,
                                            xlsxp=TRUE, csvp=FALSE, tsvp=FALSE)

meta2gois(meta.df, cnts.avg.nrm, cnts.nrm, res, gois, outdirpath=goi.outdir,
          dataname = dataname, title.leader = "GOIs-heatmap_",
          alpha = alpha, lFC = lFC)

meta2gois(meta.df, cnts.avg.nrm, cnts.nrm, res, up.genes, outdirpath=lgenes.outdir,
          dataname = dataname, title.leader = "LISTgenes-UP-heatmap_",
          alpha = alpha, lFC = lFC)

meta2gois(meta.df, cnts.avg.nrm, cnts.nrm, res, down.genes, outdirpath=lgenes.outdir,
          dataname = dataname, title.leader = "LISTgenes-DOWN-heatmap_",
          alpha = alpha, lFC = lFC)

meta2gois(meta.df, cnts.avg.nrm, cnts.nrm, res, other.gois, outdirpath=goi.outdir,
          dataname = dataname, title.leader = "Ronan.gois.heatmap_",
          alpha = alpha, lFC = lFC)

fpath <- file.path(outdirpath, paste0("run", time.now(), ".RData"))
save.image(file = fpath)



# 
# ### eoGFP ordering 180528
# 
# metapath <- "/media/josephus/archive_big/kworkspace/fr_gk1029-count-0/2018_05_eobr6CR/meta-eoGFP-comp.four.samples.txt"
# 
# ####################################################################
# # get meta information
# ####################################################################
# 
# meta.df <- read.table(metapath, sep = '\t',
#                       header = TRUE,
#                       stringsAsFactors = FALSE)
# meta.df$condition <- factor(meta.df$condition, # ensures preferred order
#                             levels = unique(meta.df$condition))
# 
# denom <- function(meta.df) as.character(meta.df$condition[meta.df$testing == "denom"][1])
# num   <- function(meta.df) as.character(meta.df$condition[meta.df$testing == "num"][1])
# core.name <- function(meta.df) paste0(num(meta.df), "-vs-", denom(meta.df), collapse = '')
# 
# new.order <- c(1, 2, c(4, 3), 7, c(6, 5), 8)
# gaps.after <- c(1, 2, 4, 5, 7, 8)
# dataname <- "eoGFP-compare-0525"
# 
# s.K.path <- "/media/josephus/Elements/DEanalysis/180525-eoGFP/eoGFPcomp0525/k2-vs-wt/heatmap/all/scaledata.Kmolten.core.score.list.eoGFP-compare-0525.180525165104.all.rds"
# 
# hm.outdir <- "/media/josephus/Elements/DEanalysis/180525-eoGFP/eoGFPcomp0525/k2-vs-wt/heatmap"
# hm.all.outdir <- file.path(hm.outdir, "all")
# hm.final.outdir <- file.path(hm.outdir, "final")
# 
# UO.heatmap.path <- file.path(hm.all.outdir, filename)
# 
# # readRDS(s.K.path)
# 
# UO2O.hm(filepath = s.K.path, res, # res is alibi
#         new.order,
#         outdirpath=hm.final.outdir,
#         sheet="all",
#         gaps.after.blocks= gaps.after,
#         alpha=0.05, lFC=1, dataname = dataname)
# 
# 
# new.order <- c(1, 2, c(4, 3), c(6, 5), 7, 8)
# gaps.after <- c(1, 2, 4, 6, 7, 8)
# 
# UO2O.hm(filepath = s.K.path, res, # res is alibi
#         new.order,
#         outdirpath=hm.final.outdir,
#         sheet="all",
#         gaps.after.blocks= gaps.after,
#         alpha=0.05, lFC=1, dataname = dataname)
# 
# # ##############################################################
