##########################################
# Generate heatmaps for Fig 1e and 2a
##########################################

rm(list=ls())


####################################################################
# required packages
####################################################################
# sudo apt install libssl-dev

require(DESeq2)        # DE analysis
require(xlsx)          # xlsx printing
require(edgeR)         # DE analysis 

require(biomaRt)       # for Annotations 

require(org.Mm.eg.db)
require(org.Hs.eg.db)

require(GO.db)         # for GO analysis
require(GOstats)       # for GO analysis

# sudo apt install libcairo2-dev libxt-dev
require(Cairo)         # pdf output
CairoFonts(
  regular="Arial:style=Regular",
  bold="Arial:style=Bold",
  italic="Arial:style=Italic",
  bolditalic="Arial:style=Bold Italic,BoldItalic",
  symbol="Symbol")

require(ggplot2)       # for plots

require(SummarizedExperiment)

require(reshape)       # for heatmap
require(gplots)        # for heatmap
require(RColorBrewer)  # for heatmap
require(pheatmap)      # for heatmap

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

#######################################################################
# get sheet names of a xlsx file
#######################################################################

getSheetNames <- function(xlsx.path) {
  require(xlsx)
  wb <- loadWorkbook(xlsx.path)
  sheets <- getSheets(wb)
  names(sheets)
}
# getSheetNames(xlsx.path)

#######################################################################
# read xlsx file into list with sheet names as names
#######################################################################

xlsx2df.list <- function(xlsx.path) {
  gc()
  require(xlsx)
  sheetNames <- getSheetNames(xlsx.path)
  res <- list()
  for (i in seq_along(sheetNames)) {
    sheet <- sheetNames[i]
    res[[i]] <- read.xlsx(xlsx.path, sheetName = sheet, row.names = 1, header = TRUE)
  }
  names(res) <- sheetNames
  res
}

####################################################################
# helper functions for printing count tables
# df or df lists are printed out into excel sheets
# sheet name are names of list dfs
#
# too complex and java problems ("xlsx" package is java-dependent)
# in later scripts replaced by "openxlsx" (much better! - no java)
####################################################################

# keep memory usage low
jgc <- function() {
  gc()
  .jcall("java/lang/System", method = "gc")
}   
# https://stackoverflow.com/questions/21937640/handling-java-lang-outofmemoryerror-when-writing-to-excel-from-r

# for tsv rownames
add_rownames_col_head <- function(df, rownames_header) {
  setNames(data.frame(rownames(df), df), c(rownames_header, colnames(df)))
}

# very complex function (in later versions replaced)
triple.write <- function(obj, core_filename, outpath = ".", sheetName="raw_counts",
                         xlsxp=TRUE, csvp=FALSE, tsvp=FALSE, sorted="padj") {
  # ensure existence of outpath
  dir.create(path=outpath, recursive = TRUE, showWarnings = FALSE)
  
  # for tsv rownames
  add_rownames_col_head <- function(df, rownames_header) {
    setNames(data.frame(rownames(df), df), c(rownames_header, colnames(df)))
  }
  
  obj2df <- function(obj, sorted) {
    df <- as.data.frame(obj)
    if (sorted != "") {
      if (sorted == "name" && !("names" %in% colnames(obj))) {
        df <- df[order(rownames(df)), , drop = FALSE]
      } else if (sorted == "name" && "names" %in% colnames(obj)) {
        df <- df[order(df[, "names"]), , drop = FALSE]
      } else if (sorted %in% colnames(obj)) {
        df <- df[order(df[, sorted]), , drop=FALSE]
      } else {
        warning(paste0('"', sorted, '" not in obj!'))
      }
    }
    df
  }
  
  if (is.matrix(obj) || is.data.frame(obj) || is.DESeq2.obj(obj) || is.DGEList.obj(obj)) {
    df <- obj2df(obj, sorted)
    if (csvp) {
      write.csv(df, file = paste0(file.path(outpath, core_filename), ".csv"))
    }
    if (tsvp) {
      write.table(as.data.frame(add_rownames_col_head(df, "")), 
                  file = paste0(file.path(outpath, core_filename), ".tsv"),
                  sep = "\t", row.names = FALSE)
    }
    if (xlsxp) {
      gc()
      jgc()
      require(xlsx)
      write.xlsx(df, file = paste0(file.path(outpath, core_filename), ".xlsx"), 
                 sheetName = sheetName)
    }
  } else if (is.matrix.list(obj) || is.data.frame.list(obj) || is.DESeq2.obj.list(obj) ||
             is.DESeq2.result.list(obj) || is.DGEList.obj.list(obj)) {
    if (csvp) {
      for (i in 1:length(obj)) {
        gc()
        df <- obj2df(obj[[i]], sorted)
        write.csv(df, file = paste0(file.path(outpath, core_filename),
                                    "_", names(obj)[i], ".csv"))
      }
    }
    if (tsvp) {
      for (i in 1:length(obj)) {
        gc()
        df <- obj2df(obj[[i]], sorted)
        write.table(as.data.frame(add_rownames_col_head(df, "")), 
                    file = paste0(file.path(outpath, core_filename),
                                  "_", names(obj)[i], ".tsv"),
                    sep = "\t", row.names = FALSE)
      }
      
    }
    if (xlsxp) {
      for (i in 1:length(obj)) {
        gc()
        jgc()
        df <- obj2df(obj[[i]], sorted)
        write.xlsx(df, file = paste0(file.path(outpath, core_filename), ".xlsx"),
                   sheetName = names(obj)[i], append = ifelse(i > 1, TRUE, FALSE))
      }
    }
  }
}

# write.table - no columnname over row.names
# https://stackoverflow.com/questions/2478352/write-table-writes-unwanted-leading-empty-column-to-header-when-has-rownames


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
                          function(idxs) rowMeans(cnts.df[, idxs])))
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
                          function(idxs) rowSds(cnts.df[, idxs])))
  colnames(cnts.df) <- new_col_names
  cnts.df[, (ncol_old + 1):ncol_limit]
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
# create a  DESeq2 obj out of
# - path to counts
# - file names in meta.txt
# - (outputpath)
#######################################################################

meta2DESeq2.obj <- function(meta.df, indirpath, normalized=FALSE) {
  # require(DESeq2)
  DESeq2.obj <- DESeqDataSetFromHTSeqCount(sampleTable = meta.df,
                                           directory = indirpath,
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
    require(xlsx)
    triple.write(cnts, 
                 paste0(ifelse(normalized, "nrm-counts-", "raw-counts-"),
                        ifelse(averaged, "avg-", ""),
                        dataname, "-", core.name(meta.df)), 
                 outdirpath, 
                 sheetName = ifelse(normalized, "nrm-counts", "raw-counts"))
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
                       alpha=0.05, lFC=1, filterp=FALSE, sorted="name",
                       xlsxp=TRUE, csvp=FALSE, tsvp=FALSE) {
  dds <- DESeq(DESeq2.obj)
  res <- results(dds, contrast=c("condition", num(meta.df), denom(meta.df)), 
                 cooksCutoff = Inf, 
                 independentFiltering = FALSE) # those two avoid NA!
  if (filterp) {
    res <- subset(subset(res, padj < alpha), abs(log2FoldChange) > lFC)
  }
  if (printp || onlynamesp) {
    up <- subset(res, res$log2FoldChange > 0)
    down <- subset(res, res$log2FoldChange < 0)
    res <- list(all=res, up=up, down=down)
  }
  if (printp) {
    filecorename <- paste0("DE_", ifelse(filterp, "sig_", ""),
                           dataname, "_", core.name(meta.df), "_", time.now())
    if (filterp) { filecorename <- paste0(filecorename, "_", alpha, "_", lFC) }
    triple.write(res, filecorename, outpath = outdirpath, 
                 xlsxp=xlsxp, csvp=csvp, tsvp=tsvp,
                 sorted=sorted)
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
                         dataname = dataname,
                         k = k, printp=FALSE,
                         alpha = alpha, lFC= lFC, filterp=FALSE,
                         xlsxp=TRUE, csvp=FALSE, tsvp=FALSE) {
  res.names <- rownames(resSig$all)
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
                  "png('", outdirpath, "/k", k, "_", core.name(meta.df), "-Cluster", i, "_", dataname, "_", alpha, "_", lFC, ".png'); ",
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
  outname <- file.path(outdirpath, paste0("k", k, "_khmap_", dataname, "_", time.now(), "_UO_", alpha, "_", lFC))
  
  # unordered clusters
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
  {
    for (i in 1:k) {
      eval(parse(text = paste0("write.xlsx(scaledata_list[[", i, "]], '", 
                               outfname, "' , sheetName = 'clust", i, 
                               "', append = ", ifelse(i <= 1, FALSE, TRUE), ")")))
    }
    write.xlsx(scaledata.k, paste0(outname, ".xlsx"), sheetName = "all", append = TRUE)
  }
  
  eval(parse(text = paste0("Kmolten.list <- list(", paste(paste("K", 1:k, "molten", sep = ''), collapse = ", "), ")")))
  names(Kmolten.list) <- paste("K", 1:k, "molten", sep = '')
  
  eval(parse(text = paste0("core.list <- list(", paste(paste("core", 1:k, sep = ''), collapse = ", "), ")")))
  names(Kmolten.list) <- paste("core", 1:k, sep = '')
  
  eval(parse(text = paste0("score.list <- list(", paste(paste("score", 1:k, sep = ''), collapse = ", "), ")")))
  names(Kmolten.list) <- paste("score", 1:k, sep = '')

  data <- list(scaledata_list, Kmolten.list, core.list, score.list)
  names(data) <- c("scaledata_list", "Kmolten.list", "core.list", "score.list")
  saveRDS(data, file = paste0(outdirpath, "/scaledata.Kmolten.core.score.list.", time.now(), ".rds"))
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
  {
    for (i in 1:k) {
      eval(parse(text = paste0("write.xlsx(scaledata_list[[", i, "]], '", outfname, "', sheetName = 'Clust", i, "', append = ", ifelse(i <= 1, FALSE, TRUE), ")")))
    }
    write.xlsx(scaledata.k.ordered, outfname, sheetName = "All", append = TRUE)
  }
  
  for (i in 1:k) {
  
  }
  
}

# reordering of heatmap clusters for nicer visualization.
# ugly function.

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
  {
    for (i in 1:k) {
      eval(parse(text = paste0("write.xlsx(scaledata_list[[", i, "]], '", outfname, "', sheetName = 'Clust", i, "', append = ", ifelse(i <= 1, FALSE, TRUE), ")")))
    }
    write.xlsx(scaledata.k.ordered, outfname, sheetName = "All", append = TRUE)
  }
  
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

####################################################################
# DE settings
####################################################################
alpha <- 0.05       # padj BH limit
FDR   <- 0.05       # false discovery rate, q-value
lFC   <- 2.5        # log2FC limit


####################################################################
# paths
####################################################################

dataname   <- "JTcellsHyNew0417"
outdirpath <- "/media/josephus/Elements/DEanalysis/JT-cells-hybrid-new"
indirpath  <- "/media/josephus/Elements/count/JT-cells-hybrid-new"
metapath   <- "/media/josephus/Elements/count/JT-cells-hybrid-new/meta-k2-vs-wtn-four.txt"

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
meta.outdir <- file.path(outdirpath, "meta")
cnts.outdir <- file.path(outdirpath, "count-table")
DE.outdir <- file.path(outdirpath, "DE-table")
plot.outdir <- file.path(outdirpath, "glimma-plots")
hm.outdir <- file.path(outdirpath, "heatmap")
hm.final.outdir <- file.path(hm.outdir, "final")
GO.outdir <- file.path(outdirpath, "GO")
pathway.outdir <- file.path(outdirpath, "pathway")
goi.outdir <- file.path(outdirpath, "goi")

#####################################################################
# ensure existence of output paths
#####################################################################

dir.create(path=outdirpath, recursive = TRUE, showWarnings = FALSE)
dir.create(path=meta.outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(path=cnts.outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(path=DE.outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(path=plot.outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(path=hm.final.outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(path=GO.outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(path=pathway.outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(path=goi.outdir, recursive = TRUE, showWarnings = FALSE)


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

DESeq2.obj <- meta2DESeq2.obj(meta.df, indirpath)
DESeq2.obj.disp <- meta2DESeq2.obj(meta.df, indirpath, normalized = TRUE)

#######################################################################
# Create the count tables
#######################################################################

cnts.raw <- meta2cnts(meta.df, DESeq2.obj, outdirpath = cnts.outdir,
                      dataname = dataname,
                      printp = FALSE, normalized = FALSE, averaged = FALSE,
                      sheetName = "raw.all")
cnts.nrm <- meta2cnts(meta.df, DESeq2.obj.disp, outdirpath = cnts.outdir,
                      dataname = dataname,
                      printp = FALSE, normalized = TRUE, averaged = FALSE,
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
                  printp = TRUE, sorted = "none",
                  xlsxp = FALSE, tsvp = TRUE)
resSig <- DEanalysis(meta.df, DESeq2.obj.disp, outdirpath = DE.outdir,
                     dataname = dataname,
                     printp = TRUE, 
                     filterp = TRUE, alpha = alpha, lFC = lFC,
                     sorted = "none",
                     xlsxp = TRUE)

# collect from counts order according to 'padj' value
cnts.norm.sig.list <- list(cnts.nrm[rownames(resSig$all), ][order(resSig$all$padj), ],
                           cnts.nrm[rownames(resSig$up), ][order(resSig$up$padj), ],
                           cnts.nrm[rownames(resSig$down), ][order(resSig$down$padj), ])
cnts.avg.norm.sig.list <- list(cnts.avg.nrm[rownames(resSig$all), ][order(resSig$all$padj), ],
                               cnts.avg.nrm[rownames(resSig$up), ][order(resSig$up$padj), ],
                               cnts.avg.nrm[rownames(resSig$down), ][order(resSig$down$padj), ])
names(cnts.norm.sig.list) <- c("all", "up", "down")
names(cnts.avg.norm.sig.list) <- c("all", "up", "down")
triple.write(cnts.norm.sig.list, 
             paste0("cnts_DE_sig_", dataname, "_", core.name(meta.df), "_",
                    time.now(), "_", alpha, "_", lFC),
             outpath = DE.outdir,
             xlsxp = TRUE, sorted = "as.is")
triple.write(cnts.avg.norm.sig.list, paste0("cnts_avg_DE_sig_", dataname, "_", core.name(meta.df), "_",
                                            time.now(), "_", alpha, "_", lFC),
             outpath = DE.outdir,
             xlsxp = TRUE, sorted = "as.is")

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
                                            outdirpath = hm.outdir, 
                                            dataname = dataname,
                                            k = 8, printp=TRUE,
                                            alpha=alpha, lFC=lFC, filterp=TRUE, 
                                            xlsxp=TRUE, csvp=FALSE, tsvp=FALSE)

meta2gois(meta.df, cnts.avg.nrm, res, gois, outdirpath=goi.outdir, 
          dataname = dataname,
          alpha = alpha, lFC= lFC)


fpath <- file.path(outdirpath, paste0("run", time.now(), ".RData"))
save.image(file = fpath)




###############################
# reordering clusters
###############################

new.order  <- c(5, 7, 4, 16, 19, 3, 13, 17, 20, 1, 10, 15, 9, 6, 18, 11, 14, 12, 2, 8)
gaps.after <- c(2, 5, 9, 11, 13, 15, 18, 20)
dataname <- "JTcellsHyNew0417-k2-option2js"
hm.outdir <- "/media/josephus/Elements/DEanalysis/JT-cells-hybrid-new/k2-vs-wtn/heatmap/"
UO.heatmap.path <- file.path(hm.outdir, ".scaledata.Kmolten.core.list.180417193446.RData")
UO2O.hm(filepath = UO.heatmap.path, res, # res is alibi
        new.order, outdirpath = file.path(hm.outdir, "final"),
        gaps.after.blocks = gaps.after,
        alpha=0.05,
        lFC=2.5,
        dataname=dataname)
# UO.heatmap.path <- file.path(hm.outdir, "k20_khmap_JTcellsHyNew0417-k2_180417193441_UO_0.05_2.5.xlsx")
# unordered2orderedheatmap(UO.heatmap.path, order.vec=new.order, outdirpath=file.path(hm.outdir, "final"),
#                          sheet="all", gaps.after.blocks=gaps.after,
#                          alpha=0.05, lFC=2.5, dataname=dataname)

## however, afterwards when creating final figure order has been changed!
