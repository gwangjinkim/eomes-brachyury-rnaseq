
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

inquire(Cairo)       # pdf output
# necessary:
# sudo apt install libcairo2-dev libxt-dev
CairoFonts(
  regular="Arial:style=Regular",
  bold="Arial:style=Bold",
  italic="Arial:style=Italic",
  bolditalic="Arial:style=Bold Italic,BoldItalic",
  symbol="Symbol")

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







#######################################################################
# is cpm-normalized count table?
#######################################################################

is.cpm <- function(cnts) length(unique(colSum(cnts))) == 1

#######################################################################
# is raw count table?
#######################################################################

is.wholenumber <- function(x, tol = .Machine$double.eps) min(abs(c(x%%1, x%%1-1))) < tol

is.raw.count.table <- function(cnts) all(cnts %% 1 == 0) && !is.cpm(cnts)

is.DESeq2.obj  <- function(obj) "DESeqDataSet" %in% class(obj)

is.DGEList.obj <- function(obj) "DGEList" %in% class(obj)

is.DESeq2.result <- function(obj) "DESeqResults" %in% class(obj)

## all this checks will be easy if attr(x, "class") <- 'whatever'
## by this, all things created using these functions will be recognized
## and I could record metainformations for logging
## so there will be a 'memory' what was done before, which parameter were
## used etc

# is.data.frame  is also is.list # is already implemeneted
is.matrix.list <- function(obj) {
  !is.matrix(obj) && is.list(obj) && all(unlist(lapply(obj, is.matrix)))
}

is.data.frame.list <- function(obj) { 
  !is.data.frame(obj) && is.list(obj) && all(unlist(lapply(obj, is.data.frame)))
}

is.DESeq2.obj.list <- function(obj) { 
  !is.DESeq2.obj(obj) && is.list(obj) && all(unlist(lapply(obj, is.DESeq2.obj)))
}

is.DGEList.obj.list <- function(obj) { 
  !is.DGEList.obj(obj) && is.list(obj) && all(unlist(lapply(obj, is.DGEList.obj)))
}

is.DESeq2.result.list <- function(obj) { 
  !is.DGEList.obj(obj) && is.list(obj) && all(unlist(lapply(obj, is.DESeq2.result)))
}

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

# TODO: should be programmed as a method with typedispatch
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

#########################################################################
# vector annotation functions hs and mm
# # usage:
# hs.annotables.ensembl(test$ensgene)$symbol
# organism      input-type            output-type
# output-type can be:
# ensgene, entrez, symbol, chr, start, end, strand, biotype, description
#########################################################################

if (!require(annotables)) {
  inquire(devtools)
  # necessary:
  # update.packages(ask = FALSE)
  # inquire(devtools)
  binquire("stephenturner/annotables")
}



inquire(dplyr)

hs.annotables.ensembl <- function(ensembl_ids) {
  require(annotables)
  eids <- data.frame(ensgene = ensembl_ids)
  symbols <- eids %>%
    dplyr::inner_join(grch38, by = c("ensgene" = "ensgene")) %>%
    dplyr::select(ensgene, entrez, symbol, chr, start, end, strand, biotype, description)
  symbols
}

hs.annotables.symbol <- function(symbols) {
  require(annotables)
  ss <- data.frame(symbol = symbols)
  df <- ss %>%
    dplyr::inner_join(grch38, by = c("symbol" = "symbol")) %>%
    dplyr::select(ensgene, entrez, symbol, chr, start, end, strand, biotype, description)
  df
}


hs.annotables.entrez <- function(entrez) {
  require(annotables)
  eid <- data.frame(entrez = entrez)
  df <- eid %>%
    dplyr::inner_join(grch38, by = c("entrez" = "entrez")) %>%
    dplyr::select(ensgene, entrez, symbol, chr, start, end, strand, biotype, description)
  df
}


mm.annotables.ensembl <- function(ensembl_ids) {
  require(annotables)
  eids <- data.frame(ensgene = ensembl_ids)
  symbols <- eids %>%
    dplyr::inner_join(grcm38, by = c("ensgene" = "ensgene")) %>%
    dplyr::select(ensgene, entrez, symbol, chr, start, end, strand, biotype, description)
  symbols
}


mm.annotables.symbol <- function(symbols) {
  require(annotables)
  ss <- data.frame(symbol = symbols)
  df <- ss %>%
    dplyr::inner_join(grcm38, by = c("symbol" = "symbol")) %>%
    dplyr::select(ensgene, entrez, symbol, chr, start, end, strand, biotype, description)
  df
}


mm.annotables.entrez <- function(entrez) {
  require(annotables)
  eid <- data.frame(entrez = entrez)
  df <- eid %>%
    dplyr::inner_join(grcm38, by = c("entrez" = "entrez")) %>%
    dplyr::select(ensgene, entrez, symbol, chr, start, end, strand, biotype, description)
  df
}

# get MGI ids out of description
# gsub("]", "", sapply(strsplit(test.d, "MGI:"), function(x) x[[2]]))

#########################################################################













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



meta2gois <- function(meta.df, cnts.avg.nrm, res, gois, outdirpath=".", 
                      dataname = dataname,
                      alpha = alpha, lFC= lFC) {
  cnts.gois <- cnts.avg.nrm[gois, ]
  
  # gene-wise normalization
  scaledata <- t(scale(t(cnts.gois)))
  scaledata <- scaledata[complete.cases(scaledata), ]
  
  # prepare heatmap
  colfunc <- colorRampPalette(c("black", "red"))

  # outpath
  outname <- file.path(outdirpath, paste0("GOIs-heatmap_", dataname, "_", time.now(), "_UO_", alpha, "_", lFC))
  
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
  
  # {
  #   setEPS()
  #   postscript(paste0(outname, ".eps"))
  #   # svg(paste0(outname, ".svg"))
  #   pheatmap(scaledata, 
  #            cluster_rows = F, 
  #            cluster_cols = F,
  #            cellwidth = 40,
  #            col = colfunc(100), 
  #            fontsize_row = 0.1,
  #            border_color = NA
  #   )
  #   dev.off()
  # }
  
  # print counts

  write.xlsx(scaledata, paste0(outname, ".xlsx"), sheetName = "goi.counts")
  
  res.all.gois <- gois[gois %in% rownames(res$all)]
  res.up.gois <- gois[gois %in% rownames(res$up)]
  res.down.gois <- gois[gois %in% rownames(res$down)]
  
  write.xlsx(res$all[res.all.gois, ], paste0(outname, ".xlsx"), sheetName = "goi.DE.all", append = TRUE)
  write.xlsx(res$up[res.up.gois, ], paste0(outname, ".xlsx"), sheetName = "goi.DE.up", append = TRUE)
  write.xlsx(res$down[res.down.gois, ], paste0(outname, ".xlsx"), sheetName = "goi.DE.down", append = TRUE)
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
# # paths
# ####################################################################
# 
# dataname <- "JTcellsOld0417"
# outdirpath <- "/media/josephus/Elements/DEanalysis/JT-cells-second-new-new"
# indirpath  <- "/media/josephus/archive_big/count/JT-cells-second-new/subread"
# gtf.path <- paste0("/media/josephus/Elements/genome/Mus_musculus/",
#                    "Ensembl/GRCm38/Annotation/Archives/",
#                    "archive-2015-07-17-14-32-40/Genes/genes.gtf")
# 
# 
# # metapath   <- "/media/josephus/archive_big/count/JT-cells-second-new/subread/meta-Eon5AD.txt"
# # metapath   <- "/media/josephus/archive_big/count/JT-cells-second-new/subread/meta-Eon5gD.txt"
# # metapath   <- "/media/josephus/archive_big/count/JT-cells-second-new/subread/meta-Br5AD.txt"
# metapath   <- "/media/josephus/archive_big/count/JT-cells-second-new/subread/meta-Br5gD.txt"



####################################################################
# paths
####################################################################



indirpath  <- "/media/josephus/archive_big/count/JT-cells-Eo56R/subread"



# dataname <- "JTcellsEo56R-0417"
# outdirpath <- "/media/josephus/Elements/DEanalysis/JT-cells-Eo56R-0417/Act"
# metapath   <- "/media/josephus/archive_big/count/JT-cells-Eo56R/subread/meta-Act-clustering.txt"

dataname <- "JTcellsEo56RNogSB-0417"
outdirpath <- "/media/josephus/Elements/DEanalysis/JT-cells-Eo56R-0417/NogSB"
metapath   <- "/media/josephus/archive_big/count/JT-cells-Eo56R/subread/meta-NogSB-clustering.txt"



# ####################################################################
# # paths
# ####################################################################
# 
# indirpath  <- "/media/josephus/archive_big/count/JT-cells-old/subread"
# outdirpath <- "/media/josephus/Elements/DEanalysis/JT-cells-old"
# 
# dataname <- "JTcellsOld"
# 
# # metapath   <- "/media/josephus/archive_big/count/JT-cells-old/subread/meta-EoRAD.txt"
# 
# # metapath   <- "/media/josephus/archive_big/count/JT-cells-old/subread/meta-EoRgD.txt"
# 
# # metapath   <- "/media/josephus/archive_big/count/JT-cells-old/subread/meta-Eo6AD.txt"
# 
# metapath   <- "/media/josephus/archive_big/count/JT-cells-old/subread/meta-Eo6gD.txt"


# ####################################################################
# # paths   create all hybrid heatmap, plots etc
# ####################################################################
# 
# dataname <- "JTcellsOldHybrid-0417"
# outdirpath <- "/media/josephus/Elements/DEanalysis/JT-cells-hybrid"
# indirpath  <- "/media/josephus/Elements/count/JT-cells-hybrid-approach"
# metapath   <- "/media/josephus/Elements/count/JT-cells-hybrid-approach/meta-total.txt"


# ####################################################################
# # paths create hybrid o-o
# ####################################################################
# dataname <- "JTcellsOldHybridoO-O-0417"
# outdirpath <- "/media/josephus/Elements/DEanalysis/JT-cells-hybrid-o-o"
# indirpath  <- "/media/josephus/Elements/count/JT-cells-hybrid-approach"
# metapath   <- "/media/josephus/Elements/count/JT-cells-hybrid-approach/meta-total-o-o.txt"



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
# scaledata.Kmolten.core.list <- meta2heatmap(meta.df, cnts.avg.nrm, resSig, outdirpath = hm.outdir, 
#                                             dataname = dataname,
#                                             k = 20, printp=TRUE,
#                                             alpha=alpha, lFC=lFC, filterp=TRUE, 
#                                             xlsxp=TRUE, csvp=FALSE, tsvp=FALSE)



































# # # reordering four stripes hybrid-new/k2-vs-wtn/heatmap/.scaledata.Kmolten.core.list.180417193446.RData
# # # from sebastian: 
# # new.order  <- c(3, 5, 7, 4, 15, 16, 19, 1, 9, 10, 13, 17, 20, 2, 6, 11, 12, 14, 8, 18)
# # gaps.after <- c(3, 7, 13, 14, 18, 20)
# # dataname <- "JTcellsHyNew0417-k2"
# # hm.outdir <- "/media/josephus/Elements/DEanalysis/JT-cells-hybrid-new/k2-vs-wtn/heatmap/"
# # UO.heatmap.path <- file.path(hm.outdir, ".scaledata.Kmolten.core.list.180417193446.RData")
# # UO2O.hm(filepath = UO.heatmap.path, res, # res is alibi
# #         new.order, outdirpath = file.path(hm.outdir, "final"),
# #         gaps.after.blocks = gaps.after,
# #         alpha=alpha,
# #         lFC=lFC,
# #         dataname=dataname)
# # UO.heatmap.path <- file.path(hm.outdir, "k20_khmap_JTcellsHyNew0417-k2_180417193441_UO_0.05_2.5.xlsx")
# # unordered2orderedheatmap(UO.heatmap.path, order.vec=new.order, outdirpath=file.path(hm.outdir, "final"), 
# #                                      sheet="all", gaps.after.blocks=gaps.after,
# #                                      alpha=0.05, lFC=2.5, dataname=dataname)
# 
# 
# 
# 
# # reordering four stripes hybrid-new/k2-vs-wtn/heatmap/.scaledata.Kmolten.core.list.180417193446.RData
# # from Jelena: 
# new.order  <- c(5, 7, 4, 16, 19, 9, 20, 1, 3, 13, 17, 10, 12, 15, 6, 18, 11, 14, 2, 8)
# gaps.after <- c(2, 5, 8, 11, 13, 14, 16, 18, 20)
# dataname <- "JTcellsHyNew0417-k2-option1"
# hm.outdir <- "/media/josephus/Elements/DEanalysis/JT-cells-hybrid-new/k2-vs-wtn/heatmap/"
# UO.heatmap.path <- file.path(hm.outdir, ".scaledata.Kmolten.core.list.180417193446.RData")
# UO2O.hm(filepath = UO.heatmap.path, res, # res is alibi
#         new.order, outdirpath = file.path(hm.outdir, "final"),
#         gaps.after.blocks = gaps.after,
#         alpha=alpha,
#         lFC=lFC,
#         dataname=dataname)
# UO.heatmap.path <- file.path(hm.outdir, "k20_khmap_JTcellsHyNew0417-k2_180417193441_UO_0.05_2.5.xlsx")
# unordered2orderedheatmap(UO.heatmap.path, order.vec=new.order, outdirpath=file.path(hm.outdir, "final"), 
#                          sheet="all", gaps.after.blocks=gaps.after,
#                          alpha=0.05, lFC=2.5, dataname=dataname)
# 
# # reordering four stripes hybrid-new/k2-vs-wtn/heatmap/.scaledata.Kmolten.core.list.180417193446.RData
# # from Jelena: 
# new.order  <- c(5, 7, 4, 16, 19, 20, 3, 13, 17, 1, 10, 15, 9, 6, 18, 11, 14, 12, 2, 8)
# gaps.after <- c(2, 5, 6, 9, 11, 13, 15, 18, 20)
# dataname <- "JTcellsHyNew0417-k2-option2"
# hm.outdir <- "/media/josephus/Elements/DEanalysis/JT-cells-hybrid-new/k2-vs-wtn/heatmap/"
# UO.heatmap.path <- file.path(hm.outdir, ".scaledata.Kmolten.core.list.180417193446.RData")
# UO2O.hm(filepath = UO.heatmap.path, res, # res is alibi
#         new.order, outdirpath = file.path(hm.outdir, "final"),
#         gaps.after.blocks = gaps.after,
#         alpha=alpha,
#         lFC=lFC,
#         dataname=dataname)
# UO.heatmap.path <- file.path(hm.outdir, "k20_khmap_JTcellsHyNew0417-k2_180417193441_UO_0.05_2.5.xlsx")
# unordered2orderedheatmap(UO.heatmap.path, order.vec=new.order, outdirpath=file.path(hm.outdir, "final"), 
#                          sheet="all", gaps.after.blocks=gaps.after,
#                          alpha=0.05, lFC=2.5, dataname=dataname)

# # reordering four stripes hybrid-new/k2-vs-wtn/heatmap/.scaledata.Kmolten.core.list.180417193446.RData
# # from Jelena with respect Sebastian:
# new.order  <- c(5, 7, 4, 16, 19, 3, 13, 17, 20, 1, 10, 15, 9, 6, 18, 11, 14, 12, 2, 8)
# gaps.after <- c(2, 5, 9, 11, 13, 15, 18, 20)
# dataname <- "JTcellsHyNew0417-k2-option2js"
# hm.outdir <- "/media/josephus/Elements/DEanalysis/JT-cells-hybrid-new/k2-vs-wtn/heatmap/"
# UO.heatmap.path <- file.path(hm.outdir, ".scaledata.Kmolten.core.list.180417193446.RData")
# UO2O.hm(filepath = UO.heatmap.path, res, # res is alibi
#         new.order, outdirpath = file.path(hm.outdir, "final"),
#         gaps.after.blocks = gaps.after,
#         alpha=0.05,
#         lFC=2.5,
#         dataname=dataname)
# # UO.heatmap.path <- file.path(hm.outdir, "k20_khmap_JTcellsHyNew0417-k2_180417193441_UO_0.05_2.5.xlsx")
# # unordered2orderedheatmap(UO.heatmap.path, order.vec=new.order, outdirpath=file.path(hm.outdir, "final"),
# #                          sheet="all", gaps.after.blocks=gaps.after,
# #                          alpha=0.05, lFC=2.5, dataname=dataname)

# reordering four stripes hybrid-new/k2-vs-wtn/heatmap/.scaledata.Kmolten.core.list.180417193446.RData
# from Jelena with respect Sebastian:
new.order  <- c(c(5, 7), c(4, 16, 19), c(3, 13, 17, 20, 9, 1), c(10, 12), c(15), c(6, 18), c(11, 14), c(2, 8))
gaps.after <- c(2, 5, 11, 13, 14, 16, 18, 20)
dataname <- "JTcellsHyNew0417-k2-option2f"
hm.outdir <- "/media/josephus/Elements/DEanalysis/JT-cells-hybrid-new/k2-vs-wtn/heatmap/"
UO.heatmap.path <- file.path(hm.outdir, ".scaledata.Kmolten.core.list.180417193446.RData")
UO2O.hm(filepath = UO.heatmap.path, res, # res is alibi
        new.order, outdirpath = file.path(hm.outdir, "final"),
        gaps.after.blocks = gaps.after,
        alpha=0.05,
        lFC=2.5,
        dataname=dataname)






# reordering four stripes hybrid-new/k2-vs-wtn/heatmap/.scaledata.Kmolten.core.list.180417193446.RData
# from Sebastian:
new.order  <- c(c(5, 7, 8, 15, 11), c(9, 10, 17, 19, 20), c(18), c(3, 4, 12), c(1, 2, 16), c(6, 13, 14))
gaps.after <- c(5, 10, 11, 14, 17, 20)
dataname <- "JTcells-k2vswtn-nen-es-v1"
filename <- ".scaledata.Kmolten.core.list.180417195704.RData"
hm.outdir <- "/media/josephus/Elements/DEanalysis/JT-cells-hybrid-new/k2-wtn-nen-es/k2-vs-wtn/heatmap"
UO.heatmap.path <- file.path(hm.outdir, filename)
UO2O.hm(filepath = UO.heatmap.path, res, # res is alibi
        new.order, outdirpath = file.path(hm.outdir, "final"),
        gaps.after.blocks = gaps.after,
        alpha=0.05,
        lFC=2.5,
        dataname=dataname)


# reordering four stripes hybrid-new/k2-vs-wtn/heatmap/.scaledata.Kmolten.core.list.180417193446.RData
# from Sebastian:
new.order  <- c(c(5, 7, 8, 15, 11), c(9, 10, 17, 19, 20), c(18), c(3, 4, 12), c(1, 2, 16), c(6, 13, 14))
gaps.after <- 1:20
dataname <- "JTcells-k2vswtn-nen-es-v1-test"
filename <- ".scaledata.Kmolten.core.list.180417195704.RData"
hm.outdir <- "/media/josephus/Elements/DEanalysis/JT-cells-hybrid-new/k2-wtn-nen-es/k2-vs-wtn/heatmap"
UO.heatmap.path <- file.path(hm.outdir, filename)
UO2O.hm(filepath = UO.heatmap.path, res, # res is alibi
        new.order, outdirpath = file.path(hm.outdir, "final"),
        gaps.after.blocks = gaps.after,
        alpha=0.05,
        lFC=2.5,
        dataname=dataname)


# reordering four stripes hybrid-new/k2-vs-wtn/heatmap/.scaledata.Kmolten.core.list.180417193446.RData
# from Sebastian:
new.order.seb  <- c(c(5, 7, 8, 15, 11), c(9, 10, 17, 19, 20), c(18), c(3, 4, 12), c(1, 2, 16), c(6, 13, 14))
names(new.order.seb) <- 1:20
new.order.corr <- new.order.seb[new.order]
new.order.corr
# 11  6 10  7  9  8  2  1  3  4  5 17 15 16 18 12 13 14 20 19 
# 18  9 20 10 19 17  7  5  8 15 11 16  1  2  6  3  4 12 14 13

# (this was correction of name mapping misunderstanding)
# new.order  <- c(c(11), c(6, 10), c(7, 9, 8), c(2, 1, 3), c(4, 5, 17), c(15, 16, 18), c(12, 13, 14), c(20, 19))
new.order <- c(c(18), c(9, 20), c(10, 19, 17), c(7, 5, 8), c(15, 11, 16), c(1, 2, 6), c(3, 4, 12), c(14, 13))
gaps.after <- c(1, 3, 6, 9, 12, 15, 18, 20)
dataname <- "JTcells-k2vswtn-nen-es-j"
filename <- ".scaledata.Kmolten.core.list.180417195704.RData"
hm.outdir <- "/media/josephus/Elements/DEanalysis/JT-cells-hybrid-new/k2-wtn-nen-es/k2-vs-wtn/heatmap"
UO.heatmap.path <- file.path(hm.outdir, filename)
UO2O.hm(filepath = UO.heatmap.path, res, # res is alibi
        new.order, outdirpath = file.path(hm.outdir, "final"),
        gaps.after.blocks = gaps.after,
        alpha=0.05,
        lFC=2.5,
        dataname=dataname)
