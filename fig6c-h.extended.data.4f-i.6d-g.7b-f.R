#######################################################################################
# Plotting Barplots with counts as dots for Figure 6c-h, Extended Data 4f-i, 6d-g, 7b-f
#######################################################################################


##########################################
# for scale normed averaged with dots
##########################################

require(xlsx2dfs)
if (!require(rlang)) {
  install.packages("rlang")
  require(tidyverse)
}
if (!require(tidyverse)) {
  install.packages("tidyverse")
  require(tidyverse)
}
require(reshape2)

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

scale.raw.counts.with.SD <- function(cnts.DE.sig.fpath, meta.fpath, out.fpath = "") {
  cnts.DE.sig <- xlsx2dfs(cnts.DE.sig.fpath)[["all.names.srt"]]
  meta.df <- read.table(meta.fpath, sep = '\t', header = T, stringsAsFactors = F)
  meta.df$condition <- factor(meta.df$condition,
                              levels = unique(meta.df$condition))
  cond <- unique(as.character(meta.df$condition))
  cond.cols <- lapply(cond,
                      function(cd) which(as.character(meta.df$condition) == cd))
  names(cond.cols) <- cond

  cnts.avg <- counts.avg.build(cnts.DE.sig, cond.cols, cond)
  std.avg  <- counts.std.build(cnts.DE.sig, cond.cols, cond)
  
  scaledata <- t(scale(t(cnts.avg)))
  
  scaled.sds <- std.avg/apply(cnts.avg, 1, sd)
  
  upper.sds.values <- scaledata + scaled.sds

  lower.sds.values <- scaledata - scaled.sds
  
  res <- list(scaled_values = scaledata,
              scaled_StdDevs = scaled.sds,
              upper_SD_values = upper.sds.values,
              lower_SD_values = lower.sds.values)
  res <- lapply(res, as.data.frame)
  if (out.fpath != "") {
    dfs2xlsx(res, out.fpath)
  }
  res
}


scale.raw.counts.with.SD.with.orig.values <- function(cnts.DE.sig.fpath, meta.fpath, out.fpath = "") {
  cnts.DE.sig <- xlsx2dfs(cnts.DE.sig.fpath)[["all.names.srt"]]
  meta.df <- read.table(meta.fpath, sep = '\t', header = T, stringsAsFactors = F)
  meta.df$condition <- factor(meta.df$condition,
                              levels = unique(meta.df$condition))
  cond <- unique(as.character(meta.df$condition))
  cond.cols <- lapply(cond,
                      function(cd) which(as.character(meta.df$condition) == cd))
  names(cond.cols) <- cond

  cnts.avg <- counts.avg.build(cnts.DE.sig, cond.cols, cond)
  std.avg  <- counts.std.build(cnts.DE.sig, cond.cols, cond)
  
  scaledata <- t(scale(t(cnts.avg)))
  
  scaled.sds <- std.avg/apply(cnts.avg, 1, sd)
  upper.sds.values <- scaledata + scaled.sds
  lower.sds.values <- scaledata - scaled.sds
  
  ## scale the counts
  avgs.avgs <- apply(cnts.avg, MARGIN=1, FUN=mean)
  sds.avgs  <- apply(cnts.avg, MARGIN=1, FUN=sd)
  
  # scaledata <- (cnts.avg/sds.avgs - avgs.avgs/sds.avgs) # exact scaledata
  scaled.cnts.DE.sig <- (cnts.DE.sig/sds.avgs - avgs.avgs/sds.avgs)
  
  res <- list(scaled_values = scaledata,
              scaled_StdDevs = scaled.sds,
              upper_SD_values = upper.sds.values,
              lower_SD_values = lower.sds.values,
              scaled_counts = scaled.cnts.DE.sig)
  res <- lapply(res, as.data.frame)
  if (out.fpath != "") {
    dfs2xlsx(res, out.fpath)
  }
  res
}

############################################
# SEM
############################################

sem <- function(x) sd(x)/sqrt(length(x))

counts.sem.build <- function(cnts.df, cols.list, groups){
  cnts.df <- as.data.frame(cnts.df)
  rowSem <- function(df) apply(df, 1, sem)
  ncol_old <- length(colnames(cnts.df))
  ncol_limit <- length(groups) + ncol_old
  new_col_names <- c(colnames(cnts.df), groups)
  cnts.df <- cbind(cnts.df,
                   lapply(cols.list,
                          function(idxs) rowSem(cnts.df[, idxs, drop = FALSE])))
  colnames(cnts.df) <- new_col_names
  cnts.df[, (ncol_old + 1):ncol_limit]
}


scale.with.SD.SEM <- function(cnts.DE.sig.fpath, meta.fpath, out.fpath = "", tabname="all.names.srt") {
  cnts.DE.sig <- xlsx2dfs(cnts.DE.sig.fpath)[[tabname]]
  meta.df <- read.table(meta.fpath, sep = '\t', header = T, stringsAsFactors = F)
  meta.df$condition <- factor(meta.df$condition,
                              levels = unique(meta.df$condition))
  cond <- unique(as.character(meta.df$condition))
  cond.cols <- lapply(cond,
                      function(cd) which(as.character(meta.df$condition) == cd))
  names(cond.cols) <- cond

  cnts.avg <- counts.avg.build(cnts.DE.sig, cond.cols, cond)
  std.avg  <- counts.std.build(cnts.DE.sig, cond.cols, cond)
  sem.avg  <- counts.sem.build(cnts.DE.sig, cond.cols, cond)
  
  scaledata <- t(scale(t(cnts.avg)))
  
  scaled.sds <- std.avg/apply(cnts.avg, 1, sd)
  scaled.sem <- sem.avg/apply(cnts.avg, 1, sem)
  
  res <- list(scaled_values = scaledata,
              scaled_StdDevs = scaled.sds,
              scaled_Sems = scaled.sem)
  res <- lapply(res, as.data.frame)
  if (out.fpath != "") {
    dfs2xlsx(res, out.fpath)
  }
  res
}

scale.centering.with.SD.SEM <- function(cnts.DE.sig.fpath, meta.fpath, out.fpath = "", tabname="all.names.srt") {
  cnts.DE.sig <- xlsx2dfs(cnts.DE.sig.fpath)[[tabname]]
  meta.df <- read.table(meta.fpath, sep = '\t', header = T, stringsAsFactors = F)
  meta.df$condition <- factor(meta.df$condition,
                              levels = unique(meta.df$condition))
  cond <- unique(as.character(meta.df$condition))
  cond.cols <- lapply(cond,
                      function(cd) which(as.character(meta.df$condition) == cd))
  names(cond.cols) <- cond

  cnts.avg <- counts.avg.build(cnts.DE.sig, cond.cols, cond)
  std.avg  <- counts.std.build(cnts.DE.sig, cond.cols, cond)
  sem.avg  <- counts.sem.build(cnts.DE.sig, cond.cols, cond)
  
  cnts.means.wt.dko <- rowMeans(cnts.avg[, c(1, 2)])
  
  scaledata <- t(scale(t(cnts.avg), center=cnts.means.wt.dko))
  
  scaled.sds <- std.avg/apply(cnts.avg, 1, sd)
  scaled.sem <- sem.avg/apply(cnts.avg, 1, sem)
  
  res <- list(scaled_values = scaledata,
              scaled_StdDevs = scaled.sds,
              scaled_Sems = scaled.sem)
  res <- lapply(res, as.data.frame)
  if (out.fpath != "") {
    dfs2xlsx(res, out.fpath)
  }
  res
}

scale.centering.with.wt.dko.SD.SEM <- function(cnts.DE.sig.fpath, meta.fpath, out.fpath = "", tabname=1) {
  cnts.DE.sig <- xlsx2dfs(cnts.DE.sig.fpath)[[tabname]]
  meta.df <- read.table(meta.fpath, sep = '\t', header = T, stringsAsFactors = F)
  meta.df$condition <- factor(meta.df$condition,
                              levels = unique(meta.df$condition))
  cond <- unique(as.character(meta.df$condition))
  cond.cols <- lapply(cond,
                      function(cd) which(as.character(meta.df$condition) == cd))
  names(cond.cols) <- cond

  cnts.avg <- counts.avg.build(cnts.DE.sig, cond.cols, cond)
  std.avg  <- counts.std.build(cnts.DE.sig, cond.cols, cond)
  sem.avg  <- counts.sem.build(cnts.DE.sig, cond.cols, cond)
  #   std.avg.wt.dko  <- counts.std.build(cnts.DE.sig[, 1:6], cond.cols[1:2], cond[1:2])
  #   sem.avg.wt.dko  <- counts.sem.build(cnts.DE.sig[, 1:6], cond.cols[1:2], cond[1:2])

  stds <- apply(cnts.DE.sig[, 1:6], 1, sd)
  
  cnts.means.wt.dko <- rowMeans(cnts.avg[, c(1, 2)])
  
  scaledata <- t(scale(t(cnts.avg), center=cnts.means.wt.dko, scale=stds))
  
  #   scaled.sds <- std.avg.wt.dko/apply(cnts.avg[, 1:6], 1, sd)
  #   scaled.sem <- sem.avg.wt.dko/apply(cnts.avg[, 1:6], 1, sem)
  scaled.sds <- std.avg/apply(cnts.avg, 1, sd)
  scaled.sem <- sem.avg/apply(cnts.avg, 1, sem)
  res <- list(scaled_values = scaledata,
              scaled_StdDevs = scaled.sds,
              scaled_Sems = scaled.sem)
  res <- lapply(res, as.data.frame)
  if (out.fpath != "") {
    dfs2xlsx(res, out.fpath)
  }
  res
}

scale.raw.counts.with.SD.SEM.with.orig.values <- function(cnts.DE.sig.fpath, meta.fpath, out.fpath = "") {
  cnts.DE.sig <- xlsx2dfs(cnts.DE.sig.fpath)[["all.names.srt"]]
  meta.df <- read.table(meta.fpath, sep = '\t', header = T, stringsAsFactors = F)
  meta.df$condition <- factor(meta.df$condition,
                              levels = unique(meta.df$condition))
  cond <- unique(as.character(meta.df$condition))
  cond.cols <- lapply(cond,
                      function(cd) which(as.character(meta.df$condition) == cd))
  names(cond.cols) <- cond

  cnts.avg <- counts.avg.build(cnts.DE.sig, cond.cols, cond)
  std.avg  <- counts.std.build(cnts.DE.sig, cond.cols, cond)
  sem.avg  <- counts.sem.build(cnts.DE.sig, cond.cols, cond)
  
  scaledata <- t(scale(t(cnts.avg)))
  
  scaled.sds <- std.avg/apply(cnts.avg, 1, sd)
  scaled.sem <- sem.avg/apply(cnts.avg, 1, sd)
  upper.sds.values <- scaledata + scaled.sds
  lower.sds.values <- scaledata - scaled.sds
  upper.sem.values <- scaledata + scaled.sem
  lower.sem.values <- scaledata - scaled.sem
  
  ## scale the counts
  avgs.avgs <- apply(cnts.avg, MARGIN=1, FUN=mean)
  sds.avgs  <- apply(cnts.avg, MARGIN=1, FUN=sd)
  
  # scaledata <- (cnts.avg/sds.avgs - avgs.avgs/sds.avgs) # exact scaledata
  scaled.cnts.DE.sig <- (cnts.DE.sig/sds.avgs - avgs.avgs/sds.avgs)
  
  res <- list(scaled_values = scaledata,
              scaled_StdDevs = scaled.sds,
              scaled_SEM = scaled.sem,
              upper_SD_values = upper.sds.values,
              lower_SD_values = lower.sds.values,
              upper_SEM_values = upper.sem.values,
              lower_SEM_values = lower.sem.values,
              scaled_counts = scaled.cnts.DE.sig)
  res <- lapply(res, as.data.frame)
  if (out.fpath != "") {
    dfs2xlsx(res, out.fpath)
  }
  res
}

###########################################
# plotting functions
###########################################

plot.scaled.avg <- function(scale.svg.xlsx.fpath, 
                            genes, 
                            error.type="sd", 
                            add.scatter=FALSE, 
                            out.svg.fpath="") {
  dfs <- xlsx2dfs(scale.svg.xlsx.fpath)
  
  meta.fpath <- dir(file.path(dirname(dirname(scale.svg.xlsx.fpath)), "meta"), pattern=".txt", full.names=T)
  meta.df <- read.table(meta.fpath, sep = '\t', header = T, stringsAsFactors = F)
  meta.df$condition <- factor(meta.df$condition,
                              levels = unique(meta.df$condition))
  cond <- unique(as.character(meta.df$condition))
  cond.cols <- lapply(cond,
                      function(cd) which(as.character(meta.df$condition) == cd))
  names(cond.cols) <- cond

  colname2cond <- as.character(meta.df$condition)
  names(colname2cond) <- meta.df$sampleName
  
  cnts.melted <- melt(dfs[["scaled_values"]], variable.name="group")
  sd.melted   <- melt(dfs[["scaled_StdDevs"]], variable.name="group")
  sem.melted  <- melt(dfs[["scaled_SEM"]], variable.name="group")
  cnts.dots.melted <- melt(dfs[["scaled_counts"]], variable.name="group")
  cnts.melted$symbol <- rownames(dfs[["scaled_values"]])
  sd.melted$symbol   <- rownames(dfs[["scaled_StdDevs"]])
  sem.melted$symbol  <- rownames(dfs[["scaled_SEM"]])
  cnts.dots.melted$symbol <- rownames(dfs[["scaled_counts"]])
  cnts.melted <- cnts.melted[, c("symbol", "group", "value")]
  sd.melted   <- sd.melted[, c("symbol", "group", "value")]
  sem.melted  <- sem.melted[, c("symbol", "group", "value")]
  cnts.dots.melted <- cnts.dots.melted[, c("symbol", "group", "value")]
  cnts.dots.melted$group <- colname2cond[cnts.dots.melted$group]
  cnts.melted.sel <- cnts.melted[cnts.melted$symbol %in% genes, ]
  sd.melted.sel   <- sd.melted[sd.melted$symbol %in% genes, ]
  sem.melted.sel  <- sem.melted[sem.melted$symbol %in% genes, ]
  cnts.dots.melted.sel <- cnts.dots.melted[cnts.dots.melted$symbol %in% genes, ]
  cnts.sd.sem.melted.sel <- cbind(cnts.melted.sel, 
                                       sd=sd.melted.sel$value, 
                                       sem=sem.melted.sel$value)
  ## err <- if (error.type == "sd") {sd} else if (error.type == "sem") {sem} # doesn't work!
  
  ## that also simply doesn't work!
  p <- ggplot(cnts.sd.sem.melted.sel, aes(x = group, 
                                          y = value, 
                                          fill=factor(symbol, 
                                                      levels=genes))) +
      theme(panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.line = element_line(colour = "black")) + # remove background and grid
       geom_bar(stat="identity", 
                width=.75,
                color="Black",
                position = position_dodge())
  if (error.type == "sd") {
    p <- p +
       geom_errorbar(aes(ymin=value-sd, 
                         ymax=value+sd), 
                     width=.6, 
                     size=.9,
                     color="Black",
                     position=position_dodge(.75))
  } else if (error.type == "sem") {
    p <- p +
       geom_errorbar(aes(ymin=value-sem, 
                         ymax=value+sem), 
                     width=.6, 
                     size=.9,
                     color="Black",
                     position=position_dodge(.75))
  }
  if (add.scatter) {
    p <- p + geom_jitter(data = cnts.dots.melted.sel,
                         mapping = aes(x = group,
                                       y = value),
                         size = 1,
                         # width=.75,
                         position = position_dodge(0.75))
  }
  if (out.svg.fpath != "") {
    ggsave(filename=out.svg.fpath, plot=p)
  }
} ## works!


create.scaled.values.sem.dots.dir.or.fpath <- function(dir.or.fpath, genes=NULL, add = NULL, error.type = "sem") {
  if (endsWith(dir.or.fpath, ".xlsx") || endsWith(dir.or.fpath, ".txt")) {
    dir.path <- dirname(dirname(dir.or.fpath))
  } else {
    dir.path <- file.path(dir.or.fpath, "k2-vs-wtn")
  }
  meta.fpath <- dir(file.path(dir.path, "meta"), pattern = ".txt", full.names = TRUE)
  cnts.DE.sig.fpath <- dir(file.path(dir.path, "DE-table"), pattern = "DE-cnts-sig-", full.names = TRUE)
  out.sig.fpath  <- gsub("-cnts-", "-scaled-avg-dots-sem-", cnts.DE.sig.fpath)
  result.1 <- scale.raw.counts.with.SD.SEM.with.orig.values(cnts.DE.sig.fpath, meta.fpath, out.sig.fpath)
  result.2 <- scale.raw.counts.with.SD.SEM.with.orig.values(cnts.DE.sig.fpath, meta.fpath, file.path(out.revision.dir, basename(out.sig.fpath)))

  if (!is.null(genes) && !is.null(add)) {
    out.svg.fpath  <- gsub(".txt", ".svg", gsub(".xlsx", ".svg", out.sig.fpath))
    out.svg.fpath  <- gsub("-sem-", paste0("-sem-", add, "-"), out.svg.fpath)
    out.svg.fpath  <- gsub("-sem-",
                           paste0("-", error.type, "-"),
                           out.svg.fpath)
    plot.scaled.avg(out.sig.fpath,
                    genes = genes,
                    error.type=error.type,
                    add.scatter=TRUE,
                    out.svg.fpath=out.svg.fpath)
    plot.scaled.avg(out.sig.fpath,
                    genes = genes,
                    error.type=error.type,
                    add.scatter=TRUE,
                    out.svg.fpath=file.path(out.revision.dir, basename(out.svg.fpath)))
  }
}



## collect all excel output files in one folder
out.revision.dir <- "/media/josephus/Elements/DEanalysis/rev2019-fig6-scaled-counts"

############################################
# genes lists
############################################

# mesp1 and msgn1 rescue
mm.mes.genes = c("Cited1", "Myo16", "Prrx2", "Nkx2-5", "Pdgfra", "Isl1", "Gata4", "Pkdcc")
mm.epi.genes = c("Nanog", "Pou5f1", "Esrrb", "Lefty2", "Wnt8a", "Nkx1-2", "Epha2")
mm.ne.genes = c("Pou3f1", "Sox1", "Sox3", "Olig3", "Msx3", "Snph", "Zic3")

# mixl1 and foxa2 rescue
mf.mes.de.genes = c("Amot", "Eya1", "Foxa1", "Isl1", "Fbn2", "Plat", "Gata5", "Sall1")
mf.epi.genes = c("Nanog", "Pou5f1", "Sox2", "L1td1", "Lefty2", "Wnt8a", "Nkx1-2", "Epha2")
mf.ne.genes = c("Pou3f1", "Sox1", "Sox3", "Olig3", "Msx3", "Snph", "Neurod1", "Zic3")

# Bra rescue
br.mes.genes = c("Foxc1", "Foxc2", "Eya1", "Prrx1", "Msgn1", "Rspo3", "Tbx6", "Twist2", "Snai1", "Tcf7l1")
br.epi.ne.genes = c("Pou5f1", "Nanog", "Wnt8a", "Pou3f1", "Utf1", "Epha2", "Sox3", "Olig3", "Efnb1", "Ncan")
# bra rescue
brbr.mes.genes = c("Foxc2", "Myocd", "Mesp2", "Prrx2", "Rspo3", "Eya1", "Snai1")
brbr.epi.ne.genes = c("Wnt8a", "Nkx1-2", "Sox1", "Sox3", "Olig3", "Neurod1", "Phox2a")

# Eo rescue
eo.me.genes = c("Foxc1", "Myh6", "Myh7", "Tnnt1", "Eya1", "Eya2", "Prrx1", "Cer1", "Hhex", "Foxa2")
eo.epi.ne.genes = c("Pou5f1", "Nanog", "Wnt8a", "Pou3f1", "Sox1", "Sox2", "Sox3", "Olig3", "Efnb1", "Ncan")

# eotbx rescure
eo.mes.de.genes = c("Foxc2", "Myh7", "Myocd", "Mesp2", "Foxa2", "Cer1", "Hhex")
eo.epi.genes = c("Nanog", "Pou5f1", "L1td1", "Wnt8a", "Nkx1-2", "Efnb1", "Epha2", "Zic3")

# mesp1 msg1 rescue
mp.me.suppl.genes = c("Myl6", "Myh6", "Myh7", "Tcf7l2", "Tnni1", "Lmo1", "Cxcr4", "Cer1", "Dkk1")
mg.me.genes = c("Osr1", "Six2", "Hmga2", "Isl1", "Eya1", "Tbx6", "Twist2", "Gata5")

# foxa2 mixl1rescue
fx.me.suppl.genes = c("Foxj1", "Krt18", "Cldn6", "Foxd4", "Tcf7l2", "Tbx6")
mx.me.genes = c("Gsc", "Lhx1", "Prrx1", "Hhex", "Sox17", "Gata6", "Cer1")

# six2 rescue
sx.me.genes = c("Six2", "Pax2", "Pax5", "Lrp2", "Podxl", "Sall1", "Wnt8b", "Tcf7l1")

paths <- c("/media/josephus/Elements/DEanalysis/rev2019-dko-vs-wt_mesp1-mgsn1-six2/dKO-vs-WT/meta/meta-rev2019-dko-vs-wt_mesp1-mgsn1-six2.txt",
"/media/josephus/Elements/DEanalysis/rev2019-dko-vs-wt_mesp1-mgsn1-six2/dKO-vs-WT/meta/meta-rev2019-dko-vs-wt_mesp1-mgsn1-six2.txt",
"/media/josephus/Elements/DEanalysis/rev2019-dko-vs-wt_mesp1-mgsn1-six2/dKO-vs-WT/meta/meta-rev2019-dko-vs-wt_mesp1-mgsn1-six2.txt",
"/media/josephus/Elements/DEanalysis/rev2019-dko-vs-wt_mixl1-foxa2-six2/dKO-vs-WT/meta/meta-rev2019-dko-vs-wt_mixl1-foxa2-six2.txt",
"/media/josephus/Elements/DEanalysis/rev2019-dko-vs-wt_mixl1-foxa2-six2/dKO-vs-WT/meta/meta-rev2019-dko-vs-wt_mixl1-foxa2-six2.txt",
"/media/josephus/Elements/DEanalysis/rev2019-dko-vs-wt_mixl1-foxa2-six2/dKO-vs-WT/meta/meta-rev2019-dko-vs-wt_mixl1-foxa2-six2.txt",
"/media/josephus/Elements/DEanalysis/wt-vs-k2-eo5Ad-br5AD-brGFPAD/k2-vs-wtn/meta/meta-wt-k2-eo-5Ad-br-5AD-br-GFPAD.txt",
"/media/josephus/Elements/DEanalysis/wt-vs-k2-eo5Ad-br5AD-brGFPAD/k2-vs-wtn/meta/meta-wt-k2-eo-5Ad-br-5AD-br-GFPAD.txt",
"/media/josephus/Elements/DEanalysis/heatmaps-181031/br-181031/k2-vs-wtn/meta/meta-table-br-181031.txt",
"/media/josephus/Elements/DEanalysis/heatmaps-181031/br-181031/k2-vs-wtn/meta/meta-table-br-181031.txt",
"/media/josephus/Elements/DEanalysis/wt-vs-k2-eo5AdD-eoGFPAD/k2-vs-wtn/meta/meta-wt-k2-eo-5AdD-eo-GFPAD.txt",
"/media/josephus/Elements/DEanalysis/wt-vs-k2-eo5AdD-eoGFPAD/k2-vs-wtn/meta/meta-wt-k2-eo-5AdD-eo-GFPAD.txt",
"/media/josephus/Elements/DEanalysis/rev2019-jitters-190807-eo/k2-vs-wtn/meta/meta-rev2019-jitter-table-eo-190807.txt",
"/media/josephus/Elements/DEanalysis/rev2019-jitters-190807-eo/k2-vs-wtn/meta/meta-rev2019-jitter-table-eo-190807.txt",
"/media/josephus/Elements/DEanalysis/rescue-0731/mesp1D-msgn1D-hm/k2-vs-wtn/meta/meta-mesp1D-mesgn1D-hm.txt",
"/media/josephus/Elements/DEanalysis/rescue-0731/mesp1D-msgn1D-hm/k2-vs-wtn/meta/meta-mesp1D-mesgn1D-hm.txt",
"/media/josephus/Elements/DEanalysis/rescue-0731/mixl1D-foxa2D-hm/k2-vs-wtn/meta/meta-mixl1D-foxa2D-hm.txt",
"/media/josephus/Elements/DEanalysis/rescue-0731/mixl1D-foxa2D-hm/k2-vs-wtn/meta/meta-mixl1D-foxa2D-hm.txt",
"/media/josephus/Elements/DEanalysis/rev2019-six2-dko-vs-wtn/dKO-vs-WT/meta/meta-rev2019-six2-dko-vs-wtn.txt")

fpaths <- unname(sapply(paths, function(p) {
          if (endsWith(p, ".xlsx") || endsWith(p, ".txt")) {
            dir.path <- dirname(dirname(p))
          } else {
            dir.path <- file.path(p, "k2-vs-wtn")
          }
          meta.fpath <- dir(file.path(dir.path, "DE-table"), pattern = "DE-scaled-avg-dots-sem-sig-", full.names = TRUE)
          meta.fpath
          }))

genes <- list(mm.mes.genes, 
           mm.epi.genes, 
           mm.ne.genes,
           mf.mes.de.genes, 
           mf.epi.genes, 
           mf.ne.genes,
           br.mes.genes, 
           br.epi.ne.genes,
           brbr.mes.genes, 
           brbr.epi.ne.genes,
           eo.me.genes, 
           eo.epi.ne.genes,
           eo.mes.de.genes, 
           eo.epi.genes,
           mp.me.suppl.genes, 
           mg.me.genes, 
           fx.me.suppl.genes, 
           mx.me.genes,
           sx.me.genes)

name.parts <- c("mes", 
                "epi", 
                "ne", 
                "mes_de", 
                "epi", 
                "ne", 
                "mes", 
                "epi_ne", 
                "mes",
                "epi_ne", 
                "me", 
                "epi_ne", 
                "mes_de", 
                "epi_ne", 
                "me_suppl", 
                "me", 
                "me_suppl", 
                "me", 
                "me")

# > fpaths
#  [1] "/media/josephus/Elements/DEanalysis/rev2019-dko-vs-wt_mesp1-mgsn1-six2/dKO-vs-WT/DE-table/DE-scaled-avg-dots-sem-sig-rev2019-dko-vs-wt_mesp1-mgsn1-six2_dKO-vs-WT_190812144822-0.05-1.xlsx"
#  [2] "/media/josephus/Elements/DEanalysis/rev2019-dko-vs-wt_mesp1-mgsn1-six2/dKO-vs-WT/DE-table/DE-scaled-avg-dots-sem-sig-rev2019-dko-vs-wt_mesp1-mgsn1-six2_dKO-vs-WT_190812144822-0.05-1.xlsx"
#  [3] "/media/josephus/Elements/DEanalysis/rev2019-dko-vs-wt_mesp1-mgsn1-six2/dKO-vs-WT/DE-table/DE-scaled-avg-dots-sem-sig-rev2019-dko-vs-wt_mesp1-mgsn1-six2_dKO-vs-WT_190812144822-0.05-1.xlsx"
#  [4] "/media/josephus/Elements/DEanalysis/rev2019-dko-vs-wt_mixl1-foxa2-six2/dKO-vs-WT/DE-table/DE-scaled-avg-dots-sem-sig-rev2019-dko-vs-wt_mixl1-foxa2-six2_dKO-vs-WT_190812184340-0.05-1.xlsx"
#  [5] "/media/josephus/Elements/DEanalysis/rev2019-dko-vs-wt_mixl1-foxa2-six2/dKO-vs-WT/DE-table/DE-scaled-avg-dots-sem-sig-rev2019-dko-vs-wt_mixl1-foxa2-six2_dKO-vs-WT_190812184340-0.05-1.xlsx"
#  [6] "/media/josephus/Elements/DEanalysis/rev2019-dko-vs-wt_mixl1-foxa2-six2/dKO-vs-WT/DE-table/DE-scaled-avg-dots-sem-sig-rev2019-dko-vs-wt_mixl1-foxa2-six2_dKO-vs-WT_190812184340-0.05-1.xlsx"
#  [7] "/media/josephus/Elements/DEanalysis/wt-vs-k2-eo5Ad-br5AD-brGFPAD/k2-vs-wtn/DE-table/DE-scaled-avg-dots-sem-sig-wt-vs-k2-eo5Ad-br5AD-brGFPAD_k2-vs-wtn_190507114338-0.05-1.xlsx"            
#  [8] "/media/josephus/Elements/DEanalysis/wt-vs-k2-eo5Ad-br5AD-brGFPAD/k2-vs-wtn/DE-table/DE-scaled-avg-dots-sem-sig-wt-vs-k2-eo5Ad-br5AD-brGFPAD_k2-vs-wtn_190507114338-0.05-1.xlsx"            
#  [9] "/media/josephus/Elements/DEanalysis/heatmaps-181031/br-181031/k2-vs-wtn/DE-table/DE-scaled-avg-dots-sem-sig-br-hm_k2-vs-wtn_181031160803-0.05-1.xlsx"                                      
# [10] "/media/josephus/Elements/DEanalysis/heatmaps-181031/br-181031/k2-vs-wtn/DE-table/DE-scaled-avg-dots-sem-sig-br-hm_k2-vs-wtn_181031160803-0.05-1.xlsx"                                      
# [11] "/media/josephus/Elements/DEanalysis/wt-vs-k2-eo5AdD-eoGFPAD/k2-vs-wtn/DE-table/DE-scaled-avg-dots-sem-sig-wt-vs-k2-eo5AdD-eoGFPAD_k2-vs-wtn_190507113057-0.05-1.xlsx"                      
# [12] "/media/josephus/Elements/DEanalysis/wt-vs-k2-eo5AdD-eoGFPAD/k2-vs-wtn/DE-table/DE-scaled-avg-dots-sem-sig-wt-vs-k2-eo5AdD-eoGFPAD_k2-vs-wtn_190507113057-0.05-1.xlsx"                      
# [13] "/media/josephus/Elements/DEanalysis/rev2019-jitters-190807-eo/k2-vs-wtn/DE-table/DE-scaled-avg-dots-sem-sig-rev2019-jitters-190807-eo_k2-vs-wtn_190807162157-0.05-1.xlsx"                  
# [14] "/media/josephus/Elements/DEanalysis/rev2019-jitters-190807-eo/k2-vs-wtn/DE-table/DE-scaled-avg-dots-sem-sig-rev2019-jitters-190807-eo_k2-vs-wtn_190807162157-0.05-1.xlsx"                  
# [15] "/media/josephus/Elements/DEanalysis/rescue-0731/mesp1D-msgn1D-hm/k2-vs-wtn/DE-table/DE-scaled-avg-dots-sem-sig-rescue-mesp1D-msgn1D-hm_k2-vs-wtn_180801182430-0.05-1.xlsx"                 
# [16] "/media/josephus/Elements/DEanalysis/rescue-0731/mesp1D-msgn1D-hm/k2-vs-wtn/DE-table/DE-scaled-avg-dots-sem-sig-rescue-mesp1D-msgn1D-hm_k2-vs-wtn_180801182430-0.05-1.xlsx"                 
# [17] "/media/josephus/Elements/DEanalysis/rescue-0731/mixl1D-foxa2D-hm/k2-vs-wtn/DE-table/DE-scaled-avg-dots-sem-sig-rescue-mixl1D-foxa2D-hm_k2-vs-wtn_180801182607-0.05-1.xlsx"                 
# [18] "/media/josephus/Elements/DEanalysis/rescue-0731/mixl1D-foxa2D-hm/k2-vs-wtn/DE-table/DE-scaled-avg-dots-sem-sig-rescue-mixl1D-foxa2D-hm_k2-vs-wtn_180801182607-0.05-1.xlsx"                 
# [19] "/media/josephus/Elements/DEanalysis/rev2019-six2-dko-vs-wtn/dKO-vs-WT/DE-table/DE-scaled-avg-dots-sem-sig-rev2019-six2-dko-vs-wtn_dKO-vs-WT_190808141531-0.05-1.xlsx"      



for (i in 1:length(fpaths)) {
  p = fpaths[i]
  g = genes[[i]]
  nm = name.parts[i]
  # print(paste(p, " - ", nm, collapse = ""))
  create.scaled.values.sem.dots.dir.or.fpath(dir.or.fpath = p, 
                                             genes = g, 
                                             add = nm,
                                             error.type = "sem")
  create.scaled.values.sem.dots.dir.or.fpath(dir.or.fpath = p, 
                                             genes = g, 
                                             add = nm,
                                             error.type = "sd")
}

fpaths <- "/media/josephus/Elements/DEanalysis/wt-vs-k2-eo5AdD-eoGFPAD/k2-vs-wtn/DE-table/DE-scaled-avg-dots-sem-sig-wt-vs-k2-eo5AdD-eoGFPAD_k2-vs-wtn_190507113057-0.05-1.xlsx" 
genes <- list(eo.epi.ne.genes = c("Pou5f1", "Nanog", "Wnt8a", "Pou3f1", "Sox1", "Sox2", "Sox3", "Olig3", "Efnb1", "Ncan"))
name.parts <- "epi_ne"

for (i in 1:length(fpaths)) {
  p = fpaths[i]
  g = genes[[i]]
  nm = name.parts[i]
  # print(paste(p, " - ", nm, collapse = ""))
  create.scaled.values.sem.dots.dir.or.fpath(dir.or.fpath = p, 
                                             genes = g, 
                                             add = nm,
                                             error.type = "sem")
  create.scaled.values.sem.dots.dir.or.fpath(dir.or.fpath = p, 
                                             genes = g, 
                                             add = nm,
                                             error.type = "sd")
}

# DE-scaled-avg-dots-sd-epi-sig-wt-vs-k2-eo5AdD-eoGFPAD_k2-vs-wtn_190507113057-0.05-1.svg
# DE-scaled-avg-dots-sem-epi-sig-wt-vs-k2-eo5AdD-eoGFPAD_k2-vs-wtn_190507113057-0.05-1.svg
# DE-scaled-avg-dots-sem-sig-wt-vs-k2-eo5AdD-eoGFPAD_k2-vs-wtn_190507113057-0.05-1.xlsx




fpaths <- "/media/josephus/Elements/DEanalysis/rev2019-jitters-190807-eo/k2-vs-wtn/DE-table/DE-scaled-avg-dots-sem-sig-rev2019-jitters-190807-eo_k2-vs-wtn_190807162157-0.05-1.xlsx" 
genes <- list(eo.epi.genes = c("Nanog", "Pou5f1", "L1td1", "Wnt8a", "Nkx1-2", "Efnb1", "Epha2", "Zic3"))
name.parts <- "epi"

for (i in 1:length(fpaths)) {
  p = fpaths[i]
  g = genes[[i]]
  nm = name.parts[i]
  # print(paste(p, " - ", nm, collapse = ""))
  create.scaled.values.sem.dots.dir.or.fpath(dir.or.fpath = p, 
                                             genes = g, 
                                             add = nm,
                                             error.type = "sem")
  create.scaled.values.sem.dots.dir.or.fpath(dir.or.fpath = p, 
                                             genes = g, 
                                             add = nm,
                                             error.type = "sd")
}


# DE-scaled-avg-dots-sd-epi-sig-rev2019-jitters-190807-eo_k2-vs-wtn_190807162157-0.05-1.svg
# DE-scaled-avg-dots-sem-epi-sig-rev2019-jitters-190807-eo_k2-vs-wtn_190807162157-0.05-1.svg
# DE-scaled-avg-dots-sem-sig-rev2019-jitters-190807-eo_k2-vs-wtn_190807162157-0.05-1.xlsx


