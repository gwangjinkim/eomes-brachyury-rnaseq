###########################################################
# Plotting Jitter plots for Figure 5 k-r (first revision)
###########################################################

######################################################
# rev2019-jitters-190807
#######################################################
# nano /media/josephus/Elements/metas/meta-rev2019-jitter-table-all-190807.txt

# #######################################################
# # rev2019-jitters-190807-br
# #######################################################
# 
# # nano /media/josephus/Elements/metas/meta-rev2019-jitter-table-br-190807-corr.txt
# 
# sampleName    fileName    condition    testing
# wtn-1    /media/josephus/archive_big/kworkspace/count/JT-cells-second-new/subread/WTnewrep1-sym-fcount.tab    wtn    denom
# wtn-2    /media/josephus/archive_big/kworkspace/count/JT-cells-second-new/subread/WTnewrep2-sym-fcount.tab    wtn    denom
# wtn-3    /media/josephus/archive_big/kworkspace/count/JT-cells-second-new/subread/WTnewrep3-sym-fcount.tab    wtn    denom
# k2-1    /media/josephus/archive_big/count/JT-cells-hybrid-new/subread/k2-1-sym-fcount.tab    k2    num
# k2-2    /media/josephus/archive_big/count/JT-cells-hybrid-new/subread/k2-2-sym-fcount.tab    k2    num
# k2-3    /media/josephus/archive_big/count/JT-cells-hybrid-new/subread/k2-3-sym-fcount.tab    k2    num
# EoV5_nodox-1    /media/josephus/archive_big/count/JT-cells-Eo56R/subread/EoV5newActrep1-sym-fcount.tab    EoV5_nodox    
# EoV5_nodox-2    /media/josephus/archive_big/count/JT-cells-Eo56R/subread/EoV5newActrep2-sym-fcount.tab    EoV5_nodox    
# EoV5_nodox-3    /media/josephus/archive_big/count/JT-cells-Eo56R/subread/EoV5newActrep3-sym-fcount.tab    EoV5_nodox    
# br-GFP-1    /media/josephus/archive_big/kworkspace/count/BraChIP2018/rsubread/grcm38/BraGFPdox1-sym-fcount.tab    brGFP    
# br-GFP-2    /media/josephus/archive_big/kworkspace/count/BraChIP2018/rsubread/grcm38/BraGFPdox2-sym-fcount.tab    brGFP    
# br-GFP-3    /media/josephus/archive_big/kworkspace/count/BraChIP2018/rsubread/grcm38/BraGFPdox3-sym-fcount.tab    brGFP    
# br-6AD-1    /media/josephus/archive_big/kworkspace/count/2018_05_eobr6CR/subread/br-6AD-1-sym-fcount.tab    br-6AD    
# br-6AD-2    /media/josephus/archive_big/kworkspace/count/2018_05_eobr6CR/subread/br-6AD-2-sym-fcount.tab    br-6AD    
# br-6AD-3    /media/josephus/archive_big/kworkspace/count/2018_05_eobr6CR/subread/br-6AD-3-sym-fcount.tab    br-6AD    
# br-RAD-1    /media/josephus/archive_big/kworkspace/fr_gk1029-count-0/JT-cells-second-new/subread/BraEnRActdoxrep1-sym-fcount.tab    br-RAD    
# br-RAD-2    /media/josephus/archive_big/kworkspace/fr_gk1029-count-0/2018_05_eobr6CR/subread/br-RAD-2-sym-fcount.tab    br-RAD    
# br-RAD-3    /media/josephus/archive_big/kworkspace/fr_gk1029-count-0/2018_05_eobr6CR/subread/br-RAD-3-sym-fcount.tab    br-RAD    
# BraTBX-1    /media/josephus/archive_big/count/eo_bra_dko_revision_2019/H3FC5BGXC_JT_29072019_19s003327-1-1_Arnold_lane1BraTBXdoxRep1-sym-fcount.tab    BraTBX    
# BraTBX-2    /media/josephus/archive_big/count/eo_bra_dko_revision_2019/H3FC5BGXC_JT_29072019_19s003327-1-1_Arnold_lane1BraTBXdoxRep2-sym-fcount.tab    BraTBX    
# BraTBX-3    /media/josephus/archive_big/count/eo_bra_dko_revision_2019/H3FC5BGXC_JT_29072019_19s003327-1-1_Arnold_lane1BraTBXdoxRep3-sym-fcount.tab    BraTBX    


# #######################################################
# # rev2019-jitters-190807-eo
# #######################################################
# 
# # nano /media/josephus/Elements/metas/meta-rev2019-jitter-table-eo-190807.txt
# 
# sampleName    fileName    condition    testing
# wtn-1    /media/josephus/archive_big/kworkspace/count/JT-cells-second-new/subread/WTnewrep1-sym-fcount.tab    wtn    denom
# wtn-2    /media/josephus/archive_big/kworkspace/count/JT-cells-second-new/subread/WTnewrep2-sym-fcount.tab    wtn    denom
# wtn-3    /media/josephus/archive_big/kworkspace/count/JT-cells-second-new/subread/WTnewrep3-sym-fcount.tab    wtn    denom
# k2-1    /media/josephus/archive_big/count/JT-cells-hybrid-new/subread/k2-1-sym-fcount.tab    k2    num
# k2-2    /media/josephus/archive_big/count/JT-cells-hybrid-new/subread/k2-2-sym-fcount.tab    k2    num
# k2-3    /media/josephus/archive_big/count/JT-cells-hybrid-new/subread/k2-3-sym-fcount.tab    k2    num
# EoV5_nodox-1    /media/josephus/archive_big/count/JT-cells-Eo56R/subread/EoV5newActrep1-sym-fcount.tab    EoV5_nodox    
# EoV5_nodox-2    /media/josephus/archive_big/count/JT-cells-Eo56R/subread/EoV5newActrep2-sym-fcount.tab    EoV5_nodox    
# EoV5_nodox-3    /media/josephus/archive_big/count/JT-cells-Eo56R/subread/EoV5newActrep3-sym-fcount.tab    EoV5_nodox    
# eo-GFP-1    /media/josephus/archive_big/kworkspace/count/BraChIP2018/rsubread/grcm38/EoGFPdox1-sym-fcount.tab    eoGFP    
# eo-GFP-2    /media/josephus/archive_big/kworkspace/count/BraChIP2018/rsubread/grcm38/EoGFPdox2-sym-fcount.tab    eoGFP    
# eo-GFP-3    /media/josephus/archive_big/kworkspace/count/BraChIP2018/rsubread/grcm38/EoGFPdox3-sym-fcount.tab    eoGFP    
# eo-6AD-1    /media/josephus/archive_big/kworkspace/fr_gk1029-count-0/JT-cells-old/subread/eo-6AD-1-sym-fcount.tab    eo-6AD    
# eo-6AD-2    /media/josephus/archive_big/kworkspace/fr_gk1029-count-0/JT-cells-old/subread/eo-6AD-2-sym-fcount.tab    eo-6AD    
# eo-6AD-3    /media/josephus/archive_big/kworkspace/fr_gk1029-count-0/JT-cells-old/subread/eo-6AD-3-sym-fcount.tab    eo-6AD    
# eo-RAD-1    /media/josephus/archive_big/kworkspace/fr_gk1029-count-0/JT-cells-old/subread/eo-RAD-1-sym-fcount.tab    eo-RAD    
# eo-RAD-2    /media/josephus/archive_big/kworkspace/fr_gk1029-count-0/JT-cells-old/subread/eo-RAD-2-sym-fcount.tab    eo-RAD    
# eo-RAD-3    /media/josephus/archive_big/kworkspace/fr_gk1029-count-0/JT-cells-old/subread/eo-RAD-3-sym-fcount.tab    eo-RAD    
# EoTBX-1    /media/josephus/archive_big/count/eo_bra_dko_revision_2019/H3FC5BGXC_JT_29072019_19s003327-1-1_Arnold_lane1EoTBXdoxRep1-sym-fcount.tab    EoTBX    
# EoTBX-2    /media/josephus/archive_big/count/eo_bra_dko_revision_2019/H3FC5BGXC_JT_29072019_19s003327-1-1_Arnold_lane1EoTBXdoxRep2-sym-fcount.tab    EoTBX    
# EoTBX-3    /media/josephus/archive_big/count/eo_bra_dko_revision_2019/H3FC5BGXC_JT_29072019_19s003327-1-1_Arnold_lane1EoTBXdoxRep3-sym-fcount.tab    EoTBX    


source(file.path("/home/josephus/Dropbox/a_myRNAseq/rdownstream", 
                 "190801.rdownstream.central.R"))
script_dir <- "/home/josephus/Dropbox/a_myRNAseq/rdownstream"
out_dir    <- "/media/josephus/Elements/DEanalysis"
meta_dir   <- "/media/josephus/Elements/metas"
dir.create(out_dir, recursive=TRUE)

####################################################################
# DE settings
####################################################################

alpha <- 0.05       # padj BH limit
FDR   <- 0.05       # false discovery rate, q-value
lFC   <- 1          # log2FC limit

data_names  <-  c("rev2019-jitters-190807-br",
                  "rev2019-jitters-190807-eo")
outdir_names <- data_names
meta_fnames <- c("meta-rev2019-jitter-table-br-190807.txt",
                 "meta-rev2019-jitter-table-eo-190807.txt")

for (i in 1:2) {

  #########################################
  # set variables here which are common to all runs
  #########################################
  
  indirpath <- ""
  
  #########################################
  # define variables specific for each run
  #########################################
  
  dataname <- data_names[i]
  outdirpath <- file.path(out_dir, outdir_names[i])
  dir.create(outdirpath, recursive=TRUE)
  metapath <-   file.path(meta_dir, meta_fnames[i])
  
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
  hm.outdir <- file.path(outdirpath, "heatmap")
  hm.all.outdir <- file.path(hm.outdir, "all")
  GO.outdir <- file.path(outdirpath, "GO")

  #####################################################################
  # ensure existence of output paths
  #####################################################################
  
  dir.create(path=outdirpath, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=meta.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=cnts.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=DE.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=hm.all.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=GO.outdir, recursive = TRUE, showWarnings = FALSE)
  
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
  
  print.cnts.DE.sortings(cnts.avg.nrm$`nrm-counts-avg`, 
                         resSig, 
                         file.path(paste0("DE-cnts-avg-sig-", dataname, "_", 
                                          core.name(meta.df), "_",
                                          time.now(), "-", alpha, "-", lFC, ".xlsx")), 
                         dir = DE.outdir)
  
  #######################################################################
  # Create heatmap
  #######################################################################
  
  set.seed(123) # for reproducibility of heatmaps set seeds to a value
  try(scaledata.Kmolten.core.list <- meta2heatmap(meta.df, cnts.avg.nrm$`nrm-counts-avg`, resSig,
                                              outdirpath = hm.all.outdir, selected.genes = NULL, name.add = "all",
                                              dataname = dataname,
                                              k = 16, printp=TRUE,
                                              alpha=alpha, lFC=lFC, filterp=TRUE,
                                              xlsxp=TRUE, csvp=FALSE, tsvp=FALSE))
  
  
  #######################################################################
  # Create GO
  #######################################################################
  
  last <- function(l) l[length(l)]
  DEfpath <- last(dir(DE.outdir, pattern = "DE_sig", full.names = TRUE))

  
  try(do_all_GO_enrichment(DEfpath, GO.outdir))

  #######################################################################
  # save session image for later investigation of the run
  #######################################################################
  fpath <- file.path(outdirpath, paste0("run", time.now(), ".RData"))
  save.image(file = fpath)
  gc()
}


# jitter code is also in rdownstream.180731.BraChIP2018RNAseqs.R!

########################################
# jitter plots
########################################

require(xlsx2dfs)
require(ggplot2)

fpath2jitter <- function(fpath, clusts, outfname) {

  dirPath <- dirname(dirname(fpath))
  
  data <- if (ensWith(x = fpath, suffix = ".rds")) {
                readRDS(fpath)
          } else {
                dfs <- xlsx2dfs(fpath)
          }
  metaPath <- dir(file.path(dirPath, "meta"), full.name=T)

  meta.df <- read.table(metaPath, sep = '\t', header = T, stringsAsFactors = F)
  meta.df$condition <- factor(meta.df$condition,
                              levels = unique(meta.df$condition))
  denom <- function(meta.df) as.character(meta.df$condition[meta.df$testing == "denom"][1])
  num <- function(meta.df) as.character(meta.df$condition[meta.df$testing == "num"][1])
  core.name <- function(meta.df) paste0(num(meta.df), "-vs-", denom(meta.df), collapse = '')

  K.mgd.molten <- Reduce(rbind, data$Kmolten.list[clusts])

  outfpath <- file.path(dirPath, "heatmap", outfname)

  min.score <- min(Reduce(rbind, data$Kmolten.list)$score)
  
  mean.values <- colMeans(Reduce(rbind, lapply(data$core.list[clusts], 
                                               function(df) df$value)))
  mean.core <- data$core.list[[clusts[1]]]
  mean.core$value <- mean.values
  
  sp1 <- ggplot(K.mgd.molten, 
                aes(x = factor(sample, levels = unique(meta.df$condition)),
                    y = value)) +
         geom_jitter(aes(colour=score, group=gene)) +
         scale_color_gradientn(colours=c('black', 'red2'),
         limits = c(min.score, 1.0)) +
         geom_jitter(data=mean.core, aes(sample,value, group=cluster),
                     color='yellow', inherit.aes=F) +
         xlab('Time') +
         ylab('Expression') +
         labs(title=paste0("Cluster ", paste(clusts, collapse = ' + '), 
                           " Expression by Group", color = 'Score'))
  svg(outfpath)
    print(sp1)
  dev.off()
}

source("~/Dropbox/R/central-scripts/utils.interleave.strings.R")

fpath2jitter.genes <- function(fpath, genes, outfname, addTitle="Selected Genes", folder="") {

  dirPath <- splitStrToVec(fpath, "/heatmap/")[1]

  data <- readRDS(fpath)
  metaPath <- dir(file.path(dirPath, "meta"), full.name=T)

  meta.df <- read.table(metaPath, sep = '\t', header = T, stringsAsFactors = F)
  meta.df$condition <- factor(meta.df$condition,
                              levels = unique(meta.df$condition))
  denom <- function(meta.df) as.character(meta.df$condition[meta.df$testing == "denom"][1])
  num <- function(meta.df) as.character(meta.df$condition[meta.df$testing == "num"][1])
  core.name <- function(meta.df) paste0(num(meta.df), "-vs-", denom(meta.df), collapse = '')

  K.mgd.molten <- Reduce(rbind, data$Kmolten.list)
  K.mgd.molten <- K.mgd.molten[K.mgd.molten$gene %in% genes, ]

  if (folder == "") {
    outfpath <- file.path(dirPath, "heatmap", outfname)
  } else {
    outfpath <- file.path(dirPath, "heatmap", folder, outfname)
  }
  dir.create(dirname(outfpath), recursive=TRUE)

  min.score <- min(Reduce(rbind, data$Kmolten.list)$score)
  
  mean.values <- sapply(split(K.mgd.molten, K.mgd.molten$sample), 
                              function(df) mean(df$value))
  mean.core <- data$core.list[[1]] # just to get the df format
  mean.core$value <- mean.values[as.character(mean.core$sample)]   # fill with real mean values
  
  sp1 <- ggplot(K.mgd.molten, 
                aes(x = factor(sample, levels = unique(meta.df$condition)),
                    y = value)) +
         geom_jitter(aes(colour=score, group=gene)) +
         scale_color_gradientn(colours=c('black', 'red2'),
         limits = c(min.score, 1.0)) +
         geom_jitter(data=mean.core, aes(sample,value, group=cluster),
                     color='yellow', inherit.aes=F) +
         xlab('Time') +
         ylab('Expression') +
         labs(title=paste0("Cluster ", addTitle, 
                           " Expression by Group", color = 'Score'))
  svg(outfpath)
    print(sp1)
  dev.off()
  
  dfs2xlsx(withNames(addTitle, K.mgd.molten), gsub(".svg", ".xlsx", outfpath))
}

# require(devtools)
# install.github("tidyverse/magrittr")

require(xlsx2dfs)

genesPath <- "/media/josephus/Elements/DEanalysis/JT-cells-Eo56R-improved-0517/06112018_jitter_plot_list.xlsx"
dfs <- xlsx2dfs(genesPath)

genes.list <- lapply(dfs, rownames)
genes.list.names <- names(genes.list)

genes.UP <- genes.list[c(1, 3, 5, 7)]
genes.DOWN <- genes.list[c(2, 4, 6, 8)]

genes.names.UP <- genes.list.names[c(1, 3, 5, 7)]
genes.names.DOWN <- genes.list.names[c(2, 4, 6, 8)]

fpath.1 <- '/media/josephus/Elements/DEanalysis/rev2019-jitters-190807-br/k2-vs-wtn/heatmap/all/scaledata.Kmolten.core.score.list.rev2019-jitters-190807-br.190807170032.all.rds'
fpath.2 <- '/media/josephus/Elements/DEanalysis/rev2019-jitters-190807-eo/k2-vs-wtn/heatmap/all/scaledata.Kmolten.core.score.list.rev2019-jitters-190807-eo.190807162235.all.rds'

for (i in 1:length(genes.UP)) {
  fpath2jitter.genes(fpath.1, genes.UP[[i]], paste0(gsub(" ", "_", genes.names.UP[i]), "_jitter.svg"), addTitle=gsub(" ", "_", genes.names.UP[i]), folder="UP")
}

for (i in 1:length(genes.DOWN)) {
  fpath2jitter.genes(fpath.1, genes.DOWN[[i]], paste0(gsub(" ", "_", genes.names.DOWN[i]), "_jitter.svg"), addTitle=gsub(" ", "_", genes.names.DOWN[i]), folder="DOWN")
}

for (i in 1:length(genes.UP)) {
  fpath2jitter.genes(fpath.2, genes.UP[[i]], paste0(gsub(" ", "_", genes.names.UP[i]), "_jitter.svg"), addTitle=gsub(" ", "_", genes.names.UP[i]), folder="UP")
}

for (i in 1:length(genes.DOWN)) {
  fpath2jitter.genes(fpath.2, genes.DOWN[[i]], paste0(gsub(" ", "_", genes.names.DOWN[i]), "_jitter.svg"), addTitle=gsub(" ", "_", genes.names.DOWN[i]), folder="DOWN")
}

#########################################################
# create scaled
#########################################################

meta.fpath <- '/media/josephus/Elements/DEanalysis/rev2019-jitters-190807-br/k2-vs-wtn/meta/meta-rev2019-jitter-table-br-190807.txt' 
DE.sig.cnts.fpath <- '/media/josephus/Elements/DEanalysis/rev2019-jitters-190807-br/k2-vs-wtn/DE-table/DE-cnts-sig-rev2019-jitters-190807-br_k2-vs-wtn_190807165956-0.05-1.xlsx' 
out.fpath         <- gsub("-counts-", "-scaled-", DE.sig.cnts.fpath)
result <- scale.raw.counts.with.SD(DE.sig.cnts.fpath, meta.fpath, out.fpath, tabname=1)

meta.fpath <- '/media/josephus/Elements/DEanalysis/rev2019-jitters-190807-eo/k2-vs-wtn/meta/meta-rev2019-jitter-table-eo-190807.txt' 
DE.sig.cnts.fpath <- '/media/josephus/Elements/DEanalysis/rev2019-jitters-190807-eo/k2-vs-wtn/DE-table/DE-cnts-sig-rev2019-jitters-190807-eo_k2-vs-wtn_190807162157-0.05-1.xlsx' 
out.fpath         <- gsub("-counts-", "-scaled-", DE.sig.cnts.fpath)
result <- scale.raw.counts.with.SD(DE.sig.cnts.fpath, meta.fpath, out.fpath, tabname=1)


##########################################################
# new jitter analysis
##########################################################
# also in main revision script 190805.rdownstream.revision.analysis.txt

require(gglpot2)
require(xlsx2dfs)

data <- readRDS(file.path(in.dir, "heatmap", "scaledata.Kmolten.core.score.list.180518084645.rds"))
out.fpath <- paste0(outdirpath, "/k24_khmap_JTcellsEo56R-0517_180518084637_UO_0.05_1_K.3.5.12.molten.svg")

rds.fpath <- '/media/josephus/Elements/DEanalysis/rev2019-jitters-190807-br/k2-vs-wtn/heatmap/all/scaledata.Kmolten.core.score.list.rev2019-jitters-190807-br.190807170032.all.rds' 
xlsx.fpath <- '/media/josephus/Elements/DEanalysis/rev2019-jitters-190807-br/k2-vs-wtn/heatmap/all/k16_khmap_rev2019-jitters-190807-br_190807170031_UO_0.05_1all.xlsx'

fpath2jitter <- function(fpath, clusts, outfname) {

    dirPath <- dirname(dirname(fpath))
  
if (endsWith(x = fpath, suffix = ".rds")) {
    data.list <- readRDS(fpath)
    K.mgd.molten <- Reduce(rbind, data.list$Kmolten.list[clusts])
} else {
    dfs <- xlsx2dfs(fpath)
    data <- Reduce(rbind, dfs)
    K.mgd.molten.1 <- melt(data)
}
  metaPath <- dir(file.path(dirname(dirPath), "meta"), full.name=T)

  meta.df <- read.table(metaPath, sep = '\t', header = T, stringsAsFactors = F)
  meta.df$condition <- factor(meta.df$condition,
                              levels = unique(meta.df$condition))
  denom <- function(meta.df) as.character(meta.df$condition[meta.df$testing == "denom"][1])
  num <- function(meta.df) as.character(meta.df$condition[meta.df$testing == "num"][1])
  core.name <- function(meta.df) paste0(num(meta.df), "-vs-", denom(meta.df), collapse = '')

  K.mgd.molten <- Reduce(rbind, data$Kmolten.list[clusts])

  outfpath <- file.path(dirPath, "heatmap", outfname)

  min.score <- min(Reduce(rbind, data$Kmolten.list)$score)
  
  mean.values <- colMeans(Reduce(rbind, lapply(data$core.list[clusts], 
                                               function(df) df$value)))
  mean.core <- data$core.list[[clusts[1]]]
  mean.core$value <- mean.values
  
  sp1 <- ggplot(K.mgd.molten, 
                aes(x = factor(sample, levels = unique(meta.df$condition)),
                    y = value)) +
         geom_jitter(aes(colour=score, group=gene)) +
         scale_color_gradientn(colours=c('black', 'red2'),
         limits = c(min.score, 1.0)) +
         geom_jitter(data=mean.core, aes(sample,value, group=cluster),
                     color='yellow', inherit.aes=F) +
         xlab('Time') +
         ylab('Expression') +
         labs(title=paste0("Cluster ", paste(clusts, collapse = ' + '), 
                           " Expression by Group", color = 'Score'))
  svg(outfpath)
    print(sp1)
  dev.off()
}


path.up.to <- function(path, dirname) {
  # Chunk away path until last folder is dirname
  tmppath <- path
  while (basename(tmppath) != dirname && tmppath != "") {
    tmppath <- dirname(tmppath)
  }
  tmppath
}


fpath2jitter.clusters <- function(fpath, clusts, outfname, addTitle="Selected Clusters", folder="", metaPath="") {

  dirPath <- path.up.to(fpath, "heatmap")
  
  data <- readRDS(fpath)
  
  tmpPath <- dirPath
  metaPath <- dir(file.path(dirname(dirPath), "meta"), pattern = "meta-", full.names=T)

  meta.df <- read.table(metaPath, sep = '\t', header = T, stringsAsFactors = F)
  meta.df$condition <- factor(meta.df$condition,
                              levels = unique(meta.df$condition))
  denom <- function(meta.df) as.character(meta.df$condition[meta.df$testing == "denom"][1])
  num <- function(meta.df) as.character(meta.df$condition[meta.df$testing == "num"][1])
  core.name <- function(meta.df) paste0(num(meta.df), "-vs-", denom(meta.df), collapse = '')

  K.mgd.molten <- if(length(clusts) > 1) {
                    Reduce(rbind, data$Kmolten.list[clusts])
                  } else {
                    data$Kmolten.list[[clusts]]
                  }

  min.score <- min(Reduce(rbind, data$Kmolten.list)$score)
  
  mean.values <- if(length(clusts) > 1) {
                    colMeans(Reduce(rbind, lapply(data$core.list[clusts], 
                                                   function(df) df$value)))
                 } else {
                   data$core.list[[clusts]]$value
                 }
  mean.core <- data$core.list[[clusts[1]]]
  mean.core$value <- mean.values
  
  
  outfpath <- file.path(dirPath, folder, outfname)
  dir.create(dirname(outfpath), recursive=TRUE)
  
  sp1 <- ggplot(K.mgd.molten, 
                aes(x = factor(sample, levels = unique(meta.df$condition)),
                    y = value)) +
         geom_jitter(aes(colour=score, group=gene)) +
         scale_color_gradientn(colours=c('black', 'red2'),
         limits = c(min.score, 1.0)) +
         geom_jitter(data=mean.core, aes(sample,value, group=cluster),
                     color='yellow', inherit.aes=F) +
         xlab('Time') +
         ylab('Expression') +
         labs(title=paste0(addTitle, " Cluster ", paste(clusts, collapse = '_'), 
                           " Expression by Group", color = 'Score'))
  svg(outfpath)
    print(sp1)
  dev.off()
  
  dfs2xlsx(withNames(addTitle, K.mgd.molten), gsub(".svg", ".xlsx", outfpath))
}

source("~/Dropbox/R/central-scripts/utils.interleave.strings.R")

fpath2jitter.genes <- function(fpath, genes, outfname, addTitle="Selected Genes", folder="") {

  dirPath <- splitStrToVec(fpath, "/heatmap/")[1]

  data <- readRDS(fpath)
  metaPath <- dir(file.path(dirPath, "meta"), full.name=T)

  meta.df <- read.table(metaPath, sep = '\t', header = T, stringsAsFactors = F)
  meta.df$condition <- factor(meta.df$condition,
                              levels = unique(meta.df$condition))
  denom <- function(meta.df) as.character(meta.df$condition[meta.df$testing == "denom"][1])
  num <- function(meta.df) as.character(meta.df$condition[meta.df$testing == "num"][1])
  core.name <- function(meta.df) paste0(num(meta.df), "-vs-", denom(meta.df), collapse = '')

  K.mgd.molten <- Reduce(rbind, data$Kmolten.list)
  K.mgd.molten <- K.mgd.molten[K.mgd.molten$gene %in% genes, ]

  if (folder == "") {
    outfpath <- file.path(dirPath, "heatmap", outfname)
  } else {
    outfpath <- file.path(dirPath, "heatmap", folder, outfname)
  }
  dir.create(dirname(outfpath), recursive=TRUE)

  min.score <- min(Reduce(rbind, data$Kmolten.list)$score)
  
  mean.values <- sapply(split(K.mgd.molten, K.mgd.molten$sample), 
                              function(df) mean(df$value))
  mean.core <- data$core.list[[1]] # just to get the df format
  mean.core$value <- mean.values[as.character(mean.core$sample)]   # fill with real mean values
  
  sp1 <- ggplot(K.mgd.molten, 
                aes(x = factor(sample, levels = unique(meta.df$condition)),
                    y = value)) +
         geom_jitter(aes(colour=score, group=gene)) +
         scale_color_gradientn(colours=c('black', 'red2'),
         limits = c(min.score, 1.0)) +
         geom_jitter(data=mean.core, aes(sample,value, group=cluster),
                     color='yellow', inherit.aes=F) +
         xlab('Time') +
         ylab('Expression') +
         labs(title=paste0("Cluster ", addTitle, 
                           " Expression by Group", color = 'Score'))
  svg(outfpath)
    print(sp1)
  dev.off()
  
  dfs2xlsx(withNames(addTitle, K.mgd.molten), gsub(".svg", ".xlsx", outfpath))
}


br.k.path <- '/media/josephus/Elements/DEanalysis/rev2019-jitters-190807-br/k2-vs-wtn/heatmap/all/scaledata.Kmolten.core.score.list.rev2019-jitters-190807-br.190807170032.all.rds'
br.k.xl.path <- '/media/josephus/Elements/DEanalysis/rev2019-jitters-190807-br/k2-vs-wtn/heatmap/all/k16_khmap_rev2019-jitters-190807-br_190807170031_UO_0.05_1all.xlsx'
eo.k.path <- '/media/josephus/Elements/DEanalysis/rev2019-jitters-190807-eo/k2-vs-wtn/heatmap/all/scaledata.Kmolten.core.score.list.rev2019-jitters-190807-eo.190807162235.all.rds'
eo.k.xl.path <- '/media/josephus/Elements/DEanalysis/rev2019-jitters-190807-eo/k2-vs-wtn/heatmap/all/k16_khmap_rev2019-jitters-190807-eo_190807162234_UO_0.05_1all.xlsx'

fpath.1 <- br.k.path
fpath.2 <- eo.k.path

##################################
# new jitter plots - new clusters              
##################################

br.UP.clusters <- c(3, 8, 14)
fpath2jitter.clusters(fpath.1, br.UP.clusters, paste0("Bra_UP_", paste(br.UP.clusters, collapse="_"), "_jitter.svg"), addTitle="Bra_UP_new", folder="UP_new")
br.UP.clusters <- c(3, 8)
fpath2jitter.clusters(fpath.1, br.UP.clusters, paste0("Bra_UP_", paste(br.UP.clusters, collapse="_"), "_jitter.svg"), addTitle="Bra_UP_new", folder="UP_new")
br.UP.clusters <- c(3, 14)
fpath2jitter.clusters(fpath.1, br.UP.clusters, paste0("Bra_UP_", paste(br.UP.clusters, collapse="_"), "_jitter.svg"), addTitle="Bra_UP_new", folder="UP_new")
br.UP.clusters <- c(8, 14)
fpath2jitter.clusters(fpath.1, br.UP.clusters, paste0("Bra_UP_", paste(br.UP.clusters, collapse="_"), "_jitter.svg"), addTitle="Bra_UP_new", folder="UP_new")
br.UP.clusters <- c(14)
fpath2jitter.clusters(fpath.1, br.UP.clusters, paste0("Bra_UP_", paste(br.UP.clusters, collapse="_"), "_jitter.svg"), addTitle="Bra_UP_new", folder="UP_new")
br.UP.clusters <- c(3)
fpath2jitter.clusters(fpath.1, br.UP.clusters, paste0("Bra_UP_", paste(br.UP.clusters, collapse="_"), "_jitter.svg"), addTitle="Bra_UP_new", folder="UP_new")

br.UP.clusters <- c(3, 6, 14)
fpath2jitter.clusters(fpath.1, br.UP.clusters, paste0("Bra_UP_", paste(br.UP.clusters, collapse="_"), "_jitter.svg"), addTitle="Bra_UP_new", folder="UP_new")
br.UP.clusters <- c(6)
fpath2jitter.clusters(fpath.1, br.UP.clusters, paste0("Bra_UP_", paste(br.UP.clusters, collapse="_"), "_jitter.svg"), addTitle="Bra_UP_new", folder="UP_new")
br.UP.clusters <- c(3, 6)
fpath2jitter.clusters(fpath.1, br.UP.clusters, paste0("Bra_UP_", paste(br.UP.clusters, collapse="_"), "_jitter.svg"), addTitle="Bra_UP_new", folder="UP_new")
br.UP.clusters <- c(3, 14)
fpath2jitter.clusters(fpath.1, br.UP.clusters, paste0("Bra_UP_", paste(br.UP.clusters, collapse="_"), "_jitter.svg"), addTitle="Bra_UP_new", folder="UP_new")
br.UP.clusters <- c(6, 14)
fpath2jitter.clusters(fpath.1, br.UP.clusters, paste0("Bra_UP_", paste(br.UP.clusters, collapse="_"), "_jitter.svg"), addTitle="Bra_UP_new", folder="UP_new")

br.DOWN.clusters <- c(10, 8, 11)
fpath2jitter.clusters(fpath.1, br.DOWN.clusters, paste0("Bra_DOWN_", paste(br.DOWN.clusters, collapse="_"), "_jitter.svg"), addTitle="Bra_DOWN_new", folder="DOWN_new")
br.DOWN.clusters <- c(10, 8)
fpath2jitter.clusters(fpath.1, br.DOWN.clusters, paste0("Bra_DOWN_", paste(br.DOWN.clusters, collapse="_"), "_jitter.svg"), addTitle="Bra_DOWN_new", folder="DOWN_new")
br.DOWN.clusters <- c(8, 11)
fpath2jitter.clusters(fpath.1, br.DOWN.clusters, paste0("Bra_DOWN_", paste(br.DOWN.clusters, collapse="_"), "_jitter.svg"), addTitle="Bra_DOWN_new", folder="DOWN_new")
br.DOWN.clusters <- c(10, 11)
fpath2jitter.clusters(fpath.1, br.DOWN.clusters, paste0("Bra_DOWN_", paste(br.DOWN.clusters, collapse="_"), "_jitter.svg"), addTitle="Bra_DOWN_new", folder="DOWN_new")

eo.UP.clusters <- c(6, 8, 9)
fpath2jitter.clusters(fpath.2, eo.UP.clusters, paste0("Eo_UP_", paste(eo.UP.clusters, collapse="_"), "_jitter.svg"), addTitle="Eo_UP_new", folder="UP_new")
eo.UP.clusters <- c(6, 8)
fpath2jitter.clusters(fpath.2, eo.UP.clusters, paste0("Eo_UP_", paste(eo.UP.clusters, collapse="_"), "_jitter.svg"), addTitle="Eo_UP_new", folder="UP_new")
eo.UP.clusters <- c(8, 9)
fpath2jitter.clusters(fpath.2, eo.UP.clusters, paste0("Eo_UP_", paste(eo.UP.clusters, collapse="_"), "_jitter.svg"), addTitle="Eo_UP_new", folder="UP_new")
eo.UP.clusters <- c(6, 9)
fpath2jitter.clusters(fpath.2, eo.UP.clusters, paste0("Eo_UP_", paste(eo.UP.clusters, collapse="_"), "_jitter.svg"), addTitle="Eo_UP_new", folder="UP_new")
eo.UP.clusters <- c(6)
fpath2jitter.clusters(fpath.2, eo.UP.clusters, paste0("Eo_UP_", paste(eo.UP.clusters, collapse="_"), "_jitter.svg"), addTitle="Eo_UP_new", folder="UP_new")
eo.UP.clusters <- c(9)
fpath2jitter.clusters(fpath.2, eo.UP.clusters, paste0("Eo_UP_", paste(eo.UP.clusters, collapse="_"), "_jitter.svg"), addTitle="Eo_UP_new", folder="UP_new")

eo.UP.clusters <- c(2, 6, 9, 11)
fpath2jitter.clusters(fpath.2, eo.UP.clusters, paste0("Eo_UP_", paste(eo.UP.clusters, collapse="_"), "_jitter.svg"), addTitle="Eo_UP_new", folder="UP_new")
eo.UP.clusters <- c(2, 6)
fpath2jitter.clusters(fpath.2, eo.UP.clusters, paste0("Eo_UP_", paste(eo.UP.clusters, collapse="_"), "_jitter.svg"), addTitle="Eo_UP_new", folder="UP_new")
eo.UP.clusters <- c(2, 9)
fpath2jitter.clusters(fpath.2, eo.UP.clusters, paste0("Eo_UP_", paste(eo.UP.clusters, collapse="_"), "_jitter.svg"), addTitle="Eo_UP_new", folder="UP_new")
eo.UP.clusters <- c(2, 11)
fpath2jitter.clusters(fpath.2, eo.UP.clusters, paste0("Eo_UP_", paste(eo.UP.clusters, collapse="_"), "_jitter.svg"), addTitle="Eo_UP_new", folder="UP_new")
eo.UP.clusters <- c(9, 11)
fpath2jitter.clusters(fpath.2, eo.UP.clusters, paste0("Eo_UP_", paste(eo.UP.clusters, collapse="_"), "_jitter.svg"), addTitle="Eo_UP_new", folder="UP_new")
eo.UP.clusters <- c(6, 11)
fpath2jitter.clusters(fpath.2, eo.UP.clusters, paste0("Eo_UP_", paste(eo.UP.clusters, collapse="_"), "_jitter.svg"), addTitle="Eo_UP_new", folder="UP_new")
eo.UP.clusters <- c(6, 9, 11)
fpath2jitter.clusters(fpath.2, eo.UP.clusters, paste0("Eo_UP_", paste(eo.UP.clusters, collapse="_"), "_jitter.svg"), addTitle="Eo_UP_new", folder="UP_new")
eo.UP.clusters <- c(2, 6, 11)
fpath2jitter.clusters(fpath.2, eo.UP.clusters, paste0("Eo_UP_", paste(eo.UP.clusters, collapse="_"), "_jitter.svg"), addTitle="Eo_UP_new", folder="UP_new")
eo.UP.clusters <- c(2, 6, 9)
fpath2jitter.clusters(fpath.2, eo.UP.clusters, paste0("Eo_UP_", paste(eo.UP.clusters, collapse="_"), "_jitter.svg"), addTitle="Eo_UP_new", folder="UP_new")
eo.UP.clusters <- c(2)
fpath2jitter.clusters(fpath.2, eo.UP.clusters, paste0("Eo_UP_", paste(eo.UP.clusters, collapse="_"), "_jitter.svg"), addTitle="Eo_UP_new", folder="UP_new")
eo.UP.clusters <- c(11)
fpath2jitter.clusters(fpath.2, eo.UP.clusters, paste0("Eo_UP_", paste(eo.UP.clusters, collapse="_"), "_jitter.svg"), addTitle="Eo_UP_new", folder="UP_new")

eo.DOWN.clusters <- c(8, 10, 13, 15)
fpath2jitter.clusters(fpath.2, eo.DOWN.clusters, paste0("Eo_DOWN_", paste(eo.DOWN.clusters, collapse="_"), "_jitter.svg"), addTitle="Eo_DOWN_new", folder="DOWN_new")
eo.DOWN.clusters <- c(8, 10, 15)
fpath2jitter.clusters(fpath.2, eo.DOWN.clusters, paste0("Eo_DOWN_", paste(eo.DOWN.clusters, collapse="_"), "_jitter.svg"), addTitle="Eo_DOWN_new", folder="DOWN_new")
eo.DOWN.clusters <- c(8, 10)
fpath2jitter.clusters(fpath.2, eo.DOWN.clusters, paste0("Eo_DOWN_", paste(eo.DOWN.clusters, collapse="_"), "_jitter.svg"), addTitle="Eo_DOWN_new", folder="DOWN_new")
eo.DOWN.clusters <- c(8, 13)
fpath2jitter.clusters(fpath.2, eo.DOWN.clusters, paste0("Eo_DOWN_", paste(eo.DOWN.clusters, collapse="_"), "_jitter.svg"), addTitle="Eo_DOWN_new", folder="DOWN_new")
eo.DOWN.clusters <- c(8, 15)
fpath2jitter.clusters(fpath.2, eo.DOWN.clusters, paste0("Eo_DOWN_", paste(eo.DOWN.clusters, collapse="_"), "_jitter.svg"), addTitle="Eo_DOWN_new", folder="DOWN_new")
eo.DOWN.clusters <- c(10, 13)
fpath2jitter.clusters(fpath.2, eo.DOWN.clusters, paste0("Eo_DOWN_", paste(eo.DOWN.clusters, collapse="_"), "_jitter.svg"), addTitle="Eo_DOWN_new", folder="DOWN_new")
eo.DOWN.clusters <- c(10, 15)
fpath2jitter.clusters(fpath.2, eo.DOWN.clusters, paste0("Eo_DOWN_", paste(eo.DOWN.clusters, collapse="_"), "_jitter.svg"), addTitle="Eo_DOWN_new", folder="DOWN_new")
eo.DOWN.clusters <- c(10, 13, 15)
fpath2jitter.clusters(fpath.2, eo.DOWN.clusters, paste0("Eo_DOWN_", paste(eo.DOWN.clusters, collapse="_"), "_jitter.svg"), addTitle="Eo_DOWN_new", folder="DOWN_new")
eo.DOWN.clusters <- c(13, 15)
fpath2jitter.clusters(fpath.2, eo.DOWN.clusters, paste0("Eo_DOWN_", paste(eo.DOWN.clusters, collapse="_"), "_jitter.svg"), addTitle="Eo_DOWN_new", folder="DOWN_new")




##################################################
# GO analysis (right panels)
##################################################


####################################################################
# source 
####################################################################


source(file.path("/home/josephus/Dropbox/a_myRNAseq/rdownstream", "190801.rdownstream.central.R"))
script_dir <- "/home/josephus/Dropbox/a_myRNAseq/rdownstream"
out_dir    <- "/media/josephus/Elements/DEanalysis"
meta_dir   <- "/media/josephus/Elements/metas"
dir.create(out_dir, recursive=TRUE)

#####################################################################


fpath.1 <- '/media/josephus/Elements/DEanalysis/rev2019-jitters-190807-br/k2-vs-wtn/heatmap/all/scaledata.Kmolten.core.score.list.rev2019-jitters-190807-br.190807141553.all.rds' 
fpath.2 <- '/media/josephus/Elements/DEanalysis/rev2019-jitters-190807-eo/k2-vs-wtn/heatmap/all/scaledata.Kmolten.core.score.list.rev2019-jitters-190807-eo.190807143339.all.rds'


#######################################
# load packages
#######################################

require(clusterProfiler)
require(org.Mm.eg.db)
require(xlsx2dfs)
require(BiocParallel)
mcp <- MulticoreParam(workers=3)

#######################################
# ALIAS to ENTREZID data frame rownames
#######################################

alias2entrezid.genenames <- function(genenames) {
  "Convert using biomart ALIAS the equivalents of ENTREZID."
  "It is not a 1:1 mapping! So number of output != number of genenames!"
  genes2eg <- bitr(genenames,
                   fromType = "ALIAS",
                   toType   = "ENTREZID",
                   OrgDb    = "org.Mm.eg.db")
  genes2eg.list <- split(genes2eg, genes2eg$ALIAS)
  genes2eg.list <- lapply(genes2eg.list, function(df) df$ENTREZID)
  # keys are ALIAS and values are vector of entrezids
  entrezids.list <- genes2eg.list[genenames]
  entrezids <- unlist(entrezids.list) # flattening
  not.found.genes <- sort(setdiff(genenames, names(genes2eg.list)))
  attr(entrezids, "not.found.genes") <- not.found.genes
  found.genes     <- names(genes2eg.list)
  attr(entrezids, "found.genes") <- found.genes
  entrezids
}

#######################################
# found and not founds
#######################################

# this works also with our list of genenames

get.founds <- function(egs) {
  data.frame(alias=attr(egs, "found.genes"))
}

get.not.founds <- function(egs) {
  data.frame(alias=attr(egs, "not.found.genes"))
}

#######################################
# GO overrepresentation
# simplify is discussed here:
#######################################
# https://github.com/GuangchuangYu/clusterProfiler/issues/28

EnrichGO <- function(genenames, ont, simplifyp=FALSE) {
  # Return enchrichment for a dataframe with entrezid rownames or a vector of entrezids
  res <- enrichGO(gene          = genenames,
                  OrgDb         = org.Mm.eg.db,
                  keyType       = "ALIAS",
                  ont           = ont,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  if (simplifyp) {
      res <- simplify(res,
                      cutoff = 0.7,
                      by = "p.adjust",
                      select_fun = min)
  }
  res
}
dfEnrichGO <- EnrichGO # x is a df

#######################################################
# GO BP,CC,MF entire analysis with or without selection
# CNTS removed, since not used in GO
#######################################################


do_GO_enrichment.genenames <- function(genenames,
                                       outfpath,
                                       simplifyp=FALSE) { # when genenames, then no GSEA
  outBase <- dirname(outfpath)
  fileName   <- basename(outfpath)
  decorateFname <- function(textPiece, fname=fileName) gsub(".xlsx", 
                                                   paste0(textPiece, "_GO.xlsx"), 
                                                   fname)
  overrepFname_ont <- decorateFname("over_GO_")
  gseaFname_ont <- decorateFname("gsea_GO_")
  notFoundFnameDE <- decorateFname("not_found_")
  FoundFnameDE <- decorateFname("found_")
  
  if (!dir.exists(outBase)) {
    dir.create(outBase, recursive = TRUE)
  }
  
  # Go enrichment from genenames
  df.bp <- EnrichGO(genenames, "BP", simplifyp=simplifyp)
  df.cc <- EnrichGO(genenames, "CC", simplifyp=simplifyp)
  df.mf <- EnrichGO(genenames, "MF", simplifyp=simplifyp)
  
  # collect the founds
  genes.found <- rownames(df.bp)
  founds.df   <- data.frame(alias=genes.found)
  # founds.df <- get.founds(egs)
  not.found.genes <- sort(setdiff(genenames, genes.found))
  not.founds.df <- data.frame(alias=not.found.genes)
  dfs2xlsx(withNames("BP", df.bp,
                     "CC", df.cc,
                     "MF", df.mf,
                     "found.genes", founds.df,
                     "not.fount.genes", not.founds.df),
           outfpath)
}

do.GO.xlsx <- function(xlsx) {
  genes <- unique(xlsx2dfs::xlsx2dfs(xlsx)[[1]]$gene)
  do_GO_enrichment.genenames(genes,
                       gsub(".xlsx", "_GO.xlsx", xlsx),
                       simplify=simplifyp)
  gc()
}


do.go.analyses <- function(Kmolten.rds.fpath, parallelize.p = F, simplifyp = T) {
  
  ####################################################################
  # inferred paths
  ####################################################################
  path <- dirname(dirname((Kmolten.rds.fpath)))
  jitter.xlsx.paths <- dir(path, pattern="_jitter.xlsx", recursive=T, full.names=T)
  
  if (parallelize.p) {
    bplapply(jitter.xlsx.paths, FUN=do.GO.xlsx,
      BPPARAM=mcp)
  } else {
    for (xlsx in jitter.xlsx.paths) {
      genes <- unique(xlsx2dfs(xlsx)[[1]]$gene)
      do_GO_enrichment.genenames(genes,
                                 gsub(".xlsx", "_GO.xlsx", xlsx),
                                 simplify=simplifyp)
    }
  }
}


do.go.analyses(fpath.1)
do.go.analyses(fpath.2)

analyze.go <- function(xlsx.fpath, parallelize.p = F, simplifyp = T) {
  
  ####################################################################
  # inferred paths
  ####################################################################
  genes <- unique(xlsx2dfs(xlsx.fpath)[[1]]$gene)
  do_GO_enrichment.genenames(genes,
                             gsub(".xlsx", "_GO_extended.xlsx", xlsx),
                             simplify=simplifyp)
}

source(file.path("/home/josephus/Dropbox/a_myRNAseq/rdownstream", "190809.rdownstream.central.R"))

script_dir <- "/home/josephus/Dropbox/a_myRNAseq/rdownstream"
out_dir    <- "/media/josephus/Elements/DEanalysis"
meta_dir   <- "/media/josephus/Elements/metas"
dir.create(out_dir, recursive=TRUE)

# GO 20190812
xlsx.paths <- c("/media/josephus/Elements/DEanalysis/rev2019-jitters-190807-eo/k2-vs-wtn/heatmap/UP_new/Eo_UP_9_11_jitter.xlsx",
"/media/josephus/Elements/DEanalysis/rev2019-jitters-190807-eo/k2-vs-wtn/heatmap/DOWN_new/Eo_DOWN_8_10_13_15_jitter.xlsx",
"/media/josephus/Elements/DEanalysis/rev2019-jitters-190807-br/k2-vs-wtn/heatmap/UP_new/Bra_UP_3_14_jitter.xlsx",
"/media/josephus/Elements/DEanalysis/rev2019-jitters-190807-br/k2-vs-wtn/heatmap/DOWN_new/Bra_DOWN_10_11_jitter.xlsx")

for (xlsx in xlsx.paths) {
  analyze.go(xlsx, simplifyp=F)
}


