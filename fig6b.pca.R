########################################################
# Preparation of normalized counts for PCA in Figure 6b
########################################################


script_dir <- "/home/josephus/Dropbox/a_myRNAseq/rdownstream"
source(file.path(script_dir, "190809.rdownstream.central.R"))

script_dir <- "/home/josephus/Dropbox/a_myRNAseq/rdownstream"
out_dir    <- "/media/josephus/Elements/DEanalysis"
meta_dir   <- "/media/josephus/Elements/metas"
dir.create(out_dir, recursive=TRUE)

####################################################################
# DE settings
####################################################################
# set here tresholds
# important are lFC (log2FoldChange for DE analysis 1 or 2.5 we use)

alpha <- 0.05       # padj BH limit
FDR   <- 0.05       # false discovery rate, q-value
lFC   <- 1          # log2FC limit

data_names <- c("rev2019-overview-dko-vs-wt_mesp1-mgsn1-mixl1-foxa2-six2")
outdir_names <- data_names
meta_fnames <- paste("meta-", data_names, ".txt", sep="")

for (i in 1:1) {
  
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
  
  #####################################################################
  # ensure existence of output paths
  #####################################################################
  
  dir.create(path=outdirpath, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=meta.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=cnts.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=DE.outdir, recursive = TRUE, showWarnings = FALSE)
  

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
  
  # rpkm counts
  raw.fpath   <- dir(cnts.outdir, pattern="raw-counts-", full.names=T)
  rpkm <- print.rpkm(raw.fpath,
                     gene.lengths.fpath)

  #######################################################################
  # Create DE table
  #######################################################################
  # prints out DE tables - no changes to make
  
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
  
  # scaled normalized counts
  
  DE.cnts.fpath <- dir(DE.outdir, pattern=paste0("DE-cnts-sig-", data_names[i], "_", 
                                                 core.name(meta.df)), full.names=T)[1]
  out.fpath <- gsub("-cnts-", "-scaled-avg-", DE.cnts.fpath)
  result <- scale.raw.counts.with.SD(DE.cnts.fpath, metapath, out.fpath)
  
  # robust scaled normalized counts
  
  out.fpath.r <- gsub("-cnts-", "-robustscaled-mad-avg-", DE.cnts.fpath)
  result.robust <- robustscale.counts.with.MAD(DE.cnts.fpath, metapath, out.fpath.r)
  # out.fpath.rsd <- gsub("-cnts-", "-robustscaled-sd-avg-", DE.cnts.fpath)
  # result.robust.sd <- robustscale.counts.with.SD(DE.cnts.fpath, metapath, out.fpath.rsd) # rownames length error
  
  # rpkm counts
  raw.fpath   <- dir(cnts.outdir, pattern="raw-counts-", full.names=T)
  rpkm <- print.rpkm(raw.fpath,
                     gene.lengths.fpath)

  #######################################################################
  # save session image for later investigation of the run
  #######################################################################
  fpath <- file.path(outdirpath, paste0("run", time.now(), ".RData"))
  save.image(file = fpath)
  gc()
}

