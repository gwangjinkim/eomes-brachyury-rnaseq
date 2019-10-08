########################################################
# Creation of heatmap for Figure 2c
########################################################

####################################################################
# load functions
####################################################################

source("/home/josephus/Dropbox/a_myRNAseq/rdownstream/181031.rdownstream.central.R")

####################################################################
# DE settings
####################################################################
# set here tresholds
# important are lFC (log2FoldChange for DE analysis 1 or 2.5 we use)

alpha <- 0.05       # padj BH limit
FDR   <- 0.05       # false discovery rate, q-value
lFC   <- 2.5          # log2FC limit

for (i in 1:1) {
  
  #########################################
  # set variables 
  #########################################
  
  indirpath <- ""
  dataname <- "k2-wtn-epi-Bao2018-up-16-no-es"
  outdirpath <- "/media/josephus/Elements/DEanalysis/JT-cells-hybrid-epi-up-16-no-es/k2-wtn-epi-Bao2018"
  metapath   <- "/media/josephus/Elements/metas/meta-k2-vs-wt-nen-episc-Bao2018.txt"
  

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
  # infer paths
  ####################################################################
  
  outdirpath <- file.path(outdirpath, core.name(meta.df))
  meta.outdir <- file.path(outdirpath, "meta")
  cnts.outdir <- file.path(outdirpath, "count-table")
  DE.outdir <- file.path(outdirpath, "DE-table")
  plot.outdir <- file.path(outdirpath, "glimma-plots")
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
  dir.create(path=plot.outdir, recursive = TRUE, showWarnings = FALSE)
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
  
  print.cnts.DE.sortings(cnts.avg.nrm, 
                         resSig, 
                         file.path(paste0("DE-cnts-avg-sig-", dataname, "_", 
                                          core.name(meta.df), "_",
                                          time.now(), "-", alpha, "-", lFC, ".xlsx")), 
                         dir = DE.outdir)
  
  #######################################################################
  # Create heatmap
  #######################################################################
  
  set.seed(123) # for reproducibility of heatmaps set seeds to a value
  scaledata.Kmolten.core.list <- meta2heatmap(meta.df, cnts.avg.nrm, resSig,
                                              outdirpath = hm.all.outdir, selected.genes = "up", name.add = "all",
                                              dataname = dataname,
                                              k = 16, printp=TRUE,
                                              alpha=alpha, lFC=lFC, filterp=TRUE,
                                              xlsxp=TRUE, csvp=FALSE, tsvp=FALSE)
  
  
  # #######################################################################
  # # Create GO
  # #######################################################################
  # # GO analysis for GO folder - actually no changes needed
  # # but if you don't want it - comment this section out
  # 
  # # inDir is outdirpath
  # # outBase is GO.outdir
  # last <- function(l) l[length(l)]
  # DEfpath <- last(dir(DE.outdir, pattern = "DE_sig", full.names = TRUE))
  # 
  # 
  # do_all_GO_enrichment(DEfpath, GO.outdir)
  
  #######################################################################
  # save session image for later investigation of the run
  #######################################################################
  fpath <- file.path(outdirpath, paste0("run", time.now(), ".RData"))
  save.image(file = fpath)
}

# reorder heatmap
# however, it is later reordered.

heatmap.path <- '/media/josephus/Elements/DEanalysis/JT-cells-hybrid-epi-up-16-no-es/k2-wtn-epi-Bao2018/k2-vs-wtn/heatmap/all/k16_khmap_k2-wtn-epi-Bao2018-up-16-no-es_181207172906_UO_0.05_2.5all.xlsx'

unordered2orderedheatmap(heatmap.path, c(c(13), c(12), c(1, 9, 5, 4, 15), c(10, 2, 7, 8, 11), c(6, 14, 16, 3)), outdirpath=dirname(heatmap.path),
                                     gaps.after.blocks=c(1, 2, 7, 12, 16), alpha=alpha, lFC=lFC, dataname="k2-wtn-epi-Bao2018_epi-up-16-no-es")





















