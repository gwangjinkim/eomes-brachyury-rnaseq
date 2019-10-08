computer <- "work"

if (computer == "ss") {
  script_dir <- "/media/daten/arnold/josephus/script"
  out_dir    <- "/media/daten/arnold/josephus/results/DEanalaysis"
  meta_dir   <- "/media/daten/arnold/josephus/metas"
} else {
  script_dir <- "/home/josephus/Dropbox/a_myRNAseq/rdownstream"
  out_dir    <- "/media/josephus/Elements/DEanalysis"
  meta_dir   <- "/media/josephus/Elements/metas"
}

source(file.path(script_dir, "190809.rdownstream.central.R"))

computer <- "work"

if (computer == "ss") {
  script_dir <- "/media/daten/arnold/josephus/script"
  out_dir    <- "/media/daten/arnold/josephus/results/DEanalaysis"
  meta_dir   <- "/media/daten/arnold/josephus/metas"
} else {
  script_dir <- "/home/josephus/Dropbox/a_myRNAseq/rdownstream"
  out_dir    <- "/media/josephus/Elements/DEanalysis"
  meta_dir   <- "/media/josephus/Elements/metas"
}
dir.create(out_dir, recursive=TRUE)








  ####################################################################
  # DE settings
  ####################################################################
  # set here tresholds
  # important are lFC (log2FoldChange for DE analysis 1 or 2.5 we use)

  alpha <- 0.05       # padj BH limit
  FDR   <- 0.05       # false discovery rate, q-value
  lFC   <- 1       # log2FC limit
  GO.alpha <- 0.01
  GO.FDR   <- 0.05
  GO.N     <- 30
  GO.cl.N  <- 10

data_names <- c("rev2019-overview-dko-vs-wt_mesp1-mgsn1-mixl1-foxa2-six2")
outdir_names <- data_names
meta_fnames <- paste("meta-", data_names, ".txt", sep="")

file.exists(file.path(meta_dir, meta_fnames))




for (i in 1:length(data_names)) {
  
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
  # this section reads in the meta file and makes it available
  # no need to change anything here
  
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
  # prepares subfolders - no changes to make
  
  outdirpath <- file.path(outdirpath, core.name(meta.df))
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
  gskb.outdir <- file.path(outdirpath, "gskb")
  PCA.outdirpath <- file.path(outdirpath, "pca")
  
  #####################################################################
  # ensure existence of output paths
  #####################################################################
  # creates  subfolders - no changes to make
  
  dir.create(path=outdirpath, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=meta.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=cnts.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=DE.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=plot.outdir, recursive = TRUE, showWarnings = FALSE)
  # dir.create(path=hm.final.outdir, recursive = TRUE, showWarnings = FALSE)
  # dir.create(path=hm.up.outdir, recursive = TRUE, showWarnings = FALSE)
  # dir.create(path=hm.down.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=hm.all.outdir, recursive = TRUE, showWarnings = FALSE)
  # dir.create(path=hm.up.list.outdir, recursive = TRUE, showWarnings = FALSE)
  # dir.create(path=hm.down.list.outdir, recursive = TRUE, showWarnings = FALSE)
  # dir.create(path=hm.all.list.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=GO.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=pathway.outdir, recursive = TRUE, showWarnings = FALSE)
  # dir.create(path=goi.outdir, recursive = TRUE, showWarnings = FALSE)
  # dir.create(path=lgenes.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=gskb.outdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(path=PCA.outdirpath, recursive = TRUE, showWarnings = FALSE)
  
  # https://guangchuangyu.github.io/2016/01/go-analysis-using-clusterprofiler/
  #######################################################################
  # Create meta table
  #######################################################################
  # creates the meta file in subfolder 'meta' - no changes to make
  
  write.table(meta.df,
              file = file.path(meta.outdir, basename(metapath)),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  #######################################################################
  # Create DESeq2.obj and DESeq2.obj.disp
  #######################################################################
  # does the DESeq2 analysis - no changes to make
  
  DESeq2.obj <- meta2DESeq2.obj(meta.df, indirpath)
  DESeq2.obj.disp <- meta2DESeq2.obj(meta.df, indirpath, normalized = TRUE)
  
  #######################################################################
  # Create the count tables
  #######################################################################
  # prints out count tables - no changes to make
  
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
  # Create plots
  #######################################################################
  # creates interactive plots in folder "glimma" - no changes to make
  
  meta2iMDSplot(meta.df, DESeq2.obj.disp, outdirpath,
                dataname, top=300, launch = FALSE)
  
  meta2iVolcano(meta.df, DESeq2.obj, DESeq2.obj.disp, outdirpath,
                dataname, alpha = alpha, lFC = lFC,
                launch = FALSE)
  
  #######################################################################
  # Create heatmap
  #######################################################################
  # creates heatmap with scaled values in 'heatmap' folder
  # you can change 'k =' to what you want. e.g. 12, 24, etc.
  # set the seed for reproducibility
  
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
  # GO analysis for GO folder - actually no changes needed
  # but if you don't want it - comment this section out
  
  # inDir is outdirpath
  # outBase is GO.outdir
  last <- function(l) l[length(l)]
  DEfpath <- last(dir(DE.outdir, pattern = "DE_sig", full.names = TRUE))
  
  try(do_all_GO_enrichment(DEfpath, GO.outdir, gsea.p=F))
  
  #######################################################################
  # KEGG
  #######################################################################
  # KEGG analysis for GO folder - actually no changes needed
  # but if you don't want it - comment this section out
  
  # DEfpath is like above
  # outBase <- pathway.outdir
  try(do_KEGG_enrichment(DEfpath, pathway.outdir, gsea.p=F))

  #######################################################################
  # gskb
  #######################################################################
  # gskb analysis for mouse - actually no changes needed
  # but if you don't want it (quicker analysis) - comment this section out
  
  # DEfpath is like above
  # outBase <- gskb.outdir
  try(do_all_gskb_enrichments(DEfpath,
                          gskb.outdir,
                          all_ez_grouped,
                          gsea.p=F))
  
  #######################################################################
  # nicer MDS/PCA plots using plotly
  #######################################################################
  # creates nicer non-interactive PCA plots (for publication)
  # - no changes needed
  # - but if PCA plot not important - comment this section out
  
  # request DE folder for significant genes path in previous results
  last <- function(l) l[length(l)]
  cntsDEfpath <- last(dir(DE.outdir, pattern = "DE-cnts-sig", full.names = TRUE))
  
#   try(meta2plyMDSplot(meta.df, 
#                   DESeq2.obj.disp, 
#                   PCA.outdirpath,
#                   dataname, 
#                   top=500, 
#                   launch = FALSE))
#   
#   try(meta2plyMDSplotDE(cntsDEfpath = cntsDEfpath, 
#                     outdirpath = PCA.outdirpath,
#                     dataname = dataname, 
#                     corename = core.name(meta.df), 
#                     top=500, 
#                     launch = FALSE))
  
  #######################################################################
  # save session image for later investigation of the run
  #######################################################################
  fpath <- file.path(outdirpath, paste0("run", time.now(), ".RData"))
  save.image(file = fpath)
  gc()
}

