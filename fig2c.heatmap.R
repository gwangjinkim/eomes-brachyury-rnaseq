

####################################################################
# source 
####################################################################
# this loads the functions necessary

source("/home/josephus/Dropbox/a_myRNAseq/rdownstream/181031.rdownstream.central.R")
# coloring changed of heatmap!

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

# ####################################################################
# # Use functions for analysis
# ####################################################################
# this big for loop handles the runs
# I change 1:2 to the n needed
# like 1:7 etc
# if you want to jump over other - just give the sequence
# for (i in c(1, 4, 7)) { ... }

for (i in 1:5) {
  
  #########################################
  # set variables here which are common to all runs
  #########################################
  
  ## i <- 2
  indirpath <- ""
  lFC <- 2.5
  
  #########################################
  # define variables specific for each run
  #########################################
  
  #   if (i == 1) {
  #     dataname <- "brGFP" # the data name is in nearly ALL output files
  #     outdirpath <- "/media/josephus/Elements/DEanalysis/eobrGFP-181102/brGFP-181102"
  #     # outdir folder into which all subfolders will be create
  #     # all folders or files will be generated automatically if not existent
  #     metapath   <- "/media/josephus/Elements/metas/meta-brGFPD-vs-d-181102.txt"
  #     # path to the meta path - nowadays I collect them all in the folder "metas"
  #     # and the meta folder will be created in the results as subfolder
  #     # for documentation reasons
  #   } # so put everything into such a if claus
  #     # and don't forget to change the `(i == n)` to the correct number!
  
  # here, I have 2 analysis
  if (i == 1) {
    dataname <- "k2-wtn-epi-Bao2018-up-16-no-es"
    outdirpath <- "/media/josephus/Elements/DEanalysis/JT-cells-hybrid-epi-up-16-no-es/k2-wtn-epi-Bao2018"
    metapath   <- "/media/josephus/Elements/metas/meta-k2-vs-wt-nen-episc-Bao2018.txt"
  }

  if (i == 2) {
    dataname <- "k2-wtn-epi-Senft2018-up-16-no-es"
    outdirpath <- "/media/josephus/Elements/DEanalysis/JT-cells-hybrid-epi-up-16-no-es/k2-wtn-epi-Senft2018"
    metapath   <- "/media/josephus/Elements/metas/meta-k2-vs-wt-nen-episc-Senft2018.txt"
  }
  
  if (i == 3) {
    dataname <- "k2-wtn-epi-Liu_epi-up-16-no-es"
    outdirpath <- "/media/josephus/Elements/DEanalysis/JT-cells-hybrid-epi-up-16-no-es/k2-wtn-epi-Liu_epi"
    metapath   <- "/media/josephus/Elements/metas/meta-k2-vs-wt-nen-episc-Liu2018_epi.txt"
  }
  
  if (i == 4) {
    dataname <- "k2-wtn-epi-Liu_epi-up-16-no-es"
    outdirpath <- "/media/josephus/Elements/DEanalysis/JT-cells-hybrid-epi-up-16-no-es/k2-wtn-epi-Liu_5_5"
    metapath   <- "/media/josephus/Elements/metas/meta-k2-vs-wt-nen-episc-Liu2018_epi5_5.txt"
  }

  if (i == 5) {
    dataname <- "k2-wtn-epi-Liu_epi-up-16-no-es"
    outdirpath <- "/media/josephus/Elements/DEanalysis/JT-cells-hybrid-epi-up-16-no-es/k2-wtn-epi-Liu_6_5"
    metapath   <- "/media/josephus/Elements/metas/meta-k2-vs-wt-nen-episc-Liu2018_epi6_5.txt"
  }

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
  
  print.cnts.DE.sortings(cnts.avg.nrm, 
                         resSig, 
                         file.path(paste0("DE-cnts-avg-sig-", dataname, "_", 
                                          core.name(meta.df), "_",
                                          time.now(), "-", alpha, "-", lFC, ".xlsx")), 
                         dir = DE.outdir)
  
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
  scaledata.Kmolten.core.list <- meta2heatmap(meta.df, cnts.avg.nrm, resSig,
                                              outdirpath = hm.all.outdir, selected.genes = "up", name.add = "all",
                                              dataname = dataname,
                                              k = 16, printp=TRUE,
                                              alpha=alpha, lFC=lFC, filterp=TRUE,
                                              xlsxp=TRUE, csvp=FALSE, tsvp=FALSE)
  
  
#   #######################################################################
#   # Create GO
#   #######################################################################
#   # GO analysis for GO folder - actually no changes needed
#   # but if you don't want it - comment this section out
#   
#   # inDir is outdirpath
#   # outBase is GO.outdir
#   last <- function(l) l[length(l)]
#   DEfpath <- last(dir(DE.outdir, pattern = "DE_sig", full.names = TRUE))
# 
#   
#   do_all_GO_enrichment(DEfpath, GO.outdir)
#   
#   #######################################################################
#   # KEGG
#   #######################################################################
#   # KEGG analysis for GO folder - actually no changes needed
#   # but if you don't want it - comment this section out
#   
#   # DEfpath is like above
#   # outBase <- pathway.outdir
#   do_KEGG_enrichment(DEfpath, pathway.outdir)
# 
#   #######################################################################
#   # gskb
#   #######################################################################
#   # gskb analysis for mouse - actually no changes needed
#   # but if you don't want it (quicker analysis) - comment this section out
#   
#   # DEfpath is like above
#   # outBase <- gskb.outdir
#   do_all_gskb_enrichments(DEfpath,
#                           gskb.outdir,
#                           all_ez_grouped)
#   
#   #######################################################################
#   # nicer MDS/PCA plots using plotly
#   #######################################################################
#   # creates nicer non-interactive PCA plots (for publication)
#   # - no changes needed
#   # - but if PCA plot not important - comment this section out

  
#   # request DE folder for significant genes path in previous results
#   cntsDEfpath <- last(dir(DE.outdir, pattern = "DE-cnts-sig", full.names = TRUE))
#   
#   meta2plyMDSplot(meta.df, 
#                   DESeq2.obj.disp, 
#                   PCA.outdirpath,
#                   dataname, 
#                   top=500, 
#                   launch = FALSE)
#   
#   meta2plyMDSplotDE(cntsDEfpath = cntsDEfpath, 
#                     outdirpath = PCA.outdirpath,
#                     dataname = dataname, 
#                     corename = core.name(meta.df), 
#                     top=500, 
#                     launch = FALSE)

  
  #######################################################################
  # save session image for later investigation of the run
  #######################################################################
  fpath <- file.path(outdirpath, paste0("run", time.now(), ".RData"))
  save.image(file = fpath)
}


heatmap.path <- '/media/josephus/Elements/DEanalysis/JT-cells-hybrid-epi-up-16-no-es/k2-wtn-epi-Liu_epi/k2-vs-wtn/heatmap/all/k16_khmap_k2-wtn-epi-Liu_epi-up-16-no-es_181207173157_UO_0.05_2.5all.xlsx'

unordered2orderedheatmap(heatmap.path, c(c(7), c(5, 1, 15, 3, 9), c(13, 10, 14, 2, 12), c(16, 11, 6, 8), c(4)), outdirpath=dirname(heatmap.path),
                                     gaps.after.blocks=c(1, 6, 11, 15, 16), alpha=alpha, lFC=lFC, dataname="k2-wtn-epi-Liu_epi-up-16-no-es")
                                     





heatmap.path <- '/media/josephus/Elements/DEanalysis/JT-cells-hybrid-epi-up-16-no-es/k2-wtn-epi-Bao2018/k2-vs-wtn/heatmap/all/k16_khmap_k2-wtn-epi-Bao2018-up-16-no-es_181207172906_UO_0.05_2.5all.xlsx'

unordered2orderedheatmap(heatmap.path, c(c(13), c(12), c(1, 9, 5, 4, 15), c(10, 2, 7, 8, 11), c(6, 14, 16, 3)), outdirpath=dirname(heatmap.path),
                                     gaps.after.blocks=c(1, 2, 7, 12, 16), alpha=alpha, lFC=lFC, dataname="k2-wtn-epi-Bao2018_epi-up-16-no-es")





















