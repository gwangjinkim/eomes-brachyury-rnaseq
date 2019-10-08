
###################################################################
# Code for preparation of normalized count tables for PCA plot 
# in extended data 4 e
###################################################################

# content of meta file: '/media/josephus/Elements/metas/meta-wt-k2-eo-5Ad-br-eo-5AD-br-GFPAD.txt' 
# 
# sampleName	fileName	condition	testing
# wtn-1	/media/josephus/archive_big/kworkspace/count/JT-cells-hybrid-new/rsubread/grcm38/WTnewrep1-sym-fcount.tab	wtn	denom
# wtn-2	/media/josephus/archive_big/kworkspace/count/JT-cells-hybrid-new/rsubread/grcm38/WTnewrep2-sym-fcount.tab	wtn	denom
# wtn-3	/media/josephus/archive_big/kworkspace/count/JT-cells-hybrid-new/rsubread/grcm38/WTnewrep3-sym-fcount.tab	wtn	denom
# br-1	/media/josephus/archive_big/kworkspace/count/JT-cells-hybrid-new/rsubread/grcm38/br-1-sym-fcount.tab	br	
# br-2	/media/josephus/archive_big/kworkspace/count/JT-cells-hybrid-new/rsubread/grcm38/br-2-sym-fcount.tab	br	
# br-3	/media/josephus/archive_big/kworkspace/count/JT-cells-hybrid-new/rsubread/grcm38/br-3-sym-fcount.tab	br	
# eon-1	/media/josephus/archive_big/kworkspace/count/JT-cells-hybrid-new/rsubread/grcm38/EomesKOnewrep1-sym-fcount.tab	eon	
# eon-2	/media/josephus/archive_big/kworkspace/count/JT-cells-hybrid-new/rsubread/grcm38/EomesKOnewrep2-sym-fcount.tab	eon	
# eon-3	/media/josephus/archive_big/kworkspace/count/JT-cells-hybrid-new/rsubread/grcm38/EomesKOnewrep3-sym-fcount.tab	eon	
# k2-1	/media/josephus/archive_big/kworkspace/count/JT-cells-hybrid-new/rsubread/grcm38/k2-1-sym-fcount.tab	k2	num
# k2-2	/media/josephus/archive_big/kworkspace/count/JT-cells-hybrid-new/rsubread/grcm38/k2-2-sym-fcount.tab	k2	num
# k2-3	/media/josephus/archive_big/kworkspace/count/JT-cells-hybrid-new/rsubread/grcm38/k2-3-sym-fcount.tab	k2	num
# br-5AD-1	/media/josephus/archive_big/kworkspace/count/JT-cells-second-new/subread/BraV5Actdoxrep1-sym-fcount.tab	br5AD	
# br-5AD-2	/media/josephus/archive_big/kworkspace/count/JT-cells-second-new/subread/BraV5Actdoxrep2-sym-fcount.tab	br5AD	
# br-5AD-3	/media/josephus/archive_big/kworkspace/count/JT-cells-second-new/subread/BraV5Actdoxrep3-sym-fcount.tab	br5AD	
# eo-5AD-1	/media/josephus/archive_big/kworkspace/count/JT-cells-second-new/subread/EoV5newActdoxrep1-sym-fcount.tab	eo5AD	
# eo-5AD-2	/media/josephus/archive_big/kworkspace/count/JT-cells-second-new/subread/EoV5newActdoxrep2-sym-fcount.tab	eo5AD	
# eo-5AD-3	/media/josephus/archive_big/kworkspace/count/JT-cells-second-new/subread/EoV5newActdoxrep3-sym-fcount.tab	eo5AD	
# brGFPD-1	/media/josephus/archive_big/kworkspace/count/BraChIP2018/rsubread/grcm38/BraGFPdox1-sym-fcount.tab	brGFPD	
# brGFPD-2	/media/josephus/archive_big/kworkspace/count/BraChIP2018/rsubread/grcm38/BraGFPdox2-sym-fcount.tab	brGFPD	
# brGFPD-3	/media/josephus/archive_big/kworkspace/count/BraChIP2018/rsubread/grcm38/BraGFPdox3-sym-fcount.tab	brGFPD	
# eoGFPD-1	/media/josephus/archive_big/kworkspace/count/BraChIP2018/rsubread/grcm38/EoGFPdox1-sym-fcount.tab	eoGFPD	
# eoGFPD-2	/media/josephus/archive_big/kworkspace/count/BraChIP2018/rsubread/grcm38/EoGFPdox2-sym-fcount.tab	eoGFPD	
# eoGFPD-3	/media/josephus/archive_big/kworkspace/count/BraChIP2018/rsubread/grcm38/EoGFPdox3-sym-fcount.tab	eoGFPD	

####################################################################
# load functions
####################################################################

source("/home/josephus/Dropbox/a_myRNAseq/rdownstream/190507.rdownstream.central.R")

####################################################################
# set thresholds
####################################################################

alpha <- 0.05       # padj BH limit
FDR   <- 0.05       # false discovery rate, q-value
lFC   <- 1          # log2FC limit

for (i in 1:1) {
  
  ####################################################################
  # set basic variables
  ####################################################################
  
  indirpath <- ""
  dataname <- "wt-vs-k2-breo5AD-breoGFPAD"
  outdirpath <- "/media/josephus/Elements/DEanalysis/wt-vs-k2-breo5AD-breoGFPAD"
  metapath   <- "/media/josephus/Elements/metas/meta-wt-k2-eo-5Ad-br-eo-5AD-br-GFPAD.txt"

  ####################################################################
  # get meta file information
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
  # Create meta table (for documentation)
  #######################################################################
  
  write.table(meta.df,
              file = file.path(meta.outdir, 
                               basename(metapath)),
              sep = "\t", 
              row.names = FALSE, 
              quote = FALSE)
  
  #######################################################################
  # Create DESeq2 objects
  #######################################################################
  
  DESeq2.obj <- meta2DESeq2.obj(meta.df, 
                                indirpath)
  DESeq2.obj.disp <- meta2DESeq2.obj(meta.df, 
                                     indirpath, 
                                     normalized = TRUE)
  
  #######################################################################
  # Create the count tables (raw counts, normalized and averaged normalized counts)
  #######################################################################
  
  cnts.raw <- meta2cnts(meta.df, 
                        DESeq2.obj, 
                        outdirpath = cnts.outdir,
                        dataname = dataname,
                        printp = TRUE, 
                        normalized = FALSE, 
                        averaged = FALSE,
                        sheetName = "raw.all")
  cnts.nrm <- meta2cnts(meta.df, 
                        DESeq2.obj.disp, 
                        outdirpath = cnts.outdir,
                        dataname = dataname,
                        printp = TRUE, 
                        normalized = TRUE, 
                        averaged = FALSE,
                        sheetName = "normalized.all")
  cnts.avg.nrm <- meta2cnts(meta.df, 
                            DESeq2.obj.disp, 
                            outdirpath = cnts.outdir,
                            dataname = dataname,
                            printp = TRUE, 
                            normalized = TRUE, 
                            averaged = TRUE,
                            sheetName = "avg.normalized.all")

  #######################################################################
  # Create DE table
  #######################################################################
  
  res <- DEanalysis(meta.df, 
                    DESeq2.obj.disp, 
                    outdirpath = DE.outdir,
                    dataname = dataname,
                    printp = TRUE)
  resSig <- DEanalysis(meta.df, 
                       DESeq2.obj.disp, 
                       outdirpath = DE.outdir,
                       dataname = dataname,
                       printp = TRUE,
                       filterp = TRUE, 
                       alpha = alpha, 
                       lFC = lFC)
  # for reference, save name 
  cnts.DE.sig.fname <- paste0("DE-cnts-sig-", 
                              dataname, 
                              "_", 
                              core.name(meta.df), 
                              "_",
                              time.now(), 
                              "-", 
                              alpha, 
                              "-", 
                              lFC, 
                              ".xlsx")
  # print sorting by name and p values
  print.cnts.DE.sortings(cnts.nrm, 
                         resSig, 
                         cnts.DE.sig.fname, 
                         dir = DE.outdir)
  
  print.cnts.DE.sortings(cnts.avg.nrm$`nrm-counts-avg`, 
                         resSig, 
                         file.path(paste0("DE-cnts-avg-sig-", 
                                          dataname, 
                                          "_", 
                                          core.name(meta.df), 
                                          "_",
                                          time.now(), 
                                          "-", 
                                          alpha, 
                                          "-", 
                                          lFC, 
                                          ".xlsx")), 
                         dir = DE.outdir)
  
  #######################################################################
  # Create Scaled Avg DE sig table (averaged, significant DE genes table)
  #######################################################################
  scaled.avg.DE.sig <- scale.raw.counts.with.SD(cnts.DE.sig.fpath = file.path(DE.outdir, 
                                                                              cnts.DE.sig.fname), 
                                                meta.fpath = metapath, 
                                                out.fpath="") 
  
  #######################################################################
  # backup session image
  #######################################################################
  fpath <- file.path(outdirpath, 
                     paste0("run", 
                            time.now(), 
                            ".RData"))
  save.image(file = fpath)

}





