#######################################################
# Script for preparing RPKM values for Extended Data 7a
#######################################################
# ./a_myRNAseq/rdownstream/190508.rpkm.heatmap181031.R

require(xlsx2dfs)
require(scater)

gene.lengths.fpath <- "/home/josephus/count/tools/sym_gene_lengths.txt"


print.rpkm <- function(cnts.fpath,
                       gene.lengths.fpath,
                       raw.count.sheet = "raw-counts",
                       out.fpath = "") {
  # Read raw counts (cnts.fpath), calculate rpkm and print to out.fpath (an excel file).
  # (If no out.fpath given, print into same folder with modified name "raw-counts" -> "rpkm-raw").
  cnts <- xlsx2dfs(cnts.fpath)[[raw.count.sheet]]
  sce <- SingleCellExperiment(assays = list(counts = as.matrix(cnts)))
  sce <- normalize(sce)
  
  read_gene_lengths <- function(gene.lengths.fpath) {
    gl <- read.delim(gene.lengths.fpath, 
                     sep = "\t", 
                     head = TRUE, 
                     comment.char="#")
    sym2length <- gl$Length
    names(sym2length) <- gl$Symbol
    sym2length
  }
  
  sym2len <- read_gene_lengths(gene.lengths.fpath)
  cnts.rpkm.len <- calculateFPKM(sce, 
                                 sym2len[rownames(cnts)], 
                                 use_size_factors=FALSE) # this calculates RPKM
  cnts.rpkm.len.cl <- na.omit(cnts.rpkm.len)
  if (out.fpath == "") {
    out.fpath <-gsub("raw-counts", "rpkm-raw", cnts.fpath)
  }
  dfs2xlsx(withNames("rpkm_raw_counts", cnts.rpkm.len.cl),
           out.fpath)
  cnts.rpkm.len.cl
}

# set paths

br.jitter.fpath <- '/media/josephus/Elements/DEanalysis/rev2019-jitters-190807-br/k2-vs-wtn/count-table/raw-counts-rev2019-jitters-190807-br-k2-vs-wtn.xlsx'
eo.jitter.fpath <- '/media/josephus/Elements/DEanalysis/rev2019-jitters-190807-eo/k2-vs-wtn/count-table/raw-counts-rev2019-jitters-190807-eo-k2-vs-wtn.xlsx'
all.jitter.fpath <- '/media/josephus/Elements/DEanalysis/rev2019-jitters-190807-all/k2-vs-wtn/count-table/raw-counts-rev2019-jitters-190807-all-k2-vs-wtn.xlsx' 
six2.snai1.fpath <- '/media/josephus/Elements/DEanalysis/rev2019-six2-snai1-dko-vs-wtn/dKO-vs-WT/count-table/raw-counts-rev2019-six2-snai1-dko-vs-wtn-dKO-vs-WT.xlsx'

# call the function for each path

rpkm1 <- print.rpkm(br.jitter.fpath,
                    gene.lengths.fpath)
rpkm2 <- print.rpkm(eo.jitter.fpath,
                    gene.lengths.fpath)
rpkm3 <- print.rpkm(six2.snai1.fpath,
                    gene.lengths.fpath)
rpkm4 <- print.rpkm(all.jitter.fpath,
                    gene.lengths.fpath)

