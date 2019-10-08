###################################################################
# Analyze GO terms for Figure 1f and 2b and 3b and 3c (embryo data)
###################################################################

rm(list=ls())

#######################################
# load packages
#######################################

source("~/Dropbox/R/central-scripts/io.xlsx.R")
require(clusterProfiler)
require(org.Mm.eg.db) 
# available keys: keytypes(org.Mm.eg.db)

#######################################
# ALIAS to ENTREZID data frame rownames
#######################################

alias2entrezid.df <- function(DE.df) {
  # converts rownames from ALIAS to ENTREZID
  genes <- rownames(DE.df)
  genes2eg <- bitr(genes,
                   fromType = "ALIAS",
                   toType   = "ENTREZID",
                   OrgDb    = "org.Mm.eg.db")
  
  # remove duplicated rows of aliases
  genes2eg.unique <- genes2eg[!duplicated(genes2eg$ALIAS), ]
  
  # remove duplicated rows of entrezids
  genes2eg.unique <- genes2eg[!duplicated(genes2eg$ENTREZID), ]
  df.new <- DE.df[genes2eg.unique$ALIAS, ]
  
  # again remove rows with duplicated rownames
  df.new <- df.new[!duplicated(rownames(df.new)), ]
  rownames(df.new) <- genes2eg.unique$ENTREZID
  
  # save/hide not found genes in the object # sorted looks nicer!
  not.found.genes <- sort(setdiff(genes, genes2eg.unique$ALIAS))
  attr(df.new, "not.found.genes") <- not.found.genes
  
  # save found genes in the object
  found.genes <- genes2eg.unique$ALIAS # don't sort to match with df.new!
  attr(df.new, "found.genes") <- found.genes
  
  df.new
}

#######################################
# GO overrepresentation
# simplify is discussed here:
#######################################
# https://github.com/GuangchuangYu/clusterProfiler/issues/28


dfEnrichGO <- function(eg.df, ont) {
  res <- enrichGO(gene          = rownames(eg.df),
                  OrgDb         = org.Mm.eg.db,
                  keyType       = "ENTREZID",
                  ont           = ont,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  res <- simplify(res,
                  cutoff = 0.7,
                  by = "p.adjust",
                  select_fun = min)
  res <- setReadable(res, org.Mm.eg.db, keytype = "ENTREZID")
  res
}

dfListEnrichGO <- function(eg.df.list, ont) {
  lapply(eg.df.list, function(df) dfEnrichGO(df, ont))
}

dfListNotFounds <- function(eg.df.list) {
  lapply(eg.df.list, function(df) data.frame(alias = attr(df, "not.found.genes")))
}

dfListFounds    <- function(eg.df.list) {
  names.list <- lapply(eg.df.list, function(df) attr(df, "found.genes"))
  Map(f = function(df, names.vec) {df$alias <- names.vec; df}, eg.df.list, names.list)
}


#######################################
# GO BP entire analysis
#######################################

do_GO_BP_enrichment <- function(fpathDExlsx, fpathCNTSxlsx, outBase) {
  
  fname <- basename(fpathDExlsx)
  overrepFname_BP <- gsub(".xlsx", "_over_GO_BP.xlsx", fname)
  gseaFname_BP <- gsub(".xlsx", "_gsea_GO_BP.xlsx", fname)
  notFoundFnameDE <- gsub(".xlsx", "_not_found.xlsx", fname)
  FoundFnameDE <- gsub(".xlsx", "_found.xlsx", fname)
  
  fnameCNTS <- basename(fpathCNTSxlsx)
  notFoundFnameCNTS <- gsub(".xlsx", "_not_found.xlsx", fnameCNTS)
  FoundFnameCNTS <- gsub(".xlsx", "_found.xlsx", fnameCNTS)
  
  # read-in data
  xlsxDE_dfList   <- xlsx2df.list(fpathDExlsx)[c(1, 4, 7)]
  xlsxCNTS_dfList <- xlsx2df.list(fpathCNTSxlsx)
  
  # translate to eg
  xlsxDE_dfListEg <- lapply(xlsxDE_dfList, alias2entrezid.df)
  xlsxCNTS_dfListEg <- lapply(xlsxCNTS_dfList, alias2entrezid.df)

  # not founds
  not_founds_list_DE       <- dfListNotFounds(xlsxDE_dfListEg)
  not_founds_list_CNTS     <- dfListNotFounds(xlsxCNTS_dfListEg)
  
  # founds
  founds_list_DE           <- dfListFounds(xlsxDE_dfListEg)
  founds_list_CNTS         <- dfListFounds(xlsxCNTS_dfListEg)
  
  # do overrep BP
  xlsxDE_dfListEg_GO_over_BP <- dfListEnrichGO(xlsxDE_dfListEg, "BP") 
  
  # do gsea BP
  xlsxDE_dfListEg_GO_gsea_BP <- dfListGseGO(xlsxDE_dfListEg, "BP") 
  
  # write results
  write.dfs(not_founds_list_DE, file.path(outBase, notFoundFnameDE))
  write.dfs(founds_list_DE, file.path(outBase, FoundFnameDE))
  write.dfs(not_founds_list_CNTS, file.path(outBase, notFoundFnameCNTS))
  write.dfs(founds_list_CNTS, file.path(outBase, FoundFnameCNTS))
  write.dfs(xlsxDE_dfListEg_GO_over_BP, file.path(outBase, overrepFname_BP))
  write.dfs(xlsxDE_dfListEg_GO_gsea_BP, file.path(outBase, gseaFname_BP))
}

DEDir <- "/media/josephus/Elements/DEanalysis"

inDir <- file.path(DEDir, "JT-cells-hybrid-new-0808/k2-vs-wtn")
DEfpath <- file.path(inDir, "DE-table/DE_sig_JT-cells-hybrid-new-0808_k2-vs-wtn_180808152049_0.05_1.xlsx")
CNTSfpath <- file.path(inDir, "count-table/nrm-counts-avg-JT-cells-hybrid-new-0808-k2-vs-wtn.xlsx")
outBase <- file.path(inDir, "GO")

do_GO_BP_enrichment(DEfpath, CNTSfpath, outBase)

inDir <- file.path(DEDir, "embr-0717/dKO-vs-wt")
DEfpath <- file.path(inDir, "DE-table/DE_sig_embr-0717_dKO-vs-wt_180717171725_0.05_1.xlsx")
CNTSfpath <- file.path(inDir, "count-table/nrm-counts-avg-embr-0717-dKO-vs-wt.xlsx")
outBase <- file.path(inDir, "GO")

do_GO_BP_enrichment(DEfpath, CNTSfpath, outBase)



#######################################################
# GO BP,CC,MF entire analysis with or without selection
# CNTS removed, since not used in GO
#######################################################

do_GO_enrichment <- function(fpathDExlsx,  # fpathCNTSxlsx, 
                             outBase, 
                             ont,
                             fpathCLUSTERxlsx = "",
                             select_clusters = c()) {
  
  fname <- basename(fpathDExlsx)
  
  if (length(select_clusters > 0)) {
    clust.numbers <- paste0("_", paste(unique(select_clusters), collapse="-"), "_")
  } else {
    clust.numbers <- "_"
  }
  
  decorateFname <- function(textPiece) gsub(".xlsx", 
                                            paste0(clust.numbers, textPiece, ont, ".xlsx"), 
                                            fname)
  overrepFname_ont <- decorateFname("over_GO_")
  gseaFname_ont <- decorateFname("gsea_GO_")
  notFoundFnameDE <- decorateFname("not_found_")
  FoundFnameDE <- decoreateFname("found_")
  
  if (!dir.exists(outBase)) {
    dir.create(outBase, recursive = TRUE)
  }
  
  # read-in data
  xlsxDE_dfList   <- xlsx2df.list(fpathDExlsx)[c(1, 4, 7)]
  
  if (fpathCLUSTERxlsx != "" && length(select_clusters) > 0) {
    # read-in cluster tables
    cluster.dfs <- xlsx2df.list(fpathCLUSTERxlsx)
    
    # extract gene names from clusters
    clusternames <- lapply(cluster.dfs, rownames)
    
    # extract genes from the clusters
    select_cluster_names <- unique(unlist(clusternames[select_clusters]))
    xlsxDE_dfList <- lapply(xlsxDE_dfList, function(df) df[intersect(rownames(df), select_cluster_names), ])
    
  }
  
  # translate to eg
  xlsxDE_dfListEg <- lapply(xlsxDE_dfList, alias2entrezid.df)
  
  # not founds
  not_founds_list_DE       <- dfListNotFounds(xlsxDE_dfListEg)
  
  # founds
  founds_list_DE           <- dfListFounds(xlsxDE_dfListEg)
  
  # do overrep ont
  xlsxDE_dfListEg_GO_over_ont <- dfListEnrichGO(xlsxDE_dfListEg, ont) 
  
  # do gsea ont
  xlsxDE_dfListEg_GO_gsea_ont <- dfListGseGO(xlsxDE_dfListEg, ont) 
  
  # write results
  write.dfs(not_founds_list_DE, file.path(outBase, notFoundFnameDE))
  write.dfs(founds_list_DE, file.path(outBase, FoundFnameDE))
  write.dfs(xlsxDE_dfListEg_GO_over_ont, file.path(outBase, overrepFname_ont))
  write.dfs(xlsxDE_dfListEg_GO_gsea_ont, file.path(outBase, gseaFname_ont))
}

do_all_GO_enrichment <- function(fpathDExlsx, 
                                 outBase,
                                 fpathCLUSTERxlsx = "",
                                 select_clusters = c())  {
  do_GO_enrichment(fpathDExlsx, outBase, "BP", fpathCLUSTERxlsx, select_clusters)
  do_GO_enrichment(fpathDExlsx, outBase, "CC", fpathCLUSTERxlsx, select_clusters)
  do_GO_enrichment(fpathDExlsx, outBase, "MF", fpathCLUSTERxlsx, select_clusters)
}


###############################################################
# pathways and calls
###############################################################

DEDir <- "/media/josephus/Elements/DEanalysis"

inDir <- file.path(DEDir, "JT-cells-hybrid-new-0808/k2-vs-wtn")
DEfpath <- file.path(inDir, "DE-table/DE_sig_JT-cells-hybrid-new-0808_k2-vs-wtn_180808152049_0.05_1.xlsx")
CNTSfpath <- file.path(inDir, "count-table/nrm-counts-avg-JT-cells-hybrid-new-0808-k2-vs-wtn.xlsx")
outBase <- file.path(inDir, "GO")

do_GO_BP_enrichment(DEfpath, CNTSfpath, outBase)

inDir <- file.path(DEDir, "embr-0717/dKO-vs-wt")
DEfpath <- file.path(inDir, "DE-table/DE_sig_embr-0717_dKO-vs-wt_180717171725_0.05_1.xlsx")
CNTSfpath <- file.path(inDir, "count-table/nrm-counts-avg-embr-0717-dKO-vs-wt.xlsx")
outBase <- file.path(inDir, "GO")

do_GO_BP_enrichment(DEfpath, CNTSfpath, outBase)

