###########################################
# Plotting of Heatmap for Figure 3a
###########################################

rm(list=ls()) 

####################################################################
# helper functions for printing count tables
####################################################################

join0 <- function(...) paste0(..., collapse = '')
csv.ending <- function(filename) join0(filename, ".csv")
xlsx.ending <- function(filename) join0(filename, ".xlsx")
tab.ending <- function(filename) join0(filename, ".tab")

triple.write <- function(obj, core_filename, outpath) {
  # write.csv(obj, file = csv.ending( file.path(outpath, core_filename)))
  write.xlsx(obj, file = xlsx.ending(file.path(outpath, core_filename)), 
             sheetName = "raw_counts")
  # write.table(obj, file = tab.ending( file.path(outpath, core_filename)))
}

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
  dfl <- lapply(seq(col.vals), function(i) select_df(df, val = col.vals[i], col.selector))
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

####################################################################
# required packages
####################################################################

require(DESeq2)      # DE analysis
require(xlsx)        # xlsx printing


####################################################################
# required paths
####################################################################

outdirpath <- "/home/josephus/Dropbox/amit_scripts/clustering_non_union/plots"
indirpath  <- "~/count/raw-subread/cleaned/"
metapath   <- "~/Dropbox/amit_scripts/alaAmit/meta.txt"
e2spath    <- "~/Dropbox/amit_scripts/alaHerrmann/ensembl2symbol/ensembleID2geneSymbol.txt"

####################################################################
# meta information
####################################################################

meta.df <- read.table(metapath, sep = '\t', 
                      header = TRUE, 
                      stringsAsFactors = FALSE)
denominator <- meta.df$condition[meta.df$testing == "denom"][1]
numerator   <- meta.df$condition[meta.df$testing == "num"][1]
core_name   <- paste0(numerator, "-vs-", denominator, collapse = '')

####################################################################
# give each condition a color and add to meta.df
####################################################################

col.condition <- c("grey", "blue", "orange", "red")
names(col.condition) <- unique(meta.df$condition) # in order given
meta.df$color <- col.condition[as.vector(meta.df$condition)]

####################################################################
# name annotation information
####################################################################
# load 'meta.txt' and 'ensembleID2geneSymbol.txt' and name symbol column 'symbol'
e2s.df <- read.table(e2spath, sep = '\t', 
                     header = FALSE, 
                     row.names = 1, 
                     stringsAsFactors = FALSE)
colnames(e2s.df) <- "symbol"

#####################################################################
# ensure existence of output paths
#####################################################################
dir.create(path=outdirpath, recursive = TRUE, showWarnings = FALSE)

#####################################################################
# create DESeqDataSet
#####################################################################

# use DESeqDataSetFromHTSeqCount() to generate count table
DESeq2.obj <- DESeqDataSetFromHTSeqCount(sampleTable = meta.df, 
                                         directory = indirpath, 
                                         design = ~ 0 + condition)
# ensure correct ordering of condition-levels
DESeq2.obj@colData$condition <- factor(DESeq2.obj@colData$condition,
                                       levels = unique(meta.df$condition))

#####################################################################
# get unnormalized and normalized count table
#####################################################################
count.table <- counts(DESeq2.obj)

# discard 'all-null-row's
count.table <- count.table[rowSums(count.table) > 0, ]

# normalize count table
DESeq2.norm <- estimateSizeFactors(DESeq2.obj)
sizeFactors(DESeq2.norm)

# modelling read counts through a neg binomial
DESeq2.disp <- estimateDispersions(DESeq2.norm)
plotDispEsts(DESeq2.disp)

count.table.norm <- counts(DESeq2.norm, normalized = TRUE)

#####################################################################
# set tresholds
#####################################################################
alpha       <- 0.05 # the padj threshold for significance
log2FClimit <- 1

# #####################################################################
# # DESeq
# #####################################################################
# 
# dds <- DESeq(DESeq2.obj)
# res <- results(dds, 
#                alpha=alpha, 
#                pAdjustMethod = "BH", 
#                contrast = c("condition", "dKO", "WT"))
# res <- na.omit(res)
# # it is exactly the same like 'wald test'

#####################################################################
# wald test
#####################################################################
 
wald.test <- nbinomWaldTest(DESeq2.disp)
res.DESeq2.dKOvsWT <- results(wald.test, alpha=alpha, 
                              pAdjustMethod = "BH", 
                              contrast = c("condition", "dKO", "WT"))

# remove NA rows
res.DESeq2.dKOvsWT <- na.omit(res.DESeq2.dKOvsWT)

head(res.DESeq2.dKOvsWT, 100)
head(res, 100)

all(rownames(res.DESeq2.dKOvsWT) == rownames(res)) # TRUE
all(dim(res.DESeq2.dKOvsWT) == dim(res)) # TRUE

#####################################################################
# Filter for settings
#####################################################################

dKOvsWT.kept <- rownames(res.DESeq2.dKOvsWT)[res.DESeq2.dKOvsWT$padj <= alpha & !is.na(res.DESeq2.dKOvsWT$padj)]
names.union <- c(dKOvsWT.kept)

#####################################################################
# select gene names based on union
# alpha <- 0.05
# gene.kept <- rownames(res.DESeq2)[res.DESeq2$padj <= alpha & !is.na(res.DESeq2$padj)]

gene.kept <- names.union


#####################################################################
# retrive normalized counts for gene of interest
#####################################################################

epsilon <- 1 # added to all values to avoid log2(0)
count.table.kept <- log2(count.table.norm + epsilon)[gene.kept, ]
dim(count.table.kept) # 705   8 ## non-union: 565 8

count.table.kept

#####################################################################
# build means
#####################################################################

avg.count.table.kept <- counts.avg.build(count.table.kept, 
                                         lapply(unique(meta.df$condition), 
                                                function(cond) which(meta.df$condition == cond)),
                                         unique(meta.df$condition))
dim(avg.count.table.kept) # 705    4 ## 565 4
head(avg.count.table.kept)

#####################################################################
# calculate 'distance between genes and samples as 1 minus correlation of one gene/sample to another)
#####################################################################

scaledata <- t(scale(t(avg.count.table.kept)))
scaledata <- scaledata[complete.cases(scaledata), ]


#####################################################################
# annotation of ensembl_id
#####################################################################
rownames(scaledata) <- e2s.df[rownames(scaledata), "symbol"]
head(scaledata)

#####################################################################
# load
#####################################################################
require(edgeR)
require(SummarizedExperiment)

require(ggplot2)               # for graphics
require(reshape)               # for graphics

require(gplots)                # for heatmap
require(RColorBrewer)          # for heatmap

require(pheatmap)             # for heatmap
require(xlsx)                  # for heatmap output

#####################################################################
# set k
#####################################################################
k <- 8

#####################################################################
# k means clustering k = 7
#####################################################################

kClust <- kmeans(scaledata, centers = k, nstart = 1000, iter.max = 30)
kClusters <- kClust$cluster

#####################################################################
# function to find centroid (cluster core) in cluster i
#####################################################################

clust.centroid <- function(i, dat, clusters) {
  ind = (clusters == i)
  colMeans(dat[ind, ])
}
kClustcentroids <- sapply(levels(factor(kClusters)), 
                          clust.centroid, scaledata, kClusters)


#####################################################################
# plot centroids to see how they behave
#####################################################################

# get in long form for plotting
Kmolten <- melt(kClustcentroids)
colnames(Kmolten) <- c("sample", "cluster", "value")

# ensure correct factorizing
Kmolten$sample <- factor(Kmolten$sample, levels = unique(meta.df$condition))

# plot
{
p1 <- ggplot(Kmolten, aes(x=factor(sample, levels = unique(meta.df$condition)), 
                          y=value, group = cluster,
                          colour=as.factor(cluster))) +
  geom_point() +
  geom_line() +
  xlab("Time") +
  ylab("Expression") +
  labs(title = "Cluster Expression by Group", color = "Cluster")
png(join0(outdirpath, "/", core_name, "-ClusterAll8.png"))
print(p1)
dev.off() 
}

#####################################################################
# check similarity of centroids
#####################################################################
cor(kClustcentroids)
# heatmap(cor(kClustcentroids))

#####################################################################
# subset the cores molten dataframe so we can plot the core
#####################################################################

for (i in 1:k) {
  assign(paste0("core", i), Kmolten[Kmolten$cluster == i, ])
  eval(parse(text = paste0("core", i, "$sample <- factor(core", i, 
                           "$sample, levels = unique(meta.df$condition))")))
}

#####################################################################
# get clusters
#####################################################################

for (i in 1:k) {
  assign(paste0("K", i), scaledata[kClusters == i, ])
}

#####################################################################
# calculate correlation with the core
#####################################################################

for (i in 1:k) {
  assign(paste0("corescore", i), 
         eval(parse(text = paste0("function(x) {cor(x, core", i, "$value)}"))))
}

for (i in 1:k) {
  assign(paste0("score", i), 
         eval(parse(text = paste0("apply(K", i, ", 1, corescore", i, ")"))))
}

#####################################################################
# get data frame into long format for plotting
#####################################################################

for (i in 1:k) {
  assign(paste0("K", i, "molten"), 
         eval(parse(text = paste0("melt(K", i, ")"))))
  eval(parse(text = paste0("colnames(K", i, "molten) <- c('gene', 'sample', 'value')")))
}

#####################################################################
# add the score
#####################################################################

for (i in 1:k) {
  eval(parse(text = paste0("K", i, "molten <- merge(K", i, "molten, score", i, ", by.x = 'gene', by.y = 'row.names', all.x = T)")))
  eval(parse(text = paste0("colnames(K", i, "molten) <- c('gene', 'sample', 'value', 'score')")))
}

#####################################################################
# order dataframe by score
#####################################################################

for (i in 1:k) {
  eval(parse(text = paste0("K", i, "molten <- K", i, "molten[order(K", i, "molten$score), ]")))
}


#####################################################################
# plotting cluster groups
#####################################################################

for (i in 1:k) {
  text = paste0("sp", i, " <- ggplot(K", i, "molten, aes(x=factor(sample, levels=unique(meta.df$condition)), y=value)) + ",
                "geom_line(aes(colour=score, group=gene)) + ",
                "scale_color_gradientn(colours=c('blue1', 'red2')) + ",
                # this adds the core
                "geom_line(data=core", i, ", aes(sample,value,group=cluster), ",
                "color='black', inherit.aes=FALSE) + ",
                "xlab('Time') + ",
                "ylab('Expression') + ",
                "labs(title='Cluster ", i, " Expression by Group', color = 'Score'); ",
                "png('", outdirpath, "/", core_name, "-Cluster", i, ".png'); ",
                "print(sp", i, "); dev.off()"
                )
  eval(parse(text = text))
}


#####################################################################
# prepare heatmap
#####################################################################

colfunc <- colorRampPalette(c("black", "red"))

eval(parse(text = paste0("scores <- c(", paste(paste("score", 1:k, sep = ''), 
                                               collapse = ", "), ")")))

# add cluster number and score for sorting the data
scaledata.k <-cbind(scaledata, 
                    cluster = kClust$cluster, 
                    score = scores[rownames(scaledata)])
scaledata.k <- scaledata.k[order(scaledata.k[, "cluster"], scaledata.k[, "score"]), ]

# unordered clusters

{
  svg(file.path(outdirpath, "kmeans_heatmap_non_union_embr_8_unordered.svg"))
  pheatmap(scaledata.k[, 1:length(unique(meta.df$condition))], 
           cluster_rows = F, 
           cluster_cols = F,
           cellwidth = 40,
           col = colfunc(100), 
           fontsize_row = .5
  )
  dev.off()
}

dfl <- df2dflist(scaledata.k, "cluster")
# within cluster sort by score
dflsc <- lapply(dfl, function(df) df[order(df[, "score"]), ])
# 1 2 3 4 5 6 7 8
# c(5, 8, 3, 7, 4, 6, 2, 1)
dflo <- list.order.manually.by.vec(dflsc, c(5, 7, 3, 8, 4, 6, 1, 2))
original.cluster.order <- names(dflo)
names(dflo) <- 1:k
dflo.corr <- repeated.values.into.df.list(dflo, "cluster", 1:k)
scaledata.k.ordered <- dflist2df(dflo.corr)

gaps_idxs <- cumsum(unlist(lapply(df2dflist(scaledata.k.ordered, "cluster"), function(df) dim(df)[1])))
gaps_idxs <- gaps_idxs[c(1, 3, 4, 5, 6, 8)]
{
  svg(file.path(outdirpath, "kmeans_heatmap_non_union_embr_8_man_ordered_prelim_final.svg"))
  pheatmap(scaledata.k.ordered[, 1:length(unique(meta.df$condition))], 
           cluster_rows = F, 
           cluster_cols = F,
           gaps_row = gaps_idxs,
           cellwidth = 40,
           border_color = NA,
           col = colfunc(100), 
           fontsize_row = .5
  )
  dev.off()
}

#            gap_rows = table(scaledata.k.corrected[, "cluster"])[1:k],
# put horizontal gaps
# https://stackoverflow.com/questions/41628450/r-pheatmap-change-annotation-colors-and-prevent-graphics-window-from-popping-up

#####################################################################
# print scaledata
#####################################################################
scaledata_list <- df2dflist(scaledata.k.ordered, "cluster")
{
  for (i in 1:k) {
    eval(parse(text = paste0("write.xlsx(scaledata_list[[", i, "]], '", outdirpath, "/kmeans", k, "Embr.xlsx', sheetName = 'Clust", i, "', append = ", ifelse(i <= 1, FALSE, TRUE), ")")))
  }
  write.xlsx(scaledata.k.ordered, paste0(outdirpath, "/kmeans", k, "Embr.xlsx"), sheetName = "All", append = TRUE)
}





