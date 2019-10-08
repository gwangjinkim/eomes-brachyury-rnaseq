#######################
# HOMER
#######################


wt_atac_peaks="/media/josephus/Elements/galaxy_atac_data/atac-wt-merge-86.bed"
dko_atac_peaks="/media/josephus/Elements/galaxy_atac_data/atac-dko-merge-92.bed"

wt_h3k27ac_peaks="/media/josephus/Elements/galaxy_atac_data/h3k27ac-wt-merge-98.bed"
dko_h3k27ac_peaks="/media/josephus/Elements/galaxy_atac_data/h3k27ac-dko-merge-104.bed"

wt_allup_atac_peaks="/media/josephus/Elements/galaxy_atac_data/atac-wt-vs-dko-all-up-178.bed"
wt_alldown_atac_peaks="/media/josephus/Elements/galaxy_atac_data/atac-wt-vs-dko-all-down-174.bed"

dko_up_atac_chip_peaks="/media/josephus/Elements/galaxy_atac_data/atac-and-eo-br-chip-intersect-dko-up-187.bed"
dko_down_atac_chip_peaks="/media/josephus/Elements/galaxy_atac_data/atac-and-eo-br-chip-intersect-dko-down-188.bed"

eogfp_chip_peaks="/media/josephus/Elements/galaxy_atac_data/eogfp-chip-154.bed"
hermann_tchip_peaks="/media/josephus/Elements/galaxy_atac_data/hermann-tchip-c16-157.bed"

all_atac_chip="/media/josephus/Elements/galaxy_atac_data/atac-chip-all-317.bed"
up_atac_chip="/media/josephus/Elements/galaxy_atac_data/atac-chip-all-up-318.bed"
down_atac_chip="/media/josephus/Elements/galaxy_atac_data/atac-chip-all-down-319.bed"

up_atac_rnaseq="/media/josephus/Elements/galaxy_atac_data/atac-rnaseq-up-415.bed"
down_atac_rnaseq="/media/josephus/Elements/galaxy_atac_data/atac-rnaseq-down-418.bed"


wt_atac_out="/media/josephus/Elements/ChIPresults/ATACseq2018/bowtie2-gx/mm10/atac-wt-merged-bg-atac-dko-0919"
dko_atac_out="/media/josephus/Elements/ChIPresults/ATACseq2018/bowtie2-gx/mm10/atac-dko-merged-bg-atac-wt-0919"

wt_h3k27ac_out="/media/josephus/Elements/ChIPresults/Chipseq2017/bowtie2-gx/mm10/h3k27ac-wt-merged"
dko_h3k27ac_out="/media/josephus/Elements/ChIPresults/Chipseq2017/bowtie2-gx/mm10/h3k27ac-dko-merged"

wt_allup_out="/media/josephus/Elements/ChIPresults/ATACseq2018/bowtie2-gx/mm10/atac-wt-allup"
wt_alldown_out="/media/josephus/Elements/ChIPresults/ATACseq2018/bowtie2-gx/mm10/atac-wt-alldown"

dko_up_atac_chip_out="/media/josephus/Elements/ChIPresults/ATACseq2018/bowtie2-gx/mm10/atac-chip-dko-up"
dko_down_atac_chip_out="/media/josephus/Elements/ChIPresults/ATACseq2018/bowtie2-gx/mm10/atac-chip-dko-down"

eogfp_chip_out="/media/josephus/Elements/ChIPresults/Chipseq2017/bowtie2-gx/mm10/eogfp-chip"
hermann_tchip_out="/media/josephus/Elements/ChIPresults/Chipseq2017/bowtie2-gx/mm10/hermann-tchip-c16"

all_atac_chip_out="/media/josephus/Elements/ChIPresults/ATACseq2018/all-atac-chip-317"
up_atac_chip_out="/media/josephus/Elements/ChIPresults/ATACseq2018/up-atac-chip-318"
down_atac_chip_out="/media/josephus/Elements/ChIPresults/ATACseq2018/down-atac-chip-319"

up_atac_rnaseq_out="/media/josephus/Elements/ChIPresults/ATACseq2018/up-atac-rnaseq-415"
down_atac_rnaseq_out="/media/josephus/Elements/ChIPresults/ATACseq2018/down-atac-rnaseq-418"


mkdir -p ${wt_atac_out}
mkdir -p ${dko_atac_out}
mkdir -p ${wt_h3k27ac_out}
mkdir -p ${dko_h3k27ac_out}
mkdir -p ${wt_allup_out}
mkdir -p ${wt_alldown_out}
mkdir -p ${dko_up_atac_chip_out}
mkdir -p ${dko_down_atac_chip_out}
mkdir -p ${eogfp_chip_out}
mkdir -p ${hermann_tchip_out}
mkdir -p ${all_atac_chip_out}
mkdir -p ${up_atac_chip_out}
mkdir -p ${down_atac_chip_out}
mkdir -p ${up_atac_rnaseq_out}
mkdir -p ${down_atac_rnaseq_out}


findMotifsGenome.pl ${wt_atac_peaks} mm10 ${wt_atac_out}/ -bg ${dko_atac_peaks} -size 200 -p 6 && \
findMotifsGenome.pl ${dko_atac_peaks} mm10 ${dko_atac_out} -bg ${wt_atac_peaks} -size 200 -p 6

findMotifsGenome.pl ${wt_atac_peaks} mm10 ${wt_atac_out}/ -size 200 -p 6 && \
findMotifsGenome.pl ${dko_atac_peaks} mm10 ${dko_atac_out} -size 200 -p 6 && \
findMotifsGenome.pl ${wt_h3k27ac_peaks} mm10 ${wt_h3k27ac_out} -size 200 -p 6 && \
findMotifsGenome.pl ${dko_h3k27ac_peaks} mm10 ${dko_h3k27ac_out} -size 200 -p 6 && \


findMotifsGenome.pl ${wt_allup_atac_peaks} mm10 ${wt_allup_out} -size 200 -p 6 && \
findMotifsGenome.pl ${wt_alldown_atac_peaks} mm10 ${wt_alldown_out} -size 200 -p 6 && \
findMotifsGenome.pl ${dko_up_atac_chip_peaks} mm10 ${dko_up_atac_chip_out} -size 200 -p 6 && \
findMotifsGenome.pl ${dko_down_atac_chip_peaks} mm10 ${dko_down_atac_chip_out} -size 200 -p 6

findMotifsGenome.pl ${dko_down_atac_chip_peaks} mm10 ${dko_down_atac_chip_out} -size 200 -p 6

findMotifsGenome.pl ${eogfp_chip_peaks} mm10 ${eogfp_chip_out} -size 200 -p 6 && \
findMotifsGenome.pl ${hermann_tchip_peaks} mm10 ${hermann_tchip_out} -size 200 -p 6

rsync -av /media/josephus/Elements/ChIPresults/ /media/josephus/24AA554AAA551A1E/ChIPresults/


findMotifsGenome.pl ${all_atac_chip} mm10 ${all_atac_chip_out} -size 200 -p 6 && \
findMotifsGenome.pl ${up_atac_chip} mm10 ${up_atac_chip_out} -size 200 -p 6 && \
findMotifsGenome.pl ${down_atac_chip} mm10 ${down_atac_chip_out} -size 200 -p 6 && \
findMotifsGenome.pl ${up_atac_rnaseq} mm10 ${up_atac_rnaseq_out} -size 200 -p 6 && \
findMotifsGenome.pl ${down_atac_rnaseq} mm10 ${down_atac_rnaseq_out} -size 200 -p 6

rsync -av /media/josephus/Elements/ChIPresults/ /media/josephus/24AA554AAA551A1E/ChIPresults/




## 181025

# homer analyses like written

mv '/media/josephus/Elements/galaxy_chip_final_bams/Galaxy340-[Galaxy174-[DEseq2_ATAC_all_peaks.txt].bed].bed' /media/josephus/Elements/galaxy_chip_final_bams/340-all-ATAC.bed
mv '/media/josephus/Elements/galaxy_chip_final_bams/Galaxy341-[Intersect_on_data_335_and_data_340_ChIPmerge_vs_ATACall].bed' /media/josephus/Elements/galaxy_chip_final_bams/341-EoTChIP-atac-all.bed


gr.eobrUP <- readBed(eobrchip.UP)
gr.eobrDOWN <- readBed(eobrchip.DOWN)
length(unique(gr.eobrUP))
length(unique(gr.eobrDOWN))

# > length(gr.eobrUP)
# [1] 4778
# > length(gr.eobrDOWN)
# [1] 6527
# > length(unique(gr.eobrUP))
# [1] 4778
# > length(unique(gr.eobrDOWN))
# [1] 6527



mergeBedUnique <- function(bed1, bed2) {

}


require(GenomicRanges)

g1 <- GRanges(data.frame(seqnames="Chr1",
                         strand="+",
                         start=c(2,2,2,4,5,6),
                         end=c(24,25,25,25,25,21),
                         m1=c("a","a","a","b","b","b"),
                         m2=c("A","B","A","B","A","B"),
                         m3=c(1,2,3,4,5,6), stringsAsFactors=F))
g2 <- GRanges(data.frame(seqnames="Chr1",
                         strand="+",
                         start=c(2,3,2,5,7,6,7),
                         end=c(24,25,25,25,25,21,25),
                         m1=c("a","a","a","b","b","b","b"),
                         m2=c("A","B","A","B","A","B","A"),
                         m3=7:1, stringsAsFactors=F))

g1 <- GRanges(seqnames="Chr1", strand="+", ranges=IRanges(start=c(2,2,2,4,5,6), end=c(24,25,25,25,25,21)))
g2 <- GRanges(seqnames="Chr1", strand="+", ranges=IRanges(start=c(2,3,2,5,7,6, 7), end=c(24,25,25,25,25,21, 25)))

g <- c(g1, g2)

# merge and overtake from the first
# merge and follow a fule function for each of the meta columns

# merge and just merge column content using a string ";"

g1r <- unique(g1)
g2r <- unique(g2)

# unique doesn't care of meta, only of the position information
# there might be a function which squishes the position but aggregates the meta information?

hit.idxs <- findOverlaps(g1r, g2r)

# g1r.starts <- start(g1r)
# g1r.ends <- end(g1r)
# g1r.chrs <- seqnames(g1r)
# g2r.starts <- start(g2r)
# g2r.ends <- end(g2r)
# g2r.chrs <- seqnames(g2r)
# 
# pos1ToPos2 <- lapply(g1r.starts, function(x) grep(x, g2r.starts))
# repeats <- sapply(pos1ToPos2, length)
# repeated.g1r.starts <- Map(rep, g1r.starts, repeats)
# unlist(pos1ToPos2)
# unlist(repeated.g1r.starts)
# 
# matchTable <- function(vec1, vec2) {
#   # Return match table of two vectors
#   idx2s <- lapply(vec1, function(x) grep(x, vec2))
#   idx1s <- Map(rep, 1:length(vec1), sapply(idx2s, length))
#   data.frame(queryHits=unlist(idx1s), subjectHits=unlist(idx2s), stringsAsFactors=F)
# }
# 

#########################
# find exact position matches - independent of meta data
#########################

hit.idxs <- findOverlaps(g1r, g2r, type='equal')

# query=       maxgap=      type=        ...=         
# subject=     minoverlap=  select=      
# > findOverlaps(g1r, g2r, type='exact')
# Error in match.arg(type) : 
#   'arg' should be one of “any”, “start”, “end”, “within”, “equal”
# 

## since made unique before, it will be 1:1 matches

# take two df rows and fuse their contents with ";"

# cave! changing dataframe row using simple as.character()
# cannot handle factors (represented then as integer)
# do sapply(row, as.character) instead!

r1vec <- sapply(as.data.frame(g1r[1]), as.character)
r2vec <- sapply(as.data.frame(g2r[2]), as.character)
paste0(r1vec, ";", r2vec) # does what it should
fused.row <- c(r1vec[1:5], paste0(r1vec[6:length(r1vec)], ";", r2vec[6:length(r2vec)]))
names(fused.row) <- names(r1vec)

GRanges(as.data.frame(t(data.frame(fused.row, stringsAsFactors=F))))

# how stupid that t() of a data.frame gives a matrix!

############################################################################
toDataFrame <- function(x) {
  as.data.frame(x, stringsAsFactors=F)
}
############################################################################
grangesRowToNamedVector <- function(x) {
  sapply(as.data.frame(x, stringsAsFactors=F), as.character)
}


# not necessary for following code
fuseTwoGRangesRows <- function(gr.row.1, gr.row.2) {
  # How to remove the warnings?
  row1.vec <- grangesRowToNamedVector(gr.row.1)
  row2.vec <- grangesRowToNamedVector(gr.row.2)
  if (length(row1.vec) > 5 && length(row2.vec) > 5) {
    fused.row <- c(r1vec[1:5], paste0(row1.vec[6:length(row1.vec)], ";",
                                      row2.vec[6:length(row2.vec)]))
    names(fused.row) <- names(row1.vec)
  } else {
    fused.row <- r1vec[1:5]
  }
  GRanges(as.data.frame(t(data.frame(fused.row, stringsAsFactors=F)), 
                        stringsAsFactors=F))
}
############################################################################

.fuseTwoGRangesTablesBySemicolon <- function(gr.1, gr.2) {
  # fuse gRangesTAbles
  df.1 <- as.data.frame(gr.1) # here, stringsAsFactors=F
  df.2 <- as.data.frame(gr.2) # produces warning - so put away.
  if(ncol(df.1) > 5 && ncol(df.2) > 5) {
    fused.df <- cbind(df.1[, 1:5], 
                      as.data.frame(lapply(6:ncol(df.1), 
                                           function(i) {
                                             paste0(df.1[, i], 
                                                    ";", 
                                                    df.2[, i])}),
                                    stringsAsFactors=F))
    colnames(fused.df) <- colnames(df.1)
  } else {
    fused.df <- df.1[, 1:5]
  }
  gr.common <- GRanges(as.data.frame(fused.df, stringsAsFactors=F))
}

fuseTwoGRangesTables <- function(gr1, gr2) {
  hit.idxs <- findOverlaps(gr1, gr2, type='equal')
  gr1.doubled.idxs <- hit.idxs@from
  gr2.doubled.idxs <- hit.idxs@to
  gr1.only <- gr1[-gr1.doubled.idxs]
  gr2.only <- gr2[-gr2.doubled.idxs]
  gr1.2.common <- .fuseTwoGRangesTablesBySemicolon(gr1[gr1.doubled.idxs], 
                                                   gr2[gr2.doubled.idxs])
  sort(c(gr1.only, gr1.2.common, gr2.only))
}

fuseTwoGRangesTables(g1r[hit.idxs@from], g2r[hit.idxs@to]) # this works
fuseTwoGRangesTables(g1r, g2r)


#########################
# fuse two bed files
#########################

eobrchip.UP <- "/media/josephus/Elements/galaxy_chip_final_bams/359-EoTChIP-UP.bed"
eobrchip.DOWN <- "/media/josephus/Elements/galaxy_chip_final_bams/358-EoTChIP-DOWN.bed"

eobrchip.atac.UP   <- '/media/josephus/Elements/galaxy_chip_final_bams/423-chip-atac-3FC-DE-up.bed'
eobrchip.atac.DOWN <- '/media/josephus/Elements/galaxy_chip_final_bams/422-chip-atac-3FC-DE-down.bed'


####################################
# load BED/Peak
####################################

require(ChIPseeker)
readBed <- function(bedPath) {
  readPeakFile(bedPath, header=F)
}

gr.eobrUP <- readBed(eobrchip.UP)
gr.eobrDOWN <- readBed(eobrchip.DOWN)

gr.eobrUPDOWN <- fuseTwoGRangesTables(gr.eobrUP, gr.eobrDOWN)
gr.eobrUPDOWN$V4 <- sapply(as.character(gr.eobrUPDOWN$V4), function(x) strsplit(x, ";")[[1]][1])

require(xlsx2dfs)
write.table(gr.eobrUPDOWN, 
            "/media/josephus/Elements/galaxy_chip_final_bams/359-358-fused-181026-EoTChIP-UPDOWN.bed",
            sep="\t", col.names=F, row.names=F, quote=F)

## success!! works!

gr.eobratacUP <- readBed(eobrchip.atac.UP)
gr.eobratacDOWN <- readBed(eobrchip.atac.DOWN)

gr.eobratacUPDOWN <- fuseTwoGRangesTables(gr.eobratacUP, gr.eobratacDOWN)
gr.eobratacUPDOWN$V4 <- sapply(as.character(gr.eobratacUPDOWN$V4), function(x) strsplit(x, ";")[[1]][1])

write.table(gr.eobratacUPDOWN, 
            "/media/josephus/Elements/galaxy_chip_final_bams/423-422-fused-181026-EoTChIP-atac-UPDOWN.bed",
            sep="\t", col.names=F, row.names=F, quote=F)

## there should be written methods to fuse BED files
## to intersect BED files in a meta-maintaining manner (using ";")
## but it should test ";" existence in each single one of them before







#####################################
# common lisp `apply` is do.call(),
# what=func or funcNameStr
# args=list of args, names are attr names
# quote=TRUE quote the arguments
# envir=environment to execute
#####################################

# match(g1r.starts, g2r.starts)









eobrchip_UP="/media/josephus/Elements/galaxy_chip_final_bams/359-EoTChIP-UP.bed"
eobrchip_DOWN="/media/josephus/Elements/galaxy_chip_final_bams/358-EoTChIP-DOWN.bed"

# 

eobrchip_UPDOWN="/media/josephus/Elements/galaxy_chip_final_bams/359-358-fused-181026-EoTChIP-UPDOWN.bed"
eobrchip_all="/media/josephus/Elements/galaxy_chip_final_bams/357-EoTChIP.bed"
eobrchip_atac_all='/media/josephus/Elements/galaxy_chip_final_bams/341-EoTChIP-atac-all.bed'
eobrchip_atac_UPDOWN='/media/josephus/Elements/galaxy_chip_final_bams/423-422-fused-181026-EoTChIP-atac-UPDOWN.bed'
atac_all='/media/josephus/Elements/galaxy_chip_final_bams/340-all-ATAC.bed'
eobrchip_atac_UP='/media/josephus/Elements/galaxy_chip_final_bams/423-chip-atac-3FC-DE-up.bed'
eobrchip_atac_DOWN='/media/josephus/Elements/galaxy_chip_final_bams/422-chip-atac-3FC-DE-down.bed'
eobrchip_atac_UPDOWN="/media/josephus/Elements/galaxy_chip_final_bams/359-358-fused-181026-EoTChIP-UPDOWN.bed"

EoTChIP_UP_bg_UPDOWN="/media/josephus/Elements/ChIPresults/homer/ChIPseq2018_ATACseq2018/EoTChIP_UP_359_bg_UPDOWN_359_358"
EoTChIP_UP_bg_EoTChIP="/media/josephus/Elements/ChIPresults/homer/ChIPseq2018_ATACseq2018/EoTChIP_UP_359_bg_EoTChIP_357"
EoTChIP_UP_bg_DOWN="/media/josephus/Elements/ChIPresults/homer/ChIPseq2018_ATACseq2018/EoTChIP_UP_359_bg_DOWN_358"

EoTChIP_DOWN_bg_UPDOWN="/media/josephus/Elements/ChIPresults/homer/ChIPseq2018_ATACseq2018/EoTChIP_DOWN_358_bg_UPDOWN_359_358"
EoTChIP_DOWN_bg_EoTChIP="/media/josephus/Elements/ChIPresults/homer/ChIPseq2018_ATACseq2018/EoTChIP_DOWN_358_bg_EoTChIP_357"
EoTChIP_DOWN_bg_UP="/media/josephus/Elements/ChIPresults/homer/ChIPseq2018_ATACseq2018/EoTChIP_DOWN_358_bg_UP_359"

EoTChIP_atac_UP_bg_EoTChIP_atac_UPDOWN_out="/media/josephus/Elements/ChIPresults/homer/ChIPseq2018_ATACseq2018/EoTChIP_atac_UP_423_bg_EoTChIP_atac_UPDOWN_423_422"
EoTChIP_atac_UP_bg_EoTChIP_out="/media/josephus/Elements/ChIPresults/homer/ChIPseq2018_ATACseq2018/EoTChIP_atac_UP_423_bg_EoTChIP_357"
EoTChIP_atac_UP_bg_EoTChIP_atac_all_out="/media/josephus/Elements/ChIPresults/homer/ChIPseq2018_ATACseq2018/EoTChIP_atac_UP_423_bg_EoTChIP_atac_341"
EoTChIP_atac_UP_bg_EoTChIP_atac_DOWN_422_out="/media/josephus/Elements/ChIPresults/homer/ChIPseq2018_ATACseq2018/EoTChIP_atac_UP_423_bg_EoTChIP_ata_DOWN_422"
EoTChIP_atac_UP_bg_atac_all_out="/media/josephus/Elements/ChIPresults/homer/ChIPseq2018_ATACseq2018/EoTChIP_atac_UP_423_bg_atac_all_340"

EoTChIP_atac_DOWN_bg_EoTChIP_atac_UPDOWN_out="/media/josephus/Elements/ChIPresults/homer/ChIPseq2018_ATACseq2018/EoTChIP_atac_DOWN_422_bg_EoTChIP_atac_UPDOWN_423_422"
EoTChIP_atac_DOWN_bg_EoTChIP_out="/media/josephus/Elements/ChIPresults/homer/ChIPseq2018_ATACseq2018/EoTChIP_atac_DOWN_422_bg_EoTChIP_357"
EoTChIP_atac_DOWN_bg_EoTChIP_atac_all_out="/media/josephus/Elements/ChIPresults/homer/ChIPseq2018_ATACseq2018/EoTChIP_atac_DOWN_422_bg_EoTChIP_atac_341"
EoTChIP_atac_DOWN_bg_EoTChIP_atac_UP_423_out="/media/josephus/Elements/ChIPresults/homer/ChIPseq2018_ATACseq2018/EoTChIP_atac_DOWN_422_bg_EoTChIP_ata_UP_423"
EoTChIP_atac_DOWN_bg_atac_all_out="/media/josephus/Elements/ChIPresults/homer/ChIPseq2018_ATACseq2018/EoTChIP_atac_DOWN_422_bg_atac_all_340"

EoTChIP_DOWN_bg_homer="/media/josephus/Elements/ChIPresults/homer/ChIPseq2018_ATACseq2018/EoTChIP_DOWN_358_bg_all"
EoTChIP_UP_bg_homer="/media/josephus/Elements/ChIPresults/homer/ChIPseq2018_ATACseq2018/EoTChIP_UP_359_bg_all"

findMotifsGenome.pl ${eobrchip_UP} mm10 ${EoTChIP_UP_bg_homer}/ -size 200 -p 6 && \
findMotifsGenome.pl ${eobrchip_DOWN} mm10 ${EoTChIP_DOWN_bg_homer}/ -size 200 -p 6

findMotifsGenome.pl ${eobrchip_UP} mm10 ${EoTChIP_UP_bg_UPDOWN}/ -bg ${eobrchip_UPDOWN} -size 200 -p 6 && \
findMotifsGenome.pl ${eobrchip_UP} mm10 ${EoTChIP_UP_bg_EoTChIP} -bg ${eobrchip_all} -size 200 -p 6 && \
findMotifsGenome.pl ${eobrchip_UP} mm10 ${EoTChIP_UP_bg_DOWN} -bg ${eobrchip_DOWN} -size 200 -p 6 && \
findMotifsGenome.pl ${eobrchip_DOWN} mm10 ${EoTChIP_DOWN_bg_UPDOWN}/ -bg ${eobrchip_UPDOWN} -size 200 -p 6 && \
findMotifsGenome.pl ${eobrchip_DOWN} mm10 ${EoTChIP_DOWN_bg_EoTChIP} -bg ${eobrchip_all} -size 200 -p 6 && \
findMotifsGenome.pl ${eobrchip_DOWN} mm10 ${EoTChIP_DOWN_bg_UP} -bg ${eobrchip_UP} -size 200 -p 6 && \
findMotifsGenome.pl ${eobrchip_atac_UP} mm10 ${EoTChIP_atac_UP_bg_EoTChIP_atac_UPDOWN_out}/ -bg ${eobrchip_atac_UPDOWN} -size 200 -p 6 && \
findMotifsGenome.pl ${eobrchip_atac_UP} mm10 ${EoTChIP_atac_UP_bg_EoTChIP_out} -bg ${eobrchip_all} -size 200 -p 6 && \
findMotifsGenome.pl ${eobrchip_atac_UP} mm10 ${EoTChIP_atac_UP_bg_EoTChIP_atac_all_out} -bg ${eobrchip_atac_al} -size 200 -p 6 && \
findMotifsGenome.pl ${eobrchip_atac_UP} mm10 ${EoTChIP_atac_UP_bg_EoTChIP_atac_DOWN_422_out} -bg ${eobrchip_atac_DOWN} -size 200 -p 6 && \
findMotifsGenome.pl ${eobrchip_atac_DOWN} mm10 ${EoTChIP_atac_DOWN_bg_EoTChIP_atac_UPDOWN_out} -bg ${eobrchip_atac_UPDOWN} -size 200 -p 6 && \
findMotifsGenome.pl ${eobrchip_atac_DOWN} mm10 ${EoTChIP_atac_DOWN_bg_EoTChIP_out} -bg ${eobrchip_all} -size 200 -p 6 && \
findMotifsGenome.pl ${eobrchip_atac_DOWN} mm10 ${EoTChIP_atac_DOWN_bg_EoTChIP_atac_all_out} -bg ${eobrchip_atac_all} -size 200 -p 6 && \
findMotifsGenome.pl ${eobrchip_atac_DOWN} mm10 ${EoTChIP_atac_DOWN_bg_EoTChIP_atac_UP_423_out} -bg ${eobrchip_atac_UP} -size 200 -p 6 && \
findMotifsGenome.pl ${eobrchip_atac_UP} mm10 ${EoTChIP_atac_UP_bg_atac_all_out} -bg ${atac_all} -size 200 -p 6 && \
findMotifsGenome.pl ${eobrchip_atac_DOWN} mm10 ${EoTChIP_atac_DOWN_bg_atac_all_out} -bg ${atac_all} -size 200 -p 6




findMotifsGenome.pl ${eobrchip_UP} mm10 ${EoTChIP_UP_bg_DOWN} -bg ${eobrchip_DOWN} -size 200 -p 6 && \
findMotifsGenome.pl ${eobrchip_DOWN} mm10 ${EoTChIP_DOWN_bg_UPDOWN}/ -bg ${eobrchip_UPDOWN} -size 200 -p 6 && \
findMotifsGenome.pl ${eobrchip_DOWN} mm10 ${EoTChIP_DOWN_bg_EoTChIP} -bg ${eobrchip_all} -size 200 -p 6 && \
findMotifsGenome.pl ${eobrchip_DOWN} mm10 ${EoTChIP_DOWN_bg_UP} -bg ${eobrchip_UP} -size 200 -p 6 





################################################################################

## 190510 homer analysis eo br final ChIP alone

path_to_tosicj190315="/media/josephus/My Book/1TB_backup/tosicj190315"
out="/media/josephus/Elements/ChIPresults/homer/ChIPseq2019"
mkdir -p ${out}

eo_final_peak_fpath="${path_to_tosicj190315}/Peaks_ChIP_EoGFP_pileup120.bed"
br_final_peak_fpath="${path_to_tosicj190315}/Peaks_ChIP_BraGFP_pileup40.bed"

out_eo="${out}/eo_tosicj190315"
out_br="${out}/br_tosicj190315"

out_eo_alone="${out_eo}/eo_final_bg_none"
out_br_alone="${out_br}/br_final_br_none"
# out_eo_bg_input="${out_eo}/eo_final_bg_input"
# out_br_bg_input="${out_br}/br_final_bg_input"
out_eo_bg_br="${out_eo}/eo_final_bg_br_final"
out_br_bg_eo="${out_br}/br_final_bg_eo_final"


mkdir -p ${out_eo_alone}
mkdir -p ${out_br_alone}
# mkdir -p ${out_eo_input}
# mkdir -p ${out_br_input)
mkdir -p ${out_eo_bg_br}
mkdir -p ${out_br_bg_eo}


findMotifsGenome.pl "${eo_final_peak_fpath}" mm10 "${out_eo_alone}"/ -size 200 -p 6 && \
findMotifsGenome.pl "${br_final_peak_fpath}" mm10 "${out_br_alone}"/ -size 200 -p 6 && \
findMotifsGenome.pl "${eo_final_peak_fpath}" mm10 "${out_eo_bg_br}"/ -bg "${br_final_peak_fpath}" -size 200 -p 6 && \
findMotifsGenome.pl "${br_final_peak_fpath}" mm10 "${out_br_bg_eo}"/ -bg "${eo_final_peak_fpath}" -size 200 -p 6


################################################################################


## 190510 in addition

eo_bra_intersect="/media/josephus/My Book/galaxy_chip_jelena/Galaxy620-[Intersect_intervals_on_data_435_and_data_617_EoBra_Overlapping].bed"

out_eo_br_intersect_bg_none="${out}/eo_intersect_br/eo_intersect_br_bg_none"
mkdir -p ${out_eo_br_intersect_bg_none}
findMotifsGenome.pl "${eo_bra_intersect}" mm10 "${out_eo_br_intersect_bg_none}"/ -size 200 -p 6



################################################################################

## 190531

## galaxy 1118 # 1kb promoters
## galaxy 1119 # 1kb non-promoters

## no background
## all background galaxy 974



promoters_1kb="/media/josephus/My Book/galaxy_chip_jelena/tbx_peaks_promoters/Galaxy1118-[Intersect_intervals_on_data_1117_and_data_974_Tbxpekas_promoters_1kb].bed"
promoters_1kb_complement="/media/josephus/My Book/galaxy_chip_jelena/tbx_peaks_promoters/Galaxy1119-[Intersect_intervals_on_data_1117_and_data_974_Tbxpeaks_more1kb].bed"

all="/media/josephus/My Book/galaxy_chip_jelena/tbx_peaks_promoters/Galaxy974-[Galaxy448-(Merged_SortBed_on_EoGFP_mergeN_PU120_AND_BraGFP_mergeN_PU40)GREAT.bed].bed"

out_1kb="/media/josephus/Elements/ChIPresults/homer/ChIPseq2019/promoter_1kb"
out_1kb_compl="/media/josephus/Elements/ChIPresults/homer/ChIPseq2019/promoter_1kb_complement"

mkdir -p "${out_1kb}"
mkdir -p "${out_1kb_compl}"

out_1kb_bg_none="${out_1kb}/1kb_bg_none"
out_1kb_bg_all="${out_1kb}/1kb_bg_Eo_T_all_peaks"

out_1kb_compl_bg_none="${out_1kb}/1kb_compl_bg_none"
out_1kb_compl_bg_all="${out_1kb}/1kb_compl_bg_Eo_T_all_peaks"

mkdir -p "${out_1kb_bg_none}"
mkdir -p "${out_1kb_bg_all}"
mkdir -p "${out_1kb_compl_bg_none}"
mkdir -p "${out_1kb_compl_bg_all}"

source activate homer

findMotifsGenome.pl "${promoters_1kb}" mm10 "${out_1kb_bg_none}"/ -size 200 -p 6 && \
findMotifsGenome.pl "${promoters_1kb_complement}" mm10 "${out_1kb_compl_bg_none}"/ -size 200 -p 6 && \
findMotifsGenome.pl "${promoters_1kb}" mm10 "${out_1kb_bg_all}"/ -bg "${all}" -size 200 -p 6 && \
findMotifsGenome.pl "${promoters_1kb_complement}" mm10 "${out_1kb_compl_bg_all}"/ -bg "${all}" -size 200 -p 6


###############################################################################









