#####################################################################
# Load openxlsx
#####################################################################

require(openxlsx)


c.dir <- "~/Dropbox/R/central-scripts"

source(file.path(c.dir, "utils.lists.R"))

#####################################################################
# read xlsx files to dfs list
#####################################################################

xlsx2df.list <- function(xlsx.path, rowNames = TRUE, colNames = TRUE, ...) {
  wb <- loadWorkbook(xlsx.path)
  sheetNames <- names(wb)
  res <- lapply(sheetNames, function(sheetName) {
    read.xlsx(wb, sheet = sheetName, rowNames = rowNames, colNames = colNames, ...)
  })
  names(res) <- sheetNames
  res
}


#####################################################################
# printing dfs to xlsx files
#####################################################################

withNames <- function(...) {
  # Returns a list constructed by using 
  # alterning name, obj, name, obj arguments
  p.l <- list(...)
  len <- length(p.l)
  if (len %% 2 == 1) {
    stop("withNames call with odd numbers of arguments")
    print()
  }
  seconds <- p.l[seq(2, len, 2)]
  firsts <- p.l[seq(1, len, 2)]
  names(seconds) <- unlist(firsts)
  seconds
}

write.dfs <- function(df.list, fpath) {
  wb <- createWorkbook()
  Map(function(data, name) {
    addWorksheet(wb, name)
    writeData(wb, name, data, rowNames = TRUE, colNames = TRUE)
  }, df.list, names(df.list))
  saveWorkbook(wb, file = fpath, overwrite = TRUE)
}


#####################################################################
# select a df by list of rownames -> df.list
#####################################################################

select.by.names.vec.list <- function(data.df, list.of.names.vecs) {
  lapply(list.of.names.vecs, function(names.vec) data.df[names.vec, ])
}

#####################################################################
# Return for a df, its sortings by l2FC, p-val and rownames
#####################################################################

df.sortings <- function(DE.res.df) {
  res <- list(DE.res.df[order(DE.res.df$log2FoldChange, decreasing = TRUE), ],
              DE.res.df[order(DE.res.df$pvalue, decreasing = FALSE), ],
              DE.res.df[order(rownames(DE.res.df), decreasing = FALSE), ])
  names(res) <- c("l2FC.srt", "p.srt", "names.srt")
  res
}

#####################################################################
# print DE sortings
#####################################################################

print.DE.sortings <- function(DE.list, fname, dir = ".") {
  DE.list.with.sortings <- lapply(once.flatten(lapply(DE.list, df.sortings)), as.data.frame)
  write.dfs(DE.list.with.sortings,
            file.path(dir, fname))
}

#####################################################################
# print cnts selections of DE 
#####################################################################

print.cnts.DE.sortings <- function(cnts, DE.list, fname, dir = ".") {
  DE.list.with.sortings <- lapply(once.flatten(lapply(DE.list, df.sortings)),
                                  as.data.frame)
  DE.list.with.sortings.names <- lapply(DE.list.with.sortings, rownames)
  DE.cnts.list.with.sortings <- select.by.names.vec.list(as.data.frame(cnts),
                                   DE.list.with.sortings.names)
  write.dfs(DE.cnts.list.with.sortings,
            file.path(dir, fname))
}
















