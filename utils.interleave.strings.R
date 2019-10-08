
# interleaveVecs <- function(vec1, vec2, acc=c()) {
#   len1 <- length(vec1)
#   len2 <- length(vec2)
#   if (len1 == 0) {
#     if (len2 == 0) {
#       return(acc)
#     } else {
#       return(interleaveVecs(vec1, vec2[2:len2], acc=c(acc, vec2[1])))
#     }
#   } else {
#     if (len2 == 0) {
#       return(interleaveVecs(vec1[2:len1], vec2, acc=c(acc, vec1[1])))
#     } else {
#       return(interleaveVecs(vec1[2:len1], vec2[2:len2], acc=c(acc, vec1[1], vec2[1])))
#     }
#   }
# } # wrong



# # 
# interleaveVecs <- function(vec1, vec2, acc=c()) {
#   if (len1 == 0) {
#     if (len2 == 0) {
#       return(acc)
#     } else {
#       return(interleaveVecs(vec1, vec2[2:len2], acc=c(acc, vec2[1])))
#     }
#   } else {
#     if (len2 == 0) {
#       return(interleaveVecs(vec1[2:len1], vec2, acc=c(acc, vec1[1])))
#     } else {
#       return(interleaveVecs(vec1[2:len1], vec2[2:len2], acc=c(acc, vec1[1], vec2[1])))
#     }
#   }
# } # wrong! 
# > interleaveVecs(c("a", "b", "c"), c(1, 2, 3 ,4, 5))
# Error: C stack usage  7971300 is too close to the limit


########################
# interleave for vecs and lists
########################
source("~/Dropbox/R/central-scripts/utils.lisp.R") # later: make package for lisp in R
.interleave <- function(l1, l2, acc = list()) {
  if (null(l1) && null(l2)) {
    reverse(acc)
  } else if (null(l1)) {
    .interleave(cdr(l1), cdr(l2), acc = cons(car(l2), acc))
  } else if (null(l2)) {
    .interleave(cdr(l1), cdr(l2), acc = cons(car(l1), acc))
  } else {
    .interleave(cdr(l1), cdr(l2), acc = cons(car(l2), cons(car(l1), acc)))
  }
}

interleaveVecs <- function(vec1, vec2) {
  l1 <- as.list(vec1)
  l2 <- as.list(vec2)
  unlist(.interleave(l1, l2))
} # works!

###############################
# vector substring regex split
###############################

# aim was to avoid `strsplit()` which is very slow!
# splitting a string should be done using substr() and its vectorized form
# which is super fast.


vecSubstr <- Vectorize(substr, USE.NAMES=F)

splitSingleStr <- function(str, pattern, keep=FALSE) {
  # Return vector of one.string split into strings vector by pattern
  first.pos <- 1
  end.pos <- nchar(str)
  pattern.starts <- gregexpr(pattern, str)[[1]]
  pattern.lengths <- attr(pattern.starts, "match.length")
  pattern.ends <- pattern.starts + pattern.lengths - 1
  rest.starts <- c(first.pos, pattern.ends + 1)
  rest.ends <- c(pattern.starts - 1, end.pos)
  if (keep) {
    interleaveVecs(vecSubstr(str, rest.starts, rest.ends), 
                   vecSubstr(str, pattern.starts, pattern.ends))
  } else {
    # key is that substr(x, 1, 0) gives ""
    vecSubstr(str, rest.starts, rest.ends)
  }
}

# called `splitStrToVec` orgiinally
splitStr <- function(str, pattern, keep=FALSE,
                          start = 1, end = NULL,
                          fromN = NULL,
                          maxN = NULL) {
  # much faster than strsplit()
  if (start > 1) {
    startPiece <- substr(str, 1, start-1)
  } else {
    startPiece <- ""
  }
  if (!is.null(end)) {
    endPiece <- substr(str, end+1, nchar(str))
  } else {
    endPiece <- ""
    end <- nchar(str)
  }
  res <- splitSingleStr(substr(str, start, end), pattern, keep)
  res[1] <- paste0(startPiece, res[1])
  res[length(res)] <- paste0(res[length(res)], endPiece)
  res
}
## this is NOT ready and inverleave_strsplit.R has better version
splitStrToVec <- splitStr


# str <- "abcdabababcdefababdgf"
# pattern <- "(ab)+"
# # if we know where all first occurrence positions are from pattern
# # and the end positions from pattern,
# then 1 to first occurrence - 1 is first part,
# c(first.pos, next.ends) to c(next.starts, end.pos)
# 
# interleaveVecs(vecSubstr(str, res.starts, rest.ends), vecSubstr(str, pattern.starts, pattern.ends))

# # this works as expected! :)
# str <- "abcdabgababcdefababdgfab"
# pattern <- "(ab(g)*ab)+"
# > splitStr(str, pattern, keep=T, start=4, end=11)
# [1] "abcd"            "abgab"           "abcdefababdgfab"
# > splitStr(str, pattern, keep=T, start=5, end=20)
# [1] "abcd"   "abgab"  "abcdef" "abab"   "dgfab" 
# > splitStr(str, pattern, keep=T, start=6, end=20)
# [1] "abcdabg" "abab"    "cdef"    "abab"    "dgfab"  
# > splitStr(str, pattern, keep=T, start=9, end=20)
# [1] "abcdabgababcdef" "abab"            "dgfab"      













