# ##############################################################################
# utility functions to imitate lisp functions in R
# ##############################################################################

######################################################
## NIL == list()
######################################################



######################################################
## predicates
######################################################

# cl's null checking for if a list is '()/NIL
# in contrast: R's `is.null` is checking whether sth is `NULL`` or not
null <- function(l) is.list(l) && length(l) == 0

# atom
# atom <- is.atom

# pair - asks for a `cons`
# corresponds to `is.list`

listp <- function(l) is.list(l) && l[[length(l)]] == list()


######################################################
## conditions testing NIL == FALSE else TRUE
######################################################

nn <- function(x) !null(x)

################
# car first
################

car <- function(l) {
  if(is.list(l)) {
    if (null(l)) {
      list()
    } else {
      l[[1]]
    }
  } else {
    stop("Not a list.")
  }
}

first <- car

################
# cdr rest
################

cdr <- function(l) {
  if (is.list(l)) {
    if (null(l) || length(l) == 1) {
      list()
    } else {
      l[2:length(l)]
    }
  } else {
    stop("Not a list.")
  }
}

rest <- cdr


################
# caar
# cadr
# cdar
# cddr
#
# caaar
# caadr
# cadar
# cdaar
# caddr
# cdadr
# cddar
# cdddr
# 
# caaaar
# caaadr
# caadar
# cadaar
# cdaaar
# caaddr
# cadadr
# cdaadr
# cdadar
# cddaar
# caddar
# cadddr
# cdaddr
# cddadr
# cdddar
# cddddr
################

caar <- function(l) car(car(l))
cadr <- function(l) car(cdr(l))
cdar <- function(l) cdr(car(l))
cddr <- function(l) cdr(cdr(l))

caaar <- function(l) car(car(car(l)))
caadr <- function(l) car(car(cdr(l)))
cadar <- function(l) car(cdr(car(l)))
cdaar <- function(l) cdr(car(car(l)))
caddr <- function(l) car(cdr(cdr(l)))
cdadr <- function(l) cdr(car(cdr(l)))
cddar <- function(l) cdr(cdr(car(l)))
cdddr <- function(l) cdr(cdr(cdr(l)))

caaaar <- function(l) car(car(car(car(l))))
caaadr <- function(l) car(car(car(cdr(l))))
caadar <- function(l) car(car(cdr(car(l))))
cadaar <- function(l) car(cdr(car(car(l))))
cdaaar <- function(l) cdr(car(car(car(l))))
caaddr <- function(l) car(car(cdr(cdr(l))))
cadadr <- function(l) car(cdr(car(cdr(l))))
cdaadr <- function(l) cdr(car(car(cdr(l))))
cdadar <- function(l) cdr(car(cdr(car(l))))
cddaar <- function(l) cdr(cdr(car(car(l))))
caddar <- function(l) car(cdr(cdr(car(l))))
cadddr <- function(l) car(cdr(cdr(cdr(l))))
cdaddr <- function(l) cdr(car(cdr(cdr(l))))
cddadr <- function(l) cdr(cdr(car(cdr(l))))
cdddar <- function(l) cdr(cdr(cdr(car(l))))
cddddr <- function(l) cdr(cdr(cdr(cdr(l))))

################
# cons
################

cons <- function(el, l) {
  if (!typeof(l) == "list") {
    stop("Cons on non-list! Please provide a list as second argument of cons!")
  } else {
    c(list(el), l)
  }
}


################
# list
################

# last element of a list is NIL


################
# list*
################


################
# reverse
################

reverse <- rev




######################################################
## flatten
######################################################

flatten <- function(nested.list, predicate.element.type = is.atom) {
  # Implemented Lisp's flatten tail call recursively -> `..flatten()`
  # `predicate.element.type`` determines when to stop flattening.
  # here: if the list element is an atom `is.atom`, the element is collected into the
  # result list, otherwise further flattened.
  ..flatten <- function(l, acc.l) { 
    if (null(l)) {
      acc.l
    } else if (predicate.element.type(l)) {   
      acc.l[[length(acc.l) + 1]] <- l
      acc.l # kind of (list* l acc.l)
    } else {
      ..flatten(car(l), ..flatten(cdr(l), acc.l))
    }
  }
  ..flatten(nested.list, list())
}

flatten.to.df.list <- function(nested.df.list) {
  # flatten nested data frame list constructs to flat
  # data frame list.
  flatten(nested.df.list, is.data.frame)
}
