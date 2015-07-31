
#PACKAGE FILE
# QR add/remove of row/column wrappers
# see http://eprints.ma.man.ac.uk/1192/01/covered/MIMS_ep2008_111.pdf
# for the theory of QR factorization update

qrupdate <- function(...) UseMethod("qrupdate")
deleterow <- function(...) UseMethod("deleterow")
addrow <- function(...) UseMethod("addrow")
deletecolumn <- function(...) UseMethod("deletecolumn")
addcolumn <- function(...) UseMethod("addcolumn")

qrupdate.default <- function(...) 
    stop("qrupdate() takes only numeric 'qr' object!")

qrupdate.qr <- function(decomp) {
    if (!is.numeric(decomp$qr))
        stop("Complex matrices are not supported now")
    structure(list(Q = qr.Q(decomp, complete = TRUE), 
                   R = qr.R(decomp, complete = TRUE)), class = "qrupdate")
}

deleterow.qrupdate <- function(prevdec, k, clean = TRUE) {
    k <- as.integer(k)
    n <- ncol(prevdec$R)
    m <- nrow(prevdec$R)
    
    if (!(k >= 1 & k <= m)) stop("K isn't in 1:nrow(X)")
    
    result <- .Call("miqrdeleterow", prevdec$R, prevdec$Q, k, 
                    PACKAGE = "qrupdate")
    if (clean) qrupdate::cleanmem()
    class(result) <- "qrupdate"
    result
}

addrow.qrupdate <- function(prevdec, k, u, clean = TRUE) {
    stopifnot(is.numeric(u))
    stopifnot(is.finite(u))
    k <- as.integer(k)
    n <- ncol(prevdec$R)
    m <- nrow(prevdec$R)
    
    if (!(k >= 1 & k <= m + 1)) stop("k isn't in 1:(nrow(X) + 1)")
    if (length(u) != n) stop("Check the length of vector u")
    
    result <- .Call("miqraddrow", prevdec$R, prevdec$Q, k, u, 
                    PACKAGE = "qrupdate")
    if (clean) qrupdate::cleanmem()
    class(result) <- "qrupdate"
    result
}

deletecolumn.qrupdate <- function(prevdec, k, clean = TRUE) {
    k <- as.integer(k)    
    n <- ncol(prevdec$R)
    m <- nrow(prevdec$R)
    
    if (!(k >= 1 & k <= n)) stop("k isn't in 1:ncol(X)")
    result <- .Call("miqrdeletecolumn", prevdec$R, prevdec$Q, k, 
                    PACKAGE = "qrupdate")
    if (clean) qrupdate::cleanmem()
    class(result) <- "qrupdate"
    result
}

addcolumn.qrupdate <- function(prevdec, k, u, clean = TRUE) {
    stopifnot(is.numeric(u))
    stopifnot(is.finite(u))
    k <- as.integer(k)    
    
    n <- ncol(prevdec$R)
    m <- nrow(prevdec$R)
    
    if (!(k >= 1 & k <= n + 1)) stop("k isn't in 1:(nrow(X) + 1)")
    if (length(u) != m) stop("Check the length of vector u")
    
    result <- .Call("miqraddcolumn", prevdec$R, prevdec$Q, k, u, 
                    PACKAGE = "qrupdate")
    if (clean) qrupdate::cleanmem()
    class(result) <- "qrupdate"
    result
}

solve.qrupdate <- function(decomp, b) { 
    stopifnot(is.numeric(b))
    stopifnot(is.finite(b))
    if (nrow(decomp$R) != ncol(decomp$R)) stop("Sorry, R must be square matrix")
    if (nrow(decomp$R) != length(b)) stop("Check the length of vector b")
    
    .Call("miqrsolve", decomp$R, decomp$Q, b, PACKAGE = "qrupdate")
}

cleanmem <- function() {
    .Call("miclean")
}
