
#PACKAGE FILE
# QR add/remove row/column functions
# see http://eprints.ma.man.ac.uk/1192/01/covered/MIMS_ep2008_111.pdf

qr.delete.row <- function(Q, R, k) {
    stopifnot(is.numeric(Q))
    stopifnot(is.numeric(R))
    stopifnot(is.finite(Q))
    stopifnot(is.finite(R))
    stopifnot(is.matrix(Q))
    stopifnot(is.matrix(R))
    stopifnot(is.integer(k))
    n <- ncol(R)
    m <- nrow(R)
    
    if (!(k >= 1 & k <= m)) stop("K isn't in 1:nrow(X)")
    if (m == 1) stop("resultant matrix is empty")
    if (nrow(R) != nrow(Q) | nrow(R) != ncol(Q)) 
        stop("check the dimensions of matrices")
    
    .Call("modernfastqrdeleterow", R, Q, k, PACKAGE = "qrupdate")
}

qr.add.row <- function(Q, R, k, u) {
    stopifnot(is.numeric(Q))
    stopifnot(is.numeric(R))
    stopifnot(is.numeric(u))
    stopifnot(is.finite(Q))
    stopifnot(is.finite(R))
    stopifnot(is.finite(u))
    stopifnot(is.matrix(Q))
    stopifnot(is.matrix(R))
    stopifnot(is.integer(k))
    n <- ncol(R)
    m <- nrow(R)
    
    if (!(k >= 1 & k <= m + 1)) stop("k isn't in 1:(nrow(X) + 1)")
    if (nrow(R) != nrow(Q) | nrow(R) != ncol(Q)) 
        stop("check the dimensions of matrices")
    if (length(u) != n) stop("check the length of vector u")
    
    .Call("modernfastqraddrow", R, Q, k, u, PACKAGE = "qrupdate")
}

qr.delete.column <- function(Q, R, k) {
    stopifnot(is.numeric(Q))
    stopifnot(is.numeric(R))
    stopifnot(is.finite(Q))
    stopifnot(is.finite(R))
    stopifnot(is.matrix(Q))
    stopifnot(is.matrix(R))
    stopifnot(is.integer(k))
    n <- ncol(R)
    m <- nrow(R)
    
    if (!(k >= 1 & k <= n)) stop("k isn't in 1:ncol(X)")
    if (n == 1) stop("resultant matrix is empty")
    if (nrow(R) != nrow(Q) | nrow(R) != ncol(Q)) 
        stop("check the dimensions of matrices")
    .Call("modernfastqrdeletecolumn", R, Q, k, PACKAGE = "qrupdate")
}

qr.add.column <- function(Q, R, k, u) {
    stopifnot(is.numeric(Q))
    stopifnot(is.numeric(R))
    stopifnot(is.numeric(u))
    stopifnot(is.finite(Q))
    stopifnot(is.finite(R))
    stopifnot(is.finite(u))
    stopifnot(is.matrix(Q))
    stopifnot(is.matrix(R))
    stopifnot(is.integer(k))
    n <- ncol(R)
    m <- nrow(R)
    
    if (!(k >= 1 & k <= n + 1)) stop("k isn't in 1:(nrow(X) + 1)")
    if (nrow(R) != nrow(Q) | nrow(R) != ncol(Q)) 
        stop("check the dimensions of matrices")
    if (length(u) != m) stop("check the length of vector u")
    
    .Call("modernfastqraddcolumn", R, Q, k, u, PACKAGE = "qrupdate")
}

solve.myqr <- function(Q, R, b) { 
    stopifnot(is.numeric(Q))
    stopifnot(is.numeric(R))
    stopifnot(is.numeric(b))
    stopifnot(is.finite(Q))
    stopifnot(is.finite(R))
    stopifnot(is.finite(b))
    stopifnot(is.matrix(Q))
    stopifnot(is.matrix(R))   
    if (nrow(R) != ncol(R)) stop("Sorry, R must be square")
    if (nrow(R) != nrow(Q) | nrow(R) != ncol(Q)) 
        stop("check the dimensions of matrices")
    if (nrow(R) != length(b)) stop("check the length of vector b")
    
    .Call("modernfastqrsolve", R, Q, b)
}