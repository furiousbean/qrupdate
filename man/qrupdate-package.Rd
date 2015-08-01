\name{qrupdate-package}
\alias{qrupdate-package}
\alias{qrupdate}
\docType{package}
\title{
Fast QR factorization updating operations
}
\description{
A tiny package which implements fast addition/deletion of columns/rows in QR factorization of a matrix.
}
\details{
The package uses object of class \emph{qrupdate} to hold full QR factorization of a real matrix. You need to call \code{\link[base:qr]{qr}} function to obtain factorization, then use \code{\link[qrupdate:qrupdate.qr]{qrupdate}} function to get \emph{qrupdate} object from a \emph{qr} object.

You can use \code{\link[qrupdate:addcolumn]{addcolumn}}, \code{\link[qrupdate:addrow]{addrow}}, \code{\link[qrupdate:deletecolumn]{deletecolumn}}, \code{\link[qrupdate:deleterow]{deleterow}} functions to modify the QR factorization, \code{\link[qrupdate:solve]{solve}} function to solve the system of linear equations given by QR factorization of square matrix.
}
\references{
\emph{Hammarling, Sven and Craig, Lucas}. Updating the QR factorization and the least squares problem. (2008). \emph{Available at} http://eprints.ma.man.ac.uk/1192/01/covered/MIMS_ep2008_111.pdf
}
\keyword{ package }
\seealso{
\code{\link[qrupdate:qrupdate.qr]{qrupdate.qr}}, \code{\link[qrupdate:addcolumn]{addcolumn}}, \code{\link[qrupdate:addrow]{addrow}}, \code{\link[qrupdate:deletecolumn]{deletecolumn}}, \code{\link[qrupdate:deleterow]{deleterow}}, \code{\link[qrupdate:solve]{solve}}, \code{\link[qrupdate:cleanmem]{cleanmem}}.
}
\examples{
#simply add one column, one row and solve a system
updateobj <- qrupdate(qr(matrix(0, nrow = 1, ncol = 1)))
updateobj <- addcolumn(updateobj, 2, 1)
updateobj <- addrow(updateobj, 2, c(1, 0))
x <- solve(updateobj, c(5, 3))
print(x) #(3, 5)
}
