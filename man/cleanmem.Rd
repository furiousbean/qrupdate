\name{cleanmem}
\alias{cleanmem}
\title{Clean inner memory}

\description{
  Clean inner memory after series of matrix manipulations.
}

\usage{
\method{cleanmem}{cleanmem}()
}

\details{
  You should use this method to clean inner extra memory if you modifed \emph{qrupdate} object with \code{clean = FALSE} parameter.

  Use GitHub manual for the example of significant speed improvement.
}

\examples{
#simply add one column, one row and solve a system
updateobj <- qrupdate(qr(matrix(0, nrow = 1, ncol = 1)))
updateobj <- addcolumn(updateobj, 2, 1, clean = FALSE)
updateobj <- addrow(updateobj, 2, c(1, 0), clean = FALSE)
x <- solve(updateobj, c(5, 3))
cleanmem()
print(x) #(3, 5)
}

\seealso{
  \code{\link[qrupdate:qrupdate]{qrupdate}}, \code{\link[qrupdate:addcolumn]{addcolumn}}, \code{\link[qrupdate:addrow]{addrow}}, \code{\link[qrupdate:solve]{solve.qrupdate}}.
}
