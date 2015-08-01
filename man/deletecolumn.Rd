\name{deletecolumn}
\alias{deletecolumn}
\alias{deletecolumn.qrupdate}
\title{Delete column from a QR factorization}

\description{
  Remove a given column from a given position in QR factorization.
}

\usage{
\method{deletecolumn}{deletecolumn}(prevdec, k, clean = TRUE)
}

\arguments{
  \item{prevdec}{object of class qrupdate containing previous QR factorization}
  \item{k}{index of a column to delete, 1 <= k <= ncol(prevdec \$ R)}
  \item{clean}{clean extra inner memory? \code{TRUE} by default}  
}

\details{
  You should use \code{clean = FALSE} and then call \code{\link[qrupdate:cleanmem]{cleanmem}} after all operations to increase calculation speed.

  Use GitHub manual for the example of significant speed improvement.
}

\value{
  An object of class \emph{qrupdate} which contains QR decomposition of an updated matrix.
}

\examples{
#simply delete one column and print a resultant matrix
updateobj <- qrupdate(qr(matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2)))
updateobj <- deletecolumn(updateobj, 1)
print(updateobj$Q \%*\% updateobj$R)
#(3)
#(4)
}

\seealso{
  \code{\link[qrupdate:addcolumn]{addcolumn}}.
}
