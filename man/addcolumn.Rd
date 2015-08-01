\name{addcolumn}
\alias{addcolumn}
\alias{addcolumn.qrupdate}
\title{Add column to a QR factorization}

\description{
  Insert a given column into a given position in QR factorization.
}

\usage{
\method{addcolumn}{addcolumn}(prevdec, k, u, clean = TRUE)
}

\arguments{
  \item{prevdec}{object of class qrupdate containing previous QR factorization}
  \item{k}{index of a new column, 1 <= k <= ncol(prevdec \$ R) + 1}
  \item{u}{vector containing new column}
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
#simply add one column and print a resultant matrix
updateobj <- qrupdate(qr(matrix(c(1, 2), nrow = 2, ncol = 1)))
updateobj <- addcolumn(updateobj, 2, c(3, 4))
print(updateobj$Q \%*\% updateobj$R)
#(1, 3)
#(2, 4)
}

\seealso{
  \code{\link[qrupdate:deletecolumn]{deletecolumn}}.
}
