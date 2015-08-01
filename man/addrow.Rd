\name{addrow}
\alias{addrow}
\alias{addrow.qrupdate}
\title{Add row to a QR factorization}

\description{
  Insert a given row into a given position in QR factorization.
}

\usage{
\method{addrow}{addrow}(prevdec, k, u, clean = TRUE)
}

\arguments{
  \item{prevdec}{object of class qrupdate containing previous QR factorization}
  \item{k}{index of a new column, 1 <= k <= nrow(prevdec \$ R) + 1}
  \item{u}{vector containing new row}
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
#simply add one row and print a resultant matrix
updateobj <- qrupdate(qr(matrix(c(1, 2), nrow = 1, ncol = 2)))
updateobj <- addrow(updateobj, 2, c(3, 4))
print(updateobj$Q \%*\% updateobj$R)
#(1, 2)
#(3, 4)
}

\seealso{
  \code{\link[qrupdate:deleterow]{deleterow}}.
}
