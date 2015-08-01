\name{solve}
\alias{solve.qrupdate}
\title{Solve a system of linear equations}

\description{
   Solve a system of linear equations given by QR factorization of a square matrix
}

\usage{
\method{solve}{solve.qrupdate}(decomp, b)
}

\arguments{
  \item{decomp}{object of class qrupdate containing QR factorization of a square matrix A}
  \item{b}{numeric vector in the right part of equation}
}

\details{
  Method solves following linear equation: \code{A \%*\% x = b}, where \code{A} is numeric square matrix of full rank, \code{b} is the vector of the same length as the length of the side of \code{A}.

  Computation time is O(N^2), where N is the length of the side of \code{A}.
}

\value{
  Numeric vector x.
}

\examples{
#simply add one column, one row and solve a system
updateobj <- qrupdate(qr(matrix(0, nrow = 1, ncol = 1)))
updateobj <- addcolumn(updateobj, 2, 1)
updateobj <- addrow(updateobj, 2, c(1, 0))
x <- solve(updateobj, c(5, 3))
print(x) #(3, 5)
}

\seealso{
  \code{\link[qrupdate:qrupdate]{qrupdate}}, \code{\link[qrupdate:addcolumn]{addcolumn}}, \code{\link[qrupdate:addrow]{addrow}}.
}
