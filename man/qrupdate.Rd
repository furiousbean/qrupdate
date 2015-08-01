\name{qrupdate}
\alias{qrupdate-class}
\alias{qrupdate.qr}
\alias{qrupdate.default}
\title{qrupdate-class constructor}

\description{
  Obtain a qrupdate-class object from a qr object.
}

\usage{
\method{qrupdate}{qrupdate}(x)
}

\arguments{
  \item{x}{object of class qr}
}

\details{
  \emph{qrupdate} retrieves QR factorization of a matrix from a \emph{qr} object and stores it in a full form.
}

\value{
  An object of class \emph{qrupdate}, which is list containing \emph{R} and \emph{Q} matrices.
}

\seealso{
  \code{\link[qrupdate:solve]{solve}}.
}
