# qrupdate
A tiny R package which implements fast update of QR decomposition of a matrix.

These operations are implemented:

    -   Add a column to QR factorization
    -   Remove a column from QR factorization
    -   Add a row to QR factorization
    -   Remove a row from QR factorization
    -   Solve a system of order n by n for given QR factorization

## How to install
```r
# install.packages("devtools")
devtools::install_github("furiousbean/qrupdate")
```

## Example of usage
### Slow method (but with input checking)
```r
updateobj <- qrupdate(qr(matrix(0, nrow = 1, ncol = 1)))
updateobj <- addcolumn(updateobj, 2, 1)
updateobj <- addrow(updateobj, 2, c(1, 0))
x <- solve(updateobj, c(5, 3))
print(x) #(3, 5)
```

### Fast method (unsafe)
```r
R <- matrix(0, nrow = 1)
Q <- matrix(1, nrow = 1)
u <- .Call("miqraddcolumn", R, Q, as.integer(2), 1, PACKAGE = "qrupdate")
u <- .Call("miqraddrow", u$R, u$Q, as.integer(2), c(1, 0), PACKAGE = "qrupdate")
x <- .Call("miqrsolve", u$R, u$Q, c(5, 3), PACKAGE = "qrupdate")
cleanmem()
print(x) #(3, 5)
```

## References
Hammarling, Sven, and Craig Lucas. "Updating the QR factorization and the least squares problem." (2008). Available at http://eprints.ma.man.ac.uk/1192/01/covered/MIMS_ep2008_111.pdf

## License
Beerware license

## TODO
    -   Manual
    -   Bug fixes
    -   Implement fast manipulation  of blocks of columns/rows
    -   Tests
    -   Checking
    -   ???
    -   profit
