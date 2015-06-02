# qrupdate
A tiny R package which implements fast update of QR decomposition of a matrix.

These operations are implemented:

    -   Add a column to QR decomposition
    -   Remove a column from QR decomposition
    -   Add a row to QR decomposition
    -   Remove a row from QR decomposition
    -   Solve a system of order n by n for given QR decomposition

## How to install
```r
# install.packages("devtools")
devtools::install_github("furiousbean/qrupdate")
```

## Example of usage
### Slow method (but with input checking)
```r
R <- matrix(0, nrow = 1)
Q <- matrix(0, nrow = 1)
u <- qr.add.column(Q, R, 2, 1)
u <- qr.add.row(Q, R, 2, c(1, 0))
x <- solve.qrupdate(Q, R, c(5, 3))
print(x) #(3, 5)
```

### Fast method (unsafe)
```r
R <- matrix(0, nrow = 1)
Q <- matrix(1, nrow = 1)
u <- .Call("modernfastqraddcolumn", R, Q, as.integer(2), 1, PACKAGE = "qrupdate")
u <- .Call("modernfastqraddrow", u$R, u$Q, as.integer(2), c(1, 0), PACKAGE = "qrupdate")
x <- .Call("modernfastqrsolve", u$R, u$Q, c(5, 3), PACKAGE = "qrupdate")
qrupdate.clean()
print(x) #(3, 5)
```

## References
Hammarling, Sven, and Craig Lucas. "Updating the QR factorization and the least squares problem." (2008). Available at http://eprints.ma.man.ac.uk/1192/01/covered/MIMS_ep2008_111.pdf

## License
Beerware license

## TODO
    -   Manual
    -   Bug fixes
    -   Working with block of columns/rows
    -   Tests
    -   Checking
    -   Maybe use "qr" object?
    -   ???
    -   Maybe profit
