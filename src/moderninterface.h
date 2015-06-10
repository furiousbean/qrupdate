#ifndef MODERNINTERFACE_H
#define MODERNINTERFACE_H

#include <R.h>
#include <Rinternals.h>

SEXP miqrsolve(SEXP Rexp, SEXP Qexp, SEXP bexp);
SEXP miqrdeleterow(SEXP Rexp, SEXP Qexp, SEXP kexp);
SEXP miqraddrow(SEXP Rexp, SEXP Qexp, SEXP kexp, SEXP uexp);
SEXP miqrdeletecolumn(SEXP Rexp, SEXP Qexp, SEXP kexp);
SEXP miqraddcolumn(SEXP Rexp, SEXP Qexp, SEXP kexp, SEXP uexp);
SEXP miclean();

#endif //MODERNINTERFACE_H