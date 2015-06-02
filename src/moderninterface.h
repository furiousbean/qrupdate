#ifndef MODERNINTERFACE_H
#define MODERNINTERFACE_H

#include <R.h>
#include <Rinternals.h>

SEXP modernfastqrsolve(SEXP Rexp, SEXP Qexp, SEXP bexp);
SEXP modernfastqrdeleterow(SEXP Rexp, SEXP Qexp, SEXP kexp);
SEXP modernfastqraddrow(SEXP Rexp, SEXP Qexp, SEXP kexp, SEXP uexp);
SEXP modernfastqrdeletecolumn(SEXP Rexp, SEXP Qexp, SEXP kexp);
SEXP modernfastqraddcolumn(SEXP Rexp, SEXP Qexp, SEXP kexp, SEXP uexp);
SEXP modernclean();

#endif //MODERNINTERFACE_H