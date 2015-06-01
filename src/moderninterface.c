#include <R.h>
#include <Rinternals.h>
#include "fastqrupdate.h"
#include "string.h"
#include "stdlib.h"

SEXP modernfastqrsolve(SEXP Rexp, SEXP Qexp, SEXP bexp) {
    int m = INTEGER(getAttrib(Qexp, R_DimSymbol))[0];
    SEXP xexp = PROTECT(allocVector(REALSXP, m));
    fastqrsolve(&m, REAL(Rexp), REAL(Qexp), REAL(bexp), REAL(xexp));
    UNPROTECT(1);
    
    return xexp;
}

SEXP modernfastqrdeleterow(SEXP Rexp, SEXP Qexp, SEXP kexp) {
    int m = INTEGER(getAttrib(Rexp, R_DimSymbol))[0];
    int n = INTEGER(getAttrib(Rexp, R_DimSymbol))[1];
    int pk = INTEGER(kexp)[0];
    double* workR = (double*) calloc(m * n, sizeof(double));
    double* workQ = (double*) calloc(m * m, sizeof(double));
    memcpy(workR, REAL(Rexp), m * n * sizeof(double));
    memcpy(workQ, REAL(Qexp), m * m * sizeof(double));
    SEXP ansR = PROTECT(allocMatrix(REALSXP, m - 1, n));
    SEXP ansQ = PROTECT(allocMatrix(REALSXP, m - 1, m - 1));
    fastqrdeleterow(&m, &n, &pk, workR, workQ, REAL(ansR), REAL(ansQ));
    free(workR);
    free(workQ);
    
    SEXP ansexp = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ansexp, 0, ansQ);
    SET_VECTOR_ELT(ansexp, 1, ansR);
    
    SEXP ansnames = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(ansnames, 0, mkChar("Q"));
    SET_STRING_ELT(ansnames, 1, mkChar("R"));
    setAttrib(ansexp, R_NamesSymbol, ansnames);

    UNPROTECT(4);
    return ansexp;
}

SEXP modernfastqraddrow(SEXP Rexp, SEXP Qexp, SEXP kexp, SEXP uexp) {
    int m = INTEGER(getAttrib(Rexp, R_DimSymbol))[0];
    int n = INTEGER(getAttrib(Rexp, R_DimSymbol))[1];
    int pk = INTEGER(kexp)[0];
    double* workR = (double*) calloc(m * n, sizeof(double));
    double* workQ = (double*) calloc(m * m, sizeof(double));
    memcpy(workR, REAL(Rexp), m * n * sizeof(double));
    memcpy(workQ, REAL(Qexp), m * m * sizeof(double));
    SEXP ansR = PROTECT(allocMatrix(REALSXP, m + 1, n));
    SEXP ansQ = PROTECT(allocMatrix(REALSXP, m + 1, m + 1));
    fastqraddrow(&m, &n, &pk, workR, workQ, REAL(uexp), REAL(ansR), REAL(ansQ));
    free(workR);
    free(workQ);
    
    SEXP ansexp = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ansexp, 0, ansQ);
    SET_VECTOR_ELT(ansexp, 1, ansR);
    
    SEXP ansnames = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(ansnames, 0, mkChar("Q"));
    SET_STRING_ELT(ansnames, 1, mkChar("R"));
    setAttrib(ansexp, R_NamesSymbol, ansnames);

    UNPROTECT(4);
    return ansexp;
}

SEXP modernfastqrdeletecolumn(SEXP Rexp, SEXP Qexp, SEXP kexp) {
    int m = INTEGER(getAttrib(Rexp, R_DimSymbol))[0];
    int n = INTEGER(getAttrib(Rexp, R_DimSymbol))[1];
    int pk = INTEGER(kexp)[0];
    double* workR = (double*) calloc(m * n, sizeof(double));

    memcpy(workR, REAL(Rexp), m * n * sizeof(double));
    SEXP ansR = PROTECT(allocMatrix(REALSXP, m, n - 1));
    SEXP ansQ = PROTECT(allocMatrix(REALSXP, m, m));
    memcpy(REAL(ansQ), REAL(Qexp), m * m * sizeof(double));
    fastqrdeletecolumn(&m, &n, &pk, workR, REAL(ansQ), REAL(ansR));
    free(workR);
    
    SEXP ansexp = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ansexp, 0, ansQ);
    SET_VECTOR_ELT(ansexp, 1, ansR);
    
    SEXP ansnames = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(ansnames, 0, mkChar("Q"));
    SET_STRING_ELT(ansnames, 1, mkChar("R"));
    setAttrib(ansexp, R_NamesSymbol, ansnames);

    UNPROTECT(4);
    return ansexp;
}

SEXP modernfastqraddcolumn(SEXP Rexp, SEXP Qexp, SEXP kexp, SEXP uexp) {
    int m = INTEGER(getAttrib(Rexp, R_DimSymbol))[0];
    int n = INTEGER(getAttrib(Rexp, R_DimSymbol))[1];
    int pk = INTEGER(kexp)[0];
    double* workR = (double*) calloc(m * n, sizeof(double));

    memcpy(workR, REAL(Rexp), m * n * sizeof(double));
    SEXP ansR = PROTECT(allocMatrix(REALSXP, m, n + 1));
    SEXP ansQ = PROTECT(allocMatrix(REALSXP, m, m));
    memcpy(REAL(ansQ), REAL(Qexp), m * m * sizeof(double));
    fastqraddcolumn(&m, &n, &pk, workR, REAL(ansQ), REAL(uexp), REAL(ansR));
    free(workR);
    
    SEXP ansexp = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ansexp, 0, ansQ);
    SET_VECTOR_ELT(ansexp, 1, ansR);
    
    SEXP ansnames = PROTECT(allocVector(STRSXP, 2));
    SET_STRING_ELT(ansnames, 0, mkChar("Q"));
    SET_STRING_ELT(ansnames, 1, mkChar("R"));
    setAttrib(ansexp, R_NamesSymbol, ansnames);

    UNPROTECT(4);
    return ansexp;
}

