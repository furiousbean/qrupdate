
#include "moderninterface.h"
#include <R.h>
#include <Rinternals.h>
#include "fastqrupdate.h"
#include "string.h"
#include "stdlib.h"

//const int THRESHOLD = 0;

static double* workR = NULL;
static double* workQ = NULL;
//static int workRsize = -1;
//static int workQsize = -1;

/*void workRrealloc(int m, int n) {
    if (m * n > workRsize) {
        workRsize = (m + THRESHOLD) * (n + THRESHOLD);
        workR = (double*) realloc(workR, workRsize * sizeof(double));
    }
}

void workQrealloc(int m) {
    if (m * m > workQsize) {
        workQsize = (m + THRESHOLD) * (m + THRESHOLD);
        workQ = (double*) realloc(workQ, workQsize * sizeof(double));
    }
}
*/
SEXP miqrsolve(SEXP Rexp, SEXP Qexp, SEXP bexp) {
    int m = INTEGER(getAttrib(Qexp, R_DimSymbol))[0];
    SEXP xexp = PROTECT(allocVector(REALSXP, m));
    fastqrsolve(&m, REAL(Rexp), REAL(Qexp), REAL(bexp), REAL(xexp));
    UNPROTECT(1);
    
    return xexp;
}

SEXP miqrdeleterow(SEXP Rexp, SEXP Qexp, SEXP kexp) {
    int m = INTEGER(getAttrib(Rexp, R_DimSymbol))[0];
    int n = INTEGER(getAttrib(Rexp, R_DimSymbol))[1];
    int pk = INTEGER(kexp)[0];

    workR = (double*) realloc(workR, m * n * sizeof(double));
    workQ = (double*) realloc(workQ, m * m * sizeof(double));
    //workRrealloc(m, n);
    //workQrealloc(m);
    
    memcpy(workR, REAL(Rexp), m * n * sizeof(double));
    memcpy(workQ, REAL(Qexp), m * m * sizeof(double));
    SEXP ansR = PROTECT(allocMatrix(REALSXP, m - 1, n));
    SEXP ansQ = PROTECT(allocMatrix(REALSXP, m - 1, m - 1));
    fastqrdeleterow(&m, &n, &pk, workR, workQ, REAL(ansR), REAL(ansQ));
    
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

SEXP miqraddrow(SEXP Rexp, SEXP Qexp, SEXP kexp, SEXP uexp) {
    int m = INTEGER(getAttrib(Rexp, R_DimSymbol))[0];
    int n = INTEGER(getAttrib(Rexp, R_DimSymbol))[1];
    int pk = INTEGER(kexp)[0];
    
    workR = (double*) realloc(workR, m * n * sizeof(double));
    workQ = (double*) realloc(workQ, m * m * sizeof(double));
    //workRrealloc(m, n);
    //workQrealloc(m);
    
    memcpy(workR, REAL(Rexp), m * n * sizeof(double));
    memcpy(workQ, REAL(Qexp), m * m * sizeof(double));
    SEXP ansR = PROTECT(allocMatrix(REALSXP, m + 1, n));
    SEXP ansQ = PROTECT(allocMatrix(REALSXP, m + 1, m + 1));
    fastqraddrow(&m, &n, &pk, workR, workQ, REAL(uexp), REAL(ansR), REAL(ansQ));
    
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

SEXP miqrdeletecolumn(SEXP Rexp, SEXP Qexp, SEXP kexp) {
    int m = INTEGER(getAttrib(Rexp, R_DimSymbol))[0];
    int n = INTEGER(getAttrib(Rexp, R_DimSymbol))[1];
    int pk = INTEGER(kexp)[0];
    workR = (double*) realloc(workR, m * n * sizeof(double));
    //workRrealloc(m, n);
    
    memcpy(workR, REAL(Rexp), m * n * sizeof(double));
    SEXP ansR = PROTECT(allocMatrix(REALSXP, m, n - 1));
    SEXP ansQ = PROTECT(allocMatrix(REALSXP, m, m));
    memcpy(REAL(ansQ), REAL(Qexp), m * m * sizeof(double));
    fastqrdeletecolumn(&m, &n, &pk, workR, REAL(ansQ), REAL(ansR));
    
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

SEXP miqraddcolumn(SEXP Rexp, SEXP Qexp, SEXP kexp, SEXP uexp) {
    int m = INTEGER(getAttrib(Rexp, R_DimSymbol))[0];
    int n = INTEGER(getAttrib(Rexp, R_DimSymbol))[1];
    int pk = INTEGER(kexp)[0];
    workR = (double*) realloc(workR, m * n * sizeof(double));
    //workRrealloc(m, n);
    
    memcpy(workR, REAL(Rexp), m * n * sizeof(double));
    SEXP ansR = PROTECT(allocMatrix(REALSXP, m, n + 1));
    SEXP ansQ = PROTECT(allocMatrix(REALSXP, m, m));
    memcpy(REAL(ansQ), REAL(Qexp), m * m * sizeof(double));
    fastqraddcolumn(&m, &n, &pk, workR, REAL(ansQ), REAL(uexp), REAL(ansR));
    
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

SEXP miclean() {
    free(workR);
    free(workQ);
    workR = NULL;
    workQ = NULL;
    //workRsize = -1;
    //workQsize = -1;
    SEXP retnull = PROTECT(allocVector(NILSXP, 1));
    UNPROTECT(1);
    return retnull;
}