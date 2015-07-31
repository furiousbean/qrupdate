
// PACKAGE FILE
// QR add/remove of row/column functions
// see http://eprints.ma.man.ac.uk/1192/01/covered/MIMS_ep2008_111.pdf
// for the theory of QR factorization update

#include <R_ext/BLAS.h>
#include <R_ext/LAPACK.h>
#include "string.h"
#include "stdlib.h"
#include "fastqrupdate.h"

void fastqrdeleterow(int *pm, int *pn, int* pk, double *R, double *Q, 
double *newR, double *newQ) {
    int n = *pn;
    int m = *pm;
    int k = *pk - 1;
        
    double *cvec = (double *) calloc(m - 1, sizeof(double));
    double *svec = (double *) calloc(m - 1, sizeof(double));
    double *q = (double *) calloc(m, sizeof(double));
    
    for (int j = 0; j < m; j++) q[j] = Q[k + j * m];
    
    for (int i = m - 2; i >= 0; --i) {
        double a = q[i];
        double b = q[i + 1];
        F77_CALL(drotg)(&a, &b, cvec + i, svec + i);
        q[i] = cvec[i] * q[i] + svec[i] * q[i + 1];
        
        double tmp;
        for (int j = i; j < n; ++j) {
            tmp = R[i + j * m];
            R[i + j * m] = cvec[i] * R[i + j * m] + svec[i] * R[i + j * m + 1];
            R[i + j * m + 1] = -svec[i] * tmp + cvec[i] * R[i + j * m + 1];
        }
    }
    
    for (int i = 0; i < n; ++i) {
        memcpy(newR + i * (m - 1), R + 1 + i * m, (m - 1) * sizeof(double));
    }
    
    F77_CALL(dlasr)("R", "V", "B", &m, &m, cvec, svec, Q, &m);
    
    for (int i = 0; i < m - 1; ++i) {
        memcpy(newQ + i * (m - 1), Q + (i + 1) * m, 
        k * sizeof(double));
            
        memcpy(newQ + k + i * (m - 1), Q + k + 1 + (i + 1) * m, 
        (m - k - 1) * sizeof(double));
    }

    
    free(cvec);
    free(svec);
    free(q);
}

void fastqraddrow(int *pm, int *pn, int* pk, double *R, double *Q, double *u,
double *newR, double *newQ) {
    int n = *pn;
    int m = *pm;
    int k = *pk - 1;
        
    double *cvec = (double *) calloc(m, sizeof(double));
    double *svec = (double *) calloc(m, sizeof(double));
    
    int lim;
    if (m < n) lim = m; else lim = n;
    
    for (int i = 0; i < n; ++i) {
        memcpy(newR + i * (m + 1), R + i * m, m * sizeof(double));
        newR[m + i * (m + 1)] = u[i];
    }

    for (int i = 0; i < m; ++i) {
        memcpy(newQ + i * (m + 1), Q + i * m, k * sizeof(double));
        newQ[k + i * (m + 1)] = 0;
        memcpy(newQ + k + 1 + i * (m + 1), Q + k + i * m, 
        (m - k) * sizeof(double));
    }
    
    memset(newQ + m * (m + 1), 0, (m + 1) * sizeof(double));
    newQ[k + m * (m + 1)] = 1;
    
    for (int i = 0; i < lim; ++i) {
        double a = newR[i + i * (m + 1)];
        double b = newR[m + i * (m + 1)];
        F77_CALL(drotg)(&a, &b, cvec + i, svec + i);
        double tmp;
        for (int j = i; j < n; ++j) {
            tmp = newR[i + j * (m + 1)];
            newR[i + j * (m + 1)] = cvec[i] * newR[i + j * (m + 1)] + 
            svec[i] * newR[m + j * (m + 1)];
            
            newR[m + j * (m + 1)] = -svec[i] * tmp + 
            cvec[i] * newR[m + j * (m + 1)];
        }
    }
    
    for (int i = lim; i < m; ++i) {
        cvec[i] = 1.0;
    }
    
    int m1 = m + 1;
    
    F77_CALL(dlasr)("R", "B", "F", &m1, &m1, cvec, svec, newQ, &m1);
    
    free(cvec);
    free(svec);
}

void fastqrdeletecolumn(int *pm, int *pn, int* pk, double *R, double *Q, 
double *newR) {
    int n = *pn;
    int m = *pm;
    int k = *pk - 1;   

    memcpy(newR, R, (k * m) * sizeof(double));
    memcpy(newR + k * m, R + (k + 1) * m, ((n - k - 1) * m) * sizeof(double));
    
    int lim;
    if (n < m) lim = n - 1; else lim = m - 1;
    if (k >= lim) return;

    double *cvec = (double *) calloc(m - 1, sizeof(double));
    double *svec = (double *) calloc(m - 1, sizeof(double)); 

    for (int i = k; i < lim; ++i) {
        double a = newR[i * m + i];
        double b = newR[i * m + i + 1];
        F77_CALL(drotg)(&a, &b, cvec + i, svec + i);
        double tmp;
        for (int j = i; j < n - 1; ++j) {
            tmp = newR[i + j * m];
            newR[i + j * m] = cvec[i] * newR[i + j * m] + 
            svec[i] * newR[i + j * m + 1];
            
            newR[i + j * m + 1] = -svec[i] * tmp + 
            cvec[i] * newR[i + j * m + 1];
        }
    }
    
    for (int i = 0; i < k; ++i) cvec[i] = 1.0;
    for (int i = lim; i < m - 1; ++i) cvec[i] = 1.0;
    
    F77_CALL(dlasr)("R", "V", "F", &m, &m, cvec, svec, Q, &m);
    
    free(cvec);
    free(svec);
}

void fastqraddcolumn(int *pm, int *pn, int* pk, double *R, double *Q, double *u, 
double *newR) {
    int n = *pn;
    int m = *pm;
    int k = *pk - 1;   
    
    double one = 1.0;
    double zero = 0.0;
    int oneint = 1;
    
    double *qtu = (double *) malloc(m * sizeof(double));
    F77_CALL(dgemv)("T", &m, &m, &one, Q, &m, u, &oneint, &zero, qtu, &oneint);
    
    memcpy(newR, R, (k * m) * sizeof(double));
    memcpy(newR + k * m, qtu, m * sizeof(double));
    memcpy(newR + (k + 1) * m, R + k * m, ((n - k) * m) * sizeof(double));
    
    free(qtu);
    
    if (k >= m - 1) return;

    double *cvec = (double *) calloc(m - 1, sizeof(double));
    double *svec = (double *) calloc(m - 1, sizeof(double)); 

    for (int i = m - 2; i >= k; --i) {
        double a = newR[i + k * m];
        double b = newR[i + 1 + k * m];
        F77_CALL(drotg)(&a, &b, cvec + i, svec + i);
        double tmp;
        for (int j = k; j < n + 1; ++j) {
            tmp = newR[i + j * m];
            newR[i + j * m] = cvec[i] * newR[i + j * m] + 
                svec[i] * newR[i + j * m + 1];
            newR[i + j * m + 1] = -svec[i] * tmp + 
                cvec[i] * newR[i + j * m + 1];
        }
    }
    
    for (int i = 0; i < k; ++i) cvec[i] = 1.0;
    
    F77_CALL(dlasr)("R", "V", "B", &m, &m, cvec, svec, Q, &m);
    
    free(cvec);
    free(svec);
}

void fastqrsolve(int *pm, double *R, double *Q, double *b, double *x) {
    int m = *pm;
    
    double one = 1.0;
    double zero = 0.0;
    int oneint = 1;
    double *qtb = (double *) malloc(m * sizeof(double));
    F77_CALL(dgemv)("T", &m, &m, &one, Q, &m, b, &oneint, &zero, qtb, &oneint);
    F77_CALL(dtrsv)("U", "N", "N", &m, R, &m, qtb, &oneint);
    memcpy(x, qtb, m * sizeof(double));
    free(qtb);
}
