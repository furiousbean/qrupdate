#ifndef FASTQRUPDATE_H
#define FASTQRUPDATE_H

void fastqrdeleterow(int *pm, int *pn, int* pk, double *R, double *Q, 
    double *newR, double *newQ);
    
void fastqraddrow(int *pm, int *pn, int* pk, double *R, double *Q, double *u,
    double *newR, double *newQ);
    
void fastqrdeletecolumn(int *pm, int *pn, int* pk, double *R, double *Q, 
    double *newR);

void fastqraddcolumn(int *pm, int *pn, int* pk, double *R, double *Q, double *u, 
    double *newR);
    
void fastqrsolve(int *pm, double *R, double *Q, double *b, double *x);

#endif //FASTQRUPDATE_H
