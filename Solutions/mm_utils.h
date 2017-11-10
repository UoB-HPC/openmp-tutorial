#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#define AVAL 3.0
#define BVAL 5.0
#define TYPE double
#define BIG 10000000.0
#define SMALL 0.00000001
#define TOL 0.001
//#define DEBUG 1

double errsqr(int Ndim, int Mdim, TYPE *C, TYPE *Cref);

void mm_clear(int Ndim, int Mdim, TYPE *C);

void mm_print(int Ndim, int Mdim, TYPE *C);

void init_const_matrix(int Ndim, int Mdim, int Pdim, TYPE *A, TYPE *B, TYPE *C);

void init_progression_matrix(int Ndim, int Mdim, int Pdim, TYPE *A, TYPE *B,
                             TYPE *C);

void output_results(int Ndim, int Mdim, int Pdim, int nerr, double ave_t,
                    double min_t, double max_t);

void mm_tst_cases(int NTRIALS, int Ndim, int Mdim, int Pdim, TYPE *A, TYPE *B,
                  TYPE *C,
                  void (*mm_func)(int, int, int, TYPE *, TYPE *, TYPE *));

void init_diag_dom_matrix(int Ndim, TYPE *A);

void init_diag_dom_near_identity_matrix(int Ndim, TYPE *A);
void init_diag_dom_near_identity_matrix_colmaj(int Ndim, TYPE *A);
