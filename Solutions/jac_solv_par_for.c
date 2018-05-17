/*
**  PROGRAM: jacobi Solver ... parallel region containing for constructs
**
**  PURPOSE: This program will explore use of a jacobi iterative
**           method to solve a system of linear equations (Ax= b).
**
**           Here is the basic idea behind the method.   Rewrite
**           the matrix A as a Lower Triangular (L), upper triangular
**           (U) and diagonal matrix (D)
**
**                Ax = (L + D + U)x = b
**
**            Carry out the multiplication and rearrange:
**
**                Dx = b - (L+U)x  -->   x = (b-(L+U)x)/D
**
**           We can do this iteratively
**
**                x_new = (b-(L+U)x_old)/D
**
**  USAGE:   Run wtihout arguments to use default SIZE.
**
**              ./jac_solv
**
**           Run with a single argument for the order of the A
**           matrix ... for example
**
**              ./jac_solv 2500
**
**  HISTORY: Written by Tim Mattson, Oct 2015
*/

#include "mm_utils.h" //a library of basic matrix utilities functions
#include <math.h>
#include <omp.h>
// and some key constants used in this program
//(such as TYPE)

#define TOLERANCE 0.001
#define DEF_SIZE 1000
#define MAX_ITERS 5000
#define LARGE 1000000.0

//#define DEBUG    1     // output a small subset of intermediate values
//#define VERBOSE  1

int main(int argc, char **argv) {
  int Ndim; // A[Ndim][Ndim]
  int i, j, iters;
  double start_time, elapsed_time;
  TYPE conv, tmp, err, chksum;
  TYPE *A, *b, *x1, *x2, *xtmp;

  // set matrix dimensions and allocate memory for matrices
  if (argc == 2) {
    Ndim = atoi(argv[1]);
  } else {
    Ndim = DEF_SIZE;
  }

  printf(" \n\n jacobi solver parallel (parallel + for version): ndim = %d\n",
         Ndim);

  A = (TYPE *)malloc(Ndim * Ndim * sizeof(TYPE));
  b = (TYPE *)malloc(Ndim * sizeof(TYPE));
  x1 = (TYPE *)malloc(Ndim * sizeof(TYPE));
  x2 = (TYPE *)malloc(Ndim * sizeof(TYPE));

  if (!A || !b || !x1 || !x2) {
    printf("\n memory allocation error\n");
    exit(-1);
  }

  // generate our diagonally dominant matrix, A
  init_diag_dom_near_identity_matrix(Ndim, A);

#ifdef VERBOSE
  mm_print(Ndim, Ndim, A);
#endif

  //
  // Initialize x and just give b some non-zero random values
  //
  for (i = 0; i < Ndim; i++) {
    x1[i] = (TYPE)0.0;
    x2[i] = (TYPE)0.0;
    b[i] = (TYPE)(rand() % 51) / 100.0;
  }

  start_time = omp_get_wtime();
  //
  // jacobi iterative solver
  //
  conv = LARGE;
  iters = 0;

#pragma omp parallel default(none) private(tmp) shared(Ndim, conv, iters, b,   \
                                                       A, x2, x1, xtmp)
  {
    // note: i am comparing against the convergence sqaured.  This saves a
    // sqrt and an extra barrier.
    while ((conv > TOLERANCE * TOLERANCE) && (iters < MAX_ITERS)) {

#ifdef DEBUG
      printf("thread %d, iters=%d conv=%f\n", omp_get_thread_num(), iters,
             (float)conv);
#endif
#pragma omp for private(i, j) nowait
      for (i = 0; i < Ndim; i++) {
        x2[i] = (TYPE)0.0;
        for (j = 0; j < Ndim; j++) {
          //    if(i!=j)
          //      x2[i]+= A[i*Ndim + j]*x1[j];
          x2[i] += A[i * Ndim + j] * x1[j] * (i != j);
        }
        x2[i] = (b[i] - x2[i]) / A[i * Ndim + i];
      }
#pragma omp single
      {
        iters++;
        conv = 0.0;
      }
//
// test convergence
//
#pragma omp for private(tmp) reduction(+ : conv)
      for (i = 0; i < Ndim; i++) {
        tmp = x2[i] - x1[i];
        conv += tmp * tmp;
      }
#ifdef DEBUG
      printf(" conv = %f \n", (float)conv);
#endif

      TYPE* tmp = x1;
      x1 = x2;
      x2 = tmp;
    }
  }
  conv = sqrt((double)conv);
  elapsed_time = omp_get_wtime() - start_time;
  printf(" Convergence = %g with %d iterations and %f seconds\n", (float)conv,
      iters, (float)elapsed_time);

  //
  // test answer by multiplying my computed value of x by
  // the input A matrix and comparing the result with the
  // input b vector.
  //
  err = (TYPE)0.0;
  chksum = (TYPE)0.0;

  for (i = 0; i < Ndim; i++) {
    x1[i] = (TYPE)0.0;
    for (j = 0; j < Ndim; j++)
      x1[i] += A[i * Ndim + j] * x2[j];
    tmp = x1[i] - b[i];
#ifdef DEBUG
    printf(" i=%d, diff = %f,  computed b = %f, input b= %f \n", i, (float)tmp,
        (float)x1[i], (float)b[i]);
#endif
    chksum += x2[i];
    err += tmp * tmp;
  }
  err = sqrt((double)err);
  printf("jacobi solver: err = %f, solution checksum = %f \n", (float)err,
      (float)chksum);
  if (err > TOLERANCE)
    printf("\nWARNING: final solution error > %g\n\n", TOLERANCE);

  free(A);
  free(b);
  free(x1);
  free(x2);
}
