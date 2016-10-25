//
// This is a set of simple utility routines and test
// generators for my matrix multiplication test bed.
//
#include "mm_utils.h"

//
// Compare two matrices ... return the sum of the squares
// of the differences of the two input matrices.
//
double errsqr(int Ndim, int Mdim, TYPE *C, TYPE *Cref) {
  int i, j;
  TYPE tmp, errsqr;
  errsqr = (TYPE)0.0;
  for (i = 0; i < Ndim; i++) {
    for (j = 0; j < Mdim; j++) {
      tmp = *(C + i * Mdim + j) - (*(Cref + i * Mdim + j));
      errsqr += tmp * tmp;
    }
  }
  return errsqr;
}

//
// Clear (i.e. set to zero) the elements of a matrix
//
void mm_clear(int Ndim, int Mdim, TYPE *C) {
  int i, j;
  for (i = 0; i < Ndim; i++)
    for (j = 0; j < Mdim; j++)
      *(C + i * Mdim + j) = (TYPE)0.0;
}

//
//  Print the elements of a matrix to standard out
//  (might be useful for debugging).
//
void mm_print(int Ndim, int Mdim, TYPE *C) {
  int i, j;
  for (i = 0; i < Ndim; i++) {
    for (j = 0; j < Mdim; j++)
      printf("[%04d][%04d] = %g   ", i, j, *(C + i * Mdim + j));
    printf("\n");
  }
}

//
//  Print error and timing results to standard out.
//
void output_results(int Ndim, int Mdim, int Pdim, int nerr, double ave_t,
                    double min_t, double max_t) {

  double dN, min_flop, max_flop, ave_flop;

  if (nerr > 0)
    printf(" %d errors\n", nerr);
  printf(" mult: ave=%f, min=%f, max=%f secs \n", ave_t, min_t, max_t);
  dN = 2.0 * (double)Ndim * (double)Mdim * (double)Pdim / (1000000.0);
  ave_flop = dN / ave_t;
  max_flop = dN / min_t;
  min_flop = dN / max_t;
  printf(" mult: ave=%f, min=%f, max=%f Mflops \n", ave_flop, min_flop,
         max_flop);
}

//=========================================================
//
// Test Matrices
//
// For each case, we have a functio to generate the test
// matrices (A and B) and the expected output matrix (C).
// This can be used to test correctness of matrix
// multiplications functions.
//
// The three matrices have dimensions ...
//           A(Ndim, Pdim), B(Pdim, Mdim), C(Ndim, Mdim)

//=========================================================
// Case one:  Constant matrices (A and B) to generat a constant
// matrix C.   This one is easy, but it is not a very strict
// test in that its easy to accidently write an erroneous
// multiplication algorithm that still passes this test
//

// Input and output matrices for constant matrices A and B
void init_const_matrix(int Ndim, int Mdim, int Pdim, TYPE *A, TYPE *B,
                       TYPE *C) {
  int i, j, k;
  TYPE Cval;

  for (i = 0; i < Ndim; i++)
    for (k = 0; k < Pdim; k++)
      *(A + i * Pdim + k) = AVAL;

  for (k = 0; k < Pdim; k++)
    for (j = 0; j < Mdim; j++)
      *(B + k * Mdim + j) = BVAL;

  Cval = (double)Pdim * (double)AVAL * (double)BVAL;
  for (i = 0; i < Ndim; i++)
    for (j = 0; j < Mdim; j++)
      *(C + i * Mdim + j) = Cval;
}

//=========================================================
// Case two: progression matrices.   The A  and B matrices
// generate finite series that when combined during the
// multiplication process produces a finite series with
// a mathematically well known, closed for answer.
//
// since the input and results matrices are not simple
// constants, it does a good job of catching errors in
// matrix multiply functions.
//
// Input matrices
//   A: elements of rows run 1 to Pdim (scaled by Aval)
//   B: elements of cols run from 1 to Pdim (scaled by Bval)
//   B: columns additionally scaled by col number (1 to Pdim)
//
void init_progression_matrix(int Ndim, int Mdim, int Pdim, TYPE *A, TYPE *B,
                             TYPE *C) {

  int i, j;
  TYPE Cval, Ctmp;

  for (i = 0; i < Ndim; i++) {
    for (j = 0; j < Pdim; j++)
      *(A + i * Pdim + j) = AVAL * (double)(j + 1);
  }

  for (i = 0; i < Pdim; i++) {
    for (j = 0; j < Mdim; j++)
      *(B + i * Mdim + j) = (j + 1) * BVAL * (double)(i + 1);
  }

  // I looked up sum of k squared for k=1 to P in
  // Gradshteyn and Ryzhik page 1.  I then scaled the C
  // matrix by the AVAL and BVAL factors and accounted for
  // the column scaling of B (thereby avoiding a constant
  // result matrix).

  Ctmp = (double)Pdim;
  Cval = Ctmp * (Ctmp + (double)1.0) * ((double)2.0 * Ctmp + (double)1.0);
  Cval = Cval * AVAL * BVAL / ((double)6.0);

  for (i = 0; i < Ndim; i++)
    for (j = 0; j < Mdim; j++)
      *(C + i * Mdim + j) = Cval * (j + 1);
}

//=========================================================
// Iteratiave solver test matrix generator
//=========================================================
void init_diag_dom_matrix(int Ndim, TYPE *A) {

  int i, j;
  TYPE sum;

  //
  // Create a random, diagonally dominant matrix.  For
  // a diagonally dominant matrix, the diagonal element
  // of each row is great than the sum of the other
  // elements in the row.
  for (i = 0; i < Ndim; i++) {
    sum = (TYPE)0.0;
    for (j = 0; j < Ndim; j++) {
      *(A + i * Ndim + j) = (rand() % 23) / (TYPE)100.0;
      sum += *(A + i * Ndim + j);
    }
    *(A + i * Ndim + i) += sum;
  }
}

//=========================================================
// Iteratiave solver test matrix generator.  This one is
// freindly to the Jacobi solver
//=========================================================
void init_colmaj_diag_dom_near_identity_matrix(int Ndim, TYPE *A) {
  int i, j;
  TYPE sum;

  //
  // Create a random, diagonally dominant matrix.  For
  // a diagonally dominant matrix, the diagonal element
  // of each row is great than the sum of the other
  // elements in the row.  Then scale the matrix so the
  // result is near the identiy matrix.
  for (i = 0; i < Ndim; i++) {
    sum = (TYPE)0.0;
    for (j = 0; j < Ndim; j++) {
      *(A + j * Ndim + i) = (rand() % 23) / (TYPE)1000.0;
      sum += *(A + j * Ndim + i);
    }
    *(A + i * Ndim + i) += sum;

    // scale the row so the final matrix is almost an identity matrix;wq
    for (j = 0; j < Ndim; j++)
      *(A + j * Ndim + i) /= sum;
  }
}
//===========================================================

void init_diag_dom_near_identity_matrix(int Ndim, TYPE *A) {

  int i, j;
  TYPE sum;

  //
  // Create a random, diagonally dominant matrix.  For
  // a diagonally dominant matrix, the diagonal element
  // of each row is great than the sum of the other
  // elements in the row.  Then scale the matrix so the
  // result is near the identiy matrix.
  for (i = 0; i < Ndim; i++) {
    sum = (TYPE)0.0;
    for (j = 0; j < Ndim; j++) {
      *(A + i * Ndim + j) = (rand() % 23) / (TYPE)1000.0;
      sum += *(A + i * Ndim + j);
    }
    *(A + i * Ndim + i) += sum;

    // scale the row so the final matrix is almost an identity matrix;wq
    for (j = 0; j < Ndim; j++)
      *(A + i * Ndim + j) /= sum;
  }
}
//===========================================================
