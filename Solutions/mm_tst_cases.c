/*
**  FUNCTION: Matrix Multiplication test cases
**
**  PURPOSE: This function runs through all the test cases
**           for a specified matrix multiplication function
**
**  HISTORY: Written by Tim Mattson, July 2012. 
*/

//#define DEBUG    1
#include "mm_utils.h"

void mm_tst_cases(int NTRIALS, int Ndim, int Mdim, int Pdim, 
              TYPE* A, TYPE* B, TYPE* C, 
              void (*mm_func)(int, int, int, TYPE *, TYPE *, TYPE *))
{
   int    nerr, itrials;
   double err,  errsq, mflops;
   double start_time, run_time;
   double min_t, max_t, ave_t;
   TYPE *Cref;

   Cref = (TYPE *) malloc (Ndim * Mdim * sizeof(TYPE));

   /* Initialize matrices */

   init_const_matrix (Ndim, Mdim, Pdim, A, B, Cref);

   printf("\n constant matrices  %d %d %d\n", Ndim, Mdim, Pdim);
   nerr = 0; min_t = BIG;  max_t = SMALL; ave_t = (double) 0.0;
   for (itrials = 0; itrials<NTRIALS; itrials++){

      mm_clear(Ndim, Mdim, C);
      start_time = omp_get_wtime(); 

      mm_func(Ndim, Mdim, Pdim, A, B, C);

      run_time = omp_get_wtime() - start_time;
  
      errsq = errsqr(Ndim, Mdim, C, Cref);
      if (errsq > TOL) nerr++;
      if(run_time < min_t) min_t = run_time;
      if(run_time > max_t) max_t = run_time;
      ave_t += run_time;
   }

   ave_t = ave_t/(double)NTRIALS;
   output_results(Ndim, Mdim, Pdim, nerr, ave_t, min_t, max_t);

   init_progression_matrix (Ndim, Mdim, Pdim, A, B, Cref);

#ifdef DEBUG
   printf(" A progression Matrix input\n");
   mm_print(Ndim, Pdim, A);

   printf(" B progression Matrix input\n");
   mm_print(Pdim, Mdim, B);

   printf(" C Reference Matrix\n");
   mm_print(Ndim, Mdim, Cref);
#endif

   printf("\n progression matrices  %d %d %d\n", Ndim, Mdim, Pdim);
   nerr = 0; min_t = BIG;  max_t = SMALL; ave_t = (double) 0.0;
   for (itrials = 0; itrials<NTRIALS; itrials++){

      mm_clear(Ndim, Mdim, C);
      start_time = omp_get_wtime(); 

      mm_func(Ndim, Mdim, Pdim, A, B, C);

      run_time = omp_get_wtime() - start_time;

#ifdef DEBUG
   printf(" C progression Matrix result\n");
   mm_print(Ndim, Mdim, C);
#endif
      errsq = errsqr(Ndim, Mdim, C, Cref);
      if (errsq > TOL) nerr++;
      if(run_time < min_t) min_t = run_time;
      if(run_time > max_t) max_t = run_time;
      ave_t += run_time;
   }

   ave_t = ave_t/(double)NTRIALS;
   output_results(Ndim, Mdim, Pdim, nerr, ave_t, min_t, max_t);
}
