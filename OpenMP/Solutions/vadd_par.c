#include <stdio.h>
#include <omp.h>
#define N 100000
#define TOL  0.0000001
//
//  This is a simple program to add two vectors
//  and verify the results.
//
//  History: Written by Tim Mattson, November 2017
//
int main()
{

    float a[N], b[N], c[N], res[N],val;
    int i,err=0;

   // fill the arrays
   #pragma omp parallel for
   for (i=0; i<N; i++){
      a[i] = (float)i;
      b[i] = 2.0*(float)i;
      c[i] = 0.0;
      res[i] = i + 2*i;
   }

   // add two vectors
   #pragma omp target
   #pragma omp teams distribute parallel for simd
   for (i=0; i<N; i++){
      c[i] = a[i] + b[i];
   }

   // test results
   #pragma omp parallel for private(val) reduction(+:err)
   for(i=0;i<N;i++){
      val = c[i] - res[i];
      val = val*val;
      if(val>TOL) err++;
   }
   printf(" vectors added with %d errors\n",err);
   return 0;
}
