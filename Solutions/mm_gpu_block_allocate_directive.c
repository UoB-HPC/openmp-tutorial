/*
**  function: Matrix Multiplication ... three loop, ikj case
**            where ijk defines the order of the loops
**
**  HISTORY: Written by Tom Deakin, January 2021.
*/
#include "mm_utils.h"

#include <assert.h>

#define Bsize 8

#ifdef USE_ALLOCATE

void mm_gpu_block_allocate(int Ndim, int Mdim, int Pdim, TYPE *A, TYPE *B, TYPE *C){
  //int i, j, k;
  //int ib, jb, kb;

  /* Block size, must evenly divide the matrices */
  assert(Ndim % Bsize == 0);
  assert(Mdim % Bsize == 0);
  assert(Pdim % Bsize == 0);

  /* Number of blocks in each dimension */
  int Nblk = Ndim / Bsize;
  int Mblk = Mdim / Bsize;
  int Pblk = Pdim / Bsize;

 #pragma omp target map(tofrom:C[0:Ndim*Mdim]) map(to:B[0:Pdim*Mdim], A[0:Ndim*Pdim]) uses_allocators(omp_pteam_mem_alloc)
 #pragma omp teams  num_teams(Nblk*Mblk) thread_limit(Bsize*Bsize)
 {

  /* Team local arrays */
  TYPE Awrk[Bsize*Bsize];
  TYPE Bwrk[Bsize*Bsize];
#pragma omp allocate(Awrk, Bwrk) allocator(omp_pteam_mem_alloc)

#pragma omp distribute collapse(2)
 for (int ib=0; ib < Nblk; ib++){ /* Loop over blocks of C. One team per block */
    for (int jb=0; jb < Mblk; jb++){


#pragma omp parallel num_threads(Bsize*Bsize)
{
      for (int kb=0; kb<Pblk; kb++){
        /* Copy block of A into pteam memory */
        #pragma omp for collapse(2) nowait
        for (int i=ib*Bsize; i<((ib+1)*Bsize); i++){
          for(int k=kb*Bsize; k<((kb+1)*Bsize); k++){
            Awrk[(i%Bsize)*Bsize + (k%Bsize)] = A[i*Pdim+k];
          }
        }

        /* Copy block of B into pteam memory */
        #pragma omp for collapse(2)
        for (int j=jb*Bsize; j<((jb+1)*Bsize); j++){
          for(int k=kb*Bsize; k<((kb+1)*Bsize); k++){
            Bwrk[(k%Bsize)*Bsize + (j%Bsize)] = B[k*Mdim+j];
          }
        }


        /* matrix multiply block */
        #pragma omp for collapse(2)
        for (int i=ib*Bsize; i<((ib+1)*Bsize); i++){
          for (int j=jb*Bsize; j<((jb+1)*Bsize); j++){
            /*for(k=kb*Bsize; k<((kb+1)*Bsize); k++){*/
            for (int k=0; k < Bsize; ++k){
              /* C(i,j) = sum(over k) A(i,k) * B(k,j) */
              *(C+(i*Mdim+j)) += Awrk[(i%Bsize)*Bsize + k] * Bwrk[k*Bsize + (j%Bsize)];
            }
          }
        }
}

      }
    }
  }
}
}

#else
void mm_gpu_block_allocate(int Ndim, int Mdim, int Pdim, TYPE *A, TYPE *B, TYPE *C){ }
#endif
