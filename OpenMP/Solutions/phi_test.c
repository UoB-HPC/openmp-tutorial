/*    Copyright (c) 2014 Alin Marin Elena <alinm.elena@gmail.com>
      The MIT License http://opensource.org/licenses/MIT*/
#include <stdio.h>
#include <omp.h>

#pragma omp declare target  
  void testThreads();
#pragma omp end declare target  

int main(int argc,char **argv)
{  
  int nCards; 
  nCards = omp_get_num_devices();
  printf("No of cards: %d\n",nCards);
  #pragma omp target  
    testThreads();
  testThreads();
  return 0;
}

void testThreads()
{
  int tCount,threadsMax,tid;   
#ifdef __MIC__
  printf("Running on Xeon Phi\n");
#else
  printf("Running on host\n");
#endif //__MIC__
  threadsMax=omp_get_max_threads();
  #pragma omp parallel default(none) private(tid) shared(tCount,threadsMax)
  {
    tid = omp_get_thread_num();
    #pragma omp critical
      printf("Hello from thread id %d\n",tid);
    #pragma omp barrier 
    #pragma omp master
    {
      tCount = omp_get_num_threads();
      printf("The number of threads is %d\n", tCount);
      printf("Maximum number of threads is %d\n", threadsMax);
      printf("================================\n");
    }
  }
}

