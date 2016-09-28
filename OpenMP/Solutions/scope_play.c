/*

This is a simple program to play around with the scope 
of variables in OpenMP programs that include target
directives. 

History: Written by Tim Mattson, 11/2015. 

*/

#include <stdio.h>
#include <omp.h>

#pragma omp declare target
   void hello(){printf(" hello\n");}
   int ifs = 42;
#pragma omp end declare target

int main()
{
  int i = 1;
  int j = 2;
  int k = 3;

  int num_devices = omp_get_num_devices();
  printf(" There are %d devices.\n ",num_devices);
  if(num_devices == 0) printf(" ... if there was a device...\n");
  #pragma omp target map(to:j) map(tofrom:k) 
  {
     printf(" Inside target ifs is %d.  It should be 42\n",  ifs);
     printf(" Inside target i   is %d.  It should be undefined \n",  i);
     printf(" Inside target j   is %d.  It should be 2\n",  j);
     printf(" Inside target k   is %d.  It should be 3\n",   k);
     ifs++; i++; j++; k++;
  }
  printf(" Outside target ifs is %d.  It should be 42\n",  ifs);
  printf(" Outside target i   is %d.  It should be 1\n",  i);
  printf(" Outside target j   is %d.  It should be 2\n",  j);
  printf(" Outside target k   is %d.  It should be 4\n",   k);

}
