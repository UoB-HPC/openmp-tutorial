/*

This is a simple program to play around with the scope 
of variables in OpenMP programs that include target
directives. 

The goal is to fill in the print statements with the value
you expect based on your understanding of the scoping rules
for variables with regard to the OpenMP target directive.

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
  if(num_devices == 0) 
         printf(" ... if there was a device, the following would hold...\n");
  #pragma omp target map(to:j,i) map(tofrom:k) 
  {
     printf(" Inside target ifs is %d.  It should be <<fill_in>>\n",  ifs);
     printf(" Inside target i   is %d.  It should be <<fill_in>>\n",  i);
     printf(" Inside target j   is %d.  It should be <<fill_in>>\n",  j);
     printf(" Inside target k   is %d.  It should be <<fill_in>>\n",   k);
     ifs++; i++; j++; k++;
  }
  printf(" Outside target ifs is %d.  It should be <<fill_in>>\n",  ifs);
  printf(" Outside target i   is %d.  It should be <<fill_in>>\n",  i);
  printf(" Outside target j   is %d.  It should be <<fill_in>>\n",  j);
  printf(" Outside target k   is %d.  It should be <<fill_in>>\n",   k);

}
