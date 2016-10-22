#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include <complex.h>


double LeJ(double sigma, double epsilon, double r);  

int main() 
{

  double r; 
  double sigma; 
  double epsilon;
  double LJ;
  int i;
  

  for (i=0;i<200;i++) { 

    r = i*0.075; 
    LJ = LeJ(1.96555,2.45963,r); 

    printf("%lf, %lf \n", r,LJ); 
  
  }
  
  return(0);
}



double LeJ(double sigma, double epsilon, double r)

  {
    double LJo; 

    LJo = (27/4)*epsilon*(pow((sigma/r),9)-pow((sigma/r),6));
     

    return LJo;

  }
