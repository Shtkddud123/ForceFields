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
  

  for (i=35;i<200;i++) { 

    r = i*0.075; 
    LJ = LeJ(4.3817,0.009709,r); 

    printf("%lf, %lf \n", r,LJ); 
  
  }
  
  return(0);
}



double LeJ(double sigma, double epsilon, double r)

  {
    double LJo; 

    LJo = (46656.0/3125.0)*epsilon*(pow((sigma/r),12)-pow((sigma/r),10));
     

    return LJo;

  }
