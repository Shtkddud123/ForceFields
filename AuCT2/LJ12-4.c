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
    LJ = LeJ(3.333,0.062305,r); 

    printf("%lf, %lf \n", r,LJ); 
  
  }
  
  return(0);
}



double LeJ(double sigma, double epsilon, double r)

  {
    double LJo; 

    LJo = ((3*sqrt(3)/2)*epsilon*(pow((sigma/r),12)-pow((sigma/r),4)));
     

    return LJo;

  }
