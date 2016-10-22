#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include <complex.h>


double LeJ(double sigma, double epsilon, double r);
double MORSE(double R, double Re, double alpha, double D);   

int main() 
{

  double r; 
  double sigma; 
  double epsilon;
  double LJ;
  double R; 
  double Re;
  double alpha;
  double D;
  int i;
  

  for (i=1.0;i<200;i++) { 

    r = i*0.075; 
    LJ = MORSE(r,2.419,0.90,2.419);
    //  LJ = LeJ(2.17921,2.41927,r); 

    printf("%lf, %lf \n", r,LJ); 
  
  }
  
  return(0);
}


/* The LJ equivalent for the Morce potential! */ 

double LeJ(double sigma, double epsilon, double r)

  {
    double LJo; 

    LJo = (3*sqrt(3)/2)*(epsilon*(pow((sigma/(r)),12)-pow((sigma/(r)),2)));
     

    return LJo;

  }

double MORSE(double R, double Re, double alpha, double D)  

{
  double E; 


  E = D*(exp(-2*alpha*(R-Re)) -2*exp(-alpha*(R-Re))); 
  //E = D*pow((1-exp(-alpha*(R-Re))),2) - D; 

  return E;

}
