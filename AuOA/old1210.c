#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include <complex.h>


double LeJ(double sigma, double epsilon, double r,double pi, double goldden,double R);  

int main() 
{

  double r; 
  double sigma; 
  double epsilon;
  double LJ;
  double i;
  double pi = 4.0*atan(1.0);
  double goldden = 9.7984*pow(10.0,1.0);
  double R = 10.0;
  

  for (i=10.1;i<25.0;i+=0.01) { 

    r = i; 
    LJ = LeJ(4.3817,0.023539,r,pi,goldden,R); 

    printf("%lf, %lf \n", r,LJ); 
  
  }
  
  return(0);
}



double LeJ(double sigma, double epsilon, double r,double pi,double goldden,double R)

  {
    double HM; 

    HM = ((186624*pi*goldden*epsilon*pow(sigma,12)*pow(R,3)/140625)*((5*pow(R,6) + 45*pow(R,4)*pow(r,2) + 63*pow(R,2)*pow(r,4) + 15*pow(r,6))/pow(pow(r,2)-pow(R,2),9)))

      - ((186626*pi*goldden*epsilon*pow(sigma,10)*pow(R,3))/65625)*((5*pow(R,4) + 14*pow(R,2)*pow(r,2) + 5*pow(r,4))/pow(pow(r,2)-pow(R,2),7));
     

    return HM;

  }
