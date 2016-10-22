#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include <complex.h>


double LeJ(double sigma, double epsilon, double r,double pi, double goldden,double R);  
double interaction(double sigma, double epsilon, double r, double pi, double goldden,double R);

int main() 
{

  double r; 
  double sigma; 
  double epsilon;
  double LJ;
  double i;
  double pi = 4.0*atan(1.0);
  double goldden = 0.059; 
  double R = 20.0;
  double inter; 
  int j = 0;
  

  for (i=20.10;i<=35.0;i+=0.01) { 

    r = i; 
    LJ = LeJ(4.3817,0.009709,r,pi,goldden,R); 
    inter = interaction(4.3817,0.009709,r,pi,goldden,R);
    j++;

    printf("%d %lf %lf %lf \n",j,r,LJ,inter); 
  
  }
  
  return(0);
}



double LeJ(double sigma, double epsilon, double r,double pi, double goldden,double R)

  {
    double HM; 

  HM = (20736*pi*goldden*epsilon*pow(sigma,12)*pow(R,3)/15625)*(5*pow(R,6) + 45*pow(R,4)*pow(r,2) + 63*pow(R,2)*pow(r,4) + 15*pow(r,6))/pow(pow(r,2) - pow(R,2),9) 
      - (62208*pi*goldden*epsilon*pow(sigma,10)*pow(R,3)/21875)*(3*pow(R,4) + 14*pow(R,2)*pow(r,2) + 7*pow(r,4))/pow(pow(r,2) - pow(R,2),7); 
     

    return HM;

  }

double interaction(double sigma, double epsilon, double r, double pi, double goldden,double R)

{

 double interaction; 
  

 interaction = (746496*pi*goldden*epsilon*pow(sigma,12)*pow(R,3)/15625)*(5*pow(R,6) + 27*pow(R,4)*pow(r,2) + 27*pow(R,2)*pow(r,4) + 5*pow(r,6))*r/pow(pow(r,2) - pow(R,2),10) 
    - (124416*pi*goldden*epsilon*pow(sigma,10)*pow(R,3)/3125)*(5*pow(R,4) + 14*pow(R,2)*pow(r,2) + 5*pow(r,4))*r/pow(pow(r,2) - pow(R,2),8); 



  return interaction;

}  
