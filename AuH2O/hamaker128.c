#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include <complex.h>


double LeJ(double sigma, double epsilon, double r, double pi, double goldden,double R);
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
  
  

  for (i=20.1;i<=35.0;i+=0.01) { 

    r = i; 
    LJ = LeJ(5.6475,0.430064,r,pi,goldden,R); 
    inter = interaction(3.2529,0.666381,r,pi,goldden,R);
    j++;

    printf("%d %lf %lf %lf \n",j,r,LJ,inter); 
  
  }
  
  return(0);
}



double LeJ(double sigma, double epsilon, double r, double pi, double goldden,double R)

  {
    double HM; 
    
    
     HM = (9*pi*goldden*epsilon*pow(sigma,12)*pow(R,3)/5)*(5*pow(R,6) + 15*pow(R,4)*pow(r,2) + 21*pow(R,2)*pow(r,4) + 5*pow(r,6))/pow(pow(r,2)-pow(R,2),9)
      - (9*pi*goldden*epsilon*pow(sigma,8)/5)*(pow(R,3)*(3*pow(R,2) + 5*pow(r,2)))/pow(pow(r,2) - pow(R,2),5);
     

    return HM;

  }

double interaction(double sigma, double epsilon, double r, double pi, double goldden,double R)

{ 
  
  double interaction; 

 interaction = (108*pi*epsilon*goldden*pow(sigma,12)*pow(R,3)/5)*(5*pow(R,6) + 27*pow(R,4)*pow(r,2) + 27*pow(R,2)*pow(r,4) + 5*pow(r,6))*r/pow(pow(r,2) - pow(R,2),10) 
    - 72*pi*goldden*pow(sigma,8)*pow(R,3)*(pow(r,2) + pow(R,2))*r/pow(pow(r,2) - pow(R,2),6); 



  return interaction;
  
}  
