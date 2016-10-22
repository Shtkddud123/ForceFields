/*Edit:The interaction potentials has been checked and it is correct!*/

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include <complex.h>
#include <string.h>

double LeJ124(double sigma, double epsilon, double r, double pi, double waterden,double R);
double LeJ96(double sigma, double epsilon, double r, double pi, double waterden,double R);
double interaction124(double sigma, double epsilon, double r, double pi, double CT2den,double R);
double interaction96(double sigma, double epsilon, double r, double pi, double CT2den,double R);

int main(int argc, char *argv[]) 
{

  double r; 
  double sigma; 
  double epsilon;
  double LJ;
  double i;
  double pi = 4.0*atan(1.0);
  double waterden = 0.0334; 
  double CT2den = 0.011348; 
  double R = 10.0;
  double inter;
  double gradient; 
  double LJ1; 
  double LJ2;
  double newLJ;
  int j = 0;
  double inter1; 
  double inter2;
  double newinter;
  double k = 1.0; /*k defines the proportion of hydrophil/hydrophobicity.*/
  

 /*defining the logarithmic part*/ 
 
 
  //epsilon=atof(argv[1]);
  //sigma=atof(argv[2]);
  

  for (i=10.10;i<=25.0;i+=0.01) { 
    
    r = i;
    /*The component for the hydrophilic nanoparticle*/
    LJ1=LeJ124(3.9500,0.7000,r,pi,waterden,R);
    /*The component for the hydrophobic nanoparticle*/
    LJ2=LeJ96(3.8400,0.3800,r,pi,CT2den,R);
    newLJ = (1.0-k)*LJ1 + k*(LJ2);
    /*The component for the hydrophilic nanoparticle*/
    inter1=interaction124(3.95000,0.7000,r,pi,waterden,R);
    /*The component for the hydrophobic nanoparticle*/ 
    inter2=interaction96(3.8400,0.3800,r,pi,CT2den,R);
    newinter = (1-k)*inter1 + k*(inter2); 

    //gradient = -(LeJ(sigma,epsilon,r+0.02,pi,goldden,R) - LeJ(sigma,epsilon,r,pi,goldden,R))/0.02;
    j++; 

    /*defining the logarithmic part*/ 


    printf("%d %lf %lf %lf  \n",j,r,newLJ+0.002689,newinter); 
  
  }
  
  return(0);
}



double LeJ124(double sigma, double epsilon, double r,double pi, double waterden, double R)

  {
    double HM; 
    double fig = (R-r)/(R+r);
    double number = log(fig);
    
    HM = (2*sqrt(3)*pi*waterden*epsilon*pow(sigma,12)/15)*(pow(R,3)*(5*pow(R,6) + 45*pow(R,4)*pow(r,2) + 63*pow(R,2)*pow(r,4) + 15*pow(r,6)))/pow(pow(r,2)-pow(R,2),9)-(3*sqrt(3)*pi*waterden*epsilon*pow(sigma,4))*((R/(pow(r,2)-pow(R,2))) + (1/(2*r))*log(-fig)); 

    return HM;

  }


double interaction124(double sigma, double epsilon, double r, double pi, double waterden,double R)

{ 

  /*This still needs editing*/
  double interaction;
  double fig = (R-r)/(R+r);
  double number = log(fig);
  

  interaction = (((72*sqrt(3))*pi*epsilon*waterden*pow(sigma,12))/15)*(((pow(R,3)*r)*(5*pow(R,6) + 27*pow(R,4)*pow(r,2) + 27*pow(R,2)*pow(r,4) + 5*pow(r,6)))/pow(pow(r,2) - pow(R,2),10)) - ((3*sqrt(3))*pi*waterden*epsilon*pow(sigma,4)/2)*((2*R*(pow(R,2)+pow(r,2)))/(pow(pow(r,2)-pow(R,2),2)*r) + (1/(pow(r,2)))*log(-fig)); 


  return interaction;

} 

 double LeJ96(double sigma, double epsilon, double r,double pi, double CT2den, double R)

  {
    double HM; 
 
    
    HM = (9*pi*CT2den*epsilon*pow(sigma,9)/35)*((pow(R,3)*(3*pow(R,4) + 42*pow(R,2)*pow(r,2) + 35*pow(r,4)))/(r*pow(pow(r,2)-pow(R,2),6))) - (9*pi*CT2den*epsilon*pow(sigma,6)*pow(R,3))/(pow(pow(r,2)-pow(R,2),3)); 

    return HM;

  }

double interaction96(double sigma, double epsilon, double r, double pi, double CT2den,double R)

{ 

  /*This still needs editing*/
  double interaction;
  

  interaction = ((27*pi*epsilon*CT2den*pow(sigma,9))/35)*((pow(R,3)*(-pow(R,6) + 27*pow(R,4)*pow(r,2) + 189*pow(R,2)*pow(r,4) + 105*pow(r,6)))/(pow(r,2)*pow(pow(r,2) - pow(R,2),7))) - (54*pi*CT2den*epsilon*pow(sigma,6))*((pow(R,3)*r)/pow(pow(r,2)-pow(R,2),4)); 


  return interaction;

}
