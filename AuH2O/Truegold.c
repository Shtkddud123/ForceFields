#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include <complex.h>
#include <string.h>

double HM(double alpha, double a, double D, double r, double pi, double goldden,double b);
double interaction(double sigma, double epsilon, double r, double pi, double goldden,double R);
double MORSE(double R, double Re, double alpha, double D);

int main(int argc, char *argv[]) 
{

  double r; 
  double sigma; 
  double epsilon;
  double LJ;
  double hamaker; 
  double i;
  double pi = 4.0*atan(1.0);
  double goldden = 0.059;
  double a = 20.0;
  double inter; 
  double gradient; 
  int j = 0;
  double alpha = 0.90; 
  double D = 2.419; 
  double b = 2.419;
 
  
  
 
  // epsilon=atof(argv[1]);
   // sigma=atof(argv[2]);
  

  for (i=20.10;i<=35.0;i+=0.01) { 
    
    r = i;
    //LJ = LeJ(3.2529,0.666381,r,pi,goldden,R);
    //inter = interaction(3.2529,0.666381,r,pi,goldden,R);
    hamaker=HM(alpha,a,D,r,pi,goldden,b);
    // inter=interaction(sigma,epsilon,r,pi,goldden,R);
    // gradient = -(LeJ(sigma,epsilon,r+0.02,pi,goldden,R) - LeJ(sigma,epsilon,r,pi,goldden,R))/0.02;
    j++; 

    printf("%d %lf %lf \n",j,r,-hamaker);
  
  }
  
  return(0);
}



double HM(double alpha, double a, double D,double r, double pi, double goldden, double b)

  {
    double HM; 


    HM = (1/24.0)*1/(pow(alpha,4)*r)*
      (
       (pi*goldden*D)*(288.0*exp(alpha*(a-r+b))*alpha*a /*ok2*/
		       -96.0*exp(alpha*(a-r+b))*pow(alpha,2.0)*pow(a,2.0) /*ok2*/ 
		       - 18.0*exp(2.0*alpha*(a-r+b))*alpha*a /*ok2*/ 
		       +12.0*exp(2.0*alpha*(a-r+b))*pow(alpha,2.0)*pow(a,2.0) /*ok*/
		       + 288.0*exp(-alpha*(a+r-b))*alpha*a /*ok2*/
		       + 96.0*exp(-alpha*(a+r-b))*pow(alpha,2.0)*pow(a,2.0) /*ok2*/
		       - 18.0*exp(-2.0*alpha*(a+r-b))*alpha*a /*ok2*/
		       -12.0*exp(-2.0*alpha*(a+r-b))*pow(alpha,2.0)*pow(a,2.0) /*ok*/
		       -96.0*exp(alpha*(a-r+b))*r*alpha /*ok2*/
		       + 96.0*exp(alpha*(a-r+b))*r*pow(alpha,2.0)*a /*ok2*/
		       + 6.0*exp(2.0*alpha*(a-r+b))*alpha*r /*ok2*/
		       - 12.0*exp(2.0*alpha*(a-r+b))*r*pow(alpha,2.0)*a /*ok2*/
		       + 96.0*exp(-alpha*(a+r-b))*r*alpha /*ok3*/
		       + 96.0*exp(-alpha*(a+r-b))*r*pow(alpha,2.0)*a /*ok2*/
		       - 6.0*exp(-2.0*alpha*(a+r-b))*r*alpha /*ok2*/ 
		       - 12.0*exp(-2.0*alpha*(a+r-b))*r*pow(alpha,2.0)*a /*ok2*/
		       - 32.0*pow(a,3.0)*pow(alpha,4.0)*r /*ok2*/
		       - 288.0*exp(alpha*(a-r+b)) /*ok2*/ 
		       + 9.0*exp(2*alpha*(a-r+b)) /*ok2*/
		       + 288.0*exp(-alpha*(a+r-b)) /*ok2*/
		       - 9.0*exp(-2*alpha*(a+r-b)))); /*ok2*/


    return HM;

  }


double interaction(double sigma, double epsilon, double r, double pi, double goldden,double R)

{ 

  double interaction; 

  interaction = (108*pi*epsilon*goldden*pow(sigma,12)*pow(R,3)/5)*(5*pow(R,6) + 27*pow(R,4)*pow(r,2) + 27*pow(R,2)*pow(r,4) + 5*pow(r,6))*r/pow(pow(r,2) - pow(R,2),10)
              - 72*pi*goldden*epsilon*pow(sigma,8)*pow(R,3)*(pow(r,2) + pow(R,2))*r/pow(pow(r,2) - pow(R,2),6); 


  return interaction;

}  

double MORSE(double R, double Re, double alpha, double D)  

{
  double E; 


  E = D*(exp(-2*alpha*(R-Re)) -2*exp(-alpha*(R-Re))); 
  //E = D*pow((1-exp(-alpha*(R-Re))),2) - D; 

  return E;

}
