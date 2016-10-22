#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include <complex.h> 


int main() 
{ 
  /*                     b (CH2)      */
  /*                     |            */
  /*   a (Au)            |            */
  /*                     |            */
  /*                     c (CH3)      */
  /*--------------------------*/
  /* This is for the CT2 bead!*/
  /*--------------------------*/

  /*Declaration of variables*/
  /*Declare variables for the doublet*/ 
  double bx,by,bz;    /*Cartesian coordinates for atom b of triad*/
  double bepsilon,bsigma; /*LJ parameters for atom b*/ 
  double cx,cy,cz; /*Cartesian coordinates for atom c of triad*/
  double cepsilon,csigma; /*LJ parameters for atom c*/
  double bondlength; /*bondlength bc and bd */
  
  /*We are now declaring the variables for gold*/ 
  double ax,ay,az; /*Cartesian coordinates for gold*/
  double aepsilon,asigma; /*LJ parameters for gold*/ 
  
  /*We are now declaring the variables for the LJ calculation*/ 
  double abdist; /*Distance between a and b*/ 
  double acdist; /*Distance between a and c*/ 
  double LJ_ab; /*LJ potential between a and b */ 
  double LJ_ac; /*LJ potential between a and c */ 
  double absigma;
  double acsigma;
  double abepsilon; 
  double acepsilon;  
  double LJtotal; 

  /*These are the parameters for the centre of mass*/
  double bmass;
  double cmass; 
  double xcentre;
  double ycentre;
  double zcentre;
  double xmass;
  double ymass;
  double zmass; 

  /*Polar coordinates and weighting from leb_weight.dat*/ 
  double theta; 
  double phi;
  double weight; /*Lebetev weights divided by 4pi*/
  double unknown; /*Need to fix this ASAP*/ 

  /*These are the new gold cartesian coordinates given by the polar coordinates, and the LJ energy from the new position to the COM of the triad*/ 
  double newgoldx; 
  double newgoldy; 
  double newgoldz; 
  double goldLJ; 
  double COMdist;
  double r; /*Arbitrary value of r for our example*/ 
  double V_ave = 0; /*The average interaction energy*/ 
  double longestlength = 16; 
        
  /*Parameters for loops and ipf*/ 
  const int coordinates = 325;
  char line[100]; 
  int index = 0; 
  int newindex=0; 
  
  /*Parameters for the new x,y,z axis for the Au particle*/ 
  /*------------------------------------------------------*/
  double axx[coordinates]; 
  double ayy[coordinates]; 
  double azz[coordinates];
  /*------------------------------------------------------*/
  double newtheta;  /*theta angle input from leb_weight.dat*/   
  double newphi; /*phi angle input from leb_weight.dat*/
  double newweight; /*weights input from leb_weight.dat*/  
  double newunknown; /*more weights from leb_weight.dat!*/
  int i;       /*integers used in FOR loops*/
  int n = 0;   /*ditto*/ 
  int j = 0;   /*ditto*/ 

 
 
  /*More test parameters for calculating the average energy from a range of values for r*/
  
  double energyarray[coordinates];
 
  /*Defining the triad positions*/ 
  bondlength = 1.53; /*length between CH2 and CH3*/ 
  bx = 0.0;
  by = 0.765; /*half of bondlength*/ 
  bz = 3.581; /*calculated using tan function*/
  
  /*Given b, we can work out the position of c */ 
  
  cx = 0.0;
  cy = -0.765;
  cz = 3.581; 
  /* printf("%lf %lf %lf \n", cx, cy, cz);  */
  
  /* Defining the triad LJ parameters*/ 
  /* Units should be changed to angstroms and kJ mol^-1 for sigma and epsilon respectively*/ 
  
  bepsilon = 0.1004; /*These are parameters for CH2*/ 
  bsigma = 4.0500; /*ditto*/ 
  cepsilon = 0.1960; /*these are parameters for CH3*/
  csigma = 3.75; /*ditto*/ 
  
  /*Define the gold positions*/ 
  /*COM of the triad*/ 
  
  theta = 1.570796;
  phi = 4.712389; 
  ax = 0.0; 
  ay = 0.0;
  az = 0.0; 
  
  /*Define the gold LJ parameters*/ 
  /*asigma is in angstoms*/ 
  
  asigma  = 2.629;
  aepsilon = 5.29; /*Units are in eV, should change to kJ mol^-1 ASAP*/ 

  /*LJ calculation*/ 
  /*LJ interaction energy between a and b*/ 
  /* We need to check whether the cross terms are calculated correctly (can depend on force field!)*/
  
  absigma = sqrt(asigma*bsigma); 
  abepsilon = sqrt(aepsilon*bepsilon); 
  abdist = sqrt((bx-ax)*(bx-ax) + (by-ay)*(by-ay) + (bz-az)*(bz-az)); 
  LJ_ab = 4*abepsilon*(pow((absigma/abdist),12)-pow((absigma/abdist),6)); 
  /* printf("%lf \n", LJ_ab); */
  
  /*LJ interaction energy between a and c*/ 
  /* We need to check whether the cross terms are calculated correctly (can depend on force field!)*/

  acsigma = sqrt(asigma*csigma); 
  acepsilon = sqrt(aepsilon*cepsilon); 
  acdist = sqrt((cx-ax)*(cx-ax) + (cy-ay)*(cy-ay) + (cz-az)*(cz-az)); 
  LJ_ac = 4*acepsilon*(pow((acsigma/acdist),12)-pow((acsigma/acdist),6)); 
  /*printf("%lf \n", LJ_ac);*/


  /*Sum of LJ energies*/ 

  LJtotal = LJ_ab + LJ_ac; 
  /*printf("%lf \n", LJtotal);*/

  /*We are now calculating the COM of the triad*/

  bmass = 14.0;   /*We are assuming here that the b atom is CH2*/ 
  cmass = 15.0;   /* We are assuming here that the c atom is CH3*/ 
   

  xmass = (bx*bmass + cx*cmass)/(bmass + cmass);
  ymass = (by*bmass + cy*cmass)/(bmass + cmass);
  zmass = (bz*bmass + cz*cmass)/(bmass + cmass);
  
  /* printf("The cartesian coordinates for the centre of mass (COM) is %lf %lf %lf \n",xmass,ymass,zmass); */

  /*We are now reading in the polar coordinates*/
  /*There are 302 datasets*/ 
  
  FILE *ipf; /*input file*/

  /*open file for reading*/
  /* printf("Opening file %s\n","leb_weights.dat");*/
  ipf = fopen("leb_weights.dat", "r");

  /* check for error opening file */
  if(ipf == NULL) {
    
    printf("Error opening file\n");
    exit(1);
  }
    
  /*get a line from the file */
  /*fgets returns NULL when it reaches an error or end of file */
  while(fgets(line,sizeof(line),ipf) !=NULL) {
       
      sscanf(line,"%lf %lf %lf",&newtheta,&newphi,&newweight);
      
      axx[index] = newtheta; 
      ayy[index] = newphi; 
      azz[index] = newweight; 
      n++; 
      index++;
     
  }
  fclose(ipf);
  
  for (i=35;i<200;i++) { 
    
    r = i*0.075; /*Sets the distance between the COM and the Au particle */ 
    V_ave = 0; 
    
    for (j=0;j<n;j++) {
      
      newgoldx = r*sin(axx[j])*cos(ayy[j]);
      newgoldy = r*sin(axx[j])*sin(ayy[j]); 
      newgoldz = r*cos(axx[j]);
      COMdist = sqrt((xmass-newgoldx)*(xmass-newgoldx) + (ymass-newgoldy)*(ymass-newgoldy) + (zmass-newgoldz)*(zmass-newgoldz));
      goldLJ = 4*abepsilon*(pow((absigma/COMdist),12)-pow((absigma/COMdist),6)); /*We are assuming for now sigma for the COM is identical to that of ab,ac,and ad*/ 
      
      /*We are now moving the gold particle around according to ONE of the given coordinates of theta and phi, then calculating the LJ energy from it!*/
      
      V_ave += goldLJ*azz[j];  /*Now we calculate the average over all the polar coordinates using the weighted values*/
      
    }

    printf("%lf %lf \n",r,V_ave);
    
  }
    
  return(0);
}


