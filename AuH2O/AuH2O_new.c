#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include <complex.h> 


int main() 
{ 
  /*                     d(H)     */
  /*                    /            */
  /*   a (Au)          b (O)       */
  /*                    \            */
  /*                     c(H)      */  
  /*---------------------------------*/
  /*This code is for the Au-CM bead  */
  /*---------------------------------*/

  /*Declaration of variables*/
  /*Declare variables for the triad*/ 
  double bx,by,bz;    /*Cartesian coordinates for atom b of triad*/
  double bepsilon,bsigma; /*LJ parameters for atom b*/ 
  double cx,cy,cz; /*Cartesian coordinates for atom c of triad*/
  double cepsilon,csigma; /*LJ parameters for atom c*/
  double dx,dy,dz; /*Cartesian coordinates for atom d of triad*/ 
  double depsilon,dsigma; /*LJ parameters for atom d*/
  double bondlength; /*bondlength bc and bd */
  double bondangle; /* bondangle c-b-d*/
  
  /*We are now declaring the variables for gold*/ 
  double ax,ay,az; /*Cartesian coordinates for gold*/
  double aepsilon,asigma; /*LJ parameters for gold*/ 
  
  /*We are now declaring the variables for the LJ calculation*/ 
  double abdist; /*Distance between a and b*/ 
  double acdist; /*Distance between a and c*/ 
  double addist; /*Distance between a and d*/ 
  double LJ_ab; /*LJ potential between a and b */ 
  double LJ_ac; /*LJ potential between a and c */ 
  double LJ_ad; /*LJ potential between a and d */ 
  double absigma;
  double acsigma;
  double adsigma; 
  double abepsilon; 
  double acepsilon; 
  double adepsilon; 
  double LJtotal; 

  /*These are the parameters for the centre of mass*/
  double bmass;
  double cmass;
  double dmass; 
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
  /* b is CH2! */ 
  bondlength = 0.9587; 
  bondangle = 104.47*(3.14159/180); 
  bx = 0.0;
  by = 0.0; 
  bz = 3.1448; /*This was calculated by dividing the R_min value of O by 2^1/6*/ 
  
  /*Given b, we can work out the position of c */ 
  
  cx = bondlength*sin(bondangle/2);
  cy = 0.0;
  cz = bz + bondlength*cos(bondangle/2); 
  /* printf("%lf %lf %lf \n", cx, cy, cz);  */
  
  /*given b, we can work out the position of d */ 
  dx = -1*bondlength*sin(bondangle/2); 
  dy = 0.0;
  dz = bz + bondlength*cos(bondangle/2);
  /* printf("%lf %lf %lf \n", dx, dy, dz); */

  /* Defining the triad LJ parameters*/ 
  /* Units should be changed to angstroms and kJ mol^-1 for sigma and epsilon respectively*/ 
  
  bepsilon = 0.1921; /* taken from klein et al, 2008, in Kcal mol^-1*/
  bsigma = 3.1448;  /*calculated from a R/2 value of 2.2730 angstoms*/
  cepsilon = 0.0; /*ditto*/ /*Corrected on the 20th July*/
  csigma = 0.0;  /*ditto*/ /*Corrected on the 20th July*/
  depsilon = 0.0; /*ditto*/ /*ditto*/
  dsigma = 0.0; /*ditto*/ /*ditto*/
  
  /*Define the gold positions*/ 
  /*COM of the triad*/ 
  
  theta = 1.570796;
  phi = 4.712389; 
  ax = 0.0; 
  ay = 0.0;
  az = 0.0; 
  
  /*Define the gold LJ parameters in kcal mol^-1*/ 
  /*asigma is in angstoms*/ 
  
  asigma  = 2.629; /*calculated from r_0 with values from heinz et al, 2008*/  
  aepsilon = 5.29; /*Units are in kcal/mol, again from heinz et al, 2008*/ 
  
  /*LJ calculation*/ 
  /*LJ interaction energy between a and b*/ 
  /* We need to check whether the cross terms are calculated correctly (can depend on force field!)*/

  absigma = (asigma + bsigma)/2; 
  abepsilon = sqrt(aepsilon*bepsilon); 
  abdist = sqrt((bx-ax)*(bx-ax) + (by-ay)*(by-ay) + (bz-az)*(bz-az)); 
  LJ_ab = 4*abepsilon*(pow((absigma/abdist),12)-pow((absigma/abdist),6)); 
  /* printf("%lf \n", LJ_ab); */
  
  /*LJ interaction energy between a and c*/ 
  /* We need to check whether the cross terms are calculated correctly (can depend on force field!)*/
  
  acsigma = (asigma + csigma)/2; 
  acepsilon = sqrt(aepsilon*cepsilon); 
  acdist = sqrt((cx-ax)*(cx-ax) + (cy-ay)*(cy-ay) + (cz-az)*(cz-az)); 
  LJ_ac = 4*acepsilon*(pow((acsigma/acdist),12)-pow((acsigma/acdist),6)); 
  /*printf("%lf \n", LJ_ac);*/
  
  /*LJ interaction energy between a and d*/ 
  /* We need to check whether the cross terms are calculated correctly (can depend on force field!)*/
  
  adsigma = (asigma + dsigma)/2; 
  adepsilon = sqrt(aepsilon*depsilon); 
  addist = sqrt((dx-ax)*(dx-ax) + (dy-ay)*(dy-ay) + (dz-az)*(dz-az)); 
  LJ_ad = 4*adepsilon*(pow((adsigma/addist),12)-pow((adsigma/addist),6)); 
  /* printf("%lf \n", LJ_ad);*/

  /*Sum of LJ energies*/ 

  LJtotal = LJ_ab + LJ_ac + LJ_ad; 
  /*printf("%lf \n", LJtotal);*/

  /*We are now calculating the COM of the triad*/

  bmass = 16.0;   /*We are assuming here that the b atom is CH2*/ 
  cmass = 1.0;   /* We are assuming here that the c atom is CH2*/ 
  dmass = 1.0;   /* We are assuming here that the d atom is CH2*/ 

  xmass = (bx*bmass + cx*cmass + dx*dmass)/(bmass + cmass + dmass);
  ymass = (by*bmass + cy*cmass + dy*dmass)/(bmass + cmass + dmass);
  zmass = (bz*bmass + cz*cmass + dz*dmass)/(bmass + cmass + dmass);

  bx-=xmass;
  by-=ymass;
  bz-=zmass;
  cx-=xmass;
  cy-=ymass;
  cz-=zmass;
  dx-=xmass;
  dy-=ymass;
  dz-=zmass;

  
  printf("# The cartesian coordinates for the centre of mass (COM) is %lf %lf %lf \n",xmass,ymass,zmass); 

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
  
  double boltzmann = 0.00198721; /* Units in kcal/K/mol*/
  double temp = 298.0; 
  double bf;
  double norm;
  double normV;
  double distb; 
  double distc; 
  double distd; 
  double bLJ;
  double cLJ;
  double dLJ;
  
  for (i=30;i<200;i++) { 
    
    r = i*0.075; /*Sets the distance between the COM and the Au particle */ 
    V_ave = 0;
    norm = 0;
    
    for (j=0;j<n;j++) {
      

      /*Taking into account the different energies of permutations of the molecules*/ 

      newgoldx = r*sin(axx[j])*cos(ayy[j]);
      newgoldy = r*sin(axx[j])*sin(ayy[j]); 
      newgoldz = r*cos(axx[j]);

      /*distances between the moving gold atom and the molecule*/

      distb = sqrt((bx-newgoldx)*(bx-newgoldx) + (by-newgoldy)*(by-newgoldy) + (bz-newgoldz)*(bz-newgoldz));
      distc = sqrt((cx-newgoldx)*(cx-newgoldx) + (cy-newgoldy)*(cy-newgoldy) + (cz-newgoldz)*(cz-newgoldz));
      distd = sqrt((dx-newgoldx)*(dx-newgoldx) + (dy-newgoldy)*(dy-newgoldy) + (dz-newgoldz)*(dz-newgoldz));

      // COMdist = sqrt((xmass-newgoldx)*(xmass-newgoldx) + (ymass-newgoldy)*(ymass-newgoldy) + (zmass-newgoldz)*(zmass-newgoldz));

      bLJ = 4*abepsilon*(pow((absigma/distb),12)-pow((absigma/distb),6));
      cLJ = 4*acepsilon*(pow((acsigma/distc),12)-pow((acsigma/distc),6));
      dLJ = 4*adepsilon*(pow((adsigma/distd),12)-pow((adsigma/distd),6));

      goldLJ = bLJ + cLJ + dLJ; 
      bf = exp(-goldLJ/(boltzmann*temp)); 
      bf=1.0;
      //printf("%lf %lf %lf %lf %lf %lf %lf %lf %lf\n",newgoldx,newgoldy,newgoldz,distb,distc,distd,bx+cx+dx,by+cy+dy,bz+cz+dz);
      /*We are now moving the gold particle around according to ONE of the given coordinates of theta and phi, then calculating the LJ energy from it!*/

      //V_ave += goldLJ*azz[j];
      
       V_ave += goldLJ*azz[j]*bf;  /*Now we calculate the average over all the polar coordinates using the weighted values*/

      /*Normalising factor*/

       norm +=azz[j]*bf;
	//  normV = V_ave/norm;
      
      
    }

    printf("%lf %lf \n",r,V_ave/norm);
    
  }
    
  return(0);
}


