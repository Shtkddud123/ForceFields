#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include <complex.h> 


int main() 
{ 
  /*                     d (CH2)    */
  /*                    /           */
  /*   a (Au)          b (O)        */
  /*                    \           */
  /*                     c (CH2)    */
  /*---------------------------------*/
  /*This code is for the Au-EO bead!  */
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
  bondlength = 1.458; 
  bondangle = 111.40*3.14159/180; 
  bx = 0.0;
  by = 0.0; 
  bz = 3.154; 
  
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
  
  bepsilon = 0.1521; /*taken from shinoda et al, 2004 */
  bsigma = 3.1537; /*taken from shinoda et al, 2004*/ 
  cepsilon = 0.1004; /* taken from klein et al, 2008 */
  csigma = 4.05;  /* taken from klein et al, 2008 */
  depsilon = 0.1004; /* taken from klein et al, 2008 */
  dsigma = 4.05; /* taken from klein et al, 2008 */
  
  /*Define the gold positions*/ 
  /*COM of the triad*/ 
  
  theta = 1.570796;
  phi = 4.712389; 
  ax = 0.0; 
  ay = 0.0;
  az = 0.0; 
  
  /*Define the gold LJ parameters*/ 
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

  bmass = 16.0;   /*We are assuming here that the b atom is O*/ 
  cmass = 14.0;   /* We are assuming here that the c atom is CH2*/ 
  dmass = 14.0;   /* We are assuming here that the d atom is CH2*/ 

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



  double bdist;
  double cdist;
  double ddist;
  
  double bLJ; 
  double cLJ;
  double dLJ; 

  double boltzmann = 0.00198721; /* Units in kcal/K/mol*/
  double temp = 298.0; 
  double bf;
  double norm;


  
  for (i=35;i<200;i++) { 
    
    r = i*0.075; /*Sets the distance between the COM and the Au particle */ 
    V_ave = 0; 
    norm = 0;
    
    
    for (j=0;j<n;j++) {
      
      newgoldx = r*sin(axx[j])*cos(ayy[j]);
      newgoldy = r*sin(axx[j])*sin(ayy[j]); 
      newgoldz = r*cos(axx[j]);


      bdist = sqrt((bx-newgoldx)*(bx-newgoldx) + (by-newgoldy)*(by-newgoldy) + (bz-newgoldz)*(bz-newgoldz));
      cdist = sqrt((cx-newgoldx)*(cx-newgoldx) + (cy-newgoldy)*(cy-newgoldy) + (cz-newgoldz)*(cz-newgoldz));
      ddist = sqrt((dx-newgoldx)*(dx-newgoldx) + (dy-newgoldy)*(dy-newgoldy) + (dz-newgoldz)*(dz-newgoldz));


      //COMdist = sqrt((xmass-newgoldx)*(xmass-newgoldx) + (ymass-newgoldy)*(ymass-newgoldy) + (zmass-newgoldz)*(zmass-newgoldz));

      bLJ = 4*abepsilon*(pow((absigma/bdist),12)-pow((absigma/bdist),6));
      cLJ = 4*acepsilon*(pow((acsigma/cdist),12)-pow((acsigma/cdist),6));
      dLJ = 4*adepsilon*(pow((adsigma/ddist),12)-pow((adsigma/ddist),6));



      goldLJ = bLJ + cLJ + dLJ; /*We are assuming for now sigma for the COM is identical to that of ab,ac,and ad*/ 
       bf = exp(-goldLJ/(boltzmann*temp)); 
       bf=1.0;
      
      /*We are now moving the gold particle around according to ONE of the given coordinates of theta and phi, then calculating the LJ energy from it!*/
      
      V_ave += goldLJ*azz[j]*bf;  /*Now we calculate the average over all the polar coordinates using the weighted values*/

      norm +=azz[j]*bf;
      
    }

    printf("%lf %lf \n",r,V_ave/norm);
    
  }
    
  return(0);
}


