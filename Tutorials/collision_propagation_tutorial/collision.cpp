#include <iostream>
#include <math.h>
#include <cassert>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
using namespace std;
#define PI 3.141592653589

#define beta 9

int LX, LY, NR_iterations, nr_samples, dt_save;

double velDiffX, velDiffY;

double kin_visc = 0.1;
double tau = 1.0;	  //kinematic viscosity
double Force[2];  // force
double Left_Wall_velocity[2]; // velocity of the walls
double Right_Wall_velocity[2];

double WLIST[] = {4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};
double cX[] = {0.0, 1.0, 0.0, -1.0,  0.0, 1.0, -1.0, -1.0,  1.0};
double cY[] = {0.0, 0.0, 1.0,  0.0, -1.0, 1.0,  1.0, -1.0, -1.0};
int INVLIST[] = {0, 3, 4, 1, 2, 7, 8, 5, 6};

double ***f, ***f_temp, **rho, **u_x, **du_x, **u_y, **du_y, *f_eq, *f_eq_new;
int **obstacle;

char dens[200], velX[200], velY[200];
FILE  *densityFile, *veloFile;

void read_parameter_file(char inputname[])
{

  ifstream infile;

  infile.open(inputname);
  if(!infile)
    {
      cout << "read_parameter_file: Could not open " << inputname << endl;
      cout << "  ==> Exitting to the system..." << endl;
      exit(13);
    }

  if (infile.is_open())
    {
      cout << "read_parameter_file: reading system parameters from " << inputname << endl;
      infile.precision(11);

      cout << "reading LX and  LY =";
      infile >> LX >> LY;
      cout << LX << " and " << LY << endl;
      
      cout << "reading NR_iterations and nr_samples =";
      infile >> NR_iterations >> nr_samples;
      cout << NR_iterations << " and " <<  nr_samples << endl;
      
      cout << "relaxation parameter =";
      infile >> tau ;
      cout << tau << endl;    
      infile >> Force[0] >>  Force[1];
      cout << " Force.x=" << Force[0] << ", Frce.y=" <<  Force[1] << endl;
      
      infile.close();
      cout << "read_parameter_file: end." << endl;
    }

}







int PeriodicBC(int a, int b)
{
   int xnew = a % b;
   if(xnew < 0)
   xnew+=b;
   return xnew;

}

void allocateMemory(void) {

  f = new double**[LX];
  f_temp = new double**[LX];

  f_eq_new = new double[beta];
  f_eq = new double[beta];

  rho = new double*[LX];
  u_x = new double*[LX];
  u_y = new double*[LX];
  du_x = new double*[LX];
  du_y = new double*[LX];
  obstacle =  new int*[LX];

  for (int i=0;  i<LX; i++)
    {
      u_x[i]= new double[LY];
      u_y[i]= new double[LY];
      rho[i]= new double[LY];
	   du_x[i] = new double[LY];
  		du_y[i] = new double[LY];
      obstacle[i] = new int[LY];

      f[i] = new double*[LY];
      f_temp[i] = new double*[LY];
      for (int j=0;j<LY;j++)
	{
	  f[i][j]=new double[beta];
	  f_temp[i][j]=new double[beta];
	}
    }
}



void computeEquilibriumDistribution(int i, int j)
{
  for (int k=0; k<beta ; k++)
    {
      double cu=u_x[i][j]*cX[k] + u_y[i][j]*cY[k];
      double u2=u_x[i][j]*u_x[i][j]+u_y[i][j]*u_y[i][j];
      f_eq[k] = WLIST[k]*rho[i][j]* ( 1.0 - 1.5*u2 + 3.0*cu + 4.5*cu*cu );

//      f_eq[k] = WLIST[k]*rho[i][j]* ( 1.0 - 1.5*(u_x[i][j]*u_x[i][j] + u_y[i][j]*u_y[i][j]) + 3.0*(u_x[i][j]*cX[k]+u_y[i][j]*cY[k]) + 4.5*(u_x[i] [j]*cX[k]+u_y[i][j]*cY[k])*(u_x[i][j]*cX[k]+u_y[i][j]*cY[k]) );

      
    }
}



void computeDensityandVelocity(void)
{
double sum = 0;
double sumX, sumY ;
for(int i=0;i<LX;i++)
	{
	for(int j=0; j<LY;j++)
		{
			sum = 0; sumX = 0; sumY = 0;
			for (int k=0; k<beta ; k++)
			{
			sum = sum + f[i][j][k];
			sumX = sumX + f[i][j][k]*cX[k];
			sumY = sumY + f[i][j][k]*cY[k];
			}

			rho[i][j] = sum;
			u_x[i][j] = (sumX )/sum;
			du_x[i][j]= Force[0]/sum;
			du_y[i][j]= Force[1]/sum;
			u_y[i][j] = (sumY )/sum;

		}
	}

}


void collision(void)
{
  for(int i=0;i<LX;i++) {
    for(int j=0; j<LY;j++) {
      if( obstacle[i][j] != 1) {//apply collision for non-solid nodes only
	computeEquilibriumDistribution(i,j);
	
	
	for (int k=0; k<beta ; k++)
	{
	  f_temp[i][j][k] = f[i][j][k]*(1.0- 1.0/tau) + (1.0/tau)*f_eq[k] ;
	  f[i][j][k] = f_temp[i][j][k];
	}
      }
    }
  }
}



void writeData(int time)
{

sprintf(dens, "output/density_field_t%d.dat",time);
sprintf(velX, "output/velocity_field_t%d.dat",time);
densityFile = fopen(dens, "w");
veloFile = fopen(velX, "w");
 for (int x = 0; x < LX; x++)
   {
     for (int y = 0; y < LY; y++)
       {
	 if(obstacle[x][y]!=1) {
	   fprintf(densityFile, "%d %d %e\n", x,  y, rho[x][y]  );
	   fprintf(veloFile, "%d %d %e %e\n", x,  y, u_x[x][y], u_y[x][y]);
	 }
       }
   }
 fclose(densityFile);
 fclose(veloFile);

 int y=LY/2; int k = 5;
  sprintf(dens, "output/density_profile_y%d_t%d.dat",  y,  time);
  sprintf(velX, "output/velocity_profile_y%d_t%d.dat",  y, time);
  densityFile = fopen(dens, "w");
  veloFile = fopen(velX, "w");
  for (int x = 0; x < LX; x++)
   {
//     for (int y = 0; y < LY; y++)
       {
	 if(obstacle[x][y]!=1) {
	 
	 
	 double cu=u_x[x][y]*cX[k] + u_y[x][y]*cY[k];
         double u2=u_x[x][y]*u_x[x][y]+u_y[x][y]*u_y[x][y];
         double fTemp = f_temp[x][y][k];// WLIST[k]*rho[x][y]* ( 1.0 - 1.5*u2 + 3.0*cu + 4.5*cu*cu );
	 
	   fprintf(densityFile, "%d %e %e\n", x, rho[x][y], fTemp  );
	   fprintf(veloFile, "%d %e %e\n", x, u_x[x][y], u_y[x][y]);
	 }
       }
   }
 fclose(densityFile);
 fclose(veloFile);
}


// MAIN FUNCTION

int main()
{

double sum = 0;

 char obstacle_filename[200];
 char parameter_filename[200];


 sprintf(obstacle_filename, "obstacles.dat");
 sprintf(parameter_filename, "parameters.dat");

// read LX, LY, LZ, simulation_time, viscosity, gravity,  etc....
 read_parameter_file(parameter_filename);

 if(nr_samples>0) {
   dt_save = (int)(NR_iterations/nr_samples);
 }else{
   cout << "nr_samples must be > 0!: nr_samples=" << nr_samples << endl;
   cout << " ==> Exitting to the system..." << endl;
   exit(13);
 }

 tau = 3.0*kin_visc + 0.5 ;

 allocateMemory();

cout << " relaxation parameter tau is " <<tau <<endl;

//initializating density and velocity
cout<< "intializing desity and velocity " << endl;
for(int i=0;i<LX;i++)
  {
    for(int j=0; j<LY;j++)
      {
	rho[i][j] = 1.0; // density 
	u_x[i][j] = 0.0; // velocity in x direction
	u_y[i][j] = 0.0; // velocity in y direction
	obstacle[i][j] = 0;
		
      }

  }
cout<< "intializing desity and velocity done. "<<endl;


 //set f to f_eq for initialization
cout<< "setting f  intialization "<<endl;
// In most cases, f is set to f_eq at the start. Here we choose different initialization
 for(int i=0;i<LX;i++)
   {
     for(int j=0; j<LY;j++)
       {
	 //computeEquilibriumDistribution(i,j);
	 for (int k=0; k<beta ; k++)
	   {             
	     f[i][j][k]      =   rho[i][j]/9.0 ;
	     f_temp[i][j][k] =   rho[i][j]/9.0 ;
	   }
       }
   }
cout<< " intialization for f_eq done. "<<endl;


 //remove already existing files in the folder in case we are running the code again
(void) system("if [ -d output ] ; then rm -rf output; fi");
(void) system("mkdir output");
writeData(0);

ofstream collisionPopulationFile;
collisionPopulationFile.open("fStar.dat");

int xIndex = 5;
int yIndex = 5;

	collisionPopulationFile << 0 <<" "<<rho[xIndex][yIndex] <<" "<<u_x[xIndex][yIndex] <<" "<<u_y[xIndex][yIndex] <<" "<< f_temp[xIndex][yIndex][0] <<" "<< f_temp[xIndex][yIndex][1] <<" "<< f_temp[xIndex][yIndex][5] <<"\n";



for(int iter=1 ; iter <NR_iterations+1; iter ++ )
  {
  	
    computeDensityandVelocity();
    collision();
    
      
    if(iter%dt_save ==0)
      {
    		cout << " main: saving results on disk for t=" << iter<<endl;   
    		writeData(iter); 
		collisionPopulationFile << iter <<" "<<rho[xIndex][yIndex] <<" "<<u_x[xIndex][yIndex] <<" "<<u_y[xIndex][yIndex] <<" "<< f_temp[xIndex][yIndex][0] <<" "<< f_temp[xIndex][yIndex][1] <<" "<< f_temp[xIndex][yIndex][5] <<"\n"; 
      }
  }
  
  collisionPopulationFile.close();

 return 0;
}
