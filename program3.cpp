/* 	Verlet algorithm for the radial wave equation in n dimensions
	01.06.2012 Adrian Grütter
*/


#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <vector>
using namespace std;

int main ()
{
	int n=1000;
	double dr=0.01; //position incrementation, so r=dr*i. in other words the position is the little dr times the index (i.e. how many times did I go dr steps away from r=0?).
	double dt=0.01; // we say that the wave propagates dr (i.e. one index of phi) in the time period dt (speed=1)
	int T=1000;
	int dimension=1; //what spatial dimension are we in

	vector<double> phi(n); //the values of phi at the positions r
	vector<double> phivel(n); //the values of the rate of change of phi wrt time at the positions r
	vector<double> phiaccel(n); //the 'wave equation', essentialy d2 / dt² phi at positions r

	for(int i=0; i <= n; i++) //make a 'clean' run-through of the vectors
	{ phi[i]=0.;
	phivel[i]=0.; }

	phi[1]=0.0; //set up boundary conditions
	phivel[10]=5.0; 

	ofstream output; //declare our output

	char filename[33];
	double height=0.0;
	output.open("output.csv");//open

	for(int t=0; t<=T; t++) //run our program through time (i.e. T*dt seconds)
	{

		for(int i=1; i<=n; i++) phiaccel[i] = (phi[i+1]+phi[i-1]-2.*phi[i]) / (dr*dr) + ((double)dimension-1.)/(dr*i) * (phi[i+1]-phi[i-1])/dr/2.;	
		for(int i=1; i<=n; i++) phivel[i]  += dt*phiaccel[i]; //prev value plus the change?
		for(int i=1; i<=n; i++) phi[i]     += dt*phivel[i];
		//for(int i=1; i<=n; i++) phiaccel[i] = (phi[i+1]+phi[i-1]-2.*phi[i]) / (dr*dr) + (dimension-1)/(dr*i) * (phi[i+1]-phi[i-1])/dr/2.;	
		//for(int i=1; i<=n; i++) phivel[i]  += 0.5*dt*phiaccel[i]; //prev value plus the change?
		if(t==200) for(int i=1; i<=n; i++) if(phi[i]>height) height=phi[i];
		if(t%200==0&&t>0){
			for(int i=1; i<=n; i++) output << dr*i <<"\t" << phi[i] + (t/200-1)*height << endl; // output the values of phi each position		}
	}

	output.close();//close



return 0;
}
	
	
