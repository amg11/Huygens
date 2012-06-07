/* 	Verlet algorithm for the radial wave equation in n dimensions
	01.06.2012 Adrian Grütter

	Need to adapt it for half-timesteps
	Need to have a closer look at boundary conditions
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
	double dr=0.01; //the physical separation between each index in the phi vector, so r=dr*i. in other words the position is the little dr times the index (i.e. how many times did I go dr away from r=0?).
	double dt=0.01; // we say that the wave propagates dr (i.e. one index of phi) in the time period dt.
	int T=1000;
	int dimension=1; //what dimension are we in

	vector<double> phi(n); //the values of phi at the positions r
	vector<double> phivel(n); //the values of the rate of change of phi wrt time at the positions r
	vector<double> phiaccel(n); //the 'wave equation', essentialy d2 / dt² phi at positions r
	vector<double> phitemp(n);

	//cout >> "Choose the value of Phi at r=0";
	//cin << phi[0]; // just as an example so far, this is how we'd read in boundary conditions. make sure they don't get 'deleted' by the loops running though!
	for(int i=0; i <= n; i++)
	{ phi[i]=0.;
	phivel[i]=0.; }
	phi[1]=0.0;
	phivel[10]=5.0; //phi[1] is r=0


	ofstream output; //how to include variable names in files? i.e. one that indicates at which timestep the snapshot was taken
	char filename[33];
	double height=0.0;
	output.open("output.csv");//open
	for(int t=0; t<=T; t++)
	{

		for(int i=1; i<=n; i++) phiaccel[i] = (phi[i+1]+phi[i-1]-2.*phi[i]) / (dr*dr) + ((double)dimension-1.)/(dr*i) * (phi[i+1]-phi[i-1])/dr/2.;	
		for(int i=1; i<=n; i++) phivel[i]  += dt*phiaccel[i]; //prev value plus the change?
		for(int i=1; i<=n; i++) phi[i]     += dt*phivel[i];
		//for(int i=1; i<=n; i++) phiaccel[i] = (phi[i+1]+phi[i-1]-2.*phi[i]) / (dr*dr) + (dimension-1)/(dr*i) * (phi[i+1]-phi[i-1])/dr/2.;	
		//for(int i=1; i<=n; i++) phivel[i]  += 0.5*dt*phiaccel[i]; //prev value plus the change?
		if(t==200) for(int i=1; i<=n; i++) if(phi[i]>height) height=phi[i];
		if(t%200==0&&t>0){
			for(int i=1; i<=n; i++) output << dr*i <<"\t" << phi[i] + (t/200-1)*height << endl; // output the positions and values of phi at those positions
		}
	}

	output.close();//close



return 0;
}
	
	
