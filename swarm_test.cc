#include "./swarm_2d/swarm_2d.hh"
#include <cstdio>
#include "omp.h"

int main() {

	double J[] = {0.1,0.1,1.,1.,1.};
	double K[] = {1.,-1.,0.,-0.1,-0.75};
	char odir[256];

	for (int i=0;i<5;i++) {
		sprintf(odir,"sw.%d",i);
		swarm_2d_fsal sw(512,J[i],K[i],odir,true,true,0.5);
		double t0=omp_get_wtime();
		sw.solve(600,1e-5,1e-5,false,1800);
		printf("Time taken: %f\n",omp_get_wtime()-t0);
		printf("Min radius: %f\n",sw.min_radius());
		printf("Max radius: %f\n",sw.max_radius());
	}
}
