#include <iostream>
#include <utility>
#include <boost/numeric/odeint.hpp>
#include "omp.h"
#include "../barnes_hut/quadtree.hh"

using namespace std;
using namespace boost::numeric::odeint;

struct swarm_barnes_hut {
    vector<double> omega;
    const size_t n;
    double J, K;

    swarm_barnes_hut(const size_t n_, double J_ = 1., double K_ = 1.)
    	: n(n_), omega(n_, 0.1), J(J_), K(K_) {}

    void operator()(const vector<double> &x, vector<double> &dxdt, double t) const {
    	printf("Iteration\n");
// #pragma omp parallel for
        QuadTree tree;
        // initialize output array
	    for(size_t i = 0; i < n; i++) {
	        dxdt[3*i] = 0;
	       	dxdt[3*i + 1] = 0;
	        dxdt[3*i + 2] = 0.1;

            tree.insert(Point())
	    }

#pragma omp parallel for reduction(+:dxdt[:3*N])
        for(size_t i = 0; i < n; i++) {
        	size_t xi = 3*i, yi = 3*i + 1, ti = 3*i + 2;
            for(size_t j = 0; j < i; j++) {
	            int xj = 3*j, yj = 3*j + 1, tj = 3*j + 2;
	            double dx = x[xj] - x[xi],
	                   dy = x[yj] - x[yi],
	                   dth = x[tj] - x[ti],
	                   distance_sq = (dx*dx+dy*dy),
	                   distance = sqrt(distance_sq),
	            	   xdot_contrib = (((1. + J*cos(dth))/distance - 1./distance_sq))/n,
	            	   tdot = K/n*sin(dth)/distance,
	            	   xdot = xdot_contrib * dx,
	            	   ydot = xdot_contrib * dy;

	        	dxdt[xi] += xdot;
	        	dxdt[yi] += ydot;
	        	dxdt[ti] += tdot;
	        	dxdt[xj] -= xdot;
	        	dxdt[yj] -= ydot;
	        	dxdt[tj] -= tdot;
	        }
	    }
	}
};

int main(int argc, char **argv) {
	srand(time(NULL));
    const size_t n = 16384;
    const double J = 1., K = 1., dt = 0.1;
    vector<double> x(3*n);

#pragma omp parallel for
    for(size_t i = 0; i < n; i++) {
    	double r = ((double) rand())/((double) RAND_MAX);
    	double theta = ((double) rand())/((double) RAND_MAX)*2.*M_PI;
        x[3*i] = r*cos(theta);
        x[3*i + 1] = r*sin(theta);
        x[3*i + 2] = 0.1;
    }

    swarm_barnes_hut group(n, J, K);
    integrate_const(runge_kutta4< vector<double> >(), boost::ref(group), x, 0., 100., dt);
}
