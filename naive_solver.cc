#include <iostream>
#include <fstream>
#include <utility>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/openmp/openmp.hpp>
#include <omp.h>

using namespace std;
using namespace boost::numeric::odeint;

#pragma omp declare reduction(vec_add : std::vector<double> : \
                              std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
                    initializer(omp_priv = omp_orig)

struct swarm {
    vector<double> omega;
    const size_t n;
    double J, K;

    swarm(const size_t n_, double J_, double K_)
        : n(n_), omega(n_, 0.1), J(J_), K(K_) {}

    // Update function
    void operator()(const vector<double> &x, vector<double> &dxdt, double t) const {
    	// Initialize position and phase velocities - omega_i are all set to 0.1
#pragma omp parallel for
        for(size_t i = 0; i < n; i++) {
            dxdt[3*i] = 0.;
            dxdt[3*i + 1] = 0.;
            dxdt[3*i + 2] = 0.1;
        }
        // Calculate position and phase velocities by iterating over all points
#pragma omp parallel for reduction(vec_add:dxdt) schedule(dynamic)
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

void print_points(const size_t n, const vector<double> &x, bool final) {
   	ofstream file;
    file.open(final ? "final.csv" : "init.csv");
    for(size_t i = 0; i < n; i++) {
        file << x[3*i] << "," << x[3*i + 1] << "," << x[3*i + 2] << endl;
    }
    file.close();
}

int main(int argc, char **argv) {
	srand(6);
	// Number of points in swarm
    const size_t n = stoi(argv[2]);

    // (J = 0.1, K = 1) uniform
    // (0.1, -1) random
    // (1, 0) continuous rainbow
    // (1, -0.1) discrete rainbow
    // (1, -0.75) mixed rainbow
    const double J = 1, K = -0.75, dt = 0.1;
    vector<double> x(3*n);

#pragma omp parallel for
    for(size_t i = 0; i < n; i++) {
        double r = ((double) rand())/((double) RAND_MAX)*1.;
        double theta = ((double) rand())/((double) RAND_MAX)*2.*M_PI;
        x[3*i] = r*cos(theta);
        x[3*i + 1] = r*sin(theta);
        x[3*i + 2] = ((double) rand())/((double) RAND_MAX)*2.*M_PI;
    }

    print_points(n, x, false);

    // Number of parallel threads
    omp_set_num_threads(stoi(argv[1]));

    swarm group(n, J, K);
    double t0 = omp_get_wtime();
    // Pass to boost library integrator
    integrate_const(runge_kutta4< vector<double> >(), boost::ref(group), x, 0., 50., dt);
    printf("Time taken: %f\n", omp_get_wtime()-t0);
    print_points(n, x, true);
}
