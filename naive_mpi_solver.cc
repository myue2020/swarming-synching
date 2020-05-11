#include <iostream>
#include <fstream>
#include <utility>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/openmp/openmp.hpp>
#include <boost/numeric/odeint/external/mpi/mpi.hpp>
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
    void operator()(const mpi_state< vector<double> > &x, mpi_state< vector<double> > &dxdt, double t) const {
        assert(x().size() % 3 == 0);
        vector<double> xx(3*n);
        size_t start = 3*n / x.world.size() * x.world.rank();
        // Get all other positions and phases
        copy(x().begin(), x().end(), xx.begin() + start);
        for(size_t i = 1; i < x.world.size(); i++) {
            int in_rank = (x.world.rank() - i) % x.world.size(), out_rank = (x.world.rank() + i) % x.world.size();
            vector<double> temp;
            x.world.send(out_rank, 0, x());
            x.world.recv(in_rank, 0, temp);
            copy(temp.begin(), temp.end(), xx.begin() + 3*n / x.world.size() * in_rank);
        }
        vector<double> dxxdt = dxdt();
        // Initialize position and phase velocities - omega_i are all set to 0.1
#pragma omp parallel for
        for(size_t i = 0; i < dxdt().size() / 3; i++) {
            dxxdt[3*i] = 0.;
            dxxdt[3*i + 1] = 0.;
            dxxdt[3*i + 2] = 0.1;
        }
        // Calculate position and phase velocities by iterating over all points
#pragma omp parallel for reduction(vec_add:dxxdt) schedule(dynamic)
        for(size_t i = 0; i < dxxdt.size() / 3; i++) {
            size_t xi = start + 3*i, yi = start + 3*i + 1, ti = start + 3*i + 2;
            for(size_t j = 0; j < n; j++) {
                if (j != i + start / 3) {
                    int xj = 3*j, yj = 3*j + 1, tj = 3*j + 2;
                    double dx = xx[xj] - xx[xi],
                           dy = xx[yj] - xx[yi],
                           dth = xx[tj] - xx[ti],
                           distance_sq = (dx*dx+dy*dy),
                           distance = sqrt(distance_sq),
                           xdot_contrib = (((1. + J*cos(dth))/distance - 1./distance_sq))/n,
                           tdot = K/n*sin(dth)/distance,
                           xdot = xdot_contrib * dx,
                           ydot = xdot_contrib * dy;

                        dxxdt[xi - start] += xdot;
                        dxxdt[yi - start] += ydot;
                        dxxdt[ti - start] += tdot;
                }
            }
        }
#pragma omp parallel for
        for(size_t i = 0; i < dxxdt.size(); i++) {
            dxdt()[i] = dxxdt[i];
        }
    }
};

// Print csv of positions and phases at either initial or final time step
void print_points(const size_t n, const vector<double> &x, bool final) {
        ofstream file;
    file.open(final ? "final.csv" : "init.csv");
    for(size_t i = 0; i < n; i++) {
        file << x[3*i] << "," << x[3*i + 1] << "," << x[3*i + 2] << endl;
    }
    file.close();
}

int main(int argc, char **argv) {
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator world;

    // Number of parallel threads
    omp_set_num_threads(stoi(argv[1]));

    srand(time(NULL));
    // Number of points in swarm
    const size_t n = stoi(argv[2]);

    // (J = 0.1, K = 1) uniform
    // (0.1, -1) random
    // (1, 0) continuous rainbow
    // (1, -0.1) discrete rainbow
    // (1, -0.75) mixed rainbow
    const double J = 1., K = -0.1, dt = 0.1;
    vector<double> x(3*n);

    // Instantiate points in a circle with random positions and phases
    if (world.rank() == 0) {
#pragma omp parallel for
        for(size_t i = 0; i < n; i++) {
            double r = ((double) rand())/((double) RAND_MAX)*1.;
            double theta = ((double) rand())/((double) RAND_MAX)*2.*M_PI;
            x[3*i] = r*cos(theta);
            x[3*i + 1] = r*sin(theta);
            x[3*i + 2] = ((double) rand())/((double) RAND_MAX)*2.*M_PI;
        }

        print_points(n, x, false);
    }

    // Each processor gets own data
    mpi_state< vector<double> > x_split(world);
    split(x, x_split);

    swarm group(n, J, K);
    double t0 = omp_get_wtime();
    // Pass to boost library integrator
    integrate_const(runge_kutta4< mpi_state< vector<double> > >(), boost::ref(group), x_split, 0., 50., dt);
    if (world.rank() == 0) {
        printf("Time taken: %f\n", omp_get_wtime()-t0);
    }
    unsplit(x_split, x);
    if (world.rank() == 0) {
        print_points(n, x, true);
    }
}
