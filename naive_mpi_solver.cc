#include <iostream>
#include <fstream>
#include <utility>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/openmp/openmp.hpp>
#include <boost/numeric/odeint/external/mpi/mpi.hpp>
#include <omp.h>
#include <mpi.h>

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

    void operator()(const mpi_state< vector<double> > &x, mpi_state< vector<double> > &dxdt, double t) const {
        vector<double> xx(n);
        size_t start = n / x.world.size() * x.world.rank();
        copy(temp.begin(), temp.end(), xx.begin() + start);
        for(size_t i = 1; i < x.world.size(); i++) {
            int in_rank = (x.world.rank() - i) % x.world.size(), out_rank = (x.world.rank() + i) % x.world.size();
            vector<double> temp;
            x.world.isend(out_rank, 0, x());
            x.world.irecv(in_rank, 0, temp).wait();
            copy(temp.begin(), temp.end(), xx.begin() + n / x.world.size() * in_rank);
        }
#pragma omp parallel for
        for(size_t i = 0; i < dxdt().size() / 3; i++) {
            dxdt[3*i] = 0.;
            dxdt[3*i + 1] = 0.;
            dxdt[3*i + 2] = 0.1;
        }

#pragma omp parallel for reduction(vec_add:dxdt()) schedule(dynamic)
        for(size_t i = 0; i < dxdt().size() / 3; i++) {
                size_t xi = start + 3*i, yi = start + 3*i + 1, ti = start + 3*i + 2;
            for(size_t j = 0; j < n; j++) {
                if (j != i) {
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

                        dxdt()[xi - start] += xdot;
                        dxdt()[yi - start] += ydot;
                        dxdt()[ti - start] += tdot;
                }
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
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator world;

    srand(time(NULL));
    const size_t n = 400;
    const double J = 1., K = -0.1, dt = 0.1;
    vector<double> x(3*n);

    printf("%d\n", world.rank());
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

    mpi_state< vector<double> > x_split(world);
    split(x, x_split);

    swarm group(n, J, K);
    double t0 = omp_get_wtime();
    integrate_const(runge_kutta4< mpi_state< vector<double> > >(), boost::ref(group), x_split, 0., 50., dt);
    if (rank == 0) {
    	printf("Time taken: %f\n", omp_get_wtime()-t0);
    }
    unsplit(x_split, x);
    if (world.rank() == 0) {
        print_points(n, x, true);
    }
}