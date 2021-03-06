#include <iostream>
#include <fstream>
#include <utility>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/mpi/mpi.hpp>
//#include <boost/serialization/access.hpp>
#include <omp.h>
#include <mpi.h>
#include "./quadtree.cc"

using namespace std;
using namespace boost::numeric::odeint;

int part_size;
int p_index;
int name_length;
char name[MPI_MAX_PROCESSOR_NAME];

#pragma omp declare reduction(vec_add : std::vector<double> : \
                              std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<double>())) \
                    initializer(omp_priv = omp_orig)

struct swarm_barnes_hut {
    vector<double> omega;
    const size_t n;
    double J, K, theta;

    swarm_barnes_hut(const size_t n_, double J_, double K_, double theta_, int rank):
        n(n_), omega(n_, 0.1), J(J_), K(K_), theta(theta_) {}

    void operator()(const vector<double> &x, vector<double> &dxdt, double t) const {
        // initialize QuadTree
        QuadTree tree;
        
        

            for(size_t i = 0; i < n; i++) {
                size_t xi = 3*i, yi = 3*i + 1, ti = 3*i + 2;
                dxdt[xi] = 0.;
                dxdt[yi] = 0.;
                dxdt[ti] = 0.1;
                // insert each point into the tree
                tree.insert(Point(x[xi], x[yi], x[ti]));
            }
        
        
        
//#pragma omp parallel for reduction(vec_add:dxdt) schedule(dynamic)
//        for(size_t i = 0; i < n; i++) {
//            size_t xi = 3*i, yi = 3*i + 1, ti = 3*i + 2;
//            std::vector<double> js = tree.get_centroids(positions[i][0], positions[i][1], theta);
//            for(size_t j = 0; j < js.size() / 4; j++) {
//                size_t xj = 4*j, yj = 4*j + 1, tj = 4*j + 2, mj = 4*j + 3;
//                double dx = js[xj] - positions[i][0],
//                       dy = js[yj] - positions[i][1],
//                       dth = js[tj] - phase[i],
//                       distance_sq = (dx*dx+dy*dy),
//                       distance = sqrt(distance_sq),
//                       xdot_contrib = (((1. + J*cos(dth))/distance - 1./distance_sq))/n,
//                       tdot = K/n*sin(dth)/distance,
//                       xdot = xdot_contrib * dx,
//                       ydot = xdot_contrib * dy;
//
//                dxdt[3*i] += xdot * js[mj];
//                dxdt[yi] += ydot * js[mj];
//                dxdt[ti] += tdot * js[mj];
//            }
//        }
        
#pragma omp parallel for reduction(vec_add:dxdt) schedule(dynamic)
        for(size_t i = 0; i < n; i++) {
            size_t xi = 3*i, yi = 3*i + 1, ti = 3*i + 2;
            std::vector<double> js = tree.get_centroids(x[xi], x[yi], theta);
            for(size_t j = 0; j < js.size() / 4; j++) {
                size_t xj = 4*j, yj = 4*j + 1, tj = 4*j + 2, mj = 4*j + 3;
                double dx = js[xj] - x[xi],
                       dy = js[yj] - x[yi],
                       dth = js[tj] - x[ti],
                       distance_sq = (dx*dx+dy*dy),
                       distance = sqrt(distance_sq),
                       xdot_contrib = (((1. + J*cos(dth))/distance - 1./distance_sq))/n,
                       tdot = K/n*sin(dth)/distance,
                       xdot = xdot_contrib * dx,
                       ydot = xdot_contrib * dy;

                dxdt[xi] += xdot * js[mj];
                dxdt[yi] += ydot * js[mj];
                dxdt[ti] += tdot * js[mj];
            }
        }
        
    
        MPI_Allgather(x + (rank * part_size), part_size, MPI_DOUBLE,
                      x, part_size, MPI_DOUBLE, MPI_COMM_WORLD);
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
    
    int rank, size;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
//    MPI_Type_contiguous(2, MPI_DOUBLE, &MPI_POSITION);
//    MPI_Type_contiguous(1, MPI_DOUBLE, &MPI_PHASE);
//    MPI_Type_commit(&MPI_POSITION);
//    MPI_Type_commit(&MPI_PHASE);
    
    MPI_Get_processor_name(name, &name_length);
    
    
    srand(time(NULL));
    const size_t n = 500;
    const double dt = 0.1;
    double positions[n][2];
    double phases[n];

    
    part_size = n / size;
    p_index = rank * part_size;
    // (0.1, 1) uniform
    // (0.1, -1) random
    // (1, 0) continuous rainbow
    // (1, -0.1) discrete rainbow
    // (1, -0.75) mixed rainbow
    const double J = 1, K = -0.1;
    const double theta_threshold = 0.5;

    vector<double> x(3*n);

    #pragma omp parallel for
    for(size_t i = 0; i < n; i++) {
        double r = ((double) rand())/((double) RAND_MAX)*1.;
        double theta = ((double) rand())/((double) RAND_MAX)*2.*M_PI;
        x[3*i] = r*cos(theta);
        x[3*i + 1] = r*sin(theta);
        x[3*i + 2] = ((double) rand())/((double) RAND_MAX)*2.*M_PI;
//
//        positions[i][0] = x[3*i];
//        positions[i][1] = x[3*i+1];
//        phases = x[3*i+2];
    }
        
        
    print_points(n, x, false);
    
//    MPI_Bcast(positions, n, MPI_POSITION, 0, MPI_COMM_WORLD);
//    MPI_Bcast(phases, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(x, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    // number of parallel threads
    omp_set_num_threads(stoi(argv[1]));

    swarm_barnes_hut group(n, J, K, theta_threshold, rank);
    
    double t0 = omp_get_wtime();
    
    integrate_const(runge_kutta4< vector<double> >(), boost::ref(group), x, 0., 50., dt);

    // if (rank == 0) {
        printf("Time taken: %f\n", omp_get_wtime()-t0);
    // }
    print_points(n, x, true);
    MPI_Finalize();
    
    return 0;
}
