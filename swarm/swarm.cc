#include "swarm_2d.hh"
#include "../barnes_hut/quadtree.hh"
#include <cstdio>
#include <cmath>
#include <ctime>
#include <sys/stat.h>
#include <sys/types.h>
#include "omp.h"
#include <vector>

swarm_2d::swarm_2d(int N_,double J_,double K_,const char* filename_,bool parallel_,bool barnes_hut_,double theta_) :
    N(N_), J(J_), K(K_), filename(filename_), parallel(parallel_), barnes_hut(barnes_hut_), theta(theta_), d_o_frame(0) {

    // Create directory for output frames
    mkdir(filename,S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH);
}

void swarm_2d::swarm_ff_barnes_hut(double t_,double *in,double *out) {
    const double iN=1./N,KiN=K*iN;
    int threads=parallel?omp_get_max_threads():1;

    // Initialize quadtree
    QuadTree tree;
    // Initialize output array and insert inputs into tree
    for(int i=0;i<N;i++) {
        int xi=3*i,yi=3*i+1,ti=3*i+2;
        out[xi]=0;
        out[yi]=0;
        out[ti]=0.1;

        tree.insert(Point(in[xi],in[yi],in[ti]));
    }
    // Calculate tree centroids
    tree.update_centroid();

    // Calculate right hand side
    #pragma omp parallel for num_threads(threads) reduction(+:out[:3*N])
    for(int i=0;i<N;i++) {
        int xi=3*i,yi=3*i+1,ti=3*i+2;

        std::vector<double> js = tree.get_centroids(in[xi],in[yi],theta);

        for(int j=0; j<js.size()/4; j++){
            int xj=4*j,yj=4*j+1,tj=4*j+2,mj=4*j+3;

            double dx=js[xj]-in[xi],
                    dy=js[yj]-in[yi],
                    dth=js[tj]-in[ti],
                    irsq=1/(dx*dx+dy*dy),
                    ir=sqrt(irsq),
                    a=iN*((1+J*cos(dth))*ir-irsq),
                    vx=dx*a,vy=dy*a,
                    vth=KiN*sin(dth)*ir;

            out[xi]+=js[mj]*vx;
            out[yi]+=js[mj]*vy;
            out[ti]+=js[mj]*vth;
        }
    }
}

void swarm_2d::swarm_ff_standard(double t_,double *in,double *out) {
    const double iN=1./N,KiN=K*iN;
    int threads=parallel?omp_get_max_threads():1;

    // Initialize output array
    #pragma omp parallel for num_threads(threads)
    for(int i=0;i<N;i++) {
        int xi=3*i,yi=3*i+1,ti=3*i+2;
        out[xi]=0;
        out[yi]=0;
        out[ti]=0.1;
    }

    // Calculate right hand side
    #pragma omp parallel for num_threads(threads) reduction(+:out[:3*N])
    for(int i=0;i<N;i++) {
        int xi=3*i,yi=3*i+1,ti=3*i+2;

        for(int j=0;j<i;j++) {
            int xj=3*j,yj=3*j+1,tj=3*j+2;
            double dx=in[xj]-in[xi],
                    dy=in[yj]-in[yi],
                    dth=in[tj]-in[ti],
                    irsq=1/(dx*dx+dy*dy),
                    ir=sqrt(irsq),
                    a=iN*((1+J*cos(dth))*ir-irsq),
                    vx=dx*a,vy=dy*a,
                    vth=KiN*sin(dth)*ir;

            out[xi]+=vx;out[xj]-=vx;
            out[yi]+=vy;out[yj]-=vy;
            out[ti]+=vth;out[tj]-=vth;
        }
    }
}

void swarm_2d::swarm_ff(double t_,double *in,double *out) {
    if (barnes_hut) {
        swarm_ff_barnes_hut(t_,in,out);
    } else {
        swarm_ff_standard(t_,in,out);
    }
}

void swarm_2d::swarm_init(double *q) {
    double x,y;
    srand(time(NULL));

    for(int i=0;i<N;i++) {
        do {
            x=rnd(-1,1);
            y=rnd(-1,1);
        } while(x*x+y*y>1);

        *(q++)=x;*(q++)=y;
        *(q++)=rnd(0,2*M_PI);
    }
}

void swarm_2d::swarm_print_dense(double t_,double *in) {
    char buf[256];

    sprintf(buf,"%s/fr.%d",filename,d_o_frame++);
    FILE *f=fopen(buf,"w");
    if(f==NULL) {
        fputs("Error opening output file\n",stderr);
        exit(1);
    }

    fprintf(f,"# Time: %.10g\n",t_);
    for(double *p=in;p<in+3*N;p+=3)
        fprintf(f,"%.8g %.8g %.8g\n",*p,p[1],p[2]);

    fclose(f);
}

double* swarm_2d::swarm_center(double *q) {
    double x_cen=0.,y_cen=0.,x,y;
    for(int i=0;i<N;i++) {
        x=q[0];
        y=q[1];
        x_cen+=x;
        y_cen+=y;
        q+=3;
    }
    double* center = new double[2];
    center[0]=x_cen/N;
    center[1]=y_cen/N;
    return center;
}

double swarm_2d::swarm_min_radius(double *q) {
    double min_rad=10.,x,y,rad,xd,yd;
    double *center = swarm_center(q);
    for(int i=0;i<N;i++) {
        x=q[0];
        y=q[1];
        xd=x-center[0];
        yd=y-center[1];
        rad=sqrt(xd*xd+yd*yd);
        if(rad<min_rad) min_rad=rad;
        q+=3;
    }
    return min_rad;
}

double swarm_2d::swarm_max_radius(double *q) {
    double max_rad=-1.,x,y,rad,xd,yd;
    double *center = swarm_center(q);
    for(int i=0;i<N;i++) {
        x=q[0];
        y=q[1];
        xd=x-center[0];
        yd=y-center[1];
        rad=sqrt(xd*xd+yd*yd);
        if(rad>max_rad) max_rad=rad;
        q+=3;
    }
    return max_rad;
}