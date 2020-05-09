#ifndef SWARM_HH
#define SWARM_HH

#include <cstdlib>

#include "../solvers/sol_fsal.hh"
#include "../solvers/sol_ckr.hh"


class swarm_2d {
    public:
        /** The number of actors. */
        const int N;
        const double J;
        const double K;
        const char* filename;
        bool parallel;
        bool barnes_hut;
        double theta;
        swarm_2d(int N_,double J_,double K_,const char* filename_,bool parallel_, bool barnes_hut_,double theta_);
        void swarm_ff(double t_,double *in,double *out);
        void swarm_init(double *q);
        void swarm_print_dense(double t_,double *in);
        inline double rnd(double l,double u) {
            return l+(u-l)/RAND_MAX*static_cast<double>(rand());
        }
        double swarm_min_radius(double *q);
        double swarm_max_radius(double *q);
    private:
        int d_o_frame;
        void swarm_ff_barnes_hut(double t_,double *in,double *out);
        void swarm_ff_standard(double t_,double *in,double *out);
        double* swarm_center(double *q);
};

class swarm_2d_fsal : public sol_fsal, public swarm_2d {
    public:
        swarm_2d_fsal(int N_,double J_,double K_,const char* filename_,bool parallel_=true,bool barnes_hut_=false,double theta_=0.5) :
            sol_fsal(3*N_), swarm_2d(N_,J_,K_,filename_,parallel_,barnes_hut_,theta_) {}
        virtual void ff(double t_,double *in,double *out) {swarm_ff(t_,in,out);}
        virtual void init() {swarm_init(q);}
        virtual void print_dense(double t_,double *in) {swarm_print_dense(t_,in);}
        virtual double min_radius() {return swarm_min_radius(q);}
        virtual double max_radius() {return swarm_max_radius(q);}
};

class swarm_2d_ckr : public sol_ckr, public swarm_2d {
    public:
        swarm_2d_ckr(int N_,double J_,double K_,const char* filename_,bool parallel_=true,bool barnes_hut_=false,double theta_=0.5) :
            sol_ckr(3*N_), swarm_2d(N_,J_,K_,filename_,parallel_,barnes_hut_,theta_) {}
        virtual void ff(double t_,double *in,double *out) {swarm_ff(t_,in,out);}
        virtual void init() {swarm_init(q);}
        virtual void print_dense(double t_,double *in) {swarm_print_dense(t_,in);}
        virtual double min_radius() {return swarm_min_radius(q);}
        virtual double max_radius() {return swarm_max_radius(q);}
};

#endif
