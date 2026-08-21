#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <stdint.h>

#define NDOF 128
#define MAX_IT 24
#define AA_M 3
#define NRII_ORDER 12
#define MAX_STRESS 121
#define PI 3.141592653589793238462643383279502884

typedef enum {
    M_NEWTON=0,
    M_LINESEARCH,
    M_AITKEN,
    M_ANDERSON,
    M_BROYDEN_NEWTON,
    M_PTC_NEWTON,
    M_NRII_NEWTON,
    M_COUNT
} Method;

static const char *METHOD_NAME[M_COUNT] = {
    "Newton",
    "Newton+LineSearch",
    "Newton+Aitken",
    "Newton+Anderson",
    "Broyden->Newton",
    "PseudoTransient->Newton",
    "NRII->Newton"
};

typedef struct {
    double dt, k, alpha, kc, beta;
    double qstar[NDOF];
    double fext[NDOF];
    double rhs_scale;
} Problem;

typedef struct {
    int success;
    int branch_ok;
    int iterations;
    int residual_evals;
    int jacobian_evals;
    int linear_solves;
    int line_search_evals;
    int accel_accepts;
    int handoff_used;
    int nrii_accepted;
    int nrii_order;
    double rel_residual;
    double rel_forward_error;
    double wall_us;
} SolveStats;

static double now_sec(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9*(double)ts.tv_nsec;
}

static double infnorm(const double *x, int n) {
    double m=0.0;
    for(int i=0;i<n;i++) { double a=fabs(x[i]); if(a>m)m=a; }
    return m;
}

static double dotv(const double *a, const double *b, int n) {
    double s=0.0;
    for(int i=0;i<n;i++) s += a[i]*b[i];
    return s;
}

static int finite_vec(const double *x, int n) {
    for(int i=0;i<n;i++) if(!isfinite(x[i])) return 0;
    return 1;
}

static void internal_force(const Problem *P, const double *q, double *f) {
    for(int i=0;i<NDOF;i++) f[i] = P->k*q[i] + P->alpha*q[i]*q[i]*q[i];
    for(int i=0;i<NDOF-1;i++) {
        double d=q[i+1]-q[i];
        double e=P->kc*d + P->beta*d*d*d;
        f[i]   -= e;
        f[i+1] += e;
    }
}

static void make_problem(Problem *P, double stress) {
    P->dt=1.0; P->k=0.10; P->alpha=1.0; P->kc=0.20; P->beta=0.02;
    for(int i=0;i<NDOF;i++) {
        double x=((double)i+0.5)/(double)NDOF;
        double shape=1.0 + 0.10*sin(2.0*PI*x) + 0.05*sin(7.0*PI*x) + 0.025*cos(11.0*PI*x);
        P->qstar[i]=stress*shape;
    }
    double fint[NDOF];
    internal_force(P,P->qstar,fint);
    for(int i=0;i<NDOF;i++) {
        /* Implicit Euler, q0=v0=0, unit mass:
           q1 + dt^2 f_int(q1) - dt^2 f_ext = 0. */
        P->fext[i]=(P->qstar[i]+P->dt*P->dt*fint[i])/(P->dt*P->dt);
    }
    P->rhs_scale = 1.0 + P->dt*P->dt*infnorm(P->fext,NDOF);
}

static void residual(const Problem *P, const double *q, double *r) {
    double fint[NDOF]; internal_force(P,q,fint);
    double d2=P->dt*P->dt;
    for(int i=0;i<NDOF;i++) r[i]=q[i] + d2*fint[i] - d2*P->fext[i];
}

static void jacobian_tridiag(const Problem *P, const double *q, double *lo, double *di, double *up, double shift) {
    double d2=P->dt*P->dt;
    for(int i=0;i<NDOF;i++) di[i]=1.0 + d2*(P->k + 3.0*P->alpha*q[i]*q[i]) + shift;
    for(int i=0;i<NDOF-1;i++) {
        double d=q[i+1]-q[i];
        double ke=d2*(P->kc + 3.0*P->beta*d*d);
        di[i]+=ke; di[i+1]+=ke; lo[i]=-ke; up[i]=-ke;
    }
}

static int thomas_solve(const double *lo0,const double *di0,const double *up0,const double *b0,double *x,int n) {
    double *di=(double*)malloc((size_t)n*sizeof(double));
    double *b=(double*)malloc((size_t)n*sizeof(double));
    double *up=(double*)malloc((size_t)(n-1)*sizeof(double));
    if(!di||!b||!up){free(di);free(b);free(up);return 0;}
    memcpy(di,di0,(size_t)n*sizeof(double)); memcpy(b,b0,(size_t)n*sizeof(double)); memcpy(up,up0,(size_t)(n-1)*sizeof(double));
    for(int i=1;i<n;i++) {
        if(fabs(di[i-1])<1e-30 || !isfinite(di[i-1])) {free(di);free(b);free(up);return 0;}
        double w=lo0[i-1]/di[i-1]; di[i]-=w*up[i-1]; b[i]-=w*b[i-1];
    }
    if(fabs(di[n-1])<1e-30 || !isfinite(di[n-1])) {free(di);free(b);free(up);return 0;}
    x[n-1]=b[n-1]/di[n-1];
    for(int i=n-2;i>=0;i--) {
