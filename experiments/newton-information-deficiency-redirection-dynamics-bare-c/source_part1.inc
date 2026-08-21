#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>

#define NDOF 128
#ifndef MAX_STAGES
#define MAX_STAGES 24
#endif
#define AA_M 3
#ifndef DISABLE_AFTER_SPECIALIST_FAILURE
#define DISABLE_AFTER_SPECIALIST_FAILURE 0
#endif
#define NRII_ORDER 12
#define MAX_STRESS 121
#define PI 3.141592653589793238462643383279502884

typedef enum {
    M_NEWTON=0,
    M_NEWTON_LINESEARCH_ALT,
    M_NEWTON_AITKEN_ALT,
    M_NEWTON_ANDERSON_ALT,
    M_BROYDEN_NEWTON_ALT,
    M_PTC_NEWTON_ALT,
    M_NRII_NEWTON_ALT,
    M_COUNT
} Method;

static const char *METHOD_NAME[M_COUNT] = {
    "Newton",
    "Newton->LineSearch->Newton",
    "Newton->Aitken->Newton",
    "Newton->Anderson->Newton",
    "Newton->Broyden->Newton",
    "Newton->PseudoTransient->Newton",
    "Newton->NRII->Newton"
};

typedef struct {
    double dt, k, alpha, kc, beta;
    double qstar[NDOF];
    double fext[NDOF];
    double rhs_scale;
} Problem;

typedef struct {
    int success, branch_ok;
    int stages, stage_a, stage_b;
    int residual_evals, jacobian_evals, linear_solves, line_search_evals;
    int accel_accepts, handoff_count, nrii_attempts, nrii_accepts;
    int redirect_triggers, redirect_successes, redirect_failures;
    int trigger_budget_deficit, trigger_nondecrease, trigger_invalid;
    int newton_steps, specialist_steps;
    double rel_residual, rel_forward_error, wall_us;
} SolveStats;

static double now_sec(void){struct timespec ts;clock_gettime(CLOCK_MONOTONIC,&ts);return (double)ts.tv_sec+1e-9*(double)ts.tv_nsec;}
static double infnorm(const double*x,int n){double m=0;for(int i=0;i<n;i++){double a=fabs(x[i]);if(a>m)m=a;}return m;}
static double dotv(const double*a,const double*b,int n){double s=0;for(int i=0;i<n;i++)s+=a[i]*b[i];return s;}
static int finite_vec(const double*x,int n){for(int i=0;i<n;i++)if(!isfinite(x[i]))return 0;return 1;}

static void internal_force(const Problem*P,const double*q,double*f){
    for(int i=0;i<NDOF;i++)f[i]=P->k*q[i]+P->alpha*q[i]*q[i]*q[i];
    for(int i=0;i<NDOF-1;i++){
        double d=q[i+1]-q[i],e=P->kc*d+P->beta*d*d*d;
        f[i]-=e;f[i+1]+=e;
    }
}

static void make_problem(Problem*P,double stress){
    P->dt=1.0;P->k=0.10;P->alpha=1.0;P->kc=0.20;P->beta=0.02;
    for(int i=0;i<NDOF;i++){
        double x=((double)i+0.5)/(double)NDOF;
        double shape=1.0+0.10*sin(2*PI*x)+0.05*sin(7*PI*x)+0.025*cos(11*PI*x);
        P->qstar[i]=stress*shape;
    }
    double fint[NDOF];internal_force(P,P->qstar,fint);
    for(int i=0;i<NDOF;i++)P->fext[i]=P->qstar[i]+fint[i];
    P->rhs_scale=1.0+infnorm(P->fext,NDOF);
}

static void residual(const Problem*P,const double*q,double*r){
    double f[NDOF];internal_force(P,q,f);
    for(int i=0;i<NDOF;i++)r[i]=q[i]+f[i]-P->fext[i];
}

static void jacobian_tridiag(const Problem*P,const double*q,double*lo,double*di,double*up,double shift){
    for(int i=0;i<NDOF;i++)di[i]=1.0+P->k+3.0*P->alpha*q[i]*q[i]+shift;
    for(int i=0;i<NDOF-1;i++){
        double d=q[i+1]-q[i],ke=P->kc+3.0*P->beta*d*d;
        di[i]+=ke;di[i+1]+=ke;lo[i]=-ke;up[i]=-ke;
    }
}

static int thomas_solve(const double*lo0,const double*di0,const double*up0,const double*b0,double*x,int n){
    double di[NDOF],b[NDOF],up[NDOF-1];
    memcpy(di,di0,(size_t)n*sizeof(double));memcpy(b,b0,(size_t)n*sizeof(double));memcpy(up,up0,(size_t)(n-1)*sizeof(double));
    for(int i=1;i<n;i++){
        if(fabs(di[i-1])<1e-30||!isfinite(di[i-1]))return 0;
        double w=lo0[i-1]/di[i-1];di[i]-=w*up[i-1];b[i]-=w*b[i-1];
    }
    if(fabs(di[n-1])<1e-30||!isfinite(di[n-1]))return 0;
    x[n-1]=b[n-1]/di[n-1];
    for(int i=n-2;i>=0;i--){if(fabs(di[i])<1e-30||!isfinite(di[i]))return 0;x[i]=(b[i]-up[i]*x[i+1])/di[i];}
    return finite_vec(x,n);
}

static int newton_direction(const Problem*P,const double*q,const double*r,double shift,double*dx,SolveStats*S){
    double lo[NDOF-1],di[NDOF],up[NDOF-1],b[NDOF];
    jacobian_tridiag(P,q,lo,di,up,shift);S->jacobian_evals++;
    for(int i=0;i<NDOF;i++)b[i]=-r[i];
    if(!thomas_solve(lo,di,up,b,dx,NDOF))return 0;
    S->linear_solves++;return 1;
}

static double eval_relres(const Problem*P,const double*q,double*r,SolveStats*S){residual(P,q,r);S->residual_evals++;return infnorm(r,NDOF)/P->rhs_scale;}
static double forward_error(const Problem*P,const double*q){double m=0,s=1.0+infnorm(P->qstar,NDOF);for(int i=0;i<NDOF;i++){double e=fabs(q[i]-P->qstar[i]);if(e>m)m=e;}return m/s;}
