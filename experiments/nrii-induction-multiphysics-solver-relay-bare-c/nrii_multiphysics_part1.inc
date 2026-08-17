#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>

/*
  Literature-grounded reduced thermo-electromagneto-hydrodynamic relay benchmark.
  ISO C17, no external numerical libraries.

  This is intentionally NOT a bitwise reproduction of Bermudez et al. (2009).
  It preserves the coupled ingredients: time-harmonic eddy-current surrogate,
  Joule heating, Lorentz forcing, nonlinear convective flow, thermal diffusion,
  and phase change. Multiple solver contracts are exercised with explicit handoff.
*/

#define NX 48
#define NY 64
#define N (NX*NY)
#define EMN (2*N)
#define MAX_GMRES_BASIS 20
#define ANDERSON_M 3
#define PI 3.14159265358979323846264338327950288

typedef struct {
    double *T, *Tprev;
    double *Ar, *Ai;
    double *u;
    double coil_I;
} State;

typedef struct {
    double beta_sigma;
    double em_coupling;
    double lorentz_gain;
    double buoy_gain;
    double conv_gain;
    double thermal_feedback;
    double phase_width;
    double latent;
    double outer_coupling;
    int force_phase_corner;
    const char *name;
} Stress;

typedef struct {
    long nr_uses;
    long gmres_uses, gmres_fail;
    long bicg_uses, bicg_fail;
    long picard_uses, picard_fail;
    long newton_krylov_uses, newton_krylov_fail;
    long thermal_newton_uses, thermal_newton_fail;
    long semismooth_uses, semismooth_fail;
    long fixed_point_uses, fixed_point_fail;
    long aitken_uses, aitken_fail;
    long anderson_uses, anderson_fail;
    long total_handoffs;
} Stats;

static inline int id(int i, int j) { return j*NX + i; }
static inline double clamp(double x,double a,double b){return x<a?a:(x>b?b:x);} 
static double norm2(const double *x, int n){ double s=0.0; for(int i=0;i<n;i++) s += x[i]*x[i]; return sqrt(s); }
static double norminf(const double *x, int n){ double m=0.0; for(int i=0;i<n;i++){ double a=fabs(x[i]); if(a>m)m=a;} return m; }
static double dot(const double *a,const double*b,int n){double s=0.0; for(int i=0;i<n;i++)s+=a[i]*b[i]; return s;}
static void axpy(double *y,double a,const double*x,int n){for(int i=0;i<n;i++)y[i]+=a*x[i];}
static void copyv(double*d,const double*s,int n){memcpy(d,s,(size_t)n*sizeof(double));}
static void setv(double*d,double v,int n){for(int i=0;i<n;i++)d[i]=v;}
static int boundary(int i,int j){return i==0||j==0||i==NX-1||j==NY-1;}

/* ---------- scalar 5-point helpers ---------- */
static double lap_at(const double *x,int i,int j){
    int k=id(i,j);
    if(boundary(i,j)) return 0.0;
    return x[id(i-1,j)] + x[id(i+1,j)] + x[id(i,j-1)] + x[id(i,j+1)] - 4.0*x[k];
}
static double dx_at(const double *x,int i,int j){
    if(i==0||i==NX-1) return 0.0;
    return 0.5*(x[id(i+1,j)]-x[id(i-1,j)]);
}
static double dy_at(const double *x,int i,int j){
    if(j==0||j==NY-1) return 0.0;
    return 0.5*(x[id(i,j+1)]-x[id(i,j-1)]);
}

/* ---------- NRII scalar RL implicit-midpoint ---------- */
static double nrii_rl_step(double I0,double dt,double R,double L,double V,int order,double *resid){
    /* I1 = I0 + dt/L * (V - R*(I0+I1)/2). */
    double a = -R/L;
    double b = V/L;
    double sum = I0;
    double c = dt*(a*I0 + b);
    sum += c;
    for(int k=2;k<=order;k++){
        c = 0.5*dt*a*c;
        sum += c;
    }
    *resid = fabs(sum - I0 - dt*(b + 0.5*a*(I0+sum)));
    return sum;
}

/* ---------- conductivity and phase enthalpy ---------- */
static double conductivity(double T,const Stress *s){
    double e = s->beta_sigma*(T-0.45);
    e = clamp(e,-8.0,8.0);
    return exp(e);
}
static double liquid_fraction(double T,const Stress*s){
    double Ts=0.60, Tl=Ts+s->phase_width;
    if(T<=Ts) return 0.0;
    if(T>=Tl) return 1.0;
    return (T-Ts)/(Tl-Ts);
}
static double enthalpy(double T,const Stress*s){return T + s->latent*liquid_fraction(T,s);} 
static double dheff(double T,const Stress*s,int semismooth){
    double Ts=0.60, Tl=Ts+s->phase_width;
    double eps=2e-10;
    if(fabs(T-Ts)<eps || fabs(T-Tl)<eps){
        if(!semismooth) return NAN;
        /* generalized derivative chosen from active interval */
        return 1.0 + 0.5*s->latent/s->phase_width;
    }
    if(T>Ts && T<Tl) return 1.0 + s->latent/s->phase_width;
    return 1.0;
}

/* ---------- EM operator ---------- */
typedef struct { const double *T; const Stress *s; double source_amp; } EMCtx;
static void em_apply(const void *ctxv,const double *x,double *y){
    const EMCtx *ctx=(const EMCtx*)ctxv;
    const double *ar=x, *ai=x+N;
    double *yr=y,*yi=y+N;
    for(int j=0;j<NY;j++) for(int i=0;i<NX;i++){
        int k=id(i,j);
        if(boundary(i,j)){yr[k]=ar[k]; yi[k]=ai[k]; continue;}
        double sig=conductivity(ctx->T[k],ctx->s);
        double c=ctx->s->em_coupling*sig;
        double nu=1.0/sqrt(0.05+sig);
        yr[k] = nu*(4.0*ar[k]-ar[id(i-1,j)]-ar[id(i+1,j)]-ar[id(i,j-1)]-ar[id(i,j+1)]) + c*ai[k];
        yi[k] = nu*(4.0*ai[k]-ai[id(i-1,j)]-ai[id(i+1,j)]-ai[id(i,j-1)]-ai[id(i,j+1)]) - c*ar[k];
    }
}
static void em_rhs(double *b,double amp){
    setv(b,0.0,EMN);
    /* localized inductor excitation near outer/right boundary */
    for(int j=NY/4;j<3*NY/4;j++){
        int i=NX-2; b[id(i,j)] = amp;
    }
}
