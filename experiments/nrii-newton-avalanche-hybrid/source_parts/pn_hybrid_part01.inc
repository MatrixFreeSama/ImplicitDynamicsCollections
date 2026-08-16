#include <bits/stdc++.h>
#include <omp.h>
using namespace std;

struct Params {
    double dx=0.25;
    double lambda2=4.0;     // Debye-like electrostatic length^2
    double Dn=1.0, Dp=0.8;
    double mun=1.0, mup=0.8;
    double krec=0.02;
    double dt=0.005;
    double Nd=8.0;          // dimensionless doping magnitude
    double pulse_amp=1e-3;
    double pulse_sigma_cells=6.0;
    // Reduced reverse-bias / impact-ionization constitutive stress model.
    // Total electrostatic potential = self-consistent phi + Vext.
    double reverse_drive=0.0;
    double bias_sigma_cells=12.0;
    double ion_A=50.0;
    double ion_B=20.0;
    double ion_gamma=0.18;
    double ion_Efloor=0.05;
    int eq_half=512;
    int local_half=768;
    int pmax=96;
    double tol=1e-10;
};

static inline double bulk_n(double C){ return 0.5*(C + sqrt(C*C+4.0)); }
static inline double bulk_p(double C){ return 1.0/bulk_n(C); }
static inline double doping_at_global(long long gi, long long N, const Params& P){
    long long c=N/2;
    return gi < c ? -P.Nd : P.Nd;
}

static inline double vext_at_global(long long gi,long long N,const Params&P){
    const double x=(double)(gi-N/2)/P.bias_sigma_cells;
    return P.reverse_drive*tanh(x);
}

static inline double avalanche_G_point(double n,double p,double phi_m,double phi_p,
                                       double vext_m,double vext_p,const Params&P,
                                       double *dGdn=nullptr,double *dGdp=nullptr,
                                       double *dGdphim=nullptr,double *dGdphip=nullptr){
    const double Etotal=-((phi_p+vext_p)-(phi_m+vext_m))/(2.0*P.dx);
    const double e=sqrt(Etotal*Etotal+P.ion_Efloor*P.ion_Efloor);
    const double alpha=P.ion_A*exp(-P.ion_B/e);
    const double H=alpha*e;
    const double V=P.mun*n+P.mup*p;
    const double G=P.ion_gamma*H*V;
    if(dGdn) *dGdn=P.ion_gamma*H*P.mun;
    if(dGdp) *dGdp=P.ion_gamma*H*P.mup;
    if(dGdphim||dGdphip){
        const double dHdE=alpha*(1.0+P.ion_B/e)*(Etotal/e);
        const double dGdE=P.ion_gamma*V*dHdE;
        if(dGdphim) *dGdphim=dGdE/(2.0*P.dx);
        if(dGdphip) *dGdphip=-dGdE/(2.0*P.dx);
    }
    return G;
}

struct EqProfile { vector<double> phi,n,p; int half; };

EqProfile build_equilibrium(const Params& P){
    int H=P.eq_half, M=2*H+1;
    vector<double> phi(M), n(M), p(M);
    double CL=-P.Nd, CR=P.Nd;
    double phiL=log(bulk_n(CL)), phiR=log(bulk_n(CR));
    for(int j=0;j<M;j++){
        double t=(double)j/(M-1);
        phi[j]=phiL*(1-t)+phiR*t;
    }
    double invh2=1.0/(P.dx*P.dx);
    // Newton scalar tridiagonal for PB: -lambda2 phi_xx -(exp(-phi)-exp(phi)+C)=0
    vector<double> a(M),b(M),c(M),d(M);
    for(int it=0;it<100;it++){
        double mr=0;
        for(int j=1;j<M-1;j++){
            double C = (j<H)?-P.Nd:P.Nd;
            double epos=exp(phi[j]), eneg=exp(-phi[j]);
            double r=-P.lambda2*(phi[j+1]-2*phi[j]+phi[j-1])*invh2 - (eneg-epos+C);
            mr=max(mr,abs(r));
            a[j]=-P.lambda2*invh2;
            c[j]=-P.lambda2*invh2;
            b[j]=2*P.lambda2*invh2 + eneg + epos;
            d[j]=-r;
        }
        if(mr<1e-13) break;
        // boundaries delta=0, solve interior Thomas
        for(int j=2;j<M-1;j++){
            double m=a[j]/b[j-1];
            b[j]-=m*c[j-1]; d[j]-=m*d[j-1];
        }
        vector<double> del(M,0.0);
        del[M-2]=d[M-2]/b[M-2];
        for(int j=M-3;j>=1;j--) del[j]=(d[j]-c[j]*del[j+1])/b[j];
        double alpha=1.0;
        // conservative damping
        double md=0; for(int j=1;j<M-1;j++) md=max(md,abs(del[j]));
        if(md>0.5) alpha=0.5/md;
        for(int j=1;j<M-1;j++) phi[j]+=alpha*del[j];
        phi[0]=phiL; phi[M-1]=phiR;
    }
    for(int j=0;j<M;j++){ n[j]=exp(phi[j]); p[j]=exp(-phi[j]); }
    return {phi,n,p,H};
}

static inline void eq_at(long long gi,long long N,const Params&P,const EqProfile&E,double&n,double&p,double&phi){
    long long c=N/2; long long off=gi-c;
    if(off < -E.half){ double C=-P.Nd; n=bulk_n(C); p=1.0/n; phi=log(n); }
    else if(off > E.half){ double C=P.Nd; n=bulk_n(C); p=1.0/n; phi=log(n); }
    else { int j=(int)(off+E.half); n=E.n[j]; p=E.p[j]; phi=E.phi[j]; }
}

// Exact discrete zero-Dirichlet Poisson perturbation on a global N grid, evaluated only on local contiguous support/window.
// Equation: -lambda2 * lap(phi) = rho. Interior indices 1..N-2.
void poisson_local_from_rho(const vector<double>& rho, vector<double>& phi, long long g0, long long N, const Params&P){
    int W=(int)rho.size();
    phi.assign(W,0.0);
    const long double L=(long double)(N-1);
    const long double scale=(long double)(P.dx*P.dx/P.lambda2);
    vector<long double> pref_j(W+1,0.0L), pref_Nj(W+1,0.0L);
    for(int t=0;t<W;t++){
        long long j=g0+t;
        long double f=(long double)rho[t]*scale;
        pref_j[t+1]=pref_j[t]+(long double)j*f;
        pref_Nj[t+1]=pref_Nj[t]+(L-(long double)j)*f;
    }
    long double total_Nj=pref_Nj[W];
    for(int t=0;t<W;t++){
        long long i=g0+t;
        // j <= i uses pref through t inclusive, j > i uses suffix
        long double left=pref_j[t+1];
        long double right=total_Nj-pref_Nj[t+1];
        long double v=((L-(long double)i)/L)*left + ((long double)i/L)*right;
        phi[t]=(double)v;
    }
