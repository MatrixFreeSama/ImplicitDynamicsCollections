#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
using namespace std;

struct NewtonRes { bool ok; double P; int it; int line_search_steps; double min_alpha; double res; string reason; };
struct PadeRes { bool ok; double value; double min_den; int p,m; double res; };

static double residual(double P,double P0,double E,double dt){
    return P-P0-dt*(P-P*P*P+E);
}

static NewtonRes newton_armijo(double P0,double E,double dt,int maxit=50){
    double P=P0, min_alpha=1.0; int total_ls=0;
    for(int it=0; it<maxit; ++it){
        double r=residual(P,P0,E,dt), nr=fabs(r);
        if(nr<1e-12) return {true,P,it,total_ls,min_alpha,nr,"pass"};
        double J=1.0-dt*(1.0-3.0*P*P);
        if(!isfinite(J) || fabs(J)<1e-14)
            return {false,P,it,total_ls,min_alpha,nr,"singular_jacobian_at_iterate"};
        double d=-r/J, alpha=1.0; bool accepted=false;
        for(int ls=0; ls<80; ++ls){
            double Pn=P+alpha*d;
            double rn=fabs(residual(Pn,P0,E,dt));
            ++total_ls;
            if(isfinite(rn) && rn <= (1.0-1e-4*alpha)*nr){
                P=Pn; min_alpha=min(min_alpha,alpha); accepted=true; break;
            }
            alpha*=0.5;
        }
        if(!accepted) return {false,P,it,total_ls,min_alpha,nr,"armijo_failure"};
    }
    double nr=fabs(residual(P,P0,E,dt));
    return {nr<1e-12,P,maxit,total_ls,min_alpha,nr,"maxit"};
}

static vector<double> nrii_coeffs(double P0,double E,double dt,int pmax){
    // Formal coefficient closure of P(h)=P0+h*dt*(P(h)-P(h)^3+E)
    vector<double> c(pmax+1,0.0); c[0]=P0;
    for(int k=0;k<pmax;++k){
        long double cubic=0.0L;
        for(int a=0;a<=k;++a) for(int b=0;b<=k-a;++b){
            int cc=k-a-b;
            cubic += (long double)c[a]*c[b]*c[cc];
        }
        long double rhs=(long double)c[k]-cubic + (k==0 ? (long double)E : 0.0L);
        c[k+1]=(double)((long double)dt*rhs);
    }
    return c;
}

static bool solve_dense(vector<vector<long double>> A, vector<long double> b, vector<long double>& x){
    const int n=(int)b.size(); x.assign(n,0.0L);
    for(int i=0;i<n;++i){
        int piv=i; long double best=fabsl(A[i][i]);
        for(int r=i+1;r<n;++r){ long double v=fabsl(A[r][i]); if(v>best){best=v;piv=r;} }
        if(best<1e-30L) return false;
        if(piv!=i){ swap(A[piv],A[i]); swap(b[piv],b[i]); }
        long double d=A[i][i];
        for(int j=i;j<n;++j) A[i][j]/=d; b[i]/=d;
        for(int r=0;r<n;++r) if(r!=i){
            long double f=A[r][i]; if(f==0.0L) continue;
            for(int j=i;j<n;++j) A[r][j]-=f*A[i][j];
            b[r]-=f*b[i];
        }
    }
    x=b; return true;
}

static PadeRes pade_eval(const vector<double>& c,int p,int m,double P0,double E,double dt){
    const int L=p-m;
    vector<vector<long double>> A(m,vector<long double>(m));
    vector<long double> b(m),q;
    for(int rr=0;rr<m;++rr){
        int k=L+1+rr;
        b[rr]=-(long double)c[k];
        for(int j=1;j<=m;++j) A[rr][j-1]=(long double)c[k-j];
    }
    if(!solve_dense(A,b,q)) return {false,NAN,0.0,p,m,INFINITY};
    long double denom=1.0L;
    for(auto z:q) denom+=z;
    if(fabsl(denom)<1e-24L) return {false,NAN,0.0,p,m,INFINITY};
    long double numerator=0.0L;
    for(int k=0;k<=L;++k){
        long double ak=c[k];
        for(int j=1;j<=min(k,m);++j) ak += q[j-1]*(long double)c[k-j];
        numerator += ak;
    }
    long double min_den=1e300L;
    for(int s=0;s<=4096;++s){
        long double z=(long double)s/4096.0L, zp=z, d=1.0L;
        for(int j=1;j<=m;++j){ d += q[j-1]*zp; zp*=z; }
        min_den=min(min_den,fabsl(d));
    }
    if(min_den<1e-10L) return {false,NAN,(double)min_den,p,m,INFINITY};
    double P=(double)(numerator/denom);
    double r=fabs(residual(P,P0,E,dt));
    return {isfinite(P)&&isfinite(r),P,(double)min_den,p,m,r};
}

static double cheb_tail_ratio(const vector<double>& c,int p){
    // Chebyshev audit only. Polynomial is sampled on Gauss-Chebyshev nodes over h in [0,1].
    int n=p+1; vector<long double> v(n),a(n);
    const long double PI=acosl(-1.0L);
    for(int j=0;j<n;++j){
        long double x=cosl(PI*(j+0.5L)/n), h=(x+1.0L)/2.0L, y=0.0L;
        for(int k=p;k>=0;--k) y=y*h+c[k];
        v[j]=y;
    }
    for(int k=0;k<n;++k){
        long double s=0.0L;
        for(int j=0;j<n;++j) s += v[j]*cosl(PI*k*(j+0.5L)/n);
        a[k]=2.0L*s/n; if(k==0) a[k]/=2.0L;
    }
    int tail_start=max(0,n-(int)ceill(sqrt((long double)n)));
    long double tail=0.0L,total=0.0L;
    for(int k=0;k<n;++k){ total+=a[k]*a[k]; if(k>=tail_start) tail+=a[k]*a[k]; }
    return total>0.0L ? sqrt((double)(tail/total)) : 0.0;
}

static double continuation_reference(double P0,double E,double dt,double& min_abs_J){
    // Independent h-continuation reference. This is not used by NRII or Newton.
    double P=P0; min_abs_J=numeric_limits<double>::infinity();
    const int steps=100000;
    for(int s=1;s<=steps;++s){
        double h=(double)s/steps;
        for(int it=0;it<20;++it){
            double r=P-P0-h*dt*(P-P*P*P+E);
            double J=1.0-h*dt*(1.0-3.0*P*P);
            min_abs_J=min(min_abs_J,fabs(J));
            if(fabs(r)<1e-14) break;
            double d=-r/J, alpha=1.0, nr=fabs(r); bool accepted=false;
            for(int ls=0;ls<24;++ls){
                double Pn=P+alpha*d;
                double rn=fabs(Pn-P0-h*dt*(Pn-Pn*Pn*Pn+E));
                if(rn<nr){P=Pn;accepted=true;break;}
                alpha*=0.5;
            }
            if(!accepted) break;
        }
    }
    return P;
}

static vector<double> final_real_roots(double P0,double E,double dt){
    // dt*P^3 + (1-dt)*P - (P0+dt*E)=0; no quadratic term.
    long double pp=(1.0L-dt)/dt;
    long double qq=-(P0+dt*E)/dt;
    long double disc=qq*qq/4.0L + pp*pp*pp/27.0L;
    vector<double> roots;
    if(disc>0){
        long double s=sqrtl(disc);
        long double u=cbrtl(-qq/2.0L+s), v=cbrtl(-qq/2.0L-s);
        roots.push_back((double)(u+v));
    }else{
        long double rho=2.0L*sqrtl(max((long double)0.0L,-pp/3.0L));
        if(rho==0){ roots.push_back(0.0); }
        else{
            long double arg=(-qq/2.0L)/sqrtl(max((long double)1e-300L,-pp*pp*pp/27.0L));
            arg=max((long double)-1.0L,min((long double)1.0L,arg));
            long double th=acosl(arg);
            const long double PI=acosl(-1.0L);
            for(int k=0;k<3;++k) roots.push_back((double)(rho*cosl((th+2.0L*PI*k)/3.0L)));
            sort(roots.begin(),roots.end());
        }
    }
    return roots;
}

int main(){
    cout<<setprecision(17);
    const double P0=-0.3, E=0.2;
    const double dtcrit=1.0/(1.0-3.0*P0*P0);
    const double tol=1e-8;
    const vector<double> factors={0.80,0.90,0.95,0.98,0.99,0.995,0.999,1.0,1.001,1.005,1.01,1.02};

    ofstream sw("LK_ROOT_BASIN_SWEEP.csv"); sw<<setprecision(17);
    sw<<"factor,dt,newton_initial_jacobian,branch_min_abs_jacobian,newton_ok,newton_iterations,newton_line_search_trials,newton_min_alpha,newton_residual,newton_P,newton_reason,nrii_ok,nrii_declared_p,nrii_pade_m,nrii_residual,nrii_P,nrii_pade_min_den,cheb_tail_ratio,branch_reference_P,branch_error_abs,root_count,root0,root1,root2,newton_matches_branch,nrii_matches_branch\n";

    ofstream det("LK_CRITICAL_POINT_P_SCAN.csv"); det<<setprecision(17);
    det<<"p,best_m,best_residual,best_P,min_den,cheb_tail_ratio,passes_1e8\n";

    for(double factor:factors){
        double dt=dtcrit*factor;
        auto nr=newton_armijo(P0,E,dt);
        auto c=nrii_coeffs(P0,E,dt,160);
        PadeRes chosen{false,NAN,0,0,0,INFINITY}; double chosen_tail=NAN;
        for(int p=8;p<=160;++p){
            PadeRes local{false,NAN,0,p,0,INFINITY};
            for(int m=min(p/2,32);m>=1;--m){
                auto pr=pade_eval(c,p,m,P0,E,dt);
                if(pr.ok && pr.res<local.res) local=pr;
            }
            if(factor==1.0 && local.ok)
                det<<p<<","<<local.m<<","<<local.res<<","<<local.value<<","<<local.min_den<<","<<cheb_tail_ratio(c,p)<<","<<(local.res<tol?1:0)<<"\n";
            if(local.ok && local.res<tol){ chosen=local; chosen_tail=cheb_tail_ratio(c,p); break; }
            if(local.ok && (!chosen.ok || local.res<chosen.res)){ chosen=local; chosen_tail=cheb_tail_ratio(c,p); }
        }
        double minJ; double ref=continuation_reference(P0,E,dt,minJ);
        auto roots=final_real_roots(P0,E,dt);
        double r0=NAN,r1=NAN,r2=NAN; if(roots.size()>0)r0=roots[0]; if(roots.size()>1)r1=roots[1]; if(roots.size()>2)r2=roots[2];
        bool nmatch=nr.ok && fabs(nr.P-ref)<1e-6;
        bool pmatch=chosen.ok && fabs(chosen.value-ref)<1e-6;
        sw<<factor<<","<<dt<<","<<(1.0-dt*(1.0-3.0*P0*P0))<<","<<minJ<<","<<(nr.ok?1:0)<<","<<nr.it<<","<<nr.line_search_steps<<","<<nr.min_alpha<<","<<nr.res<<","<<nr.P<<","<<nr.reason<<","<<(chosen.ok&&chosen.res<tol?1:0)<<","<<chosen.p<<","<<chosen.m<<","<<chosen.res<<","<<chosen.value<<","<<chosen.min_den<<","<<chosen_tail<<","<<ref<<","<<fabs(chosen.value-ref)<<","<<roots.size()<<","<<r0<<","<<r1<<","<<r2<<","<<(nmatch?1:0)<<","<<(pmatch?1:0)<<"\n";
    }

    ofstream meta("MODEL_PARAMETERS.csv"); meta<<setprecision(17);
    meta<<"parameter,value\n";
    meta<<"model,dimensionless homogeneous Landau-Khalatnikov second-order ferroelectric constitutive switching\n";
    meta<<"free_energy,F(P)=-0.5*P^2+0.25*P^4-E*P\n";
    meta<<"dynamics,dP/dt=P-P^3+E\n";
    meta<<"implicit_family,P(h)=P0+h*dt*(P(h)-P(h)^3+E)\n";
    meta<<"P0,"<<P0<<"\nE,"<<E<<"\n";
    meta<<"critical_dt,"<<dtcrit<<"\nresidual_tolerance,"<<tol<<"\n";
    meta<<"nrii_max_declared_p,160\n";
    meta<<"pade_max_m,32\n";
    meta<<"pade_denominator_sample_count,4097\n";
    meta<<"newton,analytic_jacobian_plus_Armijo_from_P0_no_homotopy\n";
    meta<<"reference,independent_100000_step_h_continuation_not_used_by_solvers\n";
    return 0;
}
