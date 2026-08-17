#include <cmath>
#include <cstdio>
#include <vector>
#include <string>
#include <algorithm>
#include <limits>
#include <fstream>
#include <iomanip>

static double F(double x,double mu){ return x + mu + std::exp(-x); }
static double Fx(double x){ return 1.0 - std::exp(-x); }
static double mu_curve(double x){ return -x - std::exp(-x); }

struct NewtonResult{ bool ok=false; double x=0,res=0; int iters=0; };
static NewtonResult newton_fixed(double x,double mu,double tol=1e-12,int maxit=50){
    NewtonResult out;
    for(int k=0;k<maxit;k++){
        double f=F(x,mu), j=Fx(x);
        if(std::abs(f)<tol){ out={true,x,std::abs(f),k}; return out; }
        if(!std::isfinite(j)||std::abs(j)<1e-14){ out={false,x,std::abs(f),k}; return out; }
        double step=-f/j, alpha=1.0, base=std::abs(f); bool acc=false;
        for(int ls=0;ls<60;ls++){
            double xn=x+alpha*step, fn=F(xn,mu);
            if(std::isfinite(fn) && std::abs(fn) <= (1.0-1e-4*alpha)*base){ x=xn; acc=true; break; }
            alpha*=0.5;
        }
        if(!acc){ out={false,x,std::abs(F(x,mu)),k}; return out; }
    }
    out={std::abs(F(x,mu))<tol,x,std::abs(F(x,mu)),maxit}; return out;
}

static bool solve_dense(std::vector<double> A,std::vector<double> b,int n,std::vector<double>& x){
    for(int k=0;k<n;k++){
        int piv=k; double best=std::abs(A[k*n+k]);
        for(int i=k+1;i<n;i++){ double v=std::abs(A[i*n+k]); if(v>best){best=v;piv=i;} }
        if(best<1e-18 || !std::isfinite(best)) return false;
        if(piv!=k){ for(int j=k;j<n;j++) std::swap(A[k*n+j],A[piv*n+j]); std::swap(b[k],b[piv]); }
        double d=A[k*n+k];
        for(int i=k+1;i<n;i++){
            double q=A[i*n+k]/d;
            A[i*n+k]=0;
            for(int j=k+1;j<n;j++) A[i*n+j]-=q*A[k*n+j];
            b[i]-=q*b[k];
        }
    }
    x.assign(n,0);
    for(int i=n-1;i>=0;i--){
        double s=b[i]; for(int j=i+1;j<n;j++) s-=A[i*n+j]*x[j];
        x[i]=s/A[i*n+i];
    }
    return true;
}

struct PadeResult{ bool ok=false; double value=0,den=0; };
static std::vector<double> nrii_coeffs(double anchor,double mu,int p){
    std::vector<double> c(p+1,0), e(p+1,0);
    c[0]=anchor; e[0]=std::exp(-anchor);
    if(p>=1) c[1]=-mu-e[0]-anchor;
    for(int n=1;n<=p;n++){
        double s=0; for(int k=1;k<=n;k++) s += k*c[k]*e[n-k];
        e[n]=-s/n;
        if(n+1<=p) c[n+1]=-e[n];
    }
    return c;
}
static PadeResult nrii_pade(double anchor,double mu,int order,double theta=1.0){
    auto c=nrii_coeffs(anchor,mu,order);
    int n=order/2, m=order-n;
    std::vector<double> A(n*n,0), b(n,0), qrest;
    for(int r=0;r<n;r++){
        int k=m+1+r;
        for(int j=1;j<=n;j++){ int idx=k-j; if(idx>=0 && idx<(int)c.size()) A[r*n+j-1]=c[idx]; }
        b[r]=-c[k];
    }
    if(n>0 && !solve_dense(A,b,n,qrest)) return {};
    std::vector<double> q(n+1,0), pc(m+1,0); q[0]=1; for(int j=1;j<=n;j++) q[j]=qrest[j-1];
    for(int k=0;k<=m;k++) for(int j=0;j<=std::min(n,k);j++) pc[k]+=q[j]*c[k-j];
    auto eval=[&](const std::vector<double>& a){ double s=0,p=1; for(double v:a){s+=v*p;p*=theta;} return s; };
    double den=eval(q); if(!std::isfinite(den)||std::abs(den)<1e-12) return {};
    double val=eval(pc)/den; if(!std::isfinite(val)) return {};
    return {true,val,den};
}

struct HandoffInterval{ bool ok=false; double center=0,radius=0,min_den=0; };
static HandoffInterval make_handoff_interval(double anchor,double mu){
    const int orders[3]={8,12,16}; std::vector<double> vals; double mind=1e300;
    for(int p:orders){ auto r=nrii_pade(anchor,mu,p); if(r.ok && std::abs(r.den)>1e-10){ vals.push_back(r.value); mind=std::min(mind,std::abs(r.den)); } }
    if(vals.size()<2) return {};
    double center=vals.back(), spread=0; for(size_t i=0;i+1<vals.size();i++) spread=std::max(spread,std::abs(vals[i]-center));
    return {true,center,std::max(1e-10,4.0*spread),mind};
}

struct Kraw{ bool ok=false; double lo=0,hi=0; };
static Kraw krawczyk(double c,double r,double mu){
    double a=c-r,b=c+r,jc=Fx(c); if(std::abs(jc)<1e-15) return {};
    double C=1.0/jc, dlo=Fx(a), dhi=Fx(b); if(dlo>dhi) std::swap(dlo,dhi);
    double A1=1-C*dlo,A2=1-C*dhi,Alo=std::min(A1,A2),Ahi=std::max(A1,A2);
    double prods[4]={Alo*(-r),Alo*r,Ahi*(-r),Ahi*r};
    double mlo=*std::min_element(prods,prods+4), mhi=*std::max_element(prods,prods+4);
    double base=c-C*F(c,mu), lo=base+mlo, hi=base+mhi;
    return {lo>a && hi<b,lo,hi};
}

struct ArcResult{ bool ok=false; double x=0,mu=0,tx=0,tmu=0; int iters=0; };
static void tangent(double x,double& tx,double& tmu){
    tx=1.0; tmu=-Fx(x); double n=std::hypot(tx,tmu); tx/=n;tmu/=n;
}
static ArcResult arc_step(double x,double mu,double tx,double tmu,double ds){
    double px=x+ds*tx, pm=mu+ds*tmu, qx=px,qm=pm;
    for(int k=0;k<15;k++){
        double g1=F(qx,qm), g2=tx*(qx-px)+tmu*(qm-pm);
        if(std::max(std::abs(g1),std::abs(g2))<1e-13){
            double ntx,ntm; tangent(qx,ntx,ntm); if(ntx*tx+ntm*tmu<0){ntx=-ntx;ntm=-ntm;}
            return {true,qx,qm,ntx,ntm,k};
        }
        double a=Fx(qx),b=1,c=tx,d=tmu,det=a*d-b*c; if(std::abs(det)<1e-15) return {};
        double dx=(-g1*d + b*g2)/det;
        double dm=(-a*g2 + c*g1)/det;
        qx+=dx; qm+=dm;
    }
    return {};
}

int main(int argc,char**argv){
    std::string outdir = argc>1?argv[1]:".";
    std::ofstream relay(outdir+"/RELAY_TRACE.csv");
    relay<<"step,method,x0,mu0,x,mu,nrii_center,radius,certified,newton_iters,arc_iters,residual\n";
    double x=-1.0,mu=mu_curve(x),tx,tmu; tangent(x,tx,tmu); const double ds=0.08; int step=0;
    while(x<1.0-1e-9 && step<100){
        double x0=x,mu0=mu, mup=mu+ds*tmu; bool used=false;
        if(mup < -1.0005){
            HandoffInterval hd=make_handoff_interval(x,mup);
            if(hd.ok){ Kraw kr=krawczyk(hd.center,hd.radius,mup); if(kr.ok){ auto nr=newton_fixed(hd.center,mup); if(nr.ok){
                x=nr.x; mu=mup; double ntx,ntm;tangent(x,ntx,ntm); if(ntx*tx+ntm*tmu<0){ntx=-ntx;ntm=-ntm;} tx=ntx;tmu=ntm;
                relay<<step<<",NRII_Krawczyk_Newton,"<<std::setprecision(17)<<x0<<","<<mu0<<","<<x<<","<<mu<<","<<hd.center<<","<<hd.radius<<",1,"<<nr.iters<<",0,"<<std::abs(F(x,mu))<<"\n"; used=true;
            } } }
        }
        if(!used){
            auto ar=arc_step(x,mu,tx,tmu,ds); if(!ar.ok) break;
            x=ar.x;mu=ar.mu;tx=ar.tx;tmu=ar.tmu; Kraw kr=krawczyk(x,1e-6,mu);
            relay<<step<<",PseudoArclength,"<<std::setprecision(17)<<x0<<","<<mu0<<","<<x<<","<<mu<<",nan,1e-6,"<<(kr.ok?1:0)<<",0,"<<ar.iters<<","<<std::abs(F(x,mu))<<"\n";
        }
        step++;
    }
    relay.close();

    std::ofstream base(outdir+"/BASELINES.csv");
    base<<"case,ok,x,residual,iters,note\n";
    auto direct=newton_fixed(-1.0,mu); base<<"Direct_Newton_FinalMu,"<<direct.ok<<","<<std::setprecision(17)<<direct.x<<","<<direct.res<<","<<direct.iters<<",converges_to_wrong_branch_if_x_negative\n";
    auto fold=newton_fixed(-1.0,-1.0,1e-14,100); base<<"Newton_Exact_Double_Root,"<<fold.ok<<","<<fold.x<<","<<fold.res<<","<<fold.iters<<",linearized_convergence_at_multiplicity_two\n";
    auto hd=make_handoff_interval(-1.0,mu); if(hd.ok){ auto kr=krawczyk(hd.center,hd.radius,mu); base<<"NRII_Handoff_FinalMu,"<<kr.ok<<","<<hd.center<<","<<std::abs(F(hd.center,mu))<<",0,"<<(kr.ok?"certified_but_lower_branch":"not_certified")<<"\n"; }
    base.close();

    std::ofstream pcont(outdir+"/NAIVE_PARAMETER_CONTINUATION.csv");
    pcont<<"step,mu,xstart,x,ok,iters,residual\n";
    std::ifstream rin(outdir+"/RELAY_TRACE.csv"); std::string line; std::getline(rin,line); double pcx=-1.0; int idx=0;
    while(std::getline(rin,line)){
        std::vector<std::string> cols; size_t s=0;
        while(true){size_t p=line.find(',',s); if(p==std::string::npos){cols.push_back(line.substr(s));break;} cols.push_back(line.substr(s,p-s));s=p+1;}
        if(cols.size()<6) continue; double pmu=std::stod(cols[5]); auto nr=newton_fixed(pcx,pmu); double xs=pcx; if(nr.ok) pcx=nr.x;
        pcont<<idx++<<","<<std::setprecision(17)<<pmu<<","<<xs<<","<<nr.x<<","<<nr.ok<<","<<nr.iters<<","<<nr.res<<"\n";
    }
    pcont.close();

    std::ofstream sum(outdir+"/SUMMARY.csv");
    sum<<"metric,value\n";
    sum<<"relay_final_x,"<<std::setprecision(17)<<x<<"\n";
    sum<<"relay_final_mu,"<<mu<<"\n";
    sum<<"relay_final_residual,"<<std::abs(F(x,mu))<<"\n";
    sum<<"relay_steps,"<<step<<"\n";
    sum<<"fold_mu,-1\nfold_x,0\n";
    sum.close();
    return 0;
}
