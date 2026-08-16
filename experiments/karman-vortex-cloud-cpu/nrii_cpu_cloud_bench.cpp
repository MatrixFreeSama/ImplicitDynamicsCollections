#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <vector>
#include <omp.h>

using Clock=std::chrono::steady_clock;
struct Params{
 int nx, steps=12; size_t nc, n;
 float dt=5e-5f, nu=2.5e-4f, beta2=1.f, sigma=.05f, tol=5e-6f, ub=1.45f;
 float invdx, invdx2; int p, cheb_tail_start, pade_cap;
 int gmres_restart, gmres_max, newton_max, ls_max;
 float armijo, ew_floor, ew_cap, cheb_tail_limit;
 int cx,cy,hw,core,hh;
};
static int next_pow2(int x){int p=1;while(p<x)p<<=1;return p;}
static int decimal_demand(float tol){int d=0;double x=tol; while(x<1.0){x*=10;d++;} return d?d:1;}
static int derive_p(const Params& b){
 double chi=(double)b.dt*(8.0*b.nu*b.invdx2 + 4.0*b.ub*b.invdx + 2.0*std::sqrt((double)b.beta2)*b.invdx + b.sigma);
 double target=b.tol/std::sqrt(3.0), e=std::exp(chi), term=1.0;
 for(int p=0;p<1000;p++){term*=chi/(p+1.0); if(e*std::abs(term)<=target) return p+1;} return -1;
}
static Params make_params(int nx){Params b; b.nx=nx;b.nc=(size_t)nx*nx;b.n=b.nc*3;b.invdx=nx-1.f;b.invdx2=b.invdx*b.invdx;b.p=derive_p(b);int root=1;while(root*root<b.p+1)root++;b.cheb_tail_start=(b.p+1)-root;b.pade_cap=b.p/2;b.cheb_tail_limit=std::sqrt(b.tol);int d=decimal_demand(b.tol);b.gmres_restart=next_pow2(d+d);b.gmres_max=b.gmres_restart*next_pow2(d);b.newton_max=next_pow2(d);b.ls_max=next_pow2(d);b.armijo=std::sqrt(1.1920928955078125e-7f);b.ew_floor=std::sqrt(b.tol);b.ew_cap=(float)d/(d+1.f);b.cx=nx*30/100;b.cy=nx/2;b.hw=nx*8/100;b.core=nx*5/100;b.hh=nx*3/100;return b;}

static void build_mask(const Params& b,std::vector<unsigned char>& m){m.assign(b.nc,0);
#pragma omp parallel for schedule(static)
 for(long long c=0;c<(long long)b.nc;c++){int x=c%b.nx,y=c/b.nx;if(x==0||y==0||x+1>=b.nx||y+1>=b.nx){m[c]=1;continue;}int dx=std::abs(x-b.cx),dy=std::abs(y-b.cy);if(dx>b.hw){m[c]=0;continue;}bool br=false;if(dx<=b.core)br=dy<=b.hh;else {int hl=b.hh*(b.hw-dx)/(b.hw-b.core);br=dy<=hl;}m[c]=br?2:0;}}
static void init_field(const Params& b,const std::vector<unsigned char>& mask,std::vector<float>& u){u.assign(b.n,0.f);
#pragma omp parallel for schedule(static)
 for(long long c=0;c<(long long)b.nc;c++){size_t i=(size_t)c*3;if(mask[c]==1){u[i]=1;continue;}if(mask[c]==2)continue;int x=c%b.nx,y=c/b.nx;float X=(float)x/(b.nx-1),Y=(float)y/(b.nx-1),uu=1,vv=0; const float vx[6]={.43f,.52f,.61f,.70f,.79f,.88f}, vy[6]={.46f,.54f,.46f,.54f,.46f,.54f}, amp[6]={6.6666665f,-6.6666665f,6.6666665f,-6.6666665f,5.8333335f,-5.f}; for(int k=0;k<6;k++){float dx=X-vx[k],dy=Y-vy[k],s=1.f-(dx*dx+dy*dy)/.0036f;if(s<0)s=0;s*=s;float a=s*amp[k];uu+=a*(-dy);vv+=a*dx;}u[i]=uu;u[i+1]=vv;}}

inline size_t idx(int x,int y,int comp,int nx){return ((size_t)y*nx+x)*3+comp;}
static void residual(const Params& b,const std::vector<unsigned char>& mask,const std::vector<float>& u,const std::vector<float>& old,std::vector<float>& r){r.resize(b.n);float h=.5f*b.invdx;
#pragma omp parallel for schedule(static)
 for(long long cc=0;cc<(long long)b.nc;cc++){size_t base=(size_t)cc*3;int x=cc%b.nx,y=cc/b.nx;unsigned char mk=mask[cc];if(mk){for(int c=0;c<3;c++){float target=(mk==1&&c==0)?1.f:0.f;r[base+c]=u[base+c]-target;}continue;}float uc=u[base],vc=u[base+1],pc=u[base+2];for(int comp=0;comp<2;comp++){float q=u[base+comp];float dx=(u[idx(x+1,y,comp,b.nx)]-u[idx(x-1,y,comp,b.nx)])*h;float dy=(u[idx(x,y+1,comp,b.nx)]-u[idx(x,y-1,comp,b.nx)])*h;float lap=(u[idx(x+1,y,comp,b.nx)]+u[idx(x-1,y,comp,b.nx)]+u[idx(x,y+1,comp,b.nx)]+u[idx(x,y-1,comp,b.nx)]-4*q)*b.invdx2;float pg=(comp==0?(u[idx(x+1,y,2,b.nx)]-u[idx(x-1,y,2,b.nx)])*h:(u[idx(x,y+1,2,b.nx)]-u[idx(x,y-1,2,b.nx)])*h);float F=-(uc*dx+vc*dy)-pg+b.nu*lap;r[base+comp]=q-old[base+comp]-b.dt*F;}float div=(u[idx(x+1,y,0,b.nx)]-u[idx(x-1,y,0,b.nx)])*h+(u[idx(x,y+1,1,b.nx)]-u[idx(x,y-1,1,b.nx)])*h;float Fp=-b.beta2*div-b.sigma*pc;r[base+2]=pc-old[base+2]-b.dt*Fp;}}
static double dot(const std::vector<float>& a,const std::vector<float>& b){double s=0;
#pragma omp parallel for reduction(+:s) schedule(static)
 for(long long i=0;i<(long long)a.size();i++)s+=(double)a[i]*b[i];return s;}
static double rms(const std::vector<float>& a){return std::sqrt(dot(a,a)/a.size());}
static double normv(const std::vector<float>& a){return std::sqrt(dot(a,a));}
static void axpy(std::vector<float>& y,const std::vector<float>& x,float a){
#pragma omp parallel for schedule(static)
 for(long long i=0;i<(long long)y.size();i++)y[i]+=a*x[i];}
static void scale_copy(std::vector<float>& y,const std::vector<float>& x,float a){y.resize(x.size());
#pragma omp parallel for schedule(static)
 for(long long i=0;i<(long long)x.size();i++)y[i]=a*x[i];}

// NRII coefficient planes: coeff[k*n+dof].
static void nrii_next_coeff(const Params& b,const std::vector<unsigned char>& mask,std::vector<float>& coeff,int k){size_t n=b.n;float h=.5f*b.invdx;float* out=coeff.data()+(size_t)(k+1)*n;
#pragma omp parallel for schedule(static)
 for(long long cc=0;cc<(long long)b.nc;cc++){size_t base=(size_t)cc*3;if(mask[cc]){out[base]=out[base+1]=out[base+2]=0;continue;}int x=cc%b.nx,y=cc/b.nx;const float* ck=coeff.data()+(size_t)k*n;for(int comp=0;comp<2;comp++){float q=ck[base+comp];float lap=(ck[idx(x+1,y,comp,b.nx)]+ck[idx(x-1,y,comp,b.nx)]+ck[idx(x,y+1,comp,b.nx)]+ck[idx(x,y-1,comp,b.nx)]-4*q)*b.invdx2;float pg=(comp==0?(ck[idx(x+1,y,2,b.nx)]-ck[idx(x-1,y,2,b.nx)])*h:(ck[idx(x,y+1,2,b.nx)]-ck[idx(x,y-1,2,b.nx)])*h);float adv=0;for(int a=0;a<=k;a++){int bb=k-a;const float* ca=coeff.data()+(size_t)a*n;const float* cb=coeff.data()+(size_t)bb*n;float ua=ca[base],va=ca[base+1];float dx=(cb[idx(x+1,y,comp,b.nx)]-cb[idx(x-1,y,comp,b.nx)])*h;float dy=(cb[idx(x,y+1,comp,b.nx)]-cb[idx(x,y-1,comp,b.nx)])*h;adv+=ua*dx+va*dy;}out[base+comp]=b.dt*(-adv-pg+b.nu*lap);}float div=(ck[idx(x+1,y,0,b.nx)]-ck[idx(x-1,y,0,b.nx)])*h+(ck[idx(x,y+1,1,b.nx)]-ck[idx(x,y-1,1,b.nx)])*h;out[base+2]=b.dt*(-b.beta2*div-b.sigma*ck[base+2]);}}

static double binom(int n,int k){if(k<0||k>n)return 0;if(k>n-k)k=n-k;double r=1;for(int i=1;i<=k;i++)r=r*(n-k+i)/i;return r;}
static std::vector<float> cheb_matrix(int p){std::vector<float> t((p+1)*(p+1));for(int j=0;j<=p;j++){double zj=std::ldexp(1.0,-j);for(int rr=0;rr<=j;rr++){double cx=binom(j,rr)*zj;if(rr==0){t[j]+=cx;continue;}for(int k=(rr&1);k<=rr;k+=2){double ck=(k==0?binom(rr,rr/2):2*binom(rr,(rr-k)/2))*std::ldexp(1.0,-rr);t[(size_t)k*(p+1)+j]+=cx*ck;}}}return t;}
static double cheb_tail_ratio(const Params& b,const std::vector<float>& coeff,const std::vector<float>& T){int p=b.p;size_t n=b.n;double tail2=0,total2=0;
#pragma omp parallel for reduction(+:tail2,total2) schedule(static)
 for(long long d=0;d<(long long)n;d++){double a2=0,tt=0;for(int k=0;k<=p;k++){double s=0;for(int j=0;j<=p;j++)s+=(double)T[(size_t)k*(p+1)+j]*coeff[(size_t)j*n+d];a2+=s*s;if(k>=b.cheb_tail_start)tt+=s*s;}total2+=a2;tail2+=tt;}return total2>0?std::sqrt(tail2/total2):0;}
static bool solve_linear(std::vector<double>A,std::vector<double>b,std::vector<double>&x,int n){x.assign(n,0);double norm=0;for(double v:A)norm=std::max(norm,std::abs(v));if(norm==0)return false;for(int k=0;k<n;k++){int piv=k;double best=0;for(int i=k;i<n;i++){double v=std::abs(A[(size_t)i*n+k]);if(v>best){best=v;piv=i;}}if(best<=std::numeric_limits<double>::epsilon()*n*norm)return false;if(piv!=k){for(int j=k;j<n;j++)std::swap(A[(size_t)k*n+j],A[(size_t)piv*n+j]);std::swap(b[k],b[piv]);}for(int i=k+1;i<n;i++){double f=A[(size_t)i*n+k]/A[(size_t)k*n+k];for(int j=k+1;j<n;j++)A[(size_t)i*n+j]-=f*A[(size_t)k*n+j];b[i]-=f*b[k];}}for(int i=n-1;i>=0;i--){double v=b[i];for(int j=i+1;j<n;j++)v-=A[(size_t)i*n+j]*x[j];x[i]=v/A[(size_t)i*n+i];}return true;}
static double comp_dot(const std::vector<float>& coeff,size_t n,int ka,int kb,int comp){double s=0;const float*a=coeff.data()+(size_t)ka*n,*b=coeff.data()+(size_t)kb*n;
#pragma omp parallel for reduction(+:s) schedule(static)
 for(long long i=comp;i<(long long)n;i+=3)s+=(double)a[i]*b[i];return s;}
static bool identify_pade(const Params& b,const std::vector<float>& coeff,std::vector<std::vector<double>>& qs,std::vector<int>& ms){qs.assign(3,{});ms.assign(3,0);for(int c=0;c<3;c++){bool ok=false;for(int m=b.p/2;m>0&&!ok;m--){int L=b.p-m;std::vector<double>A((size_t)m*m),rhs(m),x;for(int i=1;i<=m;i++){for(int j=1;j<=m;j++){double s=0;for(int k=L+1;k<=b.p;k++)s+=comp_dot(coeff,b.n,k-i,k-j,c);A[(size_t)(i-1)*m+j-1]=s;}double s=0;for(int k=L+1;k<=b.p;k++)s-=comp_dot(coeff,b.n,k-i,k,c);rhs[i-1]=s;}if(solve_linear(A,rhs,x,m)){// simple dense sampling interval pole gate for CPU bench
 bool safe=true;for(int z=0;z<=4096;z++){double zz=(double)z/4096,den=1,pw=zz;for(int j=0;j<m;j++){den+=x[j]*pw;pw*=zz;}if(std::abs(den)<std::sqrt(std::numeric_limits<float>::epsilon())*(1+std::accumulate(x.begin(),x.end(),0.0,[](double a,double v){return a+std::abs(v);}))){safe=false;break;}}if(safe){qs[c]=x;ms[c]=m;ok=true;}}}if(!ok)return false;}return true;}
static void pade_eval(const Params& b,const std::vector<float>& coeff,const std::vector<std::vector<double>>& qs,const std::vector<int>& ms,std::vector<float>& out){out.resize(b.n);
#pragma omp parallel for schedule(static)
 for(long long d=0;d<(long long)b.n;d++){int c=d%3,m=ms[c],L=b.p-m;double den=1;for(double q:qs[c])den+=q;double num=0;for(int k=0;k<=L;k++){double pk=0;for(int j=0;j<=m&&j<=k;j++){double q=(j==0?1.0:qs[c][j-1]);pk+=q*coeff[(size_t)(k-j)*b.n+d];}num+=pk;}out[d]=(float)(num/den);}}

struct NriiStats{double ms=0,res=0,tail=0;int p=0;std::vector<int> m;bool ok=false;};
static NriiStats run_nrii(const Params& b,const std::vector<unsigned char>& mask){NriiStats st;st.p=b.p;std::vector<float>u,old,coeff((size_t)(b.p+1)*b.n),tmp,r;init_field(b,mask,u);auto T=cheb_matrix(b.p);auto t0=Clock::now();for(int step=0;step<b.steps;step++){old=u;std::copy(u.begin(),u.end(),coeff.begin());for(int k=0;k<b.p;k++)nrii_next_coeff(b,mask,coeff,k);double tr=cheb_tail_ratio(b,coeff,T);st.tail=std::max(st.tail,tr);if(!(tr<=b.cheb_tail_limit)){st.ok=false;return st;}std::vector<std::vector<double>>qs;std::vector<int>ms;if(!identify_pade(b,coeff,qs,ms)){st.ok=false;return st;}pade_eval(b,coeff,qs,ms,tmp);u.swap(tmp);residual(b,mask,u,old,r);double rr=rms(r);if(!(rr<=b.tol)){st.res=rr;st.m=ms;st.ok=false;return st;}st.m=ms;st.res=rr;}auto t1=Clock::now();st.ms=std::chrono::duration<double,std::milli>(t1-t0).count();st.ok=true;return st;}

static void jv(const Params& b,const std::vector<unsigned char>& mask,const std::vector<float>& u,const std::vector<float>& v,std::vector<float>& o){o.resize(b.n);float h=.5f*b.invdx;
#pragma omp parallel for schedule(static)
 for(long long cc=0;cc<(long long)b.nc;cc++){size_t base=(size_t)cc*3;if(mask[cc]){o[base]=v[base];o[base+1]=v[base+1];o[base+2]=v[base+2];continue;}int x=cc%b.nx,y=cc/b.nx;float uc=u[base],vc=u[base+1],duc=v[base],dvc=v[base+1];for(int comp=0;comp<2;comp++){float du=v[base+comp];float gradx=(u[idx(x+1,y,comp,b.nx)]-u[idx(x-1,y,comp,b.nx)])*h;float grady=(u[idx(x,y+1,comp,b.nx)]-u[idx(x,y-1,comp,b.nx)])*h;float dgradx=(v[idx(x+1,y,comp,b.nx)]-v[idx(x-1,y,comp,b.nx)])*h;float dgrady=(v[idx(x,y+1,comp,b.nx)]-v[idx(x,y-1,comp,b.nx)])*h;float dlap=(v[idx(x+1,y,comp,b.nx)]+v[idx(x-1,y,comp,b.nx)]+v[idx(x,y+1,comp,b.nx)]+v[idx(x,y-1,comp,b.nx)]-4*du)*b.invdx2;float dpg=(comp==0?(v[idx(x+1,y,2,b.nx)]-v[idx(x-1,y,2,b.nx)])*h:(v[idx(x,y+1,2,b.nx)]-v[idx(x,y-1,2,b.nx)])*h);float dF=-(duc*gradx+uc*dgradx+dvc*grady+vc*dgrady)-dpg+b.nu*dlap;o[base+comp]=du-b.dt*dF;}float ddiv=(v[idx(x+1,y,0,b.nx)]-v[idx(x-1,y,0,b.nx)])*h+(v[idx(x,y+1,1,b.nx)]-v[idx(x,y-1,1,b.nx)])*h;float dFp=-b.beta2*ddiv-b.sigma*v[base+2];o[base+2]=v[base+2]-b.dt*dFp;}}
static void precond(const Params& b,const std::vector<unsigned char>& mask,const std::vector<float>& u,const std::vector<float>& x,std::vector<float>& o){o.resize(b.n);float h=.5f*b.invdx;
#pragma omp parallel for schedule(static)
 for(long long cc=0;cc<(long long)b.nc;cc++){size_t base=(size_t)cc*3;if(mask[cc]){o[base]=x[base];o[base+1]=x[base+1];o[base+2]=x[base+2];continue;}int xx=cc%b.nx,yy=cc/b.nx;float dudx=(u[idx(xx+1,yy,0,b.nx)]-u[idx(xx-1,yy,0,b.nx)])*h;float dvdy=(u[idx(xx,yy+1,1,b.nx)]-u[idx(xx,yy-1,1,b.nx)])*h;float du=std::max(std::abs(1+b.dt*(dudx+4*b.nu*b.invdx2)),.25f);float dv=std::max(std::abs(1+b.dt*(dvdy+4*b.nu*b.invdx2)),.25f);o[base]=x[base]/du;o[base+1]=x[base+1]/dv;o[base+2]=x[base+2]/(1+b.dt*b.sigma);}}

static int gmres(const Params& b,const std::vector<unsigned char>& mask,const std::vector<float>& u,const std::vector<float>& rhs,float eta,std::vector<float>& delta){int R=b.gmres_restart,total=0;std::vector<std::vector<float>> V(R+1,std::vector<float>(b.n));std::vector<float>pre,w,trial;delta.assign(b.n,0);precond(b,mask,u,rhs,pre);double beta=normv(pre),bnorm=beta;if(beta<=b.ew_floor)return 0;while(total<b.gmres_max){std::vector<double>H((size_t)(R+1)*R),cs(R),sn(R),g(R+1),y(R);g[0]=beta;scale_copy(V[0],pre,(float)(1/beta));int kdim=0;bool conv=false;for(int j=0;j<R&&total<b.gmres_max;j++){jv(b,mask,u,V[j],w);precond(b,mask,u,w,pre);for(int i=0;i<=j;i++){double hij=dot(pre,V[i]);H[(size_t)i*R+j]=hij;axpy(pre,V[i],(float)-hij);}double hn=normv(pre);H[(size_t)(j+1)*R+j]=hn;if(hn>b.ew_floor)scale_copy(V[j+1],pre,(float)(1/hn));for(int i=0;i<j;i++){double temp=cs[i]*H[(size_t)i*R+j]+sn[i]*H[(size_t)(i+1)*R+j];H[(size_t)(i+1)*R+j]=-sn[i]*H[(size_t)i*R+j]+cs[i]*H[(size_t)(i+1)*R+j];H[(size_t)i*R+j]=temp;}double rho=std::hypot(H[(size_t)j*R+j],H[(size_t)(j+1)*R+j]);if(rho<=b.ew_floor)return -1;cs[j]=H[(size_t)j*R+j]/rho;sn[j]=H[(size_t)(j+1)*R+j]/rho;H[(size_t)j*R+j]=rho;H[(size_t)(j+1)*R+j]=0;double temp=cs[j]*g[j]+sn[j]*g[j+1];g[j+1]=-sn[j]*g[j]+cs[j]*g[j+1];g[j]=temp;total++;kdim=j+1;if(std::abs(g[j+1])<=eta*bnorm||hn<=b.ew_floor){conv=true;break;}}for(int i=kdim-1;i>=0;i--){double acc=g[i];for(int j=i+1;j<kdim;j++)acc-=H[(size_t)i*R+j]*y[j];if(std::abs(H[(size_t)i*R+i])<=b.ew_floor)return -2;y[i]=acc/H[(size_t)i*R+i];}for(int i=0;i<kdim;i++)axpy(delta,V[i],(float)y[i]);if(conv)return total;jv(b,mask,u,delta,w);trial=rhs;axpy(trial,w,-1);precond(b,mask,u,trial,pre);beta=normv(pre);if(beta<=eta*bnorm)return total;}return -3;}
struct NkStats{double ms=0,res=0;int newton=0,gmres=0,ls=0;bool ok=false;};
static NkStats run_nk(const Params& b,const std::vector<unsigned char>& mask){NkStats st;std::vector<float>u,old,r,rhs,delta,tmp,tr;init_field(b,mask,u);auto t0=Clock::now();for(int step=0;step<b.steps;step++){old=u;double last=0;for(int ni=0;ni<b.newton_max;ni++){residual(b,mask,u,old,r);double rr=rms(r);if(rr<b.tol)break;float eta;if(ni==0||last<=0)eta=1.f/(decimal_demand(b.tol)+1.f);else{float ratio=(float)(rr/(last+1.1920928955078125e-7));eta=ratio*std::sqrt(ratio);eta=std::max(b.ew_floor,std::min(eta,b.ew_cap));}scale_copy(rhs,r,-1);int gi=gmres(b,mask,u,rhs,eta,delta);if(gi<0){st.ok=false;return st;}st.gmres+=gi;float alpha=1;bool acc=false;for(int ls=0;ls<b.ls_max;ls++){tmp=u;axpy(tmp,delta,alpha);residual(b,mask,tmp,old,tr);double trr=rms(tr);st.ls++;if(trr<b.tol||trr<=(1-b.armijo*alpha)*rr){u.swap(tmp);acc=true;break;}alpha*=.5f;}if(!acc){st.ok=false;return st;}last=rr;st.newton++;}residual(b,mask,u,old,r);double rr=rms(r);if(rr>b.tol){st.res=rr;st.ok=false;return st;}st.res=rr;}auto t1=Clock::now();st.ms=std::chrono::duration<double,std::milli>(t1-t0).count();st.ok=true;return st;}

int main(int argc,char**argv){int threads=omp_get_max_threads();std::cout<<"threads,"<<threads<<"\n";std::cout<<"grid,dof,p,nrii_ok,nrii_ms,nrii_res,nrii_tail,pade_m_u,pade_m_v,pade_m_p,nk_ok,nk_ms,nk_res,newton_iters,gmres_iters,line_search,speedup_nk_over_nrii\n";std::vector<int> grids={64,128,256};if(argc>1){grids.clear();for(int i=1;i<argc;i++)grids.push_back(std::atoi(argv[i]));}for(int nx:grids){Params b=make_params(nx);std::vector<unsigned char>mask;build_mask(b,mask); // warm operator via residual
 std::vector<float>uu,rr;init_field(b,mask,uu);residual(b,mask,uu,uu,rr);
 NriiStats ns=run_nrii(b,mask);NkStats ks=run_nk(b,mask);std::cout<<nx<<","<<b.n<<","<<b.p<<","<<ns.ok<<","<<ns.ms<<","<<ns.res<<","<<ns.tail<<","<<(ns.m.size()>0?ns.m[0]:0)<<","<<(ns.m.size()>1?ns.m[1]:0)<<","<<(ns.m.size()>2?ns.m[2]:0)<<","<<ks.ok<<","<<ks.ms<<","<<ks.res<<","<<ks.newton<<","<<ks.gmres<<","<<ks.ls<<","<<(ns.ms>0?ks.ms/ns.ms:0)<<"\n";std::cout.flush();}
}
