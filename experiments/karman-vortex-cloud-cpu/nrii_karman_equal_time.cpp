#define main legacy_bench_main
#include "nrii_cpu_cloud_bench.cpp"
#undef main

#include <iomanip>
#include <map>

struct Sample { double ms; double res; int aux1; int aux2; bool ok; };

static Sample nrii_one_order(int nx, int p) {
    Params b = make_params(nx);
    b.steps = 1;
    b.p = p;
    int root=1; while(root*root < b.p+1) root++;
    b.cheb_tail_start=(b.p+1)-root;
    b.pade_cap=b.p/2;
    b.cheb_tail_limit=1e30f; // diagnostic only, do not reject by old tolerance gate
    std::vector<unsigned char> mask; build_mask(b,mask);
    std::vector<float> u, old, coeff((size_t)(p+1)*b.n), out, r;
    init_field(b,mask,u); old=u; std::copy(u.begin(),u.end(),coeff.begin());
    auto t0=Clock::now();
    for(int k=0;k<p;k++) nrii_next_coeff(b,mask,coeff,k);
    auto T=cheb_matrix(p);
    volatile double tail=cheb_tail_ratio(b,coeff,T); (void)tail;
    std::vector<std::vector<double>> qs; std::vector<int> ms;
    bool ok=identify_pade(b,coeff,qs,ms);
    if(ok) pade_eval(b,coeff,qs,ms,out);
    auto t1=Clock::now();
    double elapsed=std::chrono::duration<double,std::milli>(t1-t0).count();
    if(!ok) return {elapsed, NAN,0,0,false};
    residual(b,mask,out,old,r);
    double rr=rms(r);
    int mm=ms.size()>=3 ? std::max(ms[0],std::max(ms[1],ms[2])) : 0;
    return {elapsed,rr,p,mm,true};
}

static std::vector<Sample> newton_trace_one_step(int nx, int max_newton=10) {
    Params b = make_params(nx);
    b.steps=1;
    // Floor-seeking traditional configuration, same residual/Jv/preconditioner family.
    b.tol=1e-10f;
    b.gmres_restart=32;
    b.gmres_max=256;
    b.newton_max=max_newton;
    b.ls_max=10;
    b.ew_floor=1e-6f;
    b.ew_cap=.9f;
    b.armijo=std::sqrt(std::numeric_limits<float>::epsilon());

    std::vector<unsigned char> mask; build_mask(b,mask);
    std::vector<float> u,old,r,rhs,delta,tmp,tr;
    init_field(b,mask,u); old=u;
    std::vector<Sample> out;
    auto t0=Clock::now();
    residual(b,mask,u,old,r);
    double rr=rms(r);
    out.push_back({0.0,rr,0,0,true});
    double last=0.0;
    int gm_total=0;
    for(int ni=0; ni<max_newton; ni++) {
        float eta;
        if(ni==0 || last<=0) eta=0.2f;
        else {
            float ratio=(float)(rr/(last + std::numeric_limits<float>::epsilon()));
            eta=ratio*std::sqrt(std::max(ratio,0.0f));
            eta=std::max(1e-6f,std::min(eta,0.9f));
        }
        scale_copy(rhs,r,-1.0f);
        int gi=gmres(b,mask,u,rhs,eta,delta);
        if(gi<0){ out.push_back({std::chrono::duration<double,std::milli>(Clock::now()-t0).count(),rr,ni,gm_total,false}); break; }
        gm_total += gi;
        float alpha=1.0f; bool acc=false; double trr=rr;
        for(int ls=0;ls<b.ls_max;ls++){
            tmp=u; axpy(tmp,delta,alpha); residual(b,mask,tmp,old,tr); trr=rms(tr);
            if(trr <= (1.0-b.armijo*alpha)*rr || trr < rr){ u.swap(tmp); r.swap(tr); acc=true; break; }
            alpha*=0.5f;
        }
        if(!acc){ out.push_back({std::chrono::duration<double,std::milli>(Clock::now()-t0).count(),rr,ni,gm_total,false}); break; }
        last=rr; rr=trr;
        double ms=std::chrono::duration<double,std::milli>(Clock::now()-t0).count();
        out.push_back({ms,rr,ni+1,gm_total,true});
        if(rr < 1e-12) break;
        if(out.size()>=3){ double prev=out[out.size()-2].res; if(rr>=prev*0.999999) break; }
    }
    return out;
}

static double median(std::vector<double> v){std::sort(v.begin(),v.end());return v[v.size()/2];}

int main(int argc,char**argv){
    std::vector<int> grids={256,512};
    if(argc>1){grids.clear(); for(int i=1;i<argc;i++) grids.push_back(std::atoi(argv[i]));}
    std::cout<<std::setprecision(10);
    std::cout<<"method,grid,work_index,wall_ms,residual,aux,ok\n";
    for(int nx:grids){
        // Warmup same operators.
        Params bw=make_params(nx); std::vector<unsigned char> mw; build_mask(bw,mw); std::vector<float> uw,rw; init_field(bw,mw,uw); residual(bw,mw,uw,uw,rw);
        int pmax = nx<=256 ? 12 : (nx<=512 ? 14 : 16);
        for(int p=2;p<=pmax;p++){
            std::vector<double> times, ress; int m=0; bool allok=true;
            for(int rep=0;rep<3;rep++){
                Sample s=nrii_one_order(nx,p); allok &= s.ok; if(s.ok){times.push_back(s.ms); ress.push_back(s.res); m=s.aux2;}
            }
            if(allok && !times.empty()) std::cout<<"NRII,"<<nx<<","<<p<<","<<median(times)<<","<<median(ress)<<","<<m<<",1\n";
            else std::cout<<"NRII,"<<nx<<","<<p<<","<<(times.empty()?NAN:median(times))<<","<<NAN<<","<<m<<",0\n";
        }
        // Newton trace: repeat full trace 3 times, use middle timing run by total duration then emit that trace.
        std::vector<std::vector<Sample>> traces;
        for(int rep=0;rep<3;rep++) traces.push_back(newton_trace_one_step(nx,10));
        std::vector<std::pair<double,int>> totals;
        for(int i=0;i<3;i++) totals.push_back({traces[i].empty()?1e300:traces[i].back().ms,i});
        std::sort(totals.begin(),totals.end()); auto &tr=traces[totals[1].second];
        for(auto &s:tr) std::cout<<"Newton,"<<nx<<","<<s.aux1<<","<<s.ms<<","<<s.res<<","<<s.aux2<<","<<(s.ok?1:0)<<"\n";
        std::cout.flush();
    }
}
