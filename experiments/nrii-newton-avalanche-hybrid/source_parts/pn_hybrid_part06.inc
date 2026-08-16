            long double s=0; for(int k=L+1;k<=p;k++)s-=G[(size_t)(k-ii)*stride+k]; rhs[ii-1]=(double)s;
        }
        if(solve_linear_vec(A,rhs,x,m)){
            double md=0; if(!denominator_gate_d(x,md)) return {};
            PadeFit f; f.ok=true; f.m=m; f.L=L; f.q=x; f.min_den=md; return f;
        }
    }
    return {};
}

static double candidate_residual(const CoeffBank&B,int p,const array<PadeFit,3>&F,double&minn,double&minp){
    int W=B.W; const Params&P=B.P; double invh2=1.0/(P.dx*P.dx); vector<double>n(W),pp(W),ph(W),vext(W); minn=1e300;minp=1e300;
    for(int i=0;i<W;i++) vext[i]=vext_at_global(B.g0+i,B.N,P);
    #pragma omp parallel for reduction(min:minn,minp) schedule(static)
    for(int i=0;i<W;i++){
        for(int comp=0;comp<3;comp++){
            auto &fit=F[comp]; long double den=1.0; for(double q:fit.q)den+=q; long double num=0;
            for(int k=0;k<=fit.L;k++){ long double pk=0; for(int j=0;j<=fit.m&&j<=k;j++){ double q=(j==0?1.0:fit.q[j-1]); double c=(comp==0?B.cn[k-j][i]:(comp==1?B.cp[k-j][i]:B.cphi[k-j][i])); pk+=(long double)q*c; } num+=pk; }
            double val=(double)(num/den); if(comp==0)n[i]=val; else if(comp==1)pp[i]=val; else ph[i]=val;
        }
        minn=min(minn,n[i]);minp=min(minp,pp[i]);
    }
    double maxr=0.0;
    #pragma omp parallel for reduction(max:maxr) schedule(static)
    for(int i=2;i<W-2;i++){
        double lapn=(n[i+1]-2*n[i]+n[i-1])*invh2, lapp=(pp[i+1]-2*pp[i]+pp[i-1])*invh2;
        double dtr=(ph[i+1]+vext[i+1])-(ph[i]+vext[i]), dtl=(ph[i]+vext[i])-(ph[i-1]+vext[i-1]);
        double drn=0.5*((n[i]+n[i+1])*dtr-(n[i-1]+n[i])*dtl)*invh2;
        double drp=0.5*((pp[i]+pp[i+1])*dtr-(pp[i-1]+pp[i])*dtl)*invh2;
        double rec=P.krec*(n[i]*pp[i]-1.0); double G=avalanche_G_point(n[i],pp[i],ph[i-1],ph[i+1],vext[i-1],vext[i+1],P);
        double Fn=P.Dn*lapn-P.mun*drn-rec+G, Fp=P.Dp*lapp+P.mup*drp-rec+G;
        double rn=n[i]-B.cn[0][i]-P.dt*Fn, rp=pp[i]-B.cp[0][i]-P.dt*Fp;
        long long gi=B.g0+i; double C=doping_at_global(gi,B.N,P); double rf=-P.lambda2*(ph[i+1]-2*ph[i]+ph[i-1])*invh2-(pp[i]-n[i]+C);
        maxr=max(maxr,max(abs(rn),max(abs(rp),abs(rf))));
    }
    return maxr;
}

struct FullBodyResult { bool ok=false,cheb_ok=false,pade_ok=false; int p=0,m0=0,m1=0,m2=0; double ms=0,res=INFINITY,tail=INFINITY,minden=0,minn=0,minp=0; };
static FullBodyResult evaluate_bank_fullbody(const CoeffBank&B,int p){
    FullBodyResult R;R.p=p;if(!B.finite||p>B.K||p<2)return R;
    R.tail=cheb_tail_ratio_bank(B,p); R.cheb_ok=isfinite(R.tail)&&R.tail<=sqrt(B.P.tol); if(!R.cheb_ok)return R;
    array<PadeFit,3>F; for(int c=0;c<3;c++){ F[c]=identify_pade_component_bank(B,p,c); if(!F[c].ok)return R; }
    R.pade_ok=true;R.m0=F[0].m;R.m1=F[1].m;R.m2=F[2].m;R.minden=min(F[0].min_den,min(F[1].min_den,F[2].min_den));
    R.res=candidate_residual(B,p,F,R.minn,R.minp); R.ok=isfinite(R.res)&&R.res<=B.P.tol; return R;
}

static FullBodyResult run_fullbody_once(long long N,const Params&P,const EqProfile&E,int p){
    auto st=chrono::steady_clock::now(); CoeffBank B=build_coeff_bank(N,P,E,p); FullBodyResult R=evaluate_bank_fullbody(B,p); auto en=chrono::steady_clock::now(); R.ms=chrono::duration<double,milli>(en-st).count(); return R;
}

struct HybridPredictor {
    bool ok=false;
    int p=0,m0=0,m1=0,m2=0;
    double tail=INFINITY,min_den=0,res=INFINITY,minn=0,minp=0;
    vector<double> n,pv,ph;
};

static HybridPredictor make_hybrid_predictor(const CoeffBank&B,int p){
    HybridPredictor H; H.p=p;
    if(!B.finite || p>B.K || p<2) return H;
    H.tail=cheb_tail_ratio_bank(B,p);
    array<PadeFit,3> F;
    for(int c=0;c<3;c++){ F[c]=identify_pade_component_bank(B,p,c); if(!F[c].ok) return H; }
    H.m0=F[0].m; H.m1=F[1].m; H.m2=F[2].m;
    H.min_den=min(F[0].min_den,min(F[1].min_den,F[2].min_den));
    int W=B.W; H.n.resize(W);H.pv.resize(W);H.ph.resize(W);
    H.minn=1e300;H.minp=1e300;
    for(int i=0;i<W;i++){
        for(int comp=0;comp<3;comp++){
            auto &fit=F[comp]; long double den=1.0; for(double q:fit.q) den+=q; long double num=0;
            for(int k=0;k<=fit.L;k++){
                long double pk=0;
                for(int j=0;j<=fit.m && j<=k;j++){
                    double q=(j==0?1.0:fit.q[j-1]);
                    double c=(comp==0?B.cn[k-j][i]:(comp==1?B.cp[k-j][i]:B.cphi[k-j][i]));
                    pk+=(long double)q*c;
                }
                num+=pk;
            }
            double val=(double)(num/den);
            if(comp==0) H.n[i]=val; else if(comp==1) H.pv[i]=val; else H.ph[i]=val;
        }
        H.minn=min(H.minn,H.n[i]); H.minp=min(H.minp,H.pv[i]);
    }
    double dummy1,dummy2;
    H.res=candidate_residual(B,p,F,dummy1,dummy2);
    H.ok=isfinite(H.res)&&isfinite(H.tail)&&H.tail<=sqrt(B.P.tol)&&H.min_den>0&&H.minn>0&&H.minp>0;
    return H;
}

int main(){
    omp_set_num_threads(5);
    Params base; base.local_half=768; base.eq_half=512; base.dt=0.005; base.tol=1e-10;
    EqProfile E=build_equilibrium(base); const long long N=1LL<<34;
    const int MAXP=192;
    const vector<int> declared_ps={16,24,32,48,64,80,96,112,128,160,192};
    const vector<double> drives={64,68,70,72,74,76,77,78,79,80,82,84,85,86,87,88};
    cout<<setprecision(14);
    cout<<"drive,pure_ok,pure_iters,pure_init_res,pure_res,pure_min_alpha,pure_ms,seed_p,seed_safe,seed_res,seed_tail,seed_min_den,seed_min_n,seed_min_p,hybrid_ok,hybrid_iters,hybrid_init_res,hybrid_res,hybrid_min_alpha,predictor_ms,hybrid_newton_ms,hybrid_total_ms\n";
    for(double drv:drives){
        Params P=base; P.reverse_drive=drv;
        auto pure=run_local_newton(N,P,E);
        CoeffBank scout=build_coeff_bank(N,P,E,MAXP);
        int chosen=-1;
        HybridPredictor scout_pred;
        for(int pp:declared_ps){
            auto pr=make_hybrid_predictor(scout,pp);
            if(pr.ok){ chosen=pp; scout_pred=std::move(pr); break; }
        }
        HybridPredictor pred; LocalNewtonResult hy{}; double pred_ms=0,total_ms=0;
        if(chosen>0){
            auto t0=chrono::steady_clock::now();
            CoeffBank B=build_coeff_bank(N,P,E,chosen);
            pred=make_hybrid_predictor(B,chosen);
            auto t1=chrono::steady_clock::now();
            pred_ms=chrono::duration<double,milli>(t1-t0).count();
            if(pred.ok) hy=run_local_newton(N,P,E,&pred.n,&pred.pv,&pred.ph);
            else { hy.ok=false; hy.iters=0; hy.init_res=INFINITY; hy.res=INFINITY; hy.min_alpha=0; hy.ms=0; }
            total_ms=pred_ms+hy.ms;
        } else {
            pred.ok=false; pred.res=INFINITY; pred.tail=INFINITY; pred.min_den=0; pred.minn=pred.minp=0;
            hy.ok=false; hy.iters=0; hy.init_res=INFINITY; hy.res=INFINITY; hy.min_alpha=0; hy.ms=0;
        }
        cout<<drv<<","<<(pure.ok?1:0)<<","<<pure.iters<<","<<pure.init_res<<","<<pure.res<<","<<pure.min_alpha<<","<<pure.ms<<","<<chosen<<","<<(pred.ok?1:0)<<","<<pred.res<<","<<pred.tail<<","<<pred.min_den<<","<<pred.minn<<","<<pred.minp<<","<<(hy.ok?1:0)<<","<<hy.iters<<","<<hy.init_res<<","<<hy.res<<","<<hy.min_alpha<<","<<pred_ms<<","<<hy.ms<<","<<total_ms<<"\n";
        cerr<<"drive="<<drv<<" pure="<<pure.ok<<" seedp="<<chosen<<" hy="<<hy.ok<<"\n";
    }
}
