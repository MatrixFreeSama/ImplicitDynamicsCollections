    auto residual_max=[&](){
        double mr=0;
        #pragma omp parallel for reduction(max:mr) schedule(static)
        for(int i=1;i<N-1;i++){
            double lapn=(n[i+1]-2*n[i]+n[i-1])*invh2;
            double lapp=(p[i+1]-2*p[i]+p[i-1])*invh2;
            double gn=0.5*((n[i]+n[i+1])*(ph[i+1]-ph[i])-(n[i-1]+n[i])*(ph[i]-ph[i-1]))*invh2;
            double gp=0.5*((p[i]+p[i+1])*(ph[i+1]-ph[i])-(p[i-1]+p[i])*(ph[i]-ph[i-1]))*invh2;
            double rec=P.krec*(n[i]*p[i]-1.0);
            double rn=n[i]-on[i]-P.dt*(P.Dn*lapn-P.mun*gn-rec);
            double rp=p[i]-op[i]-P.dt*(P.Dp*lapp+P.mup*gp-rec);
            double C=doping_at_global(i,Nll,P);
            double rf=-P.lambda2*(ph[i+1]-2*ph[i]+ph[i-1])*invh2-(p[i]-n[i]+C);
            mr=max(mr,max(abs(rn),max(abs(rp),abs(rf))));
        }
        return mr;
    };
    auto st=chrono::steady_clock::now();
    int its=0; bool ok=true; double mr=residual_max();
    for(int it=0;it<12 && mr>P.tol;it++){
        its=it+1;
        Mat3 Cprev{}; Vec3 dprev{};
        for(int ii=0;ii<(int)M;ii++){
            int i=ii+1;
            Mat3 A{},B{},C{}; Vec3 rhs{};
            double lapn=(n[i+1]-2*n[i]+n[i-1])*invh2;
            double lapp=(p[i+1]-2*p[i]+p[i-1])*invh2;
            double gn=0.5*((n[i]+n[i+1])*(ph[i+1]-ph[i])-(n[i-1]+n[i])*(ph[i]-ph[i-1]))*invh2;
            double gp=0.5*((p[i]+p[i+1])*(ph[i+1]-ph[i])-(p[i-1]+p[i])*(ph[i]-ph[i-1]))*invh2;
            double rec=P.krec*(n[i]*p[i]-1.0);
            double rn=n[i]-on[i]-P.dt*(P.Dn*lapn-P.mun*gn-rec);
            double rp=p[i]-op[i]-P.dt*(P.Dp*lapp+P.mup*gp-rec);
            double dop=doping_at_global(i,Nll,P);
            double rf=-P.lambda2*(ph[i+1]-2*ph[i]+ph[i-1])*invh2-(p[i]-n[i]+dop);
            rhs.v[0]=-rn;rhs.v[1]=-rp;rhs.v[2]=-rf;
            double dgn_nim=-0.5*(ph[i]-ph[i-1])*invh2;
            double dgn_ni=0.5*((ph[i+1]-ph[i])-(ph[i]-ph[i-1]))*invh2;
            double dgn_nip=0.5*(ph[i+1]-ph[i])*invh2;
            double dgn_fim=0.5*(n[i-1]+n[i])*invh2;
            double dgn_fi=-0.5*((n[i]+n[i+1])+(n[i-1]+n[i]))*invh2;
            double dgn_fip=0.5*(n[i]+n[i+1])*invh2;
            double dgp_pim=-0.5*(ph[i]-ph[i-1])*invh2;
            double dgp_pi=0.5*((ph[i+1]-ph[i])-(ph[i]-ph[i-1]))*invh2;
            double dgp_pip=0.5*(ph[i+1]-ph[i])*invh2;
            double dgp_fim=0.5*(p[i-1]+p[i])*invh2;
            double dgp_fi=-0.5*((p[i]+p[i+1])+(p[i-1]+p[i]))*invh2;
            double dgp_fip=0.5*(p[i]+p[i+1])*invh2;
            A.a[0]=-P.dt*(P.Dn*invh2 - P.mun*dgn_nim);
            A.a[2]=-P.dt*(-P.mun*dgn_fim);
            A.a[4]=-P.dt*(P.Dp*invh2 + P.mup*dgp_pim);
            A.a[5]=-P.dt*(P.mup*dgp_fim);
            A.a[8]=-P.lambda2*invh2;
            B.a[0]=1.0 - P.dt*(-2*P.Dn*invh2 - P.mun*dgn_ni - P.krec*p[i]);
            B.a[1]=P.dt*P.krec*n[i];
            B.a[2]=-P.dt*(-P.mun*dgn_fi);
            B.a[3]=P.dt*P.krec*p[i];
            B.a[4]=1.0 - P.dt*(-2*P.Dp*invh2 + P.mup*dgp_pi - P.krec*n[i]);
            B.a[5]=-P.dt*(P.mup*dgp_fi);
            B.a[6]=1.0; B.a[7]=-1.0; B.a[8]=2*P.lambda2*invh2;
            C.a[0]=-P.dt*(P.Dn*invh2 - P.mun*dgn_nip);
            C.a[2]=-P.dt*(-P.mun*dgn_fip);
            C.a[4]=-P.dt*(P.Dp*invh2 + P.mup*dgp_pip);
            C.a[5]=-P.dt*(P.mup*dgp_fip);
            C.a[8]=-P.lambda2*invh2;
            if(ii>0){
                Mat3 AC=mulM(A,Cprev); B=subM(B,AC);
                Vec3 Ad=mulV(A,dprev); rhs=subV(rhs,Ad);
            }
            Mat3 Binv; if(!inv3(B,Binv)){ok=false;break;}
            Ct[ii]=mulM(Binv,C);
            dtv[ii]=mulV(Binv,rhs);
            Cprev=Ct[ii]; dprev=dtv[ii];
        }
        if(!ok) break;
        for(int ii=(int)M-2;ii>=0;ii--){ Vec3 cx=mulV(Ct[ii],dtv[ii+1]); dtv[ii]=subV(dtv[ii],cx); }
        double alpha=1.0;
        for(int ii=0;ii<(int)M;ii++){
            int i=ii+1;
            if(dtv[ii].v[0]<0) alpha=min(alpha,0.8*n[i]/(-dtv[ii].v[0]));
            if(dtv[ii].v[1]<0) alpha=min(alpha,0.8*p[i]/(-dtv[ii].v[1]));
        }
        alpha=min(alpha,1.0);
        #pragma omp parallel for schedule(static)
        for(int ii=0;ii<(int)M;ii++){ int i=ii+1; n[i]+=alpha*dtv[ii].v[0]; p[i]+=alpha*dtv[ii].v[1]; ph[i]+=alpha*dtv[ii].v[2]; }
        mr=residual_max();
        if(!isfinite(mr)){ok=false;break;}
    }
    auto en=chrono::steady_clock::now();
    double ms=chrono::duration<double,milli>(en-st).count();
    mr=residual_max();
    return {ms,mr,its,ok && mr<=P.tol,gib};
}

struct LocalNewtonResult {
    double ms,res,outer_res;
    int iters;
    bool ok;
    int W;
    double edge_carrier_delta;
    double init_res;
    double state_delta_max;
    double min_alpha;
};

LocalNewtonResult run_local_newton(long long N,const Params&P,const EqProfile&E, const vector<double>* seed_n=nullptr, const vector<double>* seed_p=nullptr, const vector<double>* seed_ph=nullptr){
    int H=P.local_half, W=2*H+1;
    long long center=N/2, g0=center-H, g1=center+H;
    if(g0<=1 || g1>=N-2) return {0,INFINITY,INFINITY,0,false,W,INFINITY,INFINITY,INFINITY,0.0};
    vector<double> n(W),p(W),ph(W),on(W),op(W),eqph(W),vext(W);
    for(int t=0;t<W;t++){
        long long gi=g0+t; double ne,pe,fe; eq_at(gi,N,P,E,ne,pe,fe);
        double x=(double)(gi-center)/P.pulse_sigma_cells;
        double pulse=P.pulse_amp*exp(-0.5*x*x);
        on[t]=ne+pulse; op[t]=pe+pulse;
        n[t]=on[t]; p[t]=op[t]; ph[t]=fe; eqph[t]=fe; vext[t]=vext_at_global(gi,N,P);
    }
    if(seed_n && seed_p && seed_ph && (int)seed_n->size()==W && (int)seed_p->size()==W && (int)seed_ph->size()==W){
        n=*seed_n; p=*seed_p; ph=*seed_ph;
    }
    vector<Mat3> Ct(W); vector<Vec3> dtv(W);
    double invh2=1.0/(P.dx*P.dx);
    const double alphaL=(double)(g0-1)/(double)g0;
    const double alphaR=(double)(N-2-g1)/(double)(N-1-g1);
    auto phi_left_ghost=[&](){ return eqph[0] + alphaL*(ph[0]-eqph[0]); };
    auto phi_right_ghost=[&](){ return eqph[W-1] + alphaR*(ph[W-1]-eqph[W-1]); };
