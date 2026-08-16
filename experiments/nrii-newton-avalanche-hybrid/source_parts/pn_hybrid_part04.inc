    auto residual_max=[&](){
        double mr=0.0;
        for(int i=0;i<W;i++){
            double rn=0,rp=0;
            if(i==0 || i==W-1){ rn=n[i]-on[i]; rp=p[i]-op[i]; }
            else {
                double lapn=(n[i+1]-2*n[i]+n[i-1])*invh2;
                double lapp=(p[i+1]-2*p[i]+p[i-1])*invh2;
                double dtr=(ph[i+1]+vext[i+1])-(ph[i]+vext[i]);
                double dtl=(ph[i]+vext[i])-(ph[i-1]+vext[i-1]);
                double gn=0.5*((n[i]+n[i+1])*dtr-(n[i-1]+n[i])*dtl)*invh2;
                double gp=0.5*((p[i]+p[i+1])*dtr-(p[i-1]+p[i])*dtl)*invh2;
                double rec=P.krec*(n[i]*p[i]-1.0);
                double G=avalanche_G_point(n[i],p[i],ph[i-1],ph[i+1],vext[i-1],vext[i+1],P);
                rn=n[i]-on[i]-P.dt*(P.Dn*lapn-P.mun*gn-rec+G);
                rp=p[i]-op[i]-P.dt*(P.Dp*lapp+P.mup*gp-rec+G);
            }
            long long gi=g0+i; double dop=doping_at_global(gi,N,P);
            double phim=(i==0)?phi_left_ghost():ph[i-1];
            double phip=(i==W-1)?phi_right_ghost():ph[i+1];
            double rf=-P.lambda2*(phip-2*ph[i]+phim)*invh2-(p[i]-n[i]+dop);
            mr=max(mr,max(abs(rn),max(abs(rp),abs(rf))));
        }
        return mr;
    };

    auto st=chrono::steady_clock::now();
    int its=0; bool ok=true; double mr=residual_max();
    double init_res=mr, min_alpha_seen=1.0;
    for(int it=0;it<40 && mr>P.tol;it++){
        its=it+1; Mat3 Cprev{}; Vec3 dprev{};
        for(int i=0;i<W;i++){
            Mat3 A{},B{},C{}; Vec3 rhs{}; double rn=0,rp=0;
            if(i==0 || i==W-1){ rn=n[i]-on[i]; rp=p[i]-op[i]; B.a[0]=1.0; B.a[4]=1.0; }
            else {
                double lapn=(n[i+1]-2*n[i]+n[i-1])*invh2;
                double lapp=(p[i+1]-2*p[i]+p[i-1])*invh2;
                double dtr=(ph[i+1]+vext[i+1])-(ph[i]+vext[i]);
                double dtl=(ph[i]+vext[i])-(ph[i-1]+vext[i-1]);
                double gn=0.5*((n[i]+n[i+1])*dtr-(n[i-1]+n[i])*dtl)*invh2;
                double gp=0.5*((p[i]+p[i+1])*dtr-(p[i-1]+p[i])*dtl)*invh2;
                double rec=P.krec*(n[i]*p[i]-1.0);
                double dGdn,dGdp,dGfim,dGfip;
                double G=avalanche_G_point(n[i],p[i],ph[i-1],ph[i+1],vext[i-1],vext[i+1],P,&dGdn,&dGdp,&dGfim,&dGfip);
                rn=n[i]-on[i]-P.dt*(P.Dn*lapn-P.mun*gn-rec+G);
                rp=p[i]-op[i]-P.dt*(P.Dp*lapp+P.mup*gp-rec+G);

                double dgn_nim=-0.5*dtl*invh2;
                double dgn_ni=0.5*(dtr-dtl)*invh2;
                double dgn_nip=0.5*dtr*invh2;
                double dgn_fim=0.5*(n[i-1]+n[i])*invh2;
                double dgn_fi=-0.5*((n[i]+n[i+1])+(n[i-1]+n[i]))*invh2;
                double dgn_fip=0.5*(n[i]+n[i+1])*invh2;
                double dgp_pim=-0.5*dtl*invh2;
                double dgp_pi=0.5*(dtr-dtl)*invh2;
                double dgp_pip=0.5*dtr*invh2;
                double dgp_fim=0.5*(p[i-1]+p[i])*invh2;
                double dgp_fi=-0.5*((p[i]+p[i+1])+(p[i-1]+p[i]))*invh2;
                double dgp_fip=0.5*(p[i]+p[i+1])*invh2;

                A.a[0]=-P.dt*(P.Dn*invh2 - P.mun*dgn_nim);
                A.a[2]=-P.dt*(-P.mun*dgn_fim + dGfim);
                A.a[4]=-P.dt*(P.Dp*invh2 + P.mup*dgp_pim);
                A.a[5]=-P.dt*(P.mup*dgp_fim + dGfim);
                C.a[0]=-P.dt*(P.Dn*invh2 - P.mun*dgn_nip);
                C.a[2]=-P.dt*(-P.mun*dgn_fip + dGfip);
                C.a[4]=-P.dt*(P.Dp*invh2 + P.mup*dgp_pip);
                C.a[5]=-P.dt*(P.mup*dgp_fip + dGfip);
                B.a[0]=1.0 - P.dt*(-2*P.Dn*invh2 - P.mun*dgn_ni - P.krec*p[i] + dGdn);
                B.a[1]=P.dt*P.krec*n[i] - P.dt*dGdp;
                B.a[2]=-P.dt*(-P.mun*dgn_fi);
                B.a[3]=P.dt*P.krec*p[i] - P.dt*dGdn;
                B.a[4]=1.0 - P.dt*(-2*P.Dp*invh2 + P.mup*dgp_pi - P.krec*n[i] + dGdp);
                B.a[5]=-P.dt*(P.mup*dgp_fi);
            }
            long long gi=g0+i; double dop=doping_at_global(gi,N,P);
            double phim=(i==0)?phi_left_ghost():ph[i-1];
            double phip=(i==W-1)?phi_right_ghost():ph[i+1];
            double rf=-P.lambda2*(phip-2*ph[i]+phim)*invh2-(p[i]-n[i]+dop);
            rhs.v[0]=-rn; rhs.v[1]=-rp; rhs.v[2]=-rf;
            B.a[6]=1.0; B.a[7]=-1.0;
            if(i==0){ B.a[8]=P.lambda2*(2.0-alphaL)*invh2; C.a[8]=-P.lambda2*invh2; }
            else if(i==W-1){ A.a[8]=-P.lambda2*invh2; B.a[8]=P.lambda2*(2.0-alphaR)*invh2; }
            else { A.a[8]=-P.lambda2*invh2; B.a[8]=2*P.lambda2*invh2; C.a[8]=-P.lambda2*invh2; }
            if(i>0){ Mat3 AC=mulM(A,Cprev); B=subM(B,AC); Vec3 Ad=mulV(A,dprev); rhs=subV(rhs,Ad); }
            Mat3 Binv; if(!inv3(B,Binv)){ ok=false; break; }
            Ct[i]=mulM(Binv,C); dtv[i]=mulV(Binv,rhs); Cprev=Ct[i]; dprev=dtv[i];
        }
        if(!ok) break;
        for(int i=W-2;i>=0;i--){ Vec3 cx=mulV(Ct[i],dtv[i+1]); dtv[i]=subV(dtv[i],cx); }
        double alpha=1.0;
        for(int i=0;i<W;i++){
            if(dtv[i].v[0]<0) alpha=min(alpha,0.8*n[i]/(-dtv[i].v[0]));
            if(dtv[i].v[1]<0) alpha=min(alpha,0.8*p[i]/(-dtv[i].v[1]));
        }
        alpha=min(alpha,1.0);
        vector<double> nold=n,pold=p,phold=ph; const double oldmr=mr; bool accepted=false; double newmr=INFINITY;
        for(int ls=0;ls<14;ls++){
            for(int i=0;i<W;i++){ n[i]=nold[i]+alpha*dtv[i].v[0]; p[i]=pold[i]+alpha*dtv[i].v[1]; ph[i]=phold[i]+alpha*dtv[i].v[2]; }
            newmr=residual_max();
            if(isfinite(newmr) && newmr <= oldmr*(1.0-1e-4*alpha)){ accepted=true; break; }
            alpha*=0.5;
        }
        if(!accepted){ n.swap(nold);p.swap(pold);ph.swap(phold);ok=false;break; }
        min_alpha_seen=min(min_alpha_seen,alpha); mr=newmr;
    }
    auto en=chrono::steady_clock::now(); double ms=chrono::duration<double,milli>(en-st).count(); mr=residual_max();
    double outer=0.0;
    double edge=max({abs(n[0]-on[0]),abs(p[0]-op[0]),abs(n[W-1]-on[W-1]),abs(p[W-1]-op[W-1])});
    double delta_max=0.0; for(int i=0;i<W;i++) delta_max=max(delta_max,max(abs(n[i]-on[i]),max(abs(p[i]-op[i]),abs(ph[i]-eqph[i]))));
    return {ms,mr,outer,its,ok && mr<=P.tol,W,edge,init_res,delta_max,min_alpha_seen};
}

// ===== Full NRII body: coefficient construction + Chebyshev audit + Pade rational continuation =====
static double binom_d(int n,int k){ if(k<0||k>n) return 0.0; if(k>n-k)k=n-k; double r=1.0; for(int i=1;i<=k;i++) r=r*(double)(n-k+i)/(double)i; return r; }

static vector<double> cheb_matrix_d(int p){
    vector<double> T((size_t)(p+1)*(p+1),0.0);
    for(int j=0;j<=p;j++){
        double zj=ldexp(1.0,-j);
        for(int rr=0;rr<=j;rr++){
            double cx=binom_d(j,rr)*zj;
            if(rr==0){ T[j]+=cx; continue; }
            for(int k=(rr&1);k<=rr;k+=2){
                double ck=(k==0?binom_d(rr,rr/2):2.0*binom_d(rr,(rr-k)/2))*ldexp(1.0,-rr);
                T[(size_t)k*(p+1)+j]+=cx*ck;
            }
        }
    }
    return T;
}

static bool solve_linear_vec(vector<double>A, vector<double>b, vector<double>&x, int n){
    x.assign(n,0.0); double norm=0.0; for(double v:A) norm=max(norm,abs(v)); if(!(norm>0.0)) return false;
    for(int k=0;k<n;k++){
        int piv=k; double best=0.0; for(int i=k;i<n;i++){ double v=abs(A[(size_t)i*n+k]); if(v>best){best=v;piv=i;} }
