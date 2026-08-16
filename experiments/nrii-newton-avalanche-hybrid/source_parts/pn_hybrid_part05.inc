        if(best<=numeric_limits<double>::epsilon()*(double)n*norm) return false;
        if(piv!=k){ for(int j=k;j<n;j++) swap(A[(size_t)k*n+j],A[(size_t)piv*n+j]); swap(b[k],b[piv]); }
        for(int i=k+1;i<n;i++){
            double f=A[(size_t)i*n+k]/A[(size_t)k*n+k]; A[(size_t)i*n+k]=0.0;
            for(int j=k+1;j<n;j++) A[(size_t)i*n+j]-=f*A[(size_t)k*n+j]; b[i]-=f*b[k];
        }
    }
    for(int i=n-1;i>=0;i--){ double v=b[i]; for(int j=i+1;j<n;j++) v-=A[(size_t)i*n+j]*x[j]; x[i]=v/A[(size_t)i*n+i]; }
    return true;
}

static double bern_coeff_power_d(const vector<double>&q,int m,int k){
    double s=1.0; for(int j=1;j<=k&&j<=m;j++) s+=q[j-1]*binom_d(k,j)/binom_d(m,j); return s;
}
static bool bern_safe_rec_d(const vector<double>&b,int m,int depth,int maxdepth,double floor,double&minabs){
    bool pos=true,neg=true; double mn=1e300;
    for(int i=0;i<=m;i++){ if(b[i]<=floor)pos=false; if(b[i]>=-floor)neg=false; mn=min(mn,abs(b[i])); }
    if(pos||neg){ minabs=min(minabs,mn); return true; }
    if(depth>=maxdepth) return false;
    vector<double> work=b,left(m+1),right(m+1); left[0]=work[0]; right[m]=work[m];
    for(int r=1;r<=m;r++){ for(int j=0;j<=m-r;j++) work[j]=0.5*(work[j]+work[j+1]); left[r]=work[0]; right[m-r]=work[m-r]; }
    return bern_safe_rec_d(left,m,depth+1,maxdepth,floor,minabs)&&bern_safe_rec_d(right,m,depth+1,maxdepth,floor,minabs);
}
static bool denominator_gate_d(const vector<double>&q,double&minabs){
    int m=(int)q.size(); if(m==0){minabs=1.0; return true;} double qnorm=1.0; for(double v:q)qnorm+=abs(v);
    double floor=sqrt(numeric_limits<double>::epsilon())*qnorm;
    int maxdepth=30;
    vector<double>b(m+1); for(int i=0;i<=m;i++)b[i]=bern_coeff_power_d(q,m,i);
    minabs=1e300; return bern_safe_rec_d(b,m,0,maxdepth,floor,minabs);
}

struct CoeffBank {
    int W=0,K=0; long long N=0,g0=0; Params P;
    vector<vector<double>> cn,cp,cphi,cE,ce,cinv,cz,calpha,cae,cv,cg;
    bool finite=true;
};

static CoeffBank build_coeff_bank(long long N,const Params&P,const EqProfile&E,int K){
    CoeffBank B; B.W=2*P.local_half+1; B.K=K; B.N=N; B.g0=N/2-P.local_half; B.P=P;
    int W=B.W; long long center=N/2,g0=B.g0; double invh2=1.0/(P.dx*P.dx);
    auto alloc=[&](vector<vector<double>>&x){x.assign(K+1,vector<double>(W,0.0));};
    alloc(B.cn);alloc(B.cp);alloc(B.cphi);alloc(B.cE);alloc(B.ce);alloc(B.cinv);alloc(B.cz);alloc(B.calpha);alloc(B.cae);alloc(B.cv);alloc(B.cg);
    vector<double> vext(W);
    for(int t=0;t<W;t++){
        long long gi=g0+t; double n,p,ph; eq_at(gi,N,P,E,n,p,ph); double x=(double)(gi-center)/P.pulse_sigma_cells;
        double pulse=P.pulse_amp*exp(-0.5*x*x); B.cn[0][t]=n+pulse; B.cp[0][t]=p+pulse; B.cphi[0][t]=ph; vext[t]=vext_at_global(gi,N,P);
    }
    #pragma omp parallel for schedule(static)
    for(int i=1;i<W-1;i++){
        B.cE[0][i]=-((B.cphi[0][i+1]+vext[i+1])-(B.cphi[0][i-1]+vext[i-1]))/(2.0*P.dx);
        double s0=B.cE[0][i]*B.cE[0][i]+P.ion_Efloor*P.ion_Efloor;
        B.ce[0][i]=sqrt(s0); B.cinv[0][i]=1.0/B.ce[0][i]; B.cz[0][i]=-P.ion_B*B.cinv[0][i];
        B.calpha[0][i]=P.ion_A*exp(B.cz[0][i]); B.cae[0][i]=B.calpha[0][i]*B.ce[0][i];
        B.cv[0][i]=P.mun*B.cn[0][i]+P.mup*B.cp[0][i]; B.cg[0][i]=P.ion_gamma*B.cae[0][i]*B.cv[0][i];
    }
    for(int k=0;k<K;k++){
        auto &nn=B.cn[k+1],&pp=B.cp[k+1];
        #pragma omp parallel for schedule(static)
        for(int i=1;i<W-1;i++){
            double lapn=(B.cn[k][i+1]-2*B.cn[k][i]+B.cn[k][i-1])*invh2;
            double lapp=(B.cp[k][i+1]-2*B.cp[k][i]+B.cp[k][i-1])*invh2;
            double drift_n=0.0,drift_p=0.0,recconv=0.0;
            for(int a=0;a<=k;a++){ int bb=k-a;
                double dr=B.cphi[bb][i+1]-B.cphi[bb][i], dl=B.cphi[bb][i]-B.cphi[bb][i-1];
                if(bb==0){ dr+=vext[i+1]-vext[i]; dl+=vext[i]-vext[i-1]; }
                double frn=0.5*(B.cn[a][i]+B.cn[a][i+1])*dr, fln=0.5*(B.cn[a][i-1]+B.cn[a][i])*dl; drift_n+=(frn-fln)*invh2;
                double frp=0.5*(B.cp[a][i]+B.cp[a][i+1])*dr, flp=0.5*(B.cp[a][i-1]+B.cp[a][i])*dl; drift_p+=(frp-flp)*invh2;
                recconv+=B.cn[a][i]*B.cp[bb][i];
            }
            if(k==0)recconv-=1.0;
            nn[i]=P.dt*(P.Dn*lapn-P.mun*drift_n-P.krec*recconv+B.cg[k][i]);
            pp[i]=P.dt*(P.Dp*lapp+P.mup*drift_p-P.krec*recconv+B.cg[k][i]);
        }
        vector<double>rho(W);
        #pragma omp parallel for schedule(static)
        for(int i=0;i<W;i++)rho[i]=pp[i]-nn[i];
        poisson_local_from_rho(rho,B.cphi[k+1],g0,N,P);
        const int r=k+1;
        #pragma omp parallel for schedule(static)
        for(int i=1;i<W-1;i++){
            B.cE[r][i]=-(B.cphi[r][i+1]-B.cphi[r][i-1])/(2.0*P.dx);
            double sk=0.0; for(int a=0;a<=r;a++) sk+=B.cE[a][i]*B.cE[r-a][i];
            double yy=0.0; for(int a=1;a<r;a++) yy+=B.ce[a][i]*B.ce[r-a][i];
            B.ce[r][i]=(sk-yy)/(2.0*B.ce[0][i]);
            double sr=0.0; for(int a=1;a<=r;a++) sr+=B.ce[a][i]*B.cinv[r-a][i];
            B.cinv[r][i]=-sr/B.ce[0][i]; B.cz[r][i]=-P.ion_B*B.cinv[r][i];
            double se=0.0; for(int a=1;a<=r;a++) se+=(double)a*B.cz[a][i]*B.calpha[r-a][i];
            B.calpha[r][i]=se/(double)r;
            double sae=0.0; for(int a=0;a<=r;a++) sae+=B.calpha[a][i]*B.ce[r-a][i]; B.cae[r][i]=sae;
            B.cv[r][i]=P.mun*B.cn[r][i]+P.mup*B.cp[r][i];
            double sg=0.0; for(int a=0;a<=r;a++) sg+=B.cae[a][i]*B.cv[r-a][i]; B.cg[r][i]=P.ion_gamma*sg;
        }
        double bad=0.0;
        #pragma omp parallel for reduction(max:bad) schedule(static)
        for(int i=0;i<W;i++) if(!isfinite(nn[i])||!isfinite(pp[i])||!isfinite(B.cphi[k+1][i])||!isfinite(B.cg[k+1][i])) bad=1.0;
        if(bad){B.finite=false;B.K=k+1;break;}
    }
    return B;
}

static double cheb_tail_ratio_bank(const CoeffBank&B,int p){
    auto T=cheb_matrix_d(p); int W=B.W; int root=1; while(root*root<p+1)root++; int tailstart=(p+1)-root;
    long double tail2=0,total2=0;
    #pragma omp parallel for reduction(+:tail2,total2) schedule(static)
    for(int i=0;i<W;i++){
        long double tt=0,aa=0;
        for(int comp=0;comp<3;comp++){
            for(int k=0;k<=p;k++){
                long double s=0; for(int j=k;j<=p;j++){ double c=(comp==0?B.cn[j][i]:(comp==1?B.cp[j][i]:B.cphi[j][i])); s+=(long double)T[(size_t)k*(p+1)+j]*c; }
                aa+=s*s; if(k>=tailstart)tt+=s*s;
            }
        }
        total2+=aa; tail2+=tt;
    }
    return total2>0?sqrt((double)(tail2/total2)):0.0;
}

static vector<double> component_gram(const CoeffBank&B,int p,int comp){
    vector<double> G((size_t)(p+1)*(p+1),0.0); int W=B.W;
    for(int a=0;a<=p;a++) for(int b=a;b<=p;b++){
        long double s=0;
        #pragma omp parallel for reduction(+:s) schedule(static)
        for(int i=0;i<W;i++){ double x=(comp==0?B.cn[a][i]:(comp==1?B.cp[a][i]:B.cphi[a][i])); double y=(comp==0?B.cn[b][i]:(comp==1?B.cp[b][i]:B.cphi[b][i])); s+=(long double)x*y; }
        G[(size_t)a*(p+1)+b]=G[(size_t)b*(p+1)+a]=(double)s;
    }
    return G;
}

struct PadeFit { bool ok=false; int m=0,L=0; vector<double>q; double min_den=0; };
static PadeFit identify_pade_component_bank(const CoeffBank&B,int p,int comp){
    auto G=component_gram(B,p,comp); int stride=p+1;
    for(int m=p/2;m>0;m--){ int L=p-m; vector<double>A((size_t)m*m),rhs(m),x;
        for(int ii=1;ii<=m;ii++){
            for(int jj=1;jj<=m;jj++){ long double s=0; for(int k=L+1;k<=p;k++)s+=G[(size_t)(k-ii)*stride+(k-jj)]; A[(size_t)(ii-1)*m+jj-1]=(double)s; }
