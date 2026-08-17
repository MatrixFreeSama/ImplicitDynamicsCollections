
/* ---------- outer map ---------- */
typedef struct {
    State *base;
    const Stress*s;
    Stats*st;
    double dt;
    int gmres_used_bicg;
    int flow_solver;
    int thermal_solver;
    double em_res,flow_res,therm_res,nrii_res;
} OuterCtx;

static int outer_map(OuterCtx*oc,const double*Tguess,double*Tout){
    const Stress*s=oc->s; State*base=oc->base; Stats*st=oc->st;
    double meanT=0;for(int k=0;k<N;k++)meanT+=Tguess[k];meanT/=N;
    st->nr_uses++; double R=1.0+0.8*meanT, L=0.12, V=1.0+0.25*s->outer_coupling*meanT;
    double I=nrii_rl_step(base->coil_I,oc->dt,R,L,V,64,&oc->nrii_res);
    double *x=calloc(EMN,sizeof(double)),*b=calloc(EMN,sizeof(double)); for(int k=0;k<N;k++){x[k]=base->Ar[k];x[k+N]=base->Ai[k];}
    em_rhs(b,I*(1.0+s->thermal_feedback*meanT)); EMCtx ec={Tguess,s,I}; int gi=0;double gr=0; st->gmres_uses++;
    int gok=gmres(&ec,em_apply,b,x,EMN,20,340,5e-8,&gi,&gr); oc->gmres_used_bicg=0;
    if(!gok){st->gmres_fail++;st->total_handoffs++;st->bicg_uses++;int bi=0;double br=0;gok=bicgstab(&ec,em_apply,b,x,EMN,250,5e-8,&bi,&br);gr=br;oc->gmres_used_bicg=1;if(!gok)st->bicg_fail++;}
    oc->em_res=gr; if(!gok){free(x);free(b);return 0;}
    double *Q=calloc(N,sizeof(double)),*F=calloc(N,sizeof(double));
    for(int k=0;k<N;k++){
        double ar=x[k],ai=x[k+N],sig=conductivity(Tguess[k],s),mag=ar*ar+ai*ai;
        double fb=clamp(1.0 + s->outer_coupling*(Tguess[k]-0.25),0.15,8.0);
        Q[k]=6.0*sig*mag*(1.0+0.4*I*I)*fb;
        F[k]=s->lorentz_gain*35.0*sig*mag + s->buoy_gain*8.0*(Tguess[k]-0.48); if(boundary(k%NX,k/NX)) F[k]=0.0;
    }
    double *u=calloc(N,sizeof(double));copyv(u,base->u,N);int fw=0;double fr=0;int fok=solve_flow(u,F,s,st,&fw,&fr);oc->flow_solver=fw;oc->flow_res=fr;if(!fok){free(x);free(b);free(Q);free(F);free(u);return 0;}
    copyv(Tout,Tguess,N);int tw=0;double tr=0;int tok=solve_thermal(Tout,base->Tprev,u,Q,s,oc->dt,st,&tw,&tr);oc->thermal_solver=tw;oc->therm_res=tr;
    if(tok){for(int k=0;k<N;k++){base->Ar[k]=x[k];base->Ai[k]=x[k+N];base->u[k]=u[k];}base->coil_I=I;}
    free(x);free(b);free(Q);free(F);free(u);return tok;
}

static double rel_diff(const double*a,const double*b,int n){double s=0,d=1e-30;for(int i=0;i<n;i++){double q=a[i]-b[i];s+=q*q;d+=b[i]*b[i];}return sqrt(s/d);} 

static int fixed_point_outer(OuterCtx*oc,double*T,int maxit,double tol,double*res,int*iters){
    Stats*st=oc->st;st->fixed_point_uses++;double *next=calloc(N,sizeof(double));
    for(int k=0;k<maxit;k++){
        if(!outer_map(oc,T,next)){*iters=k+1;free(next);st->fixed_point_fail++;return 0;}
        double rr=rel_diff(next,T,N);copyv(T,next,N);if(rr<tol){*res=rr;*iters=k+1;free(next);return 1;}
        if(!isfinite(rr)||rr>10.0){*res=rr;*iters=k+1;free(next);st->fixed_point_fail++;return 0;}
    }
    *res=rel_diff(next,T,N);*iters=maxit;free(next);st->fixed_point_fail++;return 0;
}
static int aitken_outer(OuterCtx*oc,double*T,int maxit,double tol,double*res,int*iters){
    Stats*st=oc->st;st->aitken_uses++;double *f=calloc(N,sizeof(double)),*r=calloc(N,sizeof(double)),*rprev=calloc(N,sizeof(double));double omega=0.5;int have=0;
    for(int k=0;k<maxit;k++){
        if(!outer_map(oc,T,f)){*iters=k+1;free(f);free(r);free(rprev);st->aitken_fail++;return 0;}
        for(int i=0;i<N;i++)r[i]=f[i]-T[i];double rr=norm2(r,N)/(norm2(f,N)+1e-30);if(rr<tol){copyv(T,f,N);*res=rr;*iters=k+1;free(f);free(r);free(rprev);return 1;}
        if(have){double num=0,den=0;for(int i=0;i<N;i++){double dr=r[i]-rprev[i];num+=rprev[i]*dr;den+=dr*dr;}if(den>1e-30)omega=clamp(-omega*num/den,0.05,1.25);}
        copyv(rprev,r,N);have=1;for(int i=0;i<N;i++)T[i]+=omega*r[i];if(!isfinite(rr)||rr>10){*res=rr;*iters=k+1;free(f);free(r);free(rprev);st->aitken_fail++;return 0;}
    }
    *res=norm2(r,N)/(norm2(f,N)+1e-30);*iters=maxit;free(f);free(r);free(rprev);st->aitken_fail++;return 0;
}
static int solve_small(double A[ANDERSON_M][ANDERSON_M],double b[ANDERSON_M],double x[ANDERSON_M],int m){
    double M[ANDERSON_M][ANDERSON_M+1];for(int i=0;i<m;i++){for(int j=0;j<m;j++)M[i][j]=A[i][j];M[i][m]=b[i];}
    for(int k=0;k<m;k++){int p=k;for(int i=k+1;i<m;i++)if(fabs(M[i][k])>fabs(M[p][k]))p=i;if(fabs(M[p][k])<1e-20)return 0;if(p!=k)for(int j=k;j<=m;j++){double t=M[k][j];M[k][j]=M[p][j];M[p][j]=t;}double d=M[k][k];for(int j=k;j<=m;j++)M[k][j]/=d;for(int i=0;i<m;i++)if(i!=k){double f=M[i][k];for(int j=k;j<=m;j++)M[i][j]-=f*M[k][j];}}
    for(int i=0;i<m;i++)x[i]=M[i][m];return 1;
}
static int anderson_outer(OuterCtx*oc,double*T,int maxit,double tol,double*res,int*iters){
    Stats*st=oc->st;st->anderson_uses++;
    double *f=calloc(N,sizeof(double)),*r=calloc(N,sizeof(double));double *hist_r[ANDERSON_M+1],*hist_f[ANDERSON_M+1];for(int i=0;i<=ANDERSON_M;i++){hist_r[i]=calloc(N,sizeof(double));hist_f[i]=calloc(N,sizeof(double));}
    int hcount=0;
    for(int k=0;k<maxit;k++){
        if(!outer_map(oc,T,f)){*iters=k+1;st->anderson_fail++;goto fail;}
        for(int i=0;i<N;i++)r[i]=f[i]-T[i];double rr=norm2(r,N)/(norm2(f,N)+1e-30);if(rr<tol){copyv(T,f,N);*res=rr;*iters=k+1;goto success;}
        int slot=hcount%(ANDERSON_M+1);copyv(hist_r[slot],r,N);copyv(hist_f[slot],f,N);hcount++;
        int m=(hcount-1<ANDERSON_M)?hcount-1:ANDERSON_M;
        if(m<=0){for(int i=0;i<N;i++)T[i]+=0.6*r[i];continue;}
        double G[ANDERSON_M][ANDERSON_M]={{0}},bb[ANDERSON_M]={0},alpha[ANDERSON_M]={0};
        for(int a=0;a<m;a++){
            int ia=(hcount-2-a)%(ANDERSON_M+1), ib=(hcount-1-a)%(ANDERSON_M+1);
            for(int bidx=0;bidx<m;bidx++){
                int ja=(hcount-2-bidx)%(ANDERSON_M+1), jb=(hcount-1-bidx)%(ANDERSON_M+1); double s=0;
                for(int i=0;i<N;i++){double da=hist_r[ib][i]-hist_r[ia][i],db=hist_r[jb][i]-hist_r[ja][i];s+=da*db;} G[a][bidx]=s;
            }
            double s=0;for(int i=0;i<N;i++){double da=hist_r[ib][i]-hist_r[ia][i];s+=da*r[i];}bb[a]=s;
        }
        for(int i=0;i<m;i++)G[i][i]+=1e-12;
        if(!solve_small(G,bb,alpha,m)){for(int i=0;i<N;i++)T[i]+=0.5*r[i];continue;}
        double *Tnew=calloc(N,sizeof(double));copyv(Tnew,f,N);
        for(int a=0;a<m;a++){
            int ia=(hcount-2-a)%(ANDERSON_M+1), ib=(hcount-1-a)%(ANDERSON_M+1);
            for(int i=0;i<N;i++)Tnew[i]-=alpha[a]*(hist_f[ib][i]-hist_f[ia][i]);
        }
        for(int i=0;i<N;i++)T[i]=0.85*Tnew[i]+0.15*f[i];free(Tnew);
    }
    *res=norm2(r,N)/(norm2(f,N)+1e-30);*iters=maxit;st->anderson_fail++;
fail:
    for(int i=0;i<=ANDERSON_M;i++){free(hist_r[i]);free(hist_f[i]);}free(f);free(r);return 0;
success:
    for(int i=0;i<=ANDERSON_M;i++){free(hist_r[i]);free(hist_f[i]);}free(f);free(r);return 1;
}

static void init_state(State*s){
    s->T=calloc(N,sizeof(double));s->Tprev=calloc(N,sizeof(double));s->Ar=calloc(N,sizeof(double));s->Ai=calloc(N,sizeof(double));s->u=calloc(N,sizeof(double));s->coil_I=0.15;
    for(int j=0;j<NY;j++)for(int i=0;i<NX;i++){int k=id(i,j);double x=(double)i/(NX-1),y=(double)j/(NY-1);s->T[k]=s->Tprev[k]=0.48+0.09*exp(-20*((x-.55)*(x-.55)+(y-.5)*(y-.5)));}
}
static void free_state(State*s){free(s->T);free(s->Tprev);free(s->Ar);free(s->Ai);free(s->u);} 
static double mean(const double*x){double s=0;for(int i=0;i<N;i++)s+=x[i];return s/N;}
static double maxv(const double*x){double m=x[0];for(int i=1;i<N;i++)if(x[i]>m)m=x[i];return m;}
static double mean_liq(const double*T,const Stress*s){double z=0;for(int i=0;i<N;i++)z+=liquid_fraction(T[i],s);return z/N;}
static void stats_zero(Stats*s){memset(s,0,sizeof(*s));}

int main(void){
    Stress cases[] = {
        {0.6,0.10,0.03,0.04,0.08,0.08,0.10,0.55,0.15,0,"weak-coupling"},
        {18.0,0.25,0.08,0.08,0.20,0.35,0.10,0.70,0.45,0,"skin-effect-contrast"},
        {2.0,0.22,0.12,0.15,1.60,0.18,0.10,0.80,0.55,0,"strong-convection"},
        {1.2,0.18,0.05,0.06,0.25,0.15,0.025,1.80,0.50,1,"phase-boundary"},
        {1.2,0.14,0.04,0.05,0.25,0.45,0.05,1.20,1.80,0,"strong-global-coupling"}
    };
    int nc=(int)(sizeof(cases)/sizeof(cases[0]));
    FILE *fc=fopen("CASE_RESULTS.csv","w"); FILE *fsc=fopen("SOLVER_COUNTS.csv","w");
    fprintf(fc,"case,success,outer_solver,outer_iters,outer_res,em_solver,flow_solver,thermal_solver,mean_T,max_T,liquid_fraction,coil_I,nrii_res,em_res,flow_res,thermal_res,handoffs\n");
    fprintf(fsc,"case,nrii_uses,gmres_uses,gmres_fail,bicgstab_uses,bicgstab_fail,picard_uses,picard_fail,newton_krylov_uses,newton_krylov_fail,thermal_newton_uses,thermal_newton_fail,semismooth_uses,semismooth_fail,fixed_point_uses,fixed_point_fail,aitken_uses,aitken_fail,anderson_uses,anderson_fail,handoffs\n");
    printf("Bare-C multiphysics solver-relay benchmark, grid %dx%d (%d cells)\n",NX,NY,N);
    int overall=1; long grand_handoffs=0;
    for(int c=0;c<nc;c++){
        State stt;init_state(&stt);Stress *sp=&cases[c];Stats stats;stats_zero(&stats);double cdt = (c==3?0.035:(c==4?0.045:0.06)); OuterCtx oc={&stt,sp,&stats,cdt,0,0,0,0,0,0,0};
        double *Twork=calloc(N,sizeof(double));copyv(Twork,stt.T,N);double ores=0;int oit=0;int success=0;const char *outer="fixed-point";
        int fp_budget=(c==4?5:24), aitken_budget=(c==4?3:28);
        success=fixed_point_outer(&oc,Twork,fp_budget,2e-5,&ores,&oit);
        if(!success){stats.total_handoffs++;outer="Aitken";copyv(Twork,stt.T,N);success=aitken_outer(&oc,Twork,aitken_budget,2e-5,&ores,&oit);}
        if(!success){stats.total_handoffs++;outer="Anderson";copyv(Twork,stt.T,N);success=anderson_outer(&oc,Twork,32,2e-5,&ores,&oit);}
        copyv(stt.T,Twork,N);
        const char *ems=oc.gmres_used_bicg?"BiCGSTAB":"GMRES";const char *fls=oc.flow_solver==1?"Picard":"Newton-Krylov";const char *ths=oc.thermal_solver==1?"ordinary-Newton":"semismooth-active-set";
        printf("%-22s success=%d outer=%-10s it=%2d res=%.3e | EM=%-8s flow=%-13s thermal=%-21s | Tmean=%.4f Tmax=%.4f liq=%.3f handoffs=%ld\n",sp->name,success,outer,oit,ores,ems,fls,ths,mean(stt.T),maxv(stt.T),mean_liq(stt.T,sp),stats.total_handoffs);
        fprintf(fc,"%s,%d,%s,%d,%.17g,%s,%s,%s,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%ld\n",sp->name,success,outer,oit,ores,ems,fls,ths,mean(stt.T),maxv(stt.T),mean_liq(stt.T,sp),stt.coil_I,oc.nrii_res,oc.em_res,oc.flow_res,oc.therm_res,stats.total_handoffs);
        fprintf(fsc,"%s,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld\n",sp->name,stats.nr_uses,stats.gmres_uses,stats.gmres_fail,stats.bicg_uses,stats.bicg_fail,stats.picard_uses,stats.picard_fail,stats.newton_krylov_uses,stats.newton_krylov_fail,stats.thermal_newton_uses,stats.thermal_newton_fail,stats.semismooth_uses,stats.semismooth_fail,stats.fixed_point_uses,stats.fixed_point_fail,stats.aitken_uses,stats.aitken_fail,stats.anderson_uses,stats.anderson_fail,stats.total_handoffs);
        if(!success)overall=0;grand_handoffs+=stats.total_handoffs;free(Twork);free_state(&stt);
    }
    fclose(fc); fclose(fsc);
    FILE*fj=fopen("RESULTS.json","w");fprintf(fj,"{\n  \"benchmark\": \"Bare-C Reduced Thermo-Electromagneto-Hydrodynamic Solver Relay\",\n  \"grid\": [%d,%d],\n  \"cells\": %d,\n  \"cases\": %d,\n  \"all_cases_completed\": %s,\n  \"grand_handoffs\": %ld,\n  \"solver_families\": [\"NRII\",\"GMRES\",\"BiCGSTAB\",\"Picard\",\"Newton-Krylov\",\"ordinary Newton\",\"semismooth active-set\",\"fixed-point coupling\",\"Aitken relaxation\",\"Anderson acceleration\"],\n  \"claim_scope\": \"literature-grounded reduced multiphysics stress benchmark; not a bitwise reproduction of the cited induction-furnace paper\"\n}\n",NX,NY,N,nc,overall?"true":"false",grand_handoffs);fclose(fj);
    return overall?0:2;
}
