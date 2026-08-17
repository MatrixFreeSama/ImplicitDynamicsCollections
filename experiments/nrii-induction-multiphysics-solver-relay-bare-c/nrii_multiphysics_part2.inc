/* ---------- generic GMRES ---------- */
typedef void (*ApplyFn)(const void*,const double*,double*);
static int gmres(const void *ctx,ApplyFn A,const double*b,double*x,int n,int restart,int maxit,double tol,int *iters,double *final_res){
    if(restart>MAX_GMRES_BASIS) restart=MAX_GMRES_BASIS;
    double *r=calloc((size_t)n,sizeof(double)), *w=calloc((size_t)n,sizeof(double));
    double **V=calloc((size_t)(restart+1),sizeof(double*));
    for(int i=0;i<=restart;i++) V[i]=calloc((size_t)n,sizeof(double));
    double *H=calloc((size_t)(restart+1)*restart,sizeof(double));
    double *cs=calloc((size_t)restart,sizeof(double)), *sn=calloc((size_t)restart,sizeof(double)), *g=calloc((size_t)(restart+1),sizeof(double));
    double *Ax=calloc((size_t)n,sizeof(double)), *y=calloc((size_t)restart,sizeof(double));
    int total=0, ok=0;
    A(ctx,x,Ax); for(int i=0;i<n;i++) r[i]=b[i]-Ax[i];
    double bnorm=norm2(b,n); if(bnorm<1e-30)bnorm=1.0;
    double beta=norm2(r,n); if(beta/bnorm<tol){ok=1; goto done;}
    while(total<maxit){
        for(int i=0;i<n;i++) V[0][i]=r[i]/beta;
        memset(H,0,(size_t)(restart+1)*restart*sizeof(double)); memset(g,0,(size_t)(restart+1)*sizeof(double)); g[0]=beta;
        int m=0;
        for(int j=0;j<restart && total<maxit;j++,total++){
            A(ctx,V[j],w);
            for(int i=0;i<=j;i++){ double hij=dot(w,V[i],n); H[i*restart+j]=hij; for(int k=0;k<n;k++) w[k]-=hij*V[i][k]; }
            double hnext=norm2(w,n); H[(j+1)*restart+j]=hnext;
            if(hnext>1e-30) for(int k=0;k<n;k++) V[j+1][k]=w[k]/hnext;
            for(int i=0;i<j;i++){double h1=H[i*restart+j],h2=H[(i+1)*restart+j]; H[i*restart+j]=cs[i]*h1+sn[i]*h2; H[(i+1)*restart+j]=-sn[i]*h1+cs[i]*h2;}
            double h1=H[j*restart+j],h2=H[(j+1)*restart+j]; double den=hypot(h1,h2); cs[j]=den? h1/den:1.0; sn[j]=den? h2/den:0.0;
            H[j*restart+j]=cs[j]*h1+sn[j]*h2; H[(j+1)*restart+j]=0.0;
            double gj=g[j]; g[j]=cs[j]*gj; g[j+1]=-sn[j]*gj;
            m=j+1;
            if(fabs(g[j+1])/bnorm<tol) break;
        }
        for(int i=m-1;i>=0;i--){ double s=g[i]; for(int j=i+1;j<m;j++) s-=H[i*restart+j]*y[j]; y[i]=s/H[i*restart+i]; }
        for(int j=0;j<m;j++) axpy(x,y[j],V[j],n);
        A(ctx,x,Ax); for(int i=0;i<n;i++) r[i]=b[i]-Ax[i]; beta=norm2(r,n);
        if(beta/bnorm<tol){ok=1; break;}
        if(!isfinite(beta)) break;
    }
done:
    *iters=total; *final_res=beta/bnorm;
    for(int i=0;i<=restart;i++) free(V[i]); free(V); free(H); free(cs); free(sn); free(g); free(r); free(w); free(Ax); free(y);
    return ok;
}

/* ---------- BiCGSTAB ---------- */
static int bicgstab(const void *ctx,ApplyFn A,const double*b,double*x,int n,int maxit,double tol,int *iters,double *final_res){
    double *r=calloc(n,sizeof(double)),*r0=calloc(n,sizeof(double)),*p=calloc(n,sizeof(double)),*v=calloc(n,sizeof(double)),*s=calloc(n,sizeof(double)),*t=calloc(n,sizeof(double)),*Ax=calloc(n,sizeof(double));
    A(ctx,x,Ax); for(int i=0;i<n;i++)r[i]=r0[i]=b[i]-Ax[i];
    double bn=norm2(b,n); if(bn<1e-30)bn=1.0; double rho=1,alpha=1,omega=1; int ok=0,k;
    for(k=0;k<maxit;k++){
        double rho1=dot(r0,r,n); if(fabs(rho1)<1e-30)break;
        double beta=(rho1/rho)*(alpha/omega);
        for(int i=0;i<n;i++)p[i]=r[i]+beta*(p[i]-omega*v[i]);
        A(ctx,p,v); double den=dot(r0,v,n); if(fabs(den)<1e-30)break; alpha=rho1/den;
        for(int i=0;i<n;i++)s[i]=r[i]-alpha*v[i];
        if(norm2(s,n)/bn<tol){axpy(x,alpha,p,n); ok=1; k++; break;}
        A(ctx,s,t); double tt=dot(t,t,n); if(tt<1e-30)break; omega=dot(t,s,n)/tt;
        for(int i=0;i<n;i++)x[i]+=alpha*p[i]+omega*s[i];
        for(int i=0;i<n;i++)r[i]=s[i]-omega*t[i];
        double rr=norm2(r,n)/bn; if(rr<tol){ok=1;k++;break;} if(!isfinite(rr)||fabs(omega)<1e-30)break; rho=rho1;
    }
    A(ctx,x,Ax); for(int i=0;i<n;i++)r[i]=b[i]-Ax[i]; *final_res=norm2(r,n)/bn; *iters=k;
    free(r);free(r0);free(p);free(v);free(s);free(t);free(Ax); return ok;
}

/* ---------- flow nonlinear scalar surrogate ---------- */
typedef struct { const double *u0; const double *force; const Stress*s; } FlowPicCtx;
static void flow_pic_apply(const void*cv,const double*x,double*y){
    const FlowPicCtx*c=cv; double nu=0.15;
    for(int j=0;j<NY;j++)for(int i=0;i<NX;i++){
        int k=id(i,j); if(boundary(i,j)){y[k]=x[k];continue;}
        double adv=c->s->conv_gain*c->u0[k]*dx_at(x,i,j);
        y[k]=x[k]-nu*lap_at(x,i,j)+adv;
    }
}
static void flow_residual(const double*u,const double*force,const Stress*s,double*r){
    double nu=0.15;
    for(int j=0;j<NY;j++)for(int i=0;i<NX;i++){
        int k=id(i,j); if(boundary(i,j)){r[k]=u[k];continue;}
        r[k]=u[k]-nu*lap_at(u,i,j)+s->conv_gain*u[k]*dx_at(u,i,j)-force[k];
    }
}
typedef struct { const double *u; const Stress*s; } FlowJCtx;
static void flow_j_apply(const void*cv,const double*v,double*y){
    const FlowJCtx*c=cv; double nu=0.15;
    for(int j=0;j<NY;j++)for(int i=0;i<NX;i++){
        int k=id(i,j); if(boundary(i,j)){y[k]=v[k];continue;}
        y[k]=v[k]-nu*lap_at(v,i,j)+c->s->conv_gain*(v[k]*dx_at(c->u,i,j)+c->u[k]*dx_at(v,i,j));
    }
}
static int solve_flow(double*u,const double*force,const Stress*s,Stats*st,int *which,double *resout){
    /* which 1 Picard, 2 Newton-Krylov */
    st->picard_uses++;
    double *un=calloc(N,sizeof(double)),*b=calloc(N,sizeof(double)),*r=calloc(N,sizeof(double)); copyv(un,u,N); copyv(b,force,N);
    int pic_ok=0;
    for(int it=0;it<18;it++){
        FlowPicCtx ctx={un,force,s}; double *cand=calloc(N,sizeof(double)); copyv(cand,un,N);
        int lit=0; double lr=0; int ok=bicgstab(&ctx,flow_pic_apply,b,cand,N,120,1e-8,&lit,&lr);
        if(!ok){free(cand);break;}
        double diff=0,den=1e-20; for(int k=0;k<N;k++){double d=cand[k]-un[k];diff+=d*d;den+=cand[k]*cand[k];}
        copyv(un,cand,N); free(cand);
        if(sqrt(diff/den)<2e-6){pic_ok=1;break;}
    }
    if(pic_ok){copyv(u,un,N);flow_residual(u,force,s,r);*resout=norminf(r,N);*which=1;free(un);free(b);free(r);return 1;}
    st->picard_fail++; st->total_handoffs++; st->newton_krylov_uses++;
    copyv(u,un,N);
    int nkok=0;
    for(int nit=0;nit<12;nit++){
        flow_residual(u,force,s,r); double rn=norm2(r,N); if(rn<1e-8){nkok=1;break;}
        double *rhs=calloc(N,sizeof(double)),*du=calloc(N,sizeof(double)); for(int k=0;k<N;k++)rhs[k]=-r[k];
        FlowJCtx jc={u,s}; int gi=0; double gr=0; int ok=gmres(&jc,flow_j_apply,rhs,du,N,18,140,1e-7,&gi,&gr);
        if(!ok){free(rhs);free(du);break;}
        double alpha=1.0, base=rn; double *trial=calloc(N,sizeof(double)),*rt=calloc(N,sizeof(double)); int accepted=0;
        for(int ls=0;ls<12;ls++){
            for(int k=0;k<N;k++)trial[k]=u[k]+alpha*du[k]; flow_residual(trial,force,s,rt); if(norm2(rt,N)<base){copyv(u,trial,N);accepted=1;break;} alpha*=0.5;
        }
        free(rhs);free(du);free(trial);free(rt); if(!accepted)break;
    }
    flow_residual(u,force,s,r); *resout=norminf(r,N); if(*resout<2e-6)nkok=1; *which=2;
    if(!nkok)st->newton_krylov_fail++;
    free(un);free(b);free(r);return nkok;
}

/* ---------- thermal nonlinear implicit solve ---------- */
typedef struct { const double *T; const double *u; const Stress*s; double dt,kappa; int semismooth; } ThermJCtx;
static void therm_j_apply(const void*cv,const double*v,double*y){
    const ThermJCtx*c=cv;
    for(int j=0;j<NY;j++)for(int i=0;i<NX;i++){
        int k=id(i,j); if(boundary(i,j)){y[k]=v[k];continue;}
        double dh=dheff(c->T[k],c->s,c->semismooth); if(!isfinite(dh)) dh=1e12;
        y[k]=dh*v[k]-c->dt*(c->kappa*lap_at(v,i,j)-c->u[k]*dx_at(v,i,j));
    }
}
static int therm_residual(const double*T,const double*Tprev,const double*u,const double*Q,const Stress*s,double dt,double*kout){
    double kappa=0.08; int corner=0;
    for(int j=0;j<NY;j++)for(int i=0;i<NX;i++){
        int k=id(i,j); if(boundary(i,j)){kout[k]=T[k]-0.48;continue;}
        kout[k]=enthalpy(T[k],s)-enthalpy(Tprev[k],s)-dt*(kappa*lap_at(T,i,j)-u[k]*dx_at(T,i,j)+Q[k]);
        double Ts=0.60,Tl=Ts+s->phase_width; if(fabs(T[k]-Ts)<2e-10||fabs(T[k]-Tl)<2e-10)corner=1;
    }
    return corner;
}
static int thermal_newton_core(double*T,const double*Tprev,const double*u,const double*Q,const Stress*s,double dt,int semismooth,double*resout){
    double *r=calloc(N,sizeof(double)); int ok=0;
    for(int it=0;it<18;it++){
        int corner=therm_residual(T,Tprev,u,Q,s,dt,r);
        double rn=norm2(r,N); if(rn<1e-9){ok=1;break;}
        if(corner && !semismooth){break;}
        double *rhs=calloc(N,sizeof(double)),*dT=calloc(N,sizeof(double)); for(int k=0;k<N;k++)rhs[k]=-r[k];
        ThermJCtx jc={T,u,s,dt,0.08,semismooth}; int gi=0; double gr=0; int gok=gmres(&jc,therm_j_apply,rhs,dT,N,18,160,1e-8,&gi,&gr);
        if(!gok){free(rhs);free(dT);break;}
        double alpha=1.0, base=rn; int acc=0; double *trial=calloc(N,sizeof(double)),*rt=calloc(N,sizeof(double));
        for(int ls=0;ls<12;ls++){
            for(int k=0;k<N;k++)trial[k]=T[k]+alpha*dT[k]; therm_residual(trial,Tprev,u,Q,s,dt,rt); if(norm2(rt,N)<base){copyv(T,trial,N);acc=1;break;} alpha*=0.5;
        }
        free(rhs);free(dT);free(trial);free(rt); if(!acc)break;
    }
    therm_residual(T,Tprev,u,Q,s,dt,r); *resout=norminf(r,N); if(*resout<2e-6)ok=1; free(r); return ok;
}
static int solve_thermal(double*T,const double*Tprev,const double*u,const double*Q,const Stress*s,double dt,Stats*st,int*which,double*resout){
    st->thermal_newton_uses++;
    if(s->force_phase_corner){
        /* Put one interior cell exactly on the solidus to exercise the ordinary derivative contract. */
        T[id(NX/2,NY/2)]=0.60;
    }
    double *backup=calloc(N,sizeof(double));copyv(backup,T,N);
    int ok=thermal_newton_core(T,Tprev,u,Q,s,dt,0,resout);
    if(ok){*which=1;free(backup);return 1;}
    st->thermal_newton_fail++;st->total_handoffs++;st->semismooth_uses++;copyv(T,backup,N);
    if(s->force_phase_corner)T[id(NX/2,NY/2)]=0.60;
    ok=thermal_newton_core(T,Tprev,u,Q,s,dt,1,resout); *which=2; if(!ok)st->semismooth_fail++;free(backup);return ok;
}
