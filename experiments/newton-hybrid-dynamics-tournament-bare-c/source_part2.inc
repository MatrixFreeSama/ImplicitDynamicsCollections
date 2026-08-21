        if(fabs(di[i])<1e-30 || !isfinite(di[i])) {free(di);free(b);free(up);return 0;}
        x[i]=(b[i]-up[i]*x[i+1])/di[i];
    }
    free(di);free(b);free(up); return finite_vec(x,n);
}

static int newton_direction(const Problem *P,const double *q,const double *r,double shift,double *dx,SolveStats *S) {
    double lo[NDOF-1],di[NDOF],up[NDOF-1],b[NDOF];
    jacobian_tridiag(P,q,lo,di,up,shift); S->jacobian_evals++;
    for(int i=0;i<NDOF;i++) b[i]=-r[i];
    if(!thomas_solve(lo,di,up,b,dx,NDOF)) return 0;
    S->linear_solves++; return 1;
}

static double relative_forward_error(const Problem *P,const double *q) {
    double m=0.0,qs=infnorm(P->qstar,NDOF);
    for(int i=0;i<NDOF;i++){double e=fabs(q[i]-P->qstar[i]);if(e>m)m=e;}
    return m/(1.0+qs);
}

static double eval_relres(const Problem *P,const double *q,double *r,SolveStats *S){
    residual(P,q,r); S->residual_evals++; return infnorm(r,NDOF)/P->rhs_scale;
}

/* Solve tiny dense normal equations (m<=3), Gaussian elimination. */
static int solve_small(double A[AA_M][AA_M],double b[AA_M],double x[AA_M],int m){
    double M[AA_M][AA_M+1];
    for(int i=0;i<m;i++){for(int j=0;j<m;j++)M[i][j]=A[i][j];M[i][m]=b[i];}
    for(int k=0;k<m;k++){
        int p=k;double best=fabs(M[k][k]);
        for(int i=k+1;i<m;i++){double a=fabs(M[i][k]);if(a>best){best=a;p=i;}}
        if(best<1e-24||!isfinite(best))return 0;
        if(p!=k)for(int j=k;j<=m;j++){double t=M[k][j];M[k][j]=M[p][j];M[p][j]=t;}
        double piv=M[k][k];for(int j=k;j<=m;j++)M[k][j]/=piv;
        for(int i=0;i<m;i++)if(i!=k){double w=M[i][k];for(int j=k;j<=m;j++)M[i][j]-=w*M[k][j];}
    }
    for(int i=0;i<m;i++){ x[i]=M[i][m]; }
    return 1;
}

/* NRII coefficients for the full implicit-Euler first-order relation:
   q(h)=h*dt*v(h), v(h)=h*dt*(fext-f_int(q(h))), q0=v0=0. */
static void cubic_coeff_vector(double coeff[NRII_ORDER+1][NDOF],int n,double *out){
    for(int i=0;i<NDOF;i++)out[i]=0.0;
    for(int a=0;a<=n;a++)for(int b=0;b<=n-a;b++){
        int c=n-a-b;
        for(int i=0;i<NDOF;i++)out[i]+=coeff[a][i]*coeff[b][i]*coeff[c][i];
    }
}

static int nrii_candidate(const Problem *P,double *cand,int order){
    static double qc[NRII_ORDER+1][NDOF];
    static double vc[NRII_ORDER+1][NDOF];
    memset(qc,0,sizeof(qc));memset(vc,0,sizeof(vc));
    for(int kk=1;kk<=order;kk++){
        for(int i=0;i<NDOF;i++)qc[kk][i]=P->dt*vc[kk-1][i];
        int n=kk-1;
        double cub[NDOF];cubic_coeff_vector(qc,n,cub);
        double fc[NDOF];
        for(int i=0;i<NDOF;i++)fc[i]=(n==0?P->fext[i]:0.0)-P->k*qc[n][i]-P->alpha*cub[i];
        for(int i=0;i<NDOF-1;i++){
            double dn=qc[n][i+1]-qc[n][i];
            double ecub=0.0;
            for(int a=0;a<=n;a++)for(int b=0;b<=n-a;b++){
                int c=n-a-b;
                double da=qc[a][i+1]-qc[a][i];
                double db=qc[b][i+1]-qc[b][i];
                double dc=qc[c][i+1]-qc[c][i];
                ecub += da*db*dc;
            }
            double edge=P->kc*dn+P->beta*ecub;
            /* -f_int edge contribution */
            fc[i]+=edge;fc[i+1]-=edge;
        }
        for(int i=0;i<NDOF;i++)vc[kk][i]=P->dt*fc[i];
        if(!finite_vec(qc[kk],NDOF)||!finite_vec(vc[kk],NDOF))return 0;
    }
    for(int i=0;i<NDOF;i++){double s=0.0;for(int k=0;k<=order;k++)s+=qc[k][i];cand[i]=s;}
    return finite_vec(cand,NDOF);
}

static SolveStats solve_problem(const Problem *P,Method M){
    SolveStats S;memset(&S,0,sizeof(S));S.nrii_order=NRII_ORDER;
    double t0=now_sec();
    double q[NDOF]={0},r[NDOF],dx[NDOF];
    double rel=eval_relres(P,q,r,&S), r0=rel;

    if(M==M_NRII_NEWTON){
        double cand[NDOF];
        if(nrii_candidate(P,cand,NRII_ORDER)){
            double rc[NDOF];double relc=eval_relres(P,cand,rc,&S);
            /* Residual is an audit/handoff gate only. It does not correct NRII coefficients. */
            if(relc<rel && infnorm(cand,NDOF)<1e100){memcpy(q,cand,sizeof(q));memcpy(r,rc,sizeof(r));rel=relc;S.nrii_accepted=1;S.handoff_used=1;}
        }
    }

    double prev_dx[NDOF]={0};int have_prev=0;double omega=1.0;
    double Xhist[AA_M+1][NDOF],Fhist[AA_M+1][NDOF];int hcount=0;
    double *H=NULL;int broyden_phase=0;
    if(M==M_BROYDEN_NEWTON){
        H=(double*)calloc((size_t)NDOF*(size_t)NDOF,sizeof(double));
        if(H){
            /* Build an exact inverse of the initial tridiagonal Jacobian once.
               This makes the Broyden phase a fair quasi-Newton reuse experiment rather
               than an intentionally poor diagonal-preconditioner experiment. */
            double lo0[NDOF-1],di0[NDOF],up0[NDOF-1];
            jacobian_tridiag(P,q,lo0,di0,up0,0.0);S.jacobian_evals++;
            int inv_ok=1;
            for(int col=0;col<NDOF;col++){
                double e[NDOF]={0},z[NDOF];e[col]=1.0;
                if(!thomas_solve(lo0,di0,up0,e,z,NDOF)){inv_ok=0;break;}
                S.linear_solves++;
                for(int row=0;row<NDOF;row++)H[(size_t)row*NDOF+col]=z[row];
            }
            if(inv_ok)broyden_phase=1;else {free(H);H=NULL;}
        }
    }
    double mu=(M==M_PTC_NEWTON)?10.0:0.0;

    for(int it=0;it<MAX_IT;it++){
        if(rel<=1e-12){S.iterations=it;break;}
        int step_ok=1;
        if(M==M_BROYDEN_NEWTON && broyden_phase){
            for(int i=0;i<NDOF;i++){double z=0.0;for(int j=0;j<NDOF;j++)z+=H[(size_t)i*NDOF+j]*r[j];dx[i]=-z;}
            S.linear_solves++;
            double qn[NDOF],rn[NDOF];for(int i=0;i<NDOF;i++)qn[i]=q[i]+dx[i];
            double reln=eval_relres(P,qn,rn,&S);
            if(!isfinite(reln)||infnorm(qn,NDOF)>1e100){step_ok=0;}
            else {
                double svec[NDOF],yvec[NDOF],Hy[NDOF];
                for(int i=0;i<NDOF;i++){svec[i]=dx[i];yvec[i]=rn[i]-r[i];}
                for(int i=0;i<NDOF;i++){double z=0.0;for(int j=0;j<NDOF;j++)z+=H[(size_t)i*NDOF+j]*yvec[j];Hy[i]=z;}
                double yy=dotv(yvec,yvec,NDOF);
                if(yy>1e-30){for(int i=0;i<NDOF;i++){double u=svec[i]-Hy[i];for(int j=0;j<NDOF;j++)H[(size_t)i*NDOF+j]+=u*yvec[j]/yy;}}
                memcpy(q,qn,sizeof(q));memcpy(r,rn,sizeof(r));rel=reln;
                if(rel<1e-1*r0 || it>=3){broyden_phase=0;S.handoff_used=1;}
            }
        } else if(M==M_PTC_NEWTON){
            if(mu<1e-10) mu=0.0;
            if(!newton_direction(P,q,r,mu,dx,&S)){step_ok=0;}
