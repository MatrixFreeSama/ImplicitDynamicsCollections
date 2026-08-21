
static int solve_small(double A[AA_M][AA_M],double b[AA_M],double x[AA_M],int m){
    double M[AA_M][AA_M+1];
    for(int i=0;i<m;i++){for(int j=0;j<m;j++)M[i][j]=A[i][j];M[i][m]=b[i];}
    for(int k=0;k<m;k++){
        int p=k;double best=fabs(M[k][k]);for(int i=k+1;i<m;i++){double a=fabs(M[i][k]);if(a>best){best=a;p=i;}}
        if(best<1e-24||!isfinite(best))return 0;
        if(p!=k)for(int j=k;j<=m;j++){double t=M[k][j];M[k][j]=M[p][j];M[p][j]=t;}
        double piv=M[k][k];for(int j=k;j<=m;j++)M[k][j]/=piv;
        for(int i=0;i<m;i++)if(i!=k){double w=M[i][k];for(int j=k;j<=m;j++)M[i][j]-=w*M[k][j];}
    }
    for(int i=0;i<m;i++){x[i]=M[i][m];}
    return 1;
}

/* coefficient of componentwise q^3 at order n */
static void cubic_coeff(double c[NRII_ORDER+1][NDOF],int n,double*out){
    for(int i=0;i<NDOF;i++)out[i]=0.0;
    for(int a=0;a<=n;a++)for(int b=0;b<=n-a;b++){int cc=n-a-b;for(int i=0;i<NDOF;i++)out[i]+=c[a][i]*c[b][i]*c[cc][i];}
}

/* Re-anchored NRII for R(q)=q+f_int(q)-fext=0.
   q(h)=S+h[(fext-S)-f_int(q(h))]. At h=1 the original target is recovered.
   Residual is used only after construction as an accept/reject handoff audit. */
static int nrii_reanchor_candidate(const Problem*P,const double*Sanchor,double*cand){
    static double c[NRII_ORDER+1][NDOF];memset(c,0,sizeof(c));memcpy(c[0],Sanchor,sizeof(c[0]));
    for(int k=1;k<=NRII_ORDER;k++){
        int n=k-1;double cub[NDOF];cubic_coeff(c,n,cub);double fn[NDOF];
        for(int i=0;i<NDOF;i++)fn[i]=P->k*c[n][i]+P->alpha*cub[i];
        for(int i=0;i<NDOF-1;i++){
            double dn=c[n][i+1]-c[n][i],ecub=0.0;
            for(int a=0;a<=n;a++)for(int b=0;b<=n-a;b++){
                int cc=n-a-b;
                double da=c[a][i+1]-c[a][i],db=c[b][i+1]-c[b][i],dc=c[cc][i+1]-c[cc][i];
                ecub+=da*db*dc;
            }
            double edge=P->kc*dn+P->beta*ecub;fn[i]-=edge;fn[i+1]+=edge;
        }
        for(int i=0;i<NDOF;i++)c[k][i]=(n==0?(P->fext[i]-Sanchor[i]):0.0)-fn[i];
        if(!finite_vec(c[k],NDOF))return 0;
    }
    for(int i=0;i<NDOF;i++){double s=0;for(int k=0;k<=NRII_ORDER;k++)s+=c[k][i];cand[i]=s;}
    return finite_vec(cand,NDOF)&&infnorm(cand,NDOF)<1e100;
}

static int build_exact_inverse(const Problem*P,const double*q,double*H,SolveStats*S){
    double lo[NDOF-1],di[NDOF],up[NDOF-1];jacobian_tridiag(P,q,lo,di,up,0.0);S->jacobian_evals++;
    for(int col=0;col<NDOF;col++){
        double e[NDOF]={0},z[NDOF];e[col]=1.0;if(!thomas_solve(lo,di,up,e,z,NDOF))return 0;S->linear_solves++;
        for(int row=0;row<NDOF;row++)H[(size_t)row*NDOF+col]=z[row];
    }
    return 1;
}

static int update_inverse_broyden(double *H,const double*svec,const double*y){
    double Hy[NDOF],v[NDOF];
    for(int i=0;i<NDOF;i++){double z=0;for(int j=0;j<NDOF;j++)z+=H[(size_t)i*NDOF+j]*y[j];Hy[i]=z;}
    for(int j=0;j<NDOF;j++){double z=0;for(int k=0;k<NDOF;k++)z+=svec[k]*H[(size_t)k*NDOF+j];v[j]=z;}
    double den=dotv(v,y,NDOF);
    if(fabs(den)<1e-24||!isfinite(den))return 0;
    for(int i=0;i<NDOF;i++){double ui=svec[i]-Hy[i];for(int j=0;j<NDOF;j++)H[(size_t)i*NDOF+j]+=ui*v[j]/den;}
    return 1;
}

static double projected_newton_steps(double rel,double ratio){
    if(rel<=1e-12)return 0.0;
    if(!(ratio>0.0&&ratio<1.0)||!isfinite(ratio))return INFINITY;
    double n=log(1e-12/rel)/log(ratio);
    return n<0.0?0.0:n;
}

/*
   Information-deficiency gate.
   Newton remains primary. After at least two Newton stages since the previous redirect,
   a specialist is requested only when the most recent Newton progress is nondecreasing
   or its observed geometric rate predicts that the declared remaining dispatch budget
   is insufficient to reach the residual tolerance. This is an operational gate, not a
   theorem about the mathematical Newton basin.
*/
static int gate_requests_redirect(double rel_before,double rel_after,int newton_since_redirect,int remaining,SolveStats*S){
    if(newton_since_redirect<2||rel_after<=1e-12||remaining<=0)return 0;
    if(!isfinite(rel_after)) {S->trigger_invalid++;return 1;}
    double ratio=rel_after/fmax(rel_before,DBL_MIN);
    if(!(ratio<1.0)){S->trigger_nondecrease++;return 1;}
    double need=projected_newton_steps(rel_after,ratio);
    if(need>(double)remaining){S->trigger_budget_deficit++;return 1;}
    return 0;
}

static int specialist_linesearch(const Problem*P,double*q,double*r,double*rel,SolveStats*S){
    double dx[NDOF];if(!newton_direction(P,q,r,0.0,dx,S))return 0;
    double phi0=dotv(r,r,NDOF),lam=1.0;
    for(int ls=0;ls<80;ls++){
        double qt[NDOF],rt[NDOF];for(int i=0;i<NDOF;i++)qt[i]=q[i]+lam*dx[i];
        double rr=eval_relres(P,qt,rt,S);S->line_search_evals++;
        double phi=dotv(rt,rt,NDOF);
        if(isfinite(rr)&&phi<=(1.0-1e-4*lam)*(1.0-1e-4*lam)*phi0&&rr<*rel){memcpy(q,qt,sizeof(qt));memcpy(r,rt,sizeof(rt));*rel=rr;return 1;}
        lam*=0.5;
    }
    return 0;
}

static int specialist_aitken(const Problem*P,double*q,double*r,double*rel,const double*prev_dx,int have_prev,double*omega,SolveStats*S){
    double dx[NDOF];if(!newton_direction(P,q,r,0.0,dx,S))return 0;
    double qr[NDOF],rrv[NDOF];for(int i=0;i<NDOF;i++)qr[i]=q[i]+dx[i];double rr=eval_relres(P,qr,rrv,S);
    double best=rr;double qb[NDOF],rb[NDOF];memcpy(qb,qr,sizeof(qb));memcpy(rb,rrv,sizeof(rb));
    if(have_prev){double df[NDOF];for(int i=0;i<NDOF;i++)df[i]=dx[i]-prev_dx[i];double den=dotv(df,df,NDOF);if(den>1e-30){double w=-(*omega)*dotv(prev_dx,df,NDOF)/den;if(!isfinite(w))w=1.0;if(w<0.05)w=0.05;if(w>1.5)w=1.5;double qa[NDOF],rav[NDOF];for(int i=0;i<NDOF;i++)qa[i]=q[i]+w*dx[i];double ra=eval_relres(P,qa,rav,S);if(ra<best){best=ra;memcpy(qb,qa,sizeof(qb));memcpy(rb,rav,sizeof(rb));*omega=w;S->accel_accepts++;}}}
    if(isfinite(best)&&best<*rel){memcpy(q,qb,sizeof(qb));memcpy(r,rb,sizeof(rb));*rel=best;return 1;}return 0;
}

static int specialist_anderson(const Problem*P,double*q,double*r,double*rel,double Xhist[AA_M+1][NDOF],double Fhist[AA_M+1][NDOF],int*hcount,SolveStats*S){
    double dx[NDOF];if(!newton_direction(P,q,r,0.0,dx,S))return 0;
    if(*hcount<AA_M+1){memcpy(Xhist[*hcount],q,sizeof(Xhist[0]));memcpy(Fhist[*hcount],dx,sizeof(Fhist[0]));(*hcount)++;}
    else{for(int h=0;h<AA_M;h++){memcpy(Xhist[h],Xhist[h+1],sizeof(Xhist[h]));memcpy(Fhist[h],Fhist[h+1],sizeof(Fhist[h]));}memcpy(Xhist[AA_M],q,sizeof(Xhist[AA_M]));memcpy(Fhist[AA_M],dx,sizeof(Fhist[AA_M]));}
    double qr[NDOF],rrv[NDOF];for(int i=0;i<NDOF;i++)qr[i]=q[i]+dx[i];double rr=eval_relres(P,qr,rrv,S);
    double best=rr,qb[NDOF],rb[NDOF];memcpy(qb,qr,sizeof(qb));memcpy(rb,rrv,sizeof(rb));
    int m=(*hcount-1<AA_M)?*hcount-1:AA_M;
    if(m>0){double DF[AA_M][NDOF],DX[AA_M][NDOF];int start=*hcount-m-1;if(start<0)start=0;
        for(int c=0;c<m;c++){int j=start+c;for(int i=0;i<NDOF;i++){DF[c][i]=Fhist[j+1][i]-Fhist[j][i];DX[c][i]=Xhist[j+1][i]-Xhist[j][i];}}
        double A[AA_M][AA_M]={{0}},b[AA_M]={0},gam[AA_M]={0};for(int a=0;a<m;a++){b[a]=dotv(DF[a],dx,NDOF);for(int c=0;c<m;c++)A[a][c]=dotv(DF[a],DF[c],NDOF);A[a][a]+=1e-14*(1.0+A[a][a]);}
