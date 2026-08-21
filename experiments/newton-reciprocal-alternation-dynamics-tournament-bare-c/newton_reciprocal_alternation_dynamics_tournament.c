#define _POSIX_C_SOURCE 200809L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>

#define NDOF 128
#ifndef MAX_STAGES
#define MAX_STAGES 24
#endif
#define AA_M 3
#define NRII_ORDER 12
#define MAX_STRESS 121
#define PI 3.141592653589793238462643383279502884

typedef enum {
    M_NEWTON=0,
    M_NEWTON_LINESEARCH_ALT,
    M_NEWTON_AITKEN_ALT,
    M_NEWTON_ANDERSON_ALT,
    M_BROYDEN_NEWTON_ALT,
    M_PTC_NEWTON_ALT,
    M_NRII_NEWTON_ALT,
    M_COUNT
} Method;

static const char *METHOD_NAME[M_COUNT] = {
    "Newton",
    "Newton<->LineSearch",
    "Newton<->Aitken",
    "Newton<->Anderson",
    "Broyden<->Newton",
    "PseudoTransient<->Newton",
    "NRII<->Newton"
};

typedef struct {
    double dt, k, alpha, kc, beta;
    double qstar[NDOF];
    double fext[NDOF];
    double rhs_scale;
} Problem;

typedef struct {
    int success, branch_ok;
    int stages, stage_a, stage_b;
    int residual_evals, jacobian_evals, linear_solves, line_search_evals;
    int accel_accepts, handoff_count, nrii_attempts, nrii_accepts;
    double rel_residual, rel_forward_error, wall_us;
} SolveStats;

static double now_sec(void){struct timespec ts;clock_gettime(CLOCK_MONOTONIC,&ts);return (double)ts.tv_sec+1e-9*(double)ts.tv_nsec;}
static double infnorm(const double*x,int n){double m=0;for(int i=0;i<n;i++){double a=fabs(x[i]);if(a>m)m=a;}return m;}
static double dotv(const double*a,const double*b,int n){double s=0;for(int i=0;i<n;i++)s+=a[i]*b[i];return s;}
static int finite_vec(const double*x,int n){for(int i=0;i<n;i++)if(!isfinite(x[i]))return 0;return 1;}

static void internal_force(const Problem*P,const double*q,double*f){
    for(int i=0;i<NDOF;i++)f[i]=P->k*q[i]+P->alpha*q[i]*q[i]*q[i];
    for(int i=0;i<NDOF-1;i++){
        double d=q[i+1]-q[i],e=P->kc*d+P->beta*d*d*d;
        f[i]-=e;f[i+1]+=e;
    }
}

static void make_problem(Problem*P,double stress){
    P->dt=1.0;P->k=0.10;P->alpha=1.0;P->kc=0.20;P->beta=0.02;
    for(int i=0;i<NDOF;i++){
        double x=((double)i+0.5)/(double)NDOF;
        double shape=1.0+0.10*sin(2*PI*x)+0.05*sin(7*PI*x)+0.025*cos(11*PI*x);
        P->qstar[i]=stress*shape;
    }
    double fint[NDOF];internal_force(P,P->qstar,fint);
    for(int i=0;i<NDOF;i++)P->fext[i]=P->qstar[i]+fint[i];
    P->rhs_scale=1.0+infnorm(P->fext,NDOF);
}

static void residual(const Problem*P,const double*q,double*r){
    double f[NDOF];internal_force(P,q,f);
    for(int i=0;i<NDOF;i++)r[i]=q[i]+f[i]-P->fext[i];
}

static void jacobian_tridiag(const Problem*P,const double*q,double*lo,double*di,double*up,double shift){
    for(int i=0;i<NDOF;i++)di[i]=1.0+P->k+3.0*P->alpha*q[i]*q[i]+shift;
    for(int i=0;i<NDOF-1;i++){
        double d=q[i+1]-q[i],ke=P->kc+3.0*P->beta*d*d;
        di[i]+=ke;di[i+1]+=ke;lo[i]=-ke;up[i]=-ke;
    }
}

static int thomas_solve(const double*lo0,const double*di0,const double*up0,const double*b0,double*x,int n){
    double di[NDOF],b[NDOF],up[NDOF-1];
    memcpy(di,di0,(size_t)n*sizeof(double));memcpy(b,b0,(size_t)n*sizeof(double));memcpy(up,up0,(size_t)(n-1)*sizeof(double));
    for(int i=1;i<n;i++){
        if(fabs(di[i-1])<1e-30||!isfinite(di[i-1]))return 0;
        double w=lo0[i-1]/di[i-1];di[i]-=w*up[i-1];b[i]-=w*b[i-1];
    }
    if(fabs(di[n-1])<1e-30||!isfinite(di[n-1]))return 0;
    x[n-1]=b[n-1]/di[n-1];
    for(int i=n-2;i>=0;i--){if(fabs(di[i])<1e-30||!isfinite(di[i]))return 0;x[i]=(b[i]-up[i]*x[i+1])/di[i];}
    return finite_vec(x,n);
}

static int newton_direction(const Problem*P,const double*q,const double*r,double shift,double*dx,SolveStats*S){
    double lo[NDOF-1],di[NDOF],up[NDOF-1],b[NDOF];
    jacobian_tridiag(P,q,lo,di,up,shift);S->jacobian_evals++;
    for(int i=0;i<NDOF;i++)b[i]=-r[i];
    if(!thomas_solve(lo,di,up,b,dx,NDOF))return 0;
    S->linear_solves++;return 1;
}

static double eval_relres(const Problem*P,const double*q,double*r,SolveStats*S){residual(P,q,r);S->residual_evals++;return infnorm(r,NDOF)/P->rhs_scale;}
static double forward_error(const Problem*P,const double*q){double m=0,s=1.0+infnorm(P->qstar,NDOF);for(int i=0;i<NDOF;i++){double e=fabs(q[i]-P->qstar[i]);if(e>m)m=e;}return m/s;}

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

static SolveStats solve_problem(const Problem*P,Method M){
    SolveStats S;memset(&S,0,sizeof(S));double t0=now_sec();
    double q[NDOF]={0},r[NDOF],dx[NDOF];double rel=eval_relres(P,q,r,&S);
    double prev_dx[NDOF]={0};int have_prev=0;double omega=1.0;
    double Xhist[AA_M+1][NDOF],Fhist[AA_M+1][NDOF];int hcount=0;
    double *H=NULL;if(M==M_BROYDEN_NEWTON_ALT){H=(double*)calloc((size_t)NDOF*NDOF,sizeof(double));if(!H||!build_exact_inverse(P,q,H,&S)){free(H);H=NULL;}}
    double mu=10.0;

    for(int st=0;st<MAX_STAGES && rel>1e-12;st++){
        int phase=(M==M_NEWTON)?0:(st&1);
        int ok=1;
        if(M==M_NEWTON){S.stage_a++;if(!newton_direction(P,q,r,0,dx,&S))ok=0;else{for(int i=0;i<NDOF;i++)q[i]+=dx[i];rel=eval_relres(P,q,r,&S);}}
        else if(M==M_NEWTON_LINESEARCH_ALT){
            if(phase==0){S.stage_a++;if(!newton_direction(P,q,r,0,dx,&S))ok=0;else{for(int i=0;i<NDOF;i++)q[i]+=dx[i];rel=eval_relres(P,q,r,&S);}}
            else {S.stage_b++;S.handoff_count++;if(!newton_direction(P,q,r,0,dx,&S))ok=0;else{
                double phi0=dotv(r,r,NDOF),lam=1.0;int acc=0;
                for(int ls=0;ls<80;ls++){double qt[NDOF],rt[NDOF];for(int i=0;i<NDOF;i++)qt[i]=q[i]+lam*dx[i];double rr=eval_relres(P,qt,rt,&S);S.line_search_evals++;double phi=dotv(rt,rt,NDOF);
                    if(isfinite(rr)&&phi<=(1.0-1e-4*lam)*(1.0-1e-4*lam)*phi0){memcpy(q,qt,sizeof(q));memcpy(r,rt,sizeof(r));rel=rr;acc=1;break;}lam*=0.5;}
                if(!acc)ok=0;
            }}
        } else if(M==M_NEWTON_AITKEN_ALT){
            if(!newton_direction(P,q,r,0,dx,&S))ok=0;
            else if(phase==0){S.stage_a++;for(int i=0;i<NDOF;i++)q[i]+=dx[i];rel=eval_relres(P,q,r,&S);memcpy(prev_dx,dx,sizeof(dx));have_prev=1;}
            else {S.stage_b++;S.handoff_count++;double w=1.0;if(have_prev){double df[NDOF];for(int i=0;i<NDOF;i++)df[i]=dx[i]-prev_dx[i];double den=dotv(df,df,NDOF);if(den>1e-30){w=-omega*dotv(prev_dx,df,NDOF)/den;if(!isfinite(w))w=1.0;if(w<0.05)w=0.05;if(w>1.5)w=1.5;}}
                double qr[NDOF],rrv[NDOF],qa[NDOF],rav[NDOF];for(int i=0;i<NDOF;i++){qr[i]=q[i]+dx[i];qa[i]=q[i]+w*dx[i];}
                double rr=eval_relres(P,qr,rrv,&S),ra=eval_relres(P,qa,rav,&S);if(ra<rr){memcpy(q,qa,sizeof(q));memcpy(r,rav,sizeof(r));rel=ra;S.accel_accepts++;omega=w;}else{memcpy(q,qr,sizeof(q));memcpy(r,rrv,sizeof(r));rel=rr;omega=1.0;}memcpy(prev_dx,dx,sizeof(dx));have_prev=1;}
        } else if(M==M_NEWTON_ANDERSON_ALT){
            if(!newton_direction(P,q,r,0,dx,&S))ok=0;else{
                if(hcount<AA_M+1){memcpy(Xhist[hcount],q,sizeof(q));memcpy(Fhist[hcount],dx,sizeof(dx));hcount++;}
                else{for(int h=0;h<AA_M;h++){memcpy(Xhist[h],Xhist[h+1],sizeof(Xhist[h]));memcpy(Fhist[h],Fhist[h+1],sizeof(Fhist[h]));}memcpy(Xhist[AA_M],q,sizeof(q));memcpy(Fhist[AA_M],dx,sizeof(dx));}
                if(phase==0){S.stage_a++;for(int i=0;i<NDOF;i++)q[i]+=dx[i];rel=eval_relres(P,q,r,&S);}
                else {S.stage_b++;S.handoff_count++;double qr[NDOF],rrv[NDOF];for(int i=0;i<NDOF;i++)qr[i]=q[i]+dx[i];double rr=eval_relres(P,qr,rrv,&S);int used=0;int m=(hcount-1<AA_M)?hcount-1:AA_M;
                    if(m>0){double DF[AA_M][NDOF],DX[AA_M][NDOF];int start=hcount-m-1;if(start<0)start=0;for(int c=0;c<m;c++){int j=start+c;for(int i=0;i<NDOF;i++){DF[c][i]=Fhist[j+1][i]-Fhist[j][i];DX[c][i]=Xhist[j+1][i]-Xhist[j][i];}}
                        double A[AA_M][AA_M]={{0}},b[AA_M]={0},gam[AA_M]={0};for(int a=0;a<m;a++){b[a]=dotv(DF[a],dx,NDOF);for(int c=0;c<m;c++)A[a][c]=dotv(DF[a],DF[c],NDOF);A[a][a]+=1e-14*(1.0+A[a][a]);}
                        if(solve_small(A,b,gam,m)){double qa[NDOF];for(int i=0;i<NDOF;i++){double z=qr[i];for(int c=0;c<m;c++)z-=(DX[c][i]+DF[c][i])*gam[c];qa[i]=z;}if(finite_vec(qa,NDOF)&&infnorm(qa,NDOF)<1e100){double rav[NDOF];double ra=eval_relres(P,qa,rav,&S);if(ra<rr){memcpy(q,qa,sizeof(q));memcpy(r,rav,sizeof(r));rel=ra;used=1;S.accel_accepts++;}}}
                    }if(!used){memcpy(q,qr,sizeof(q));memcpy(r,rrv,sizeof(r));rel=rr;}}
            }
        } else if(M==M_BROYDEN_NEWTON_ALT){
            if(phase==0){S.stage_a++;if(!H){ok=0;}else{for(int i=0;i<NDOF;i++){double z=0;for(int j=0;j<NDOF;j++)z+=H[(size_t)i*NDOF+j]*r[j];dx[i]=-z;}S.linear_solves++;
                double qn[NDOF],rn[NDOF];for(int i=0;i<NDOF;i++)qn[i]=q[i]+dx[i];double rnrel=eval_relres(P,qn,rn,&S);if(!isfinite(rnrel)||infnorm(qn,NDOF)>1e100)ok=0;else{
                    double y[NDOF],Hy[NDOF],svec[NDOF];for(int i=0;i<NDOF;i++){svec[i]=dx[i];y[i]=rn[i]-r[i];}
                    for(int i=0;i<NDOF;i++){double z=0;for(int j=0;j<NDOF;j++)z+=H[(size_t)i*NDOF+j]*y[j];Hy[i]=z;}
                    double vTHy=dotv(svec,Hy,NDOF);if(fabs(vTHy)>1e-24&&isfinite(vTHy)){for(int i=0;i<NDOF;i++){double ui=svec[i]-Hy[i];for(int j=0;j<NDOF;j++){double vj=0;for(int k=0;k<NDOF;k++)vj+=svec[k]*H[(size_t)k*NDOF+j];H[(size_t)i*NDOF+j]+=ui*vj/vTHy;}}}
                    memcpy(q,qn,sizeof(q));memcpy(r,rn,sizeof(r));rel=rnrel;}}
            } else {S.stage_b++;S.handoff_count++;if(!newton_direction(P,q,r,0,dx,&S))ok=0;else{for(int i=0;i<NDOF;i++)q[i]+=dx[i];rel=eval_relres(P,q,r,&S);if(H&&!build_exact_inverse(P,q,H,&S))ok=0;}}
        } else if(M==M_PTC_NEWTON_ALT){
            if(phase==0){S.stage_a++;if(!newton_direction(P,q,r,mu,dx,&S))ok=0;else{double qn[NDOF],rn[NDOF];for(int i=0;i<NDOF;i++)qn[i]=q[i]+dx[i];double rr=eval_relres(P,qn,rn,&S);if(isfinite(rr)&&rr<rel){memcpy(q,qn,sizeof(q));memcpy(r,rn,sizeof(r));rel=rr;mu=fmax(1e-12,mu*0.25);}else mu=fmin(1e12,fmax(1.0,mu*10.0));}}
            else {S.stage_b++;S.handoff_count++;if(!newton_direction(P,q,r,0,dx,&S))ok=0;else{for(int i=0;i<NDOF;i++)q[i]+=dx[i];rel=eval_relres(P,q,r,&S);}}
        } else if(M==M_NRII_NEWTON_ALT){
            if(phase==0){S.stage_a++;S.nrii_attempts++;double cand[NDOF];if(nrii_reanchor_candidate(P,q,cand)){double rc[NDOF];double rr=eval_relres(P,cand,rc,&S);if(rr<rel){memcpy(q,cand,sizeof(q));memcpy(r,rc,sizeof(r));rel=rr;S.nrii_accepts++;}}}
            else {S.stage_b++;S.handoff_count++;if(!newton_direction(P,q,r,0,dx,&S))ok=0;else{for(int i=0;i<NDOF;i++)q[i]+=dx[i];rel=eval_relres(P,q,r,&S);}}
        }
        S.stages=st+1;
        if(!ok||!finite_vec(q,NDOF)||infnorm(q,NDOF)>1e100)break;
    }
    free(H);
    if(rel>1e-12)rel=eval_relres(P,q,r,&S);
    S.rel_residual=rel;S.rel_forward_error=forward_error(P,q);S.branch_ok=(S.rel_forward_error<=1e-8);S.success=(rel<=1e-12&&S.branch_ok);S.wall_us=(now_sec()-t0)*1e6;return S;
}

static void write_methods(void){FILE*f=fopen("METHODS.csv","w");fprintf(f,"method,alternation,contract\n");
    fprintf(f,"Newton,single solver,24 full Newton stages\n");
    fprintf(f,"Newton<->LineSearch,1:1,raw full Newton stage alternates with Armijo line-search Newton stage\n");
    fprintf(f,"Newton<->Aitken,1:1,raw Newton stage alternates with Aitken-relaxed Newton-map stage; accelerated candidate kept only when better than raw candidate\n");
    fprintf(f,"Newton<->Anderson,1:1,raw Newton stage alternates with depth-3 Anderson-combined Newton-map stage; accelerated candidate kept only when better than raw candidate\n");
    fprintf(f,"Broyden<->Newton,1:1,inverse-Broyden stage alternates with exact-Jacobian Newton stage; Newton stage rebuilds inverse Jacobian for the next Broyden stage\n");
    fprintf(f,"PseudoTransient<->Newton,1:1,shifted-Jacobian pseudo-transient stage alternates with exact Newton stage\n");
    fprintf(f,"NRII<->Newton,1:1,re-anchored fixed-order NRII construction alternates with exact Newton; residual audits only accept/reject NRII handoff and never correct NRII coefficients\n");fclose(f);}

int main(void){
    write_methods();FILE*fs=fopen("SPEED_SWEEP.csv","w"),*ft=fopen("SELECTED_TRACES.csv","w");if(!fs||!ft)return 2;
    fprintf(fs,"stress,method,success,branch_ok,stages,stage_a,stage_b,residual_evals,jacobian_evals,linear_solves,line_search_evals,accel_accepts,handoff_count,nrii_attempts,nrii_accepts,rel_residual,rel_forward_error,wall_us\n");
    fprintf(ft,"stress,method,success,stages,stage_a,stage_b,residual_evals,linear_solves,handoff_count,nrii_accepts,rel_residual,rel_forward_error\n");
    const int K=MAX_STRESS;double stress[K],loglo=log10(0.05),loghi=log10(1e12);for(int j=0;j<K;j++)stress[j]=pow(10.0,loglo+(loghi-loglo)*(double)j/(K-1));
    int first_fail[M_COUNT],succ[M_COUNT];double contiguous[M_COUNT],farthest[M_COUNT];for(int m=0;m<M_COUNT;m++){first_fail[m]=-1;succ[m]=0;contiguous[m]=0;farthest[m]=0;}
    for(int j=0;j<K;j++){Problem P;make_problem(&P,stress[j]);for(int m=0;m<M_COUNT;m++){SolveStats S=solve_problem(&P,(Method)m);
        fprintf(fs,"%.17g,%s,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%.17g,%.17g,%.6f\n",stress[j],METHOD_NAME[m],S.success,S.branch_ok,S.stages,S.stage_a,S.stage_b,S.residual_evals,S.jacobian_evals,S.linear_solves,S.line_search_evals,S.accel_accepts,S.handoff_count,S.nrii_attempts,S.nrii_accepts,S.rel_residual,S.rel_forward_error,S.wall_us);
        if(S.success){succ[m]++;farthest[m]=stress[j];if(first_fail[m]<0)contiguous[m]=stress[j];}else if(first_fail[m]<0)first_fail[m]=j;
        if(j==0||j==15||j==30||j==45||j==60||j==75||j==90||j==105||j==120)fprintf(ft,"%.17g,%s,%d,%d,%d,%d,%d,%d,%d,%d,%.17g,%.17g\n",stress[j],METHOD_NAME[m],S.success,S.stages,S.stage_a,S.stage_b,S.residual_evals,S.linear_solves,S.handoff_count,S.nrii_accepts,S.rel_residual,S.rel_forward_error);
    }}fclose(fs);fclose(ft);
    FILE*fe=fopen("CONVERGENCE_ENVELOPE.csv","w");fprintf(fe,"method,contiguous_success_limit,first_failed_stress,farthest_isolated_success,successful_grid_points,total_grid_points\n");for(int m=0;m<M_COUNT;m++)fprintf(fe,"%s,%.17g,%.17g,%.17g,%d,%d\n",METHOD_NAME[m],contiguous[m],first_fail[m]>=0?stress[first_fail[m]]:NAN,farthest[m],succ[m],K);fclose(fe);
    FILE*fp=fopen("MODEL_PARAMETERS.csv","w");fprintf(fp,"parameter,value\nNDOF,%d\nimplicit_integrator,backward_Euler\ndt,1.0\nq0,zero\nv0,zero\nk,0.10\nalpha,1.0\nkc,0.20\nbeta,0.02\nNRII_order,%d\nmax_solver_stages,%d\nalternation_ratio,1:1\nrelative_residual_tolerance,1e-12\nrelative_branch_error_tolerance,1e-8\nstress_min,0.05\nstress_max,1e12\nstress_grid_points,%d\n",NDOF,NRII_ORDER,MAX_STAGES,K);fclose(fp);
#ifndef NO_TIMING
    FILE*fti=fopen("TIMING.csv","w");fprintf(fti,"stress,method,repeats,successes,mean_wall_us,mean_stages,mean_residual_evals,mean_linear_solves\n");double ts[]={0.1,0.5,1,2,5};int nt=5,reps=100;for(int a=0;a<nt;a++)for(int m=0;m<M_COUNT;m++){Problem P;make_problem(&P,ts[a]);double sw=0,ss=0,sr=0,sl=0;int ok=0;for(int z=0;z<reps;z++){SolveStats S=solve_problem(&P,(Method)m);sw+=S.wall_us;ss+=S.stages;sr+=S.residual_evals;sl+=S.linear_solves;ok+=S.success;}fprintf(fti,"%.17g,%s,%d,%d,%.6f,%.6f,%.6f,%.6f\n",ts[a],METHOD_NAME[m],reps,ok,sw/reps,ss/reps,sr/reps,sl/reps);}fclose(fti);
#endif
    FILE*fj=fopen("RESULTS.json","w");fprintf(fj,"{\n  \"benchmark\": \"General Nonlinear Dynamics Newton Reciprocal Alternation Tournament\",\n  \"NDOF\": %d,\n  \"stress_points\": %d,\n  \"max_solver_stages\": %d,\n  \"alternation_ratio\": \"1:1\",\n  \"success_definition\": {\"relative_residual\": 1e-12, \"relative_branch_error\": 1e-8},\n  \"methods\": [\n",NDOF,K,MAX_STAGES);for(int m=0;m<M_COUNT;m++)fprintf(fj,"    {\"name\": \"%s\", \"contiguous_success_limit\": %.17g, \"first_failed_stress\": %.17g, \"farthest_isolated_success\": %.17g, \"successful_grid_points\": %d}%s\n",METHOD_NAME[m],contiguous[m],first_fail[m]>=0?stress[first_fail[m]]:NAN,farthest[m],succ[m],m==M_COUNT-1?"":",");fprintf(fj,"  ]\n}\n");fclose(fj);
    printf("General Nonlinear Dynamics Newton Reciprocal Alternation Tournament\n");printf("NDOF=%d stress=[0.05,1e12] points=%d max_stages=%d\n",NDOF,K,MAX_STAGES);for(int m=0;m<M_COUNT;m++)printf("%-28s contiguous<=%.9g first_fail=%.9g farthest=%.9g successes=%d/%d\n",METHOD_NAME[m],contiguous[m],first_fail[m]>=0?stress[first_fail[m]]:NAN,farthest[m],succ[m],K);return 0;
}
