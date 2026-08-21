            else {
                double qn[NDOF],rn[NDOF];for(int i=0;i<NDOF;i++)qn[i]=q[i]+dx[i];
                double reln=eval_relres(P,qn,rn,&S);
                if(isfinite(reln)&&reln<rel){memcpy(q,qn,sizeof(q));memcpy(r,rn,sizeof(r));rel=reln;mu*=0.25;if(mu<1e-10){mu=0.0;S.handoff_used=1;}}
                else {mu=fmin(1e12,fmax(1.0,mu*10.0));}
            }
        } else {
            if(!newton_direction(P,q,r,0.0,dx,&S)){step_ok=0;}
            else if(M==M_NEWTON || M==M_NRII_NEWTON || M==M_BROYDEN_NEWTON){
                for(int i=0;i<NDOF;i++){ q[i]+=dx[i]; }
                rel=eval_relres(P,q,r,&S);
            } else if(M==M_LINESEARCH){
                double phi0=dotv(r,r,NDOF);double lam=1.0;int accepted=0;
                for(int ls=0;ls<80;ls++){
                    double qt[NDOF],rt[NDOF];for(int i=0;i<NDOF;i++)qt[i]=q[i]+lam*dx[i];
                    double relt=eval_relres(P,qt,rt,&S);S.line_search_evals++;
                    double phit=dotv(rt,rt,NDOF);
                    if(isfinite(relt)&&phit <= (1.0-1e-4*lam)*(1.0-1e-4*lam)*phi0){memcpy(q,qt,sizeof(q));memcpy(r,rt,sizeof(r));rel=relt;accepted=1;if(lam<0.999999)S.handoff_used=1;break;}
                    lam*=0.5;
                }
                if(!accepted)step_ok=0;
            } else if(M==M_AITKEN){
                double w=1.0;
                if(have_prev){double df[NDOF];for(int i=0;i<NDOF;i++)df[i]=dx[i]-prev_dx[i];double den=dotv(df,df,NDOF);if(den>1e-30){w=-omega*dotv(prev_dx,df,NDOF)/den;if(!isfinite(w))w=1.0;if(w<0.05)w=0.05;if(w>1.5)w=1.5;}}
                double qraw[NDOF],rraw[NDOF],qacc[NDOF],racc[NDOF];
                for(int i=0;i<NDOF;i++){qraw[i]=q[i]+dx[i];qacc[i]=q[i]+w*dx[i];}
                double rr=eval_relres(P,qraw,rraw,&S),ra=eval_relres(P,qacc,racc,&S);
                if(ra<rr){memcpy(q,qacc,sizeof(q));memcpy(r,racc,sizeof(r));rel=ra;omega=w;S.accel_accepts++;S.handoff_used=1;}
                else {memcpy(q,qraw,sizeof(q));memcpy(r,rraw,sizeof(r));rel=rr;omega=1.0;}
                memcpy(prev_dx,dx,sizeof(dx));have_prev=1;
            } else if(M==M_ANDERSON){
                double qraw[NDOF],rraw[NDOF];for(int i=0;i<NDOF;i++)qraw[i]=q[i]+dx[i];double rr=eval_relres(P,qraw,rraw,&S);
                /* append current x and Newton fixed-point residual f=dx */
                if(hcount<AA_M+1){memcpy(Xhist[hcount],q,sizeof(q));memcpy(Fhist[hcount],dx,sizeof(dx));hcount++;}
                else {for(int h=0;h<AA_M;h++){memcpy(Xhist[h],Xhist[h+1],sizeof(Xhist[h]));memcpy(Fhist[h],Fhist[h+1],sizeof(Fhist[h]));}memcpy(Xhist[AA_M],q,sizeof(q));memcpy(Fhist[AA_M],dx,sizeof(dx));}
                int m=(hcount-1<AA_M)?hcount-1:AA_M;int used=0;
                if(m>0){
                    double DF[AA_M][NDOF],DX[AA_M][NDOF];
                    int start=hcount-m-1;if(start<0)start=0;
                    for(int c=0;c<m;c++){int j=start+c;for(int i=0;i<NDOF;i++){DF[c][i]=Fhist[j+1][i]-Fhist[j][i];DX[c][i]=Xhist[j+1][i]-Xhist[j][i];}}
                    double A[AA_M][AA_M]={{0}},b[AA_M]={0},gam[AA_M]={0};
                    for(int a=0;a<m;a++){b[a]=dotv(DF[a],dx,NDOF);for(int c=0;c<m;c++)A[a][c]=dotv(DF[a],DF[c],NDOF);A[a][a]+=1e-14*(1.0+A[a][a]);}
                    if(solve_small(A,b,gam,m)){
                        double qa[NDOF];for(int i=0;i<NDOF;i++){double z=qraw[i];for(int c=0;c<m;c++)z-=(DX[c][i]+DF[c][i])*gam[c];qa[i]=z;}
                        if(finite_vec(qa,NDOF)&&infnorm(qa,NDOF)<1e100){double ra_vec[NDOF];double ra=eval_relres(P,qa,ra_vec,&S);if(ra<rr){memcpy(q,qa,sizeof(q));memcpy(r,ra_vec,sizeof(r));rel=ra;used=1;S.accel_accepts++;S.handoff_used=1;}}
                    }
                }
                if(!used){memcpy(q,qraw,sizeof(q));memcpy(r,rraw,sizeof(r));rel=rr;}
            }
        }
        if(!step_ok||!finite_vec(q,NDOF)||infnorm(q,NDOF)>1e100){S.iterations=it+1;break;}
        S.iterations=it+1;
    }
    free(H);
    if(rel>1e-12) rel=eval_relres(P,q,r,&S);
    S.rel_residual=rel;S.rel_forward_error=relative_forward_error(P,q);
    S.branch_ok=(S.rel_forward_error<=1e-8);
    S.success=(S.rel_residual<=1e-12 && S.branch_ok);
    S.wall_us=(now_sec()-t0)*1e6;
    return S;
}

static void write_methods(void){
    FILE*f=fopen("METHODS.csv","w");fprintf(f,"method,role,contract\n");
    fprintf(f,"Newton,baseline,full-step exact-Jacobian Newton\n");
    fprintf(f,"Newton+LineSearch,globalized Newton,Armijo residual line search\n");
    fprintf(f,"Newton+Aitken,accelerated Newton,Aitken dynamic relaxation candidate accepted only when residual improves over raw Newton candidate\n");
    fprintf(f,"Newton+Anderson,accelerated Newton,depth-3 Anderson candidate accepted only when residual improves over raw Newton candidate\n");
    fprintf(f,"Broyden->Newton,quasi-Newton relay,dense inverse-Broyden phase then exact Newton terminal phase\n");
    fprintf(f,"PseudoTransient->Newton,continuation relay,regularized pseudo-transient Newton phase then exact Newton when regularization vanishes\n");
    fprintf(f,"NRII->Newton,independent constructive relay,fixed-order NRII implicit-Euler coefficient candidate; residual audit gates handoff; Newton is terminal corrector\n");
    fclose(f);
}

int main(void){
    write_methods();
    FILE *fs=fopen("SPEED_SWEEP.csv","w");
    FILE *ft=fopen("SELECTED_TRACES.csv","w");
    if(!fs||!ft){fprintf(stderr,"file open failed\n");return 2;}
    fprintf(fs,"stress,method,success,branch_ok,iterations,residual_evals,jacobian_evals,linear_solves,line_search_evals,accel_accepts,handoff_used,nrii_accepted,rel_residual,rel_forward_error,wall_us\n");
    fprintf(ft,"stress,method,success,iterations,residual_evals,linear_solves,rel_residual,rel_forward_error,handoff_used\n");
    const int K=MAX_STRESS; double stress[K];
    double loglo=log10(0.05),loghi=log10(1e12);
    for(int j=0;j<K;j++)stress[j]=pow(10.0,loglo+(loghi-loglo)*(double)j/(double)(K-1));
    int first_fail[M_COUNT];double contiguous[M_COUNT],farthest[M_COUNT];int succ_count[M_COUNT];
    for(int m=0;m<M_COUNT;m++){first_fail[m]=-1;contiguous[m]=0.0;farthest[m]=0.0;succ_count[m]=0;}
    for(int j=0;j<K;j++){
        Problem P;make_problem(&P,stress[j]);
        for(int m=0;m<M_COUNT;m++){
            SolveStats S=solve_problem(&P,(Method)m);
            fprintf(fs,"%.17g,%s,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%.17g,%.17g,%.6f\n",stress[j],METHOD_NAME[m],S.success,S.branch_ok,S.iterations,S.residual_evals,S.jacobian_evals,S.linear_solves,S.line_search_evals,S.accel_accepts,S.handoff_used,S.nrii_accepted,S.rel_residual,S.rel_forward_error,S.wall_us);
            if(S.success){succ_count[m]++;farthest[m]=stress[j];if(first_fail[m]<0)contiguous[m]=stress[j];}
            else if(first_fail[m]<0){first_fail[m]=j;}
            if(j==0||j==15||j==30||j==45||j==60||j==75||j==90||j==105||j==120)fprintf(ft,"%.17g,%s,%d,%d,%d,%d,%.17g,%.17g,%d\n",stress[j],METHOD_NAME[m],S.success,S.iterations,S.residual_evals,S.linear_solves,S.rel_residual,S.rel_forward_error,S.handoff_used);
        }
    }
    fclose(fs);fclose(ft);
    FILE *fe=fopen("CONVERGENCE_ENVELOPE.csv","w");
    fprintf(fe,"method,contiguous_success_limit,first_failed_stress,farthest_isolated_success,successful_grid_points,total_grid_points\n");
    for(int m=0;m<M_COUNT;m++){
        double ff=(first_fail[m]>=0)?stress[first_fail[m]]:NAN;
        fprintf(fe,"%s,%.17g,%.17g,%.17g,%d,%d\n",METHOD_NAME[m],contiguous[m],ff,farthest[m],succ_count[m],K);
    }
    fclose(fe);

    /* Timing on common easy/moderate points. Counts are primary; this is secondary same-process timing. */
    FILE *fti=fopen("TIMING.csv","w");fprintf(fti,"stress,method,repeats,successes,mean_wall_us,mean_iterations,mean_residual_evals,mean_linear_solves\n");
    double tstress[]={0.10,0.50,1.0,2.0,5.0};int nt=5,reps=100;
    for(int a=0;a<nt;a++)for(int m=0;m<M_COUNT;m++){
        Problem P;make_problem(&P,tstress[a]);double sw=0,si=0,sr=0,sl=0;int oks=0;
        for(int r=0;r<reps;r++){SolveStats S=solve_problem(&P,(Method)m);sw+=S.wall_us;si+=S.iterations;sr+=S.residual_evals;sl+=S.linear_solves;oks+=S.success;}
        fprintf(fti,"%.17g,%s,%d,%d,%.6f,%.6f,%.6f,%.6f\n",tstress[a],METHOD_NAME[m],reps,oks,sw/reps,si/reps,sr/reps,sl/reps);
    }fclose(fti);

    FILE *fp=fopen("MODEL_PARAMETERS.csv","w");
    fprintf(fp,"parameter,value\nNDOF,%d\nimplicit_integrator,backward_Euler\ndt,1.0\nq0,zero\nv0,zero\ninternal_foundation,k*q+alpha*q^3\nk,0.10\nalpha,1.0\nedge_force,kc*d+beta*d^3\nkc,0.20\nbeta,0.02\nNRII_order,%d\nmax_nonlinear_iterations,%d\nrelative_residual_tolerance,1e-12\nrelative_branch_error_tolerance,1e-8\nstress_min,0.05\nstress_max,1e12\nstress_grid_points,%d\n",NDOF,NRII_ORDER,MAX_IT,K);fclose(fp);

    FILE *fj=fopen("RESULTS.json","w");
    fprintf(fj,"{\n  \"benchmark\": \"General Nonlinear Dynamics Newton Hybrid Tournament\",\n  \"NDOF\": %d,\n  \"stress_points\": %d,\n  \"max_iterations\": %d,\n  \"success_definition\": {\"relative_residual\": 1e-12, \"relative_branch_error\": 1e-8},\n  \"methods\": [\n",NDOF,K,MAX_IT);
    for(int m=0;m<M_COUNT;m++)fprintf(fj,"    {\"name\": \"%s\", \"contiguous_success_limit\": %.17g, \"first_failed_stress\": %.17g, \"farthest_isolated_success\": %.17g, \"successful_grid_points\": %d}%s\n",METHOD_NAME[m],contiguous[m],(first_fail[m]>=0?stress[first_fail[m]]:NAN),farthest[m],succ_count[m],m==M_COUNT-1?"":",");
    fprintf(fj,"  ]\n}\n");fclose(fj);

    printf("General Nonlinear Dynamics Newton Hybrid Tournament\n");
    printf("NDOF=%d, stress=[0.05,1e12], points=%d, max_it=%d\n",NDOF,K,MAX_IT);
    for(int m=0;m<M_COUNT;m++)printf("%-24s contiguous<=%.6g first_fail=%.6g farthest_success=%.6g success_points=%d/%d\n",METHOD_NAME[m],contiguous[m],first_fail[m]>=0?stress[first_fail[m]]:NAN,farthest[m],succ_count[m],K);
    return 0;
}
