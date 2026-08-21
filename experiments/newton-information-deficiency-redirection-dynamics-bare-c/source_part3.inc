        if(solve_small(A,b,gam,m)){double qa[NDOF];for(int i=0;i<NDOF;i++){double z=qr[i];for(int c=0;c<m;c++)z-=(DX[c][i]+DF[c][i])*gam[c];qa[i]=z;}if(finite_vec(qa,NDOF)&&infnorm(qa,NDOF)<1e100){double rav[NDOF];double ra=eval_relres(P,qa,rav,S);if(ra<best){best=ra;memcpy(qb,qa,sizeof(qb));memcpy(rb,rav,sizeof(rb));S->accel_accepts++;}}}
    }
    if(isfinite(best)&&best<*rel){memcpy(q,qb,sizeof(qb));memcpy(r,rb,sizeof(rb));*rel=best;return 1;}return 0;
}

static int specialist_broyden(const Problem*P,double*q,double*r,double*rel,double*H,SolveStats*S){
    if(!H)return 0;
    double dx[NDOF];
    for(int i=0;i<NDOF;i++){double z=0;for(int j=0;j<NDOF;j++)z+=H[(size_t)i*NDOF+j]*r[j];dx[i]=-z;}
    S->linear_solves++;
    double qn[NDOF],rn[NDOF];for(int i=0;i<NDOF;i++)qn[i]=q[i]+dx[i];double rr=eval_relres(P,qn,rn,S);if(!isfinite(rr)||rr>=*rel)return 0;
    double y[NDOF];for(int i=0;i<NDOF;i++)y[i]=rn[i]-r[i];update_inverse_broyden(H,dx,y);memcpy(q,qn,sizeof(qn));memcpy(r,rn,sizeof(rn));*rel=rr;return 1;
}

static int specialist_ptc(const Problem*P,double*q,double*r,double*rel,double*mu,SolveStats*S){
    double local=*mu;
    for(int trial=0;trial<8;trial++){
        double dx[NDOF];if(!newton_direction(P,q,r,local,dx,S))return 0;
        double qn[NDOF],rn[NDOF];for(int i=0;i<NDOF;i++)qn[i]=q[i]+dx[i];double rr=eval_relres(P,qn,rn,S);
        if(isfinite(rr)&&rr<*rel){memcpy(q,qn,sizeof(qn));memcpy(r,rn,sizeof(rn));*rel=rr;*mu=fmax(1e-12,local*0.25);return 1;}
        local=fmin(1e12,fmax(1.0,local*10.0));
    }
    *mu=local;return 0;
}

static int specialist_nrii(const Problem*P,double*q,double*r,double*rel,SolveStats*S){
    S->nrii_attempts++;double cand[NDOF];if(!nrii_reanchor_candidate(P,q,cand))return 0;double rc[NDOF];double rr=eval_relres(P,cand,rc,S);if(isfinite(rr)&&rr<*rel){memcpy(q,cand,sizeof(cand));memcpy(r,rc,sizeof(rc));*rel=rr;S->nrii_accepts++;return 1;}return 0;
}

static SolveStats solve_problem(const Problem*P,Method M){
    SolveStats S;memset(&S,0,sizeof(S));double t0=now_sec();
    double q[NDOF]={0},r[NDOF];double rel=eval_relres(P,q,r,&S);
    double prev_dx[NDOF]={0};int have_prev=0,newton_since_redirect=0;double omega=1.0,mu=10.0;
    double Xhist[AA_M+1][NDOF],Fhist[AA_M+1][NDOF];int hcount=0;
    double *H=NULL;if(M==M_BROYDEN_NEWTON_ALT){H=(double*)calloc((size_t)NDOF*NDOF,sizeof(double));if(H&&!build_exact_inverse(P,q,H,&S)){free(H);H=NULL;}}
    int pending_redirect=0,redirect_disabled=0;
    for(int st=0;st<MAX_STAGES && rel>1e-12;st++){
        int ok=1;
        if(M==M_NEWTON||!pending_redirect){
            S.stage_a++;S.newton_steps++;double qold[NDOF],rold[NDOF];memcpy(qold,q,sizeof(qold));memcpy(rold,r,sizeof(rold));double rel_before=rel;
            double dx[NDOF];if(!newton_direction(P,q,r,0.0,dx,&S))ok=0;else{
                if(hcount<AA_M+1){memcpy(Xhist[hcount],q,sizeof(Xhist[0]));memcpy(Fhist[hcount],dx,sizeof(Fhist[0]));hcount++;}
                else{for(int h=0;h<AA_M;h++){memcpy(Xhist[h],Xhist[h+1],sizeof(Xhist[h]));memcpy(Fhist[h],Fhist[h+1],sizeof(Fhist[h]));}memcpy(Xhist[AA_M],q,sizeof(Xhist[AA_M]));memcpy(Fhist[AA_M],dx,sizeof(Fhist[AA_M]));}
                for(int i=0;i<NDOF;i++)q[i]+=dx[i];
                rel=eval_relres(P,q,r,&S);
                if(H){double y[NDOF];for(int i=0;i<NDOF;i++)y[i]=r[i]-rold[i];update_inverse_broyden(H,dx,y);}
                memcpy(prev_dx,dx,sizeof(dx));have_prev=1;newton_since_redirect++;
                if(M!=M_NEWTON&&!redirect_disabled){int remaining=MAX_STAGES-(st+1);if(gate_requests_redirect(rel_before,rel,newton_since_redirect,remaining,&S)){pending_redirect=1;S.redirect_triggers++;}}
            }
        } else {
            S.stage_b++;S.specialist_steps++;S.handoff_count++;int accepted=0;
            if(M==M_NEWTON_LINESEARCH_ALT)accepted=specialist_linesearch(P,q,r,&rel,&S);
            else if(M==M_NEWTON_AITKEN_ALT)accepted=specialist_aitken(P,q,r,&rel,prev_dx,have_prev,&omega,&S);
            else if(M==M_NEWTON_ANDERSON_ALT)accepted=specialist_anderson(P,q,r,&rel,Xhist,Fhist,&hcount,&S);
            else if(M==M_BROYDEN_NEWTON_ALT)accepted=specialist_broyden(P,q,r,&rel,H,&S);
            else if(M==M_PTC_NEWTON_ALT)accepted=specialist_ptc(P,q,r,&rel,&mu,&S);
            else if(M==M_NRII_NEWTON_ALT)accepted=specialist_nrii(P,q,r,&rel,&S);
            if(accepted){S.redirect_successes++;S.accel_accepts++;}
            else {S.redirect_failures++;if(DISABLE_AFTER_SPECIALIST_FAILURE)redirect_disabled=1;}
            pending_redirect=0;newton_since_redirect=0;
        }
        S.stages=st+1;if(!ok||!finite_vec(q,NDOF)||infnorm(q,NDOF)>1e100)break;
    }
    free(H);if(rel>1e-12)rel=eval_relres(P,q,r,&S);S.rel_residual=rel;S.rel_forward_error=forward_error(P,q);S.branch_ok=(S.rel_forward_error<=1e-8);S.success=(rel<=1e-12&&S.branch_ok);S.wall_us=(now_sec()-t0)*1e6;return S;
}

static void write_methods(void){FILE*f=fopen("METHODS.csv","w");fprintf(f,"method,composition,gate,contract\n");
    fprintf(f,"Newton,standalone,none,raw Newton baseline\n");
    fprintf(f,"Newton->LineSearch->Newton,triggered redirect,budget-projected Newton information gate,Armijo line-search Newton specialist; returns to Newton after one specialist dispatch\n");
    fprintf(f,"Newton->Aitken->Newton,triggered redirect,budget-projected Newton information gate,Aitken-relaxed Newton-map specialist; better of raw and relaxed candidate is retained\n");
    fprintf(f,"Newton->Anderson->Newton,triggered redirect,budget-projected Newton information gate,depth-3 Anderson Newton-map specialist; better of raw and accelerated candidate is retained\n");
    fprintf(f,"Newton->Broyden->Newton,triggered redirect,budget-projected Newton information gate,inverse-Broyden specialist using secant-updated inverse information accumulated along Newton states\n");
    fprintf(f,"Newton->PseudoTransient->Newton,triggered redirect,budget-projected Newton information gate,shifted-Jacobian pseudo-transient specialist with bounded shift search\n");
    fprintf(f,"Newton->NRII->Newton,triggered redirect,budget-projected Newton information gate,re-anchored order-12 NRII specialist; residual is accept/reject audit only and never corrects NRII coefficients\n");fclose(f);}

int main(void){
    write_methods();FILE*fs=fopen("SPEED_SWEEP.csv","w"),*ft=fopen("SELECTED_TRACES.csv","w");if(!fs||!ft)return 2;
    fprintf(fs,"stress,method,success,branch_ok,stages,newton_steps,specialist_steps,redirect_triggers,redirect_successes,redirect_failures,trigger_budget_deficit,trigger_nondecrease,trigger_invalid,residual_evals,jacobian_evals,linear_solves,line_search_evals,accel_accepts,nrii_attempts,nrii_accepts,rel_residual,rel_forward_error,wall_us\n");
    fprintf(ft,"stress,method,success,stages,newton_steps,specialist_steps,redirect_triggers,redirect_successes,redirect_failures,rel_residual,rel_forward_error\n");
    const int K=MAX_STRESS;double stress[K],loglo=log10(0.05),loghi=log10(1e12);for(int j=0;j<K;j++)stress[j]=pow(10.0,loglo+(loghi-loglo)*(double)j/(K-1));
    int first_fail[M_COUNT],succ[M_COUNT],red[M_COUNT],redok[M_COUNT];double contiguous[M_COUNT],farthest[M_COUNT];for(int m=0;m<M_COUNT;m++){first_fail[m]=-1;succ[m]=0;red[m]=redok[m]=0;contiguous[m]=0;farthest[m]=0;}
    for(int j=0;j<K;j++){Problem P;make_problem(&P,stress[j]);for(int m=0;m<M_COUNT;m++){SolveStats S=solve_problem(&P,(Method)m);red[m]+=S.redirect_triggers;redok[m]+=S.redirect_successes;
        fprintf(fs,"%.17g,%s,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%.17g,%.17g,%.6f\n",stress[j],METHOD_NAME[m],S.success,S.branch_ok,S.stages,S.newton_steps,S.specialist_steps,S.redirect_triggers,S.redirect_successes,S.redirect_failures,S.trigger_budget_deficit,S.trigger_nondecrease,S.trigger_invalid,S.residual_evals,S.jacobian_evals,S.linear_solves,S.line_search_evals,S.accel_accepts,S.nrii_attempts,S.nrii_accepts,S.rel_residual,S.rel_forward_error,S.wall_us);
        if(S.success){succ[m]++;farthest[m]=stress[j];if(first_fail[m]<0)contiguous[m]=stress[j];}else if(first_fail[m]<0)first_fail[m]=j;
        if(j==0||j==15||j==30||j==45||j==60||j==75||j==90||j==105||j==120)fprintf(ft,"%.17g,%s,%d,%d,%d,%d,%d,%d,%d,%.17g,%.17g\n",stress[j],METHOD_NAME[m],S.success,S.stages,S.newton_steps,S.specialist_steps,S.redirect_triggers,S.redirect_successes,S.redirect_failures,S.rel_residual,S.rel_forward_error);
    }}fclose(fs);fclose(ft);
    FILE*fe=fopen("CONVERGENCE_ENVELOPE.csv","w");fprintf(fe,"method,contiguous_success_limit,first_failed_stress,farthest_isolated_success,successful_grid_points,total_grid_points,total_redirects,successful_redirects,redirect_success_fraction\n");for(int m=0;m<M_COUNT;m++)fprintf(fe,"%s,%.17g,%.17g,%.17g,%d,%d,%d,%d,%.17g\n",METHOD_NAME[m],contiguous[m],first_fail[m]>=0?stress[first_fail[m]]:NAN,farthest[m],succ[m],K,red[m],redok[m],red[m]?((double)redok[m]/red[m]):0.0);fclose(fe);
    FILE*fp=fopen("MODEL_PARAMETERS.csv","w");fprintf(fp,"parameter,value\nNDOF,%d\nimplicit_integrator,backward_Euler\ndt,1.0\nq0,zero\nv0,zero\nk,0.10\nalpha,1.0\nkc,0.20\nbeta,0.02\nNRII_order,%d\nmax_dispatch_stages,%d\nminimum_Newton_stages_between_redirects,2\nredirect_gate,projected_remaining_Newton_steps_exceed_remaining_dispatch_budget_or_nondecreasing_progress\nrelative_residual_tolerance,1e-12\nrelative_branch_error_tolerance,1e-8\nstress_min,0.05\nstress_max,1e12\nstress_grid_points,%d\n",NDOF,NRII_ORDER,MAX_STAGES,K);fclose(fp);
#ifndef NO_TIMING
    FILE*fti=fopen("TIMING.csv","w");fprintf(fti,"stress,method,repeats,successes,mean_wall_us,mean_stages,mean_redirects,mean_residual_evals,mean_linear_solves\n");double ts[]={0.1,0.5,1,2,5};int nt=5,reps=100;for(int a=0;a<nt;a++)for(int m=0;m<M_COUNT;m++){Problem P;make_problem(&P,ts[a]);double sw=0,ss=0,sd=0,sr=0,sl=0;int ok=0;for(int z=0;z<reps;z++){SolveStats S=solve_problem(&P,(Method)m);sw+=S.wall_us;ss+=S.stages;sd+=S.redirect_triggers;sr+=S.residual_evals;sl+=S.linear_solves;ok+=S.success;}fprintf(fti,"%.17g,%s,%d,%d,%.6f,%.6f,%.6f,%.6f,%.6f\n",ts[a],METHOD_NAME[m],reps,ok,sw/reps,ss/reps,sd/reps,sr/reps,sl/reps);}fclose(fti);
#endif
    FILE*fj=fopen("RESULTS.json","w");fprintf(fj,"{\n  \"benchmark\": \"General Nonlinear Dynamics Newton Information-Deficiency Redirect Tournament\",\n  \"NDOF\": %d,\n  \"stress_points\": %d,\n  \"max_dispatch_stages\": %d,\n  \"gate\": \"after at least two Newton stages, redirect if the latest Newton contraction rate predicts more Newton stages than the declared remaining dispatch budget, or if progress is nondecreasing\",\n  \"success_definition\": {\"relative_residual\": 1e-12, \"relative_branch_error\": 1e-8},\n  \"methods\": [\n",NDOF,K,MAX_STAGES);for(int m=0;m<M_COUNT;m++)fprintf(fj,"    {\"name\": \"%s\", \"contiguous_success_limit\": %.17g, \"first_failed_stress\": %.17g, \"farthest_isolated_success\": %.17g, \"successful_grid_points\": %d, \"redirects\": %d, \"successful_redirects\": %d}%s\n",METHOD_NAME[m],contiguous[m],first_fail[m]>=0?stress[first_fail[m]]:NAN,farthest[m],succ[m],red[m],redok[m],m==M_COUNT-1?"":",");fprintf(fj,"  ]\n}\n");fclose(fj);
    printf("General Nonlinear Dynamics Newton Information-Deficiency Redirect Tournament\n");printf("NDOF=%d stress=[0.05,1e12] points=%d max_dispatch_stages=%d\n",NDOF,K,MAX_STAGES);for(int m=0;m<M_COUNT;m++)printf("%-34s contiguous<=%.9g first_fail=%.9g farthest=%.9g successes=%d/%d redirects=%d success_redirects=%d\n",METHOD_NAME[m],contiguous[m],first_fail[m]>=0?stress[first_fail[m]]:NAN,farthest[m],succ[m],K,red[m],redok[m]);return 0;
}
