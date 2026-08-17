#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884
#endif

#define NV 6
#define ORDER 16

typedef struct {
    double fB_hz;
    double fA_hz;
    double fN_hz;
    double wB, wA, wN;
    double B0_T;
    double A0_m2;
    double epsA;
    double dt;
    double T;
    int order;
} Params;

typedef struct {
    double max_midpoint_resid;
    double max_state_err_midpoint_oracle;
    double max_state_err_continuous;
    double max_flux_err_continuous;
    double max_emf_err_continuous;
    double max_area_rel_err_continuous;
    double max_B_vec_err_continuous;
    double max_S_vec_err_continuous;
    double max_flux_abs;
    double max_emf_abs;
    long long steps;
    long long bx_zero_crossings;
    long long area_dev_zero_crossings;
    long long flux_zero_crossings;
    long long emf_zero_crossings;
    double elapsed_s;
    double dft_flux_amp[3];
    double dft_emf_amp[3];
    double dft_flux_phase[3];
    double dft_emf_phase[3];
} Metrics;

static inline void apply_A(const Params *p, const double x[NV], double y[NV]) {
    y[0] = -p->wB * x[1];
    y[1] =  p->wB * x[0];
    y[2] = -p->wA * x[3];
    y[3] =  p->wA * x[2];
    y[4] = -p->wN * x[5];
    y[5] =  p->wN * x[4];
}

static void nrii_midpoint_step(const Params *p, const double y0[NV], double y1[NV]) {
    double term[NV], next[NV], tmp[NV];
    apply_A(p, y0, tmp);
    for (int i=0;i<NV;i++) {
        term[i] = p->dt * tmp[i];
        y1[i] = y0[i] + term[i];
    }
    for (int k=2;k<=p->order;k++) {
        apply_A(p, term, tmp);
        for (int i=0;i<NV;i++) {
            next[i] = 0.5 * p->dt * tmp[i];
            y1[i] += next[i];
        }
        memcpy(term, next, sizeof(term));
    }
}

static inline void midpoint_pair(double w, double dt, double c0, double s0, double *c1, double *s1) {
    const double a = 0.5*w*dt;
    const double den = 1.0 + a*a;
    const double d = 1.0 - a*a;
    *c1 = (d*c0 - 2.0*a*s0)/den;
    *s1 = (2.0*a*c0 + d*s0)/den;
}

static void midpoint_oracle_step(const Params *p, const double y0[NV], double y1[NV]) {
    midpoint_pair(p->wB,p->dt,y0[0],y0[1],&y1[0],&y1[1]);
    midpoint_pair(p->wA,p->dt,y0[2],y0[3],&y1[2],&y1[3]);
    midpoint_pair(p->wN,p->dt,y0[4],y0[5],&y1[4],&y1[5]);
}

static double midpoint_residual_inf(const Params *p, const double y0[NV], const double y1[NV]) {
    double mid[NV], Am[NV], m=0.0;
    for (int i=0;i<NV;i++) mid[i]=0.5*(y0[i]+y1[i]);
    apply_A(p,mid,Am);
    for (int i=0;i<NV;i++) {
        double r = y1[i]-y0[i]-p->dt*Am[i];
        double a=fabs(r); if(a>m)m=a;
    }
    return m;
}

static inline double maxabsdiff6(const double a[NV], const double b[NV]) {
    double m=0.0;
    for(int i=0;i<NV;i++){double d=fabs(a[i]-b[i]); if(d>m)m=d;}
    return m;
}

static void exact_state(const Params *p, double t, double y[NV]) {
    y[0]=cos(p->wB*t); y[1]=sin(p->wB*t);
    y[2]=cos(p->wA*t); y[3]=sin(p->wA*t);
    y[4]=cos(p->wN*t); y[5]=sin(p->wN*t);
}

static inline void observables(const Params *p, const double y[NV],
                               double *area, double Bv[2], double Sv[2],
                               double *flux, double *emf) {
    const double af = 1.0 + p->epsA*y[2];
    const double A = p->A0_m2*af;
    const double dot = y[0]*y[4] + y[1]*y[5];
    const double ddot = (p->wB-p->wN)*(y[0]*y[5]-y[1]*y[4]);
    const double afd = -p->epsA*p->wA*y[3];
    Bv[0]=p->B0_T*y[0]; Bv[1]=p->B0_T*y[1];
    Sv[0]=A*y[4]; Sv[1]=A*y[5];
    *area=A;
    *flux=p->B0_T*A*dot;
    *emf=-p->B0_T*p->A0_m2*(afd*dot + af*ddot);
}

static inline void exact_observables(const Params *p,double t,double *area,double Bv[2],double Sv[2],double *flux,double *emf){
    const double cB=cos(p->wB*t), sB=sin(p->wB*t);
    const double cA=cos(p->wA*t), sA=sin(p->wA*t);
    const double cN=cos(p->wN*t), sN=sin(p->wN*t);
    const double af=1.0+p->epsA*cA;
    const double A=p->A0_m2*af;
    const double d=(p->wB-p->wN)*t;
    Bv[0]=p->B0_T*cB; Bv[1]=p->B0_T*sB;
    Sv[0]=A*cN; Sv[1]=A*sN;
    *area=A;
    *flux=p->B0_T*p->A0_m2*af*cos(d);
    *emf=p->B0_T*p->A0_m2*(p->epsA*p->wA*sA*cos(d) + (p->wB-p->wN)*af*sin(d));
}

static inline int crossed(double a,double b){ return (a<0.0 && b>=0.0)||(a>0.0 && b<=0.0); }
static inline double hypot2diff(const double a[2],const double b[2]){return hypot(a[0]-b[0],a[1]-b[1]);}

static int run_case(const Params *p, Metrics *m, const char *trace_path, int trace_stride) {
    memset(m,0,sizeof(*m));
    const long long steps=(long long)llround(p->T/p->dt);
    m->steps=steps;
    double y[NV]={1,0,1,0,1,0};
    double ym[NV]={1,0,1,0,1,0};
    double yn[NV], ymn[NV], ye[NV];
    FILE *tf=NULL;
    if(trace_path){
        tf=fopen(trace_path,"w"); if(!tf)return 0;
        fprintf(tf,"step,time_s,Bx_T,By_T,area_m2,Sx_m2,Sy_m2,flux_Wb,emf_V,midpoint_resid,state_err_midpoint_oracle,state_err_continuous,flux_err_continuous,emf_err_continuous\n");
    }
    const double freqs[3]={fabs(p->fB_hz-p->fN_hz), fabs(p->fA_hz-fabs(p->fB_hz-p->fN_hz)), p->fA_hz+fabs(p->fB_hz-p->fN_hz)};
    double Fr[3]={0},Fi[3]={0},Er[3]={0},Ei[3]={0};
    double area0,B0v[2],S0v[2],flux0,emf0;
    observables(p,y,&area0,B0v,S0v,&flux0,&emf0);
    double prev_bx=B0v[0], prev_adev=area0-p->A0_m2, prev_flux=flux0, prev_emf=emf0;

    clock_t c0=clock();
    for(long long step=1;step<=steps;step++){
        nrii_midpoint_step(p,y,yn);
        midpoint_oracle_step(p,ym,ymn);
        const double resid=midpoint_residual_inf(p,y,yn);
        const double e_mid=maxabsdiff6(yn,ymn);
        const double t=(double)step*p->dt;
        exact_state(p,t,ye);
        const double e_cont=maxabsdiff6(yn,ye);
        double area,Bv[2],Sv[2],flux,emf;
        double areae,Be[2],Se[2],fluxe,emfe;
        observables(p,yn,&area,Bv,Sv,&flux,&emf);
        exact_observables(p,t,&areae,Be,Se,&fluxe,&emfe);
        const double fe=fabs(flux-fluxe), ee=fabs(emf-emfe);
        const double are=fabs(area-areae)/p->A0_m2;
        const double bve=hypot2diff(Bv,Be);
        const double sve=hypot2diff(Sv,Se);
        if(resid>m->max_midpoint_resid)m->max_midpoint_resid=resid;
        if(e_mid>m->max_state_err_midpoint_oracle)m->max_state_err_midpoint_oracle=e_mid;
        if(e_cont>m->max_state_err_continuous)m->max_state_err_continuous=e_cont;
        if(fe>m->max_flux_err_continuous)m->max_flux_err_continuous=fe;
        if(ee>m->max_emf_err_continuous)m->max_emf_err_continuous=ee;
        if(are>m->max_area_rel_err_continuous)m->max_area_rel_err_continuous=are;
        if(bve>m->max_B_vec_err_continuous)m->max_B_vec_err_continuous=bve;
        if(sve>m->max_S_vec_err_continuous)m->max_S_vec_err_continuous=sve;
        if(fabs(flux)>m->max_flux_abs)m->max_flux_abs=fabs(flux);
        if(fabs(emf)>m->max_emf_abs)m->max_emf_abs=fabs(emf);
        if(crossed(prev_bx,Bv[0]))m->bx_zero_crossings++;
        if(crossed(prev_adev,area-p->A0_m2))m->area_dev_zero_crossings++;
        if(crossed(prev_flux,flux))m->flux_zero_crossings++;
        if(crossed(prev_emf,emf))m->emf_zero_crossings++;
        prev_bx=Bv[0]; prev_adev=area-p->A0_m2; prev_flux=flux; prev_emf=emf;
        for(int j=0;j<3;j++){
            const double a=2.0*M_PI*freqs[j]*t;
            const double c=cos(a), s=sin(a);
            Fr[j]+=flux*c; Fi[j]-=flux*s;
            Er[j]+=emf*c; Ei[j]-=emf*s;
        }
        if(tf && (step==1 || step==steps || step%trace_stride==0)){
            fprintf(tf,"%lld,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                step,t,Bv[0],Bv[1],area,Sv[0],Sv[1],flux,emf,resid,e_mid,e_cont,fe,ee);
        }
        memcpy(y,yn,sizeof(y)); memcpy(ym,ymn,sizeof(ym));
    }
    clock_t c1=clock();
    m->elapsed_s=(double)(c1-c0)/(double)CLOCKS_PER_SEC;
    const double scale=2.0/(double)steps;
    for(int j=0;j<3;j++){
        m->dft_flux_amp[j]=scale*hypot(Fr[j],Fi[j]);
        m->dft_emf_amp[j]=scale*hypot(Er[j],Ei[j]);
        m->dft_flux_phase[j]=atan2(Fi[j],Fr[j]);
        m->dft_emf_phase[j]=atan2(Ei[j],Er[j]);
    }
    if(tf)fclose(tf);
    return 1;
}

static void write_model(const Params *p,const char *path){
    FILE *f=fopen(path,"w"); if(!f)return;
    fprintf(f,"parameter,value,unit\n");
    fprintf(f,"field_vector_frequency,%.17g,Hz\n",p->fB_hz);
    fprintf(f,"area_oscillation_frequency,%.17g,Hz\n",p->fA_hz);
    fprintf(f,"surface_normal_rotation_frequency,%.17g,Hz\n",p->fN_hz);
    fprintf(f,"magnetic_flux_density_amplitude,%.17g,T\n",p->B0_T);
    fprintf(f,"mean_surface_area,%.17g,m^2\n",p->A0_m2);
    fprintf(f,"relative_area_modulation,%.17g,1\n",p->epsA);
    fprintf(f,"time_step,%.17g,s\n",p->dt);
    fprintf(f,"duration,%.17g,s\n",p->T);
    fprintf(f,"NRII_order,%d,1\n",p->order);
    fclose(f);
}

static void write_spectrum(const Params *p,const Metrics *m,const char *path){
    const double fd=fabs(p->fB_hz-p->fN_hz);
    const double f[3]={fd,fabs(p->fA_hz-fd),p->fA_hz+fd};
    const double flux_theory[3]={p->B0_T*p->A0_m2,0.5*p->epsA*p->B0_T*p->A0_m2,0.5*p->epsA*p->B0_T*p->A0_m2};
    FILE *o=fopen(path,"w"); if(!o)return;
    fprintf(o,"component,frequency_Hz,NRII_flux_amp_Wb,analytic_flux_amp_Wb,flux_rel_amp_error,NRII_emf_amp_V,analytic_emf_amp_V,emf_rel_amp_error\n");
    const char *name[3]={"relative_rotation_carrier","lower_area_sideband","upper_area_sideband"};
    for(int j=0;j<3;j++){
        const double ea=2.0*M_PI*f[j]*flux_theory[j];
        fprintf(o,"%s,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",name[j],f[j],m->dft_flux_amp[j],flux_theory[j],fabs(m->dft_flux_amp[j]-flux_theory[j])/flux_theory[j],m->dft_emf_amp[j],ea,fabs(m->dft_emf_amp[j]-ea)/ea);
    }
    fclose(o);
}

static void write_results(const Params *p,const Metrics *m,const char *path){
    const double fd=fabs(p->fB_hz-p->fN_hz);
    FILE *o=fopen(path,"w"); if(!o)return;
    fprintf(o,"{\n");
    fprintf(o,"  \"benchmark\": \"NRII oscillating-field moving-area flux benchmark, bare C17\",\n");
    fprintf(o,"  \"status\": \"literature-grounded combined stress benchmark; not a reproduction of a single published numerical case\",\n");
    fprintf(o,"  \"parameters\": {\"fB_Hz\": %.17g, \"fA_Hz\": %.17g, \"fN_Hz\": %.17g, \"B0_T\": %.17g, \"A0_m2\": %.17g, \"epsA\": %.17g, \"dt_s\": %.17g, \"T_s\": %.17g, \"NRII_order\": %d},\n",p->fB_hz,p->fA_hz,p->fN_hz,p->B0_T,p->A0_m2,p->epsA,p->dt,p->T,p->order);
    fprintf(o,"  \"derived_frequencies_Hz\": {\"relative_rotation\": %.17g, \"lower_area_sideband\": %.17g, \"upper_area_sideband\": %.17g},\n",fd,fabs(p->fA_hz-fd),p->fA_hz+fd);
    fprintf(o,"  \"metrics\": {\n");
    fprintf(o,"    \"steps\": %lld,\n",m->steps);
    fprintf(o,"    \"elapsed_s\": %.17g,\n",m->elapsed_s);
    fprintf(o,"    \"steps_per_second\": %.17g,\n",m->steps/m->elapsed_s);
    fprintf(o,"    \"max_midpoint_residual_inf\": %.17g,\n",m->max_midpoint_resid);
    fprintf(o,"    \"max_state_error_vs_closed_form_midpoint_oracle\": %.17g,\n",m->max_state_err_midpoint_oracle);
    fprintf(o,"    \"max_state_error_vs_continuous_exact\": %.17g,\n",m->max_state_err_continuous);
    fprintf(o,"    \"max_flux_error_vs_continuous_exact_Wb\": %.17g,\n",m->max_flux_err_continuous);
    fprintf(o,"    \"max_emf_error_vs_continuous_exact_V\": %.17g,\n",m->max_emf_err_continuous);
    fprintf(o,"    \"max_area_relative_error_vs_continuous_exact\": %.17g,\n",m->max_area_rel_err_continuous);
    fprintf(o,"    \"max_B_vector_error_T\": %.17g,\n",m->max_B_vec_err_continuous);
    fprintf(o,"    \"max_area_vector_error_m2\": %.17g,\n",m->max_S_vec_err_continuous);
    fprintf(o,"    \"max_abs_flux_Wb\": %.17g,\n",m->max_flux_abs);
    fprintf(o,"    \"max_abs_emf_V\": %.17g,\n",m->max_emf_abs);
    fprintf(o,"    \"Bx_zero_crossings\": %lld,\n",m->bx_zero_crossings);
    fprintf(o,"    \"area_deviation_zero_crossings\": %lld,\n",m->area_dev_zero_crossings);
    fprintf(o,"    \"flux_zero_crossings\": %lld,\n",m->flux_zero_crossings);
    fprintf(o,"    \"emf_zero_crossings\": %lld\n",m->emf_zero_crossings);
    fprintf(o,"  },\n");
    fprintf(o,"  \"solver_contract\": {\"state_advance\": \"fixed-order NRII implicit-midpoint coefficient cascade\", \"global_matrix_assembly\": false, \"global_linear_solve\": false, \"residual_role\": \"audit only\", \"oracle\": \"closed-form 2x2 implicit-midpoint Cayley update for each oscillator block\"}\n");
    fprintf(o,"}\n");
    fclose(o);
}

static void write_dt_sweep(const Params *base,const char *path){
    const double dts[]={2e-7,1e-7,5e-8,2e-8,1e-8};
    FILE *o=fopen(path,"w"); if(!o)return;
    fprintf(o,"dt_s,steps,max_midpoint_resid,max_state_err_midpoint_oracle,max_state_err_continuous,max_flux_err_Wb,max_emf_err_V,elapsed_s\n");
    for(size_t i=0;i<sizeof(dts)/sizeof(dts[0]);i++){
        Params p=*base; p.dt=dts[i]; Metrics m;
        if(run_case(&p,&m,NULL,0)){
            fprintf(o,"%.17g,%lld,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",p.dt,m.steps,m.max_midpoint_resid,m.max_state_err_midpoint_oracle,m.max_state_err_continuous,m.max_flux_err_continuous,m.max_emf_err_continuous,m.elapsed_s);
        }
    }
    fclose(o);
}

static void write_area_control(const Params *base,const char *path){
    const double eps_values[]={0.0,base->epsA};
    FILE *o=fopen(path,"w"); if(!o)return;
    fprintf(o,"epsA,carrier_flux_amp_Wb,lower_sideband_flux_amp_Wb,upper_sideband_flux_amp_Wb,carrier_emf_amp_V,lower_sideband_emf_amp_V,upper_sideband_emf_amp_V,flux_zero_crossings,emf_zero_crossings\n");
    for(int k=0;k<2;k++){
        Params p=*base; p.epsA=eps_values[k]; Metrics m;
        if(run_case(&p,&m,NULL,0)){
            fprintf(o,"%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%lld,%lld\n",
                p.epsA,m.dft_flux_amp[0],m.dft_flux_amp[1],m.dft_flux_amp[2],
                m.dft_emf_amp[0],m.dft_emf_amp[1],m.dft_emf_amp[2],m.flux_zero_crossings,m.emf_zero_crossings);
        }
    }
    fclose(o);
}

int main(int argc,char **argv){
    const char *outdir=(argc>1)?argv[1]:".";
    Params p;
    p.fB_hz=40000.0;
    p.fA_hz=37600.0;
    p.fN_hz=39100.0;
    p.wB=2.0*M_PI*p.fB_hz; p.wA=2.0*M_PI*p.fA_hz; p.wN=2.0*M_PI*p.fN_hz;
    p.B0_T=1.0e-3;
    p.A0_m2=1.0e-6;
    p.epsA=0.15;
    p.dt=1.0e-8;
    p.T=1.0e-2;
    p.order=ORDER;

    char path[1024]; Metrics m;
    snprintf(path,sizeof(path),"%s/TRACE_DECIMATED.csv",outdir);
    if(!run_case(&p,&m,path,100)){fprintf(stderr,"run failed\n");return 2;}
    snprintf(path,sizeof(path),"%s/MODEL_PARAMETERS.csv",outdir); write_model(&p,path);
    snprintf(path,sizeof(path),"%s/SPECTRUM.csv",outdir); write_spectrum(&p,&m,path);
    snprintf(path,sizeof(path),"%s/RESULTS.json",outdir); write_results(&p,&m,path);
    snprintf(path,sizeof(path),"%s/DT_SWEEP.csv",outdir); write_dt_sweep(&p,path);
    snprintf(path,sizeof(path),"%s/AREA_MODULATION_CONTROL.csv",outdir); write_area_control(&p,path);

    printf("steps=%lld\n",m.steps);
    printf("elapsed_s=%.9g\n",m.elapsed_s);
    printf("steps_per_second=%.9g\n",m.steps/m.elapsed_s);
    printf("max_midpoint_residual_inf=%.17g\n",m.max_midpoint_resid);
    printf("max_state_error_vs_midpoint_oracle=%.17g\n",m.max_state_err_midpoint_oracle);
    printf("max_state_error_vs_continuous_exact=%.17g\n",m.max_state_err_continuous);
    printf("max_flux_error_vs_continuous_exact_Wb=%.17g\n",m.max_flux_err_continuous);
    printf("max_emf_error_vs_continuous_exact_V=%.17g\n",m.max_emf_err_continuous);
    printf("zero_crossings Bx=%lld area_dev=%lld flux=%lld emf=%lld\n",m.bx_zero_crossings,m.area_dev_zero_crossings,m.flux_zero_crossings,m.emf_zero_crossings);
    return 0;
}
