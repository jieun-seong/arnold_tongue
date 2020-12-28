#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "define.h"
#include "goodsum.c"

#ifdef QUAD
#include <quadmath.h>
#endif
#include <assert.h>

long int m, m_t, N, N_t;
REAL lambda, omega, initx, target, rN, rN_t;
REAL r1, r0, r2, r2pi, r3, r5, r10, rhalf;
REAL *orbit_t, *tot_weight_t;
REAL **weightarray_t;
char buf[60];

void PRINT(REAL x){
#ifdef FLOAT
  printf(" %f ", x);
#endif
#ifdef DOUBLE
  printf(" %20.15f ", (double)x);
#endif
#ifdef LONG
  printf(" Lf ", (long double)x);
#endif
#ifdef QUAD
  quadmath_snprintf(buf, sizeof buf, "%+-#*.33Qe", 33, x);
#endif
}

void PRINTTOCHECK(char *s, REAL x){
#ifdef FLOAT
  printf("\n"); printf("%s", s); printf(" = ");
  printf(" %f ", x);
  printf("\n\n");
#endif
#ifdef DOUBLE
  printf("\n"); printf("%s", s); printf(" = ");
  printf(" %20.15f ", (double)x);
  printf("\n\n");
#endif
#ifdef LONG
  printf("\n"); printf("%s", s); printf(" = ");
  printf(" Lf ", (long double)x);
  printf("\n\n");
#endif
#ifdef QUAD
  printf("\n"); printf("%s", s); printf(" = ");
  quadmath_snprintf(buf, sizeof buf, "%+-#*.33Qe", 33, x);
  printf("\n\n");
#endif
}

REAL step(REAL x){
  REAL dummy;
  dummy = lambda*SIN(r2pi*x)+omega;
  return dummy;
}

REAL weight(long int i){
  REAL frac, dummy;
  frac = ((REAL)i)/rN;
  if((r0<i)&&(i<N)){
    #ifdef BIRKHOFF
      return r1;
    #else // DSY
      dummy = EXP(-r1/(frac*(r1-frac)));
      return dummy;
    #endif
  }
  return r0;
}

REAL map(REAL x){
  REAL dummy;
  dummy = (REAL)(x+lambda*SIN(r2pi*x)+omega);
  dummy = dummy-(REAL)(floor(dummy));
  return dummy;
}

REAL rotnum(REAL omegahere){
  REAL x = initx;
  omega = omegahere;
  for(long int i=0; i<N_t; i++){
    orbit_t[i] = step(x);
    x = map(x);
  }
  return (GOODPROD(orbit_t, weightarray_t[m_t], N_t)/tot_weight_t[m_t]-target);
}

REAL dfdx(REAL x) {return (REAL) (r1+r2pi*lambda*COS(r2pi*x));}

REAL ddfdxx(REAL x) {return (REAL) (-r2pi*r2pi*lambda*SIN(r2pi*x));}

REAL dfdl(REAL x) {return (REAL) (SIN(r2pi*x));}

REAL ddfdll(REAL x) {return r0;}

REAL ddfdldx(REAL x) {return (REAL) (r2pi*COS(r2pi*x));}

REAL dfdw(REAL x) {return (REAL) (r1);}

REAL ddfdww(REAL x) {return r0;}

REAL ddfdwdx(REAL x) {return r0;}

REAL ddfdwdl(REAL x) {return r0;}

#define ITMAX 100
#ifdef FLOAT
#define EPS (1.0e-7)
#endif
#ifdef DOUBLE
#define EPS (1.0e-14)
#endif
#ifdef LONG
#define EPS (1.0e-17)
#endif
#ifdef QUAD
#define EPS (1.0e-32)
#endif

void nerror(char *s){
  printf("%s\n",s);
  abort();
}

REAL zbrent(REAL (*func) (REAL), REAL x1, REAL x2, REAL tol){
  int iter;
  REAL a=x1,b=x2,c,d,e,min1,min2,fa=(*func)(a),fb=(*func)(b),fc,p,q,r,s,tol1,xm;
  void nerror();
  if (fb*fa>r0) nerror("Root must be bracketed in ZBRENT\n");
  fc=fb;
  for(iter=1;iter<=ITMAX;iter++){
    if(fb*fc>r0){c=a;fc=fa;e=d=b-a;}
    if(ABS(fc)<ABS(fb)){a=b;b=c;c=a;fa=fb;fb=fc;fc=fa;}
    tol1=r2*EPS*ABS(b)+(tol*rhalf);xm=(c-b)*rhalf;
    if(ABS(xm)<=tol1||fb==r0){
      #ifdef VERBOSE
        printf("\neps = ");
        PRINT(ABS(xm));
        printf("\n");
      #endif
      return b;
    }
    if(ABS(e)>=tol1&&ABS(fa)>ABS(fb)){s=fb/fa;
      if(a==c){p=r2*xm*s;q=r1-s;}
      else{q=fa/fc;r=fb/fc;p=s*(r2*xm*q*(q-r)-(b-a)*(r-r1));q=(q-r1)*(r-r1)*(s-r1);}
      if(p>r0) q=-q;
      p=ABS(p);min1=r3*xm*q-ABS(tol1*q);min2=ABS(e*q);
      if((r2*p)<(min1<min2?min1:min2)){e=d;d=p/q;}
      else{d=xm;e=d;}
    }
    else{d=xm;e=d;}
    a=b;fa=fb;
    if(ABS(d)>tol1) b=b+d;
    else b=b+(xm>r0?ABS(tol1):-ABS(tol1));
    fb=(*func)(b);
  }
  nerror("Maximum number of iterations exceeded in ZBRENT");
  return r0;
}

#undef ITMAX

int main(int argc, char** argv){
  long int i;
  long int NMIN_der, NMAX_der, NSTEP_der, NMIN_t, NMAX_t;
  REAL lambdamin, lambdamax, lambdastep, lambda_critical;
  REAL **weightarray_der;
  REAL *tot_weight_der, *orbit_der, *d1_lambda_array, *d2_lambda_array, *d1_omega_array, *d2_omega_array, *d2_omega_lambda_array;
  REAL x, y_lambda, y_omega, d1_lambda, d2_lambda, d1_omega, d2_omega, dw, ddw, d_omega_lambda;
  int lambda_index;
  REAL rotnum_curr;

  if(argc!=9){
    printf("usage: %s N_tongue, NMIN_deriv, NMAX_deriv, NSTEP_deriv, lambdamin, lambdamax, lambdastep, initx \n", argv[0]);
    printf("critical lambda = 0.15915494309189533576888376337251436203445964574... \n");
    abort();
  }

  sscanf(argv[1], "%ld", &N_t);
  sscanf(argv[2], "%ld", &NMIN_der);
  sscanf(argv[3], "%ld", &NMAX_der);
  sscanf(argv[4], "%ld", &NSTEP_der);
	
  #ifdef FLOAT
    sscanf(argv[5], "%f", &lambdamin);
    sscanf(argv[6], "%f", &lambdamax);
    sscanf(argv[7], "%f", &lambdastep);
    sscanf(argv[8], "%f", &initx);
  #endif
  #ifdef DOUBLE
    sscanf(argv[5], "%lf", &lambdamin);
    sscanf(argv[6], "%lf", &lambdamax);
    sscanf(argv[7], "%lf", &lambdastep);
    sscanf(argv[8], "%lf", &initx);
  #endif
  #ifdef LONG
    sscanf(argv[5], "%Lf", &lambdamin);
    sscanf(argv[6], "%Lf", &lambdamax);
    sscanf(argv[7], "%Lf", &lambdastep);
    sscanf(argv[8], "%Lf", &initx);
  #endif
  #ifdef QUAD
    lambdamin = strtoflt128(argv[5], NULL);
    lambdamax = strtoflt128(argv[6], NULL);
    lambdastep = strtoflt128(argv[7], NULL);
    initx = strtoflt128(argv[8], NULL);
  #endif
	
  printf("# N_t %ld\n", N_t);
  printf("# NMIN_deriv %ld\n", NMIN_der);
  printf("# NMAX_deriv %ld\n", NMAX_der);
  printf("# NSTEP_deriv %ld\n", NSTEP_der);
  printf("# lambdamin "); PRINT(lambdamin); printf("%s \n", buf);
  printf("# lambdamax "); PRINT(lambdamax); printf("%s \n", buf);
  printf("# lambdastep "); PRINT(lambdastep); printf("%s \n", buf);
  printf("# initx "); PRINT(initx); printf("%s \n\n", buf);

  NMIN_t = N_t;
  NMAX_t = N_t;
    
  r0 = (REAL)0; 
  r1 = (REAL)1;
  r2 = (REAL)2;
  r3 = (REAL)3;
  r5 = (REAL)5;
  r10 = (REAL)10;
  r2pi = r2*PI;
  rhalf = r1/r2;

  target = (REAL)(SQRT(r5)-r1)/r2;
  lambda_critical = (REAL)(r1/r2pi);
  
  printf("#1:N  #2:lambda_index  #3:omega  #4:lambda  #5:|lambda_critical-lambda| #6:|rotnum_curr-target|  #7:d1_lambda  #8:d1_omega  #9:d2_lambda  #10:d2_omega  #11:d_omega_lambda  #12:dw  #13:ddw \n\n");
    
  // construct the weight array for the tongue
  weightarray_t = calloc(sizeof(REAL*), NMAX_t-NMIN_t+1);
  tot_weight_t = calloc(sizeof(REAL), NMAX_t-NMIN_t+1);
  orbit_t = calloc(sizeof(REAL), NMAX_t);
  for(m_t=0; m_t<NMAX_t-NMIN_t+1; m_t++){
    N_t = m_t+NMIN_t;
    rN_t = (REAL) N_t;
    weightarray_t[m_t] = calloc(sizeof(REAL), N_t);
    N = N_t;
    rN = rN_t;
    for(i=0; i<N_t; i++){
      weightarray_t[m_t][i] = weight(i);
    }
    tot_weight_t[m_t] = GOODSUM(weightarray_t[m_t], N_t);
  }
  m_t = 0; // since for now I only look at the case NMAX_t = NMIN_t
  // need to loop through m_t if NMAX_t != NMIN_t
  // `` for derivatives
  weightarray_der = calloc(sizeof(REAL *), NMAX_der-NMIN_der+1);
  tot_weight_der = calloc(sizeof(REAL), NMAX_der-NMIN_der+1);
  //orbit_der = calloc(sizeof(REAL), NMAX_der);
  for(m=0; m<NMAX_der-NMIN_der+1; m++){
    N = m+NMIN_der;
    rN = (REAL) N;
    weightarray_der[m] = calloc(sizeof(REAL), N);
    orbit_der = calloc(sizeof(REAL), N);
    for(i=0; i<N; i++){
      weightarray_der[m][i] = weight(i);
    }
    tot_weight_der[m] = GOODSUM(weightarray_der[m], N);
    d1_lambda_array = calloc(sizeof(REAL), N);
    d1_omega_array = calloc(sizeof(REAL), N);
    d2_lambda_array = calloc(sizeof(REAL), N);
    d2_omega_array = calloc(sizeof(REAL), N);
    d2_omega_lambda_array = calloc(sizeof(REAL), N);
    //calculate derivatives
    lambda_index = 0; // for convenience in printing
    for(lambda=lambdamin; lambda<lambdamax; lambda+=lambdastep){
      omega = zbrent(rotnum, r0+r10*EPS, r1-r10*EPS, EPS*r10);
      orbit_der[0] = initx;
      d1_lambda_array[0] = r0;
      d2_lambda_array[0] = r0;
      d1_omega_array[0] = r0;
      d2_omega_array[0] = r0;
      d2_omega_lambda_array[0] = r0;
      for(i=0; i<N-1; i++){
        x = orbit_der[i];
        y_lambda = d1_lambda_array[i];
        y_omega = d1_omega_array[i];
        orbit_der[i+1] = map(x);
        d1_lambda_array[i+1] = dfdl(x)+dfdx(x)*y_lambda;
        d1_omega_array[i+1] = dfdw(x)+dfdx(x)*y_omega;
        d2_lambda_array[i+1] = ddfdll(x)+ddfdldx(x)*y_lambda*r2+ddfdxx(x)*y_lambda*y_lambda+dfdx(x)*d2_lambda_array[i];
        d2_omega_array[i+1] = ddfdww(x)+ddfdwdx(x)*y_omega*r2+ddfdxx(x)*y_omega*y_omega+dfdx(x)*d2_omega_array[i];
        d2_omega_lambda_array[i+1] = ddfdwdl(x)+ddfdldx(x)*y_omega+(ddfdwdx(x)+ddfdxx(x)*y_omega)*y_lambda+dfdx(x)*d2_omega_lambda_array[i];
      }
      rotnum_curr = GOODPROD(orbit_der, weightarray_der[m], N)/tot_weight_der[m];
      d1_lambda = GOODPROD(d1_lambda_array, weightarray_der[m], N)/tot_weight_der[m];
      d1_omega = GOODPROD(d1_omega_array, weightarray_der[m], N)/tot_weight_der[m];
      d2_lambda = GOODPROD(d2_lambda_array, weightarray_der[m], N)/tot_weight_der[m];
      d2_omega = GOODPROD(d2_omega_array, weightarray_der[m], N)/tot_weight_der[m];
      d_omega_lambda = GOODPROD(d2_omega_lambda_array, weightarray_der[m], N)/tot_weight_der[m];
      dw = -r1*d1_lambda/d1_omega;
      ddw = (d2_lambda-r2*d_omega_lambda*dw-d2_omega*dw)/d1_omega;
      // print the result
      printf("%ld   ", N);
      printf("%d   ", lambda_index); lambda_index++;
      PRINT(omega); printf("%s ", buf);
      PRINT(lambda); printf("%s ", buf);
      PRINT(lambda_critical-lambda); printf("%s ", buf);
      PRINT(ABS(rotnum_curr-target)); printf("%s ", buf);
      PRINT(d1_lambda); printf("%s ", buf);
      PRINT(d1_omega); printf("%s ", buf);
      PRINT(d2_lambda); printf("%s ", buf);
      PRINT(d2_omega); printf("%s ", buf);
      PRINT(d_omega_lambda); printf("%s ", buf);
      PRINT(dw); printf("%s ", buf);
      PRINT(ddw); printf("%s \n", buf);
    }
  }
}
