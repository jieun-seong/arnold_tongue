#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "define.h"
#include "goodsum.c"

#ifdef QUAD
#include <quadmath.h>
#endif
#include <assert.h>

long int m, N;
REAL lambda, omega, initx, target, rN;
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
  if((0<i)&&(i<N)){
    dummy = EXP(-r1/(frac*(r1-frac)));
    return dummy;
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
  for(long int i=0; i<N; i++){
    orbit_t[i] = step(x);
    x = map(x);
  }
  return (GOODPROD(orbit_t, weightarray_t[m], N)/tot_weight_t[m]-target);
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
  long int N_t, N_der;
  REAL lambdamin, lambdamax, lambdastep, lambda_critical;
  REAL **weightarray_der;
  REAL *tot_weight_der, *orbit_der, *d1_lambda_array, *d2_lambda_array, *d1_omega_array, *d2_omega_array, *d2_omega_lambda_array;
  REAL x, y_lambda, y_omega, d1_lambda, d2_lambda, d1_omega, d2_omega, dw, ddw, d_omega_lambda;

  if(argc!=7){
    printf("usage: %s N_tongue, N_deriv, lambdamin, lambdamax, lambdastep, initx \n", argv[0]);
    printf("critical lambda = 0.15915494309189533576888376337251436203445964574... \n");
    abort();
  }

  sscanf(argv[1], "%ld", &N_t);
  sscanf(argv[2], "%ld", &N_der);
	
  #ifdef FLOAT
    sscanf(argv[3], "%f", &lambdamin);
    sscanf(argv[4], "%f", &lambdamax);
    sscanf(argv[5], "%f", &lambdastep);
    sscanf(argv[6], "%f", &initx);
  #endif
  #ifdef DOUBLE
    sscanf(argv[3], "%lf", &lambdamin);
    sscanf(argv[4], "%lf", &lambdamax);
    sscanf(argv[5], "%lf", &lambdastep);
    sscanf(argv[6], "%lf", &initx);
  #endif
  #ifdef LONG
    sscanf(argv[3], "%Lf", &lambdamin);
    sscanf(argv[4], "%Lf", &lambdamax);
    sscanf(argv[5], "%Lf", &lambdastep);
    sscanf(argv[6], "%Lf", &initx);
  #endif
  #ifdef QUAD
    lambdamin = strtoflt128(argv[3], NULL);
    lambdamax = strtoflt128(argv[4], NULL);
    lambdastep = strtoflt128(argv[5], NULL);
    initx = strtoflt128(argv[6], NULL);
  #endif
	
  printf("# N_t %ld\n", N_t);
  printf("# N_der %ld\n", N_der);
  printf("# lambdamin "); PRINT(lambdamin); printf("%s \n", buf);
  printf("# lambdamax "); PRINT(lambdamax); printf("%s \n", buf);
  printf("# lambdastep "); PRINT(lambdastep); printf("%s \n", buf);
  printf("# initx "); PRINT(initx); printf("%s \n", buf);

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

  // construct the weight array for the tongue

  weightarray_t = calloc(sizeof(REAL*), N_t);
  tot_weight_t = calloc(sizeof(REAL), N_t);
  orbit_t = calloc(sizeof(REAL), N_t);

  for(m=0; m<N_t; m++){
    N = m+1;
    rN = (REAL) N;
    weightarray_t[m] = calloc(sizeof(REAL), N);
    /*#ifdef VERBOSE
      printf("\n=====================================\n");
      printf("=====================================\n");
      PRINTTOCHECK("m", m);
    #endif*/
    for(i=0; i<N; i++){
      weightarray_t[m][i] = weight(i);
    }
    tot_weight_t[m] = GOODSUM(weightarray_t[m], N);
    /*#ifdef VERBOSE
      for(i=0; i<N; i++){
        printf("\nw_t[%ld][%ld] = ", m, i); PRINT(weight(i));
      }
      PRINTTOCHECK("\ntot_w_t[m]", tot_weight_t[m]);
    #endif*/
  }
    
  // construct the weight array for the derivatives

  weightarray_der = calloc(sizeof(REAL *), N_der);
  tot_weight_der = calloc(sizeof(REAL), N_der);
  orbit_der = calloc(sizeof(REAL), N_der);

  //for(m=0; m<N_der; m++){
  m = N_der-1;
  N = m+1;
  rN = (REAL) N;
  weightarray_der[m] = calloc(sizeof(REAL), N);
  #ifdef VERBOSE
    printf("\n=====================================\n");
    printf("=====================================\n");
    PRINTTOCHECK("m_d", m);
  #endif
  for(i=0; i<N; i++){
    weightarray_der[m][i] = weight(i);
  }
  tot_weight_der[m] = GOODSUM(weightarray_der[m], N);
  #ifdef VERBOSE
    for(i=0; i<N; i++){
      printf("\nw_d[%ld][%ld] = ", m, i);
      PRINT(weightarray_der[m][i]);
    }
    PRINTTOCHECK("\ntot_w_d[m]", tot_weight_t[m]);
  #endif
  //}

  d1_lambda_array = calloc(sizeof(REAL), N_der);
  d1_omega_array = calloc(sizeof(REAL), N_der);
  #ifndef TEST1
    d2_lambda_array = calloc(sizeof(REAL), N_der);
    d2_omega_array = calloc(sizeof(REAL), N_der);
    d2_omega_lambda_array = calloc(sizeof(REAL), N_der);
  #endif

  m = N_t-1;
  //for(m=N_t-1; m<N_t; m++){
  N = m+1;
  for(lambda=(REAL)lambdamin; lambda<(REAL)lambdamax; lambda+=(REAL)lambdastep){
    omega = zbrent(rotnum, r0+r10*EPS, r1-r10*EPS, EPS*r10);
    #ifdef VERBOSE
      PRINTTOCHECK("lambda = ", lambda);
      PRINTTOCHECK("omega = ", omega);
    #endif
    orbit_der[0] = initx;
    d1_lambda_array[0] = r0;
    d1_omega_array[0] = r0;
    #ifndef TEST1
      d2_lambda_array[0] = r0;
      d2_omega_array[0] = r0;
      d2_omega_lambda_array[0] = r0;
    #endif
    for(i=0; i<N_der-1; i++){
      x = orbit_der[i];
      y_lambda = d1_lambda_array[i];
      y_omega = d1_omega_array[i];
      orbit_der[i+1] = map(x);
      d1_lambda_array[i+1] = dfdl(x)+dfdx(x)*y_lambda;
      d1_omega_array[i+1] = dfdw(x)+dfdx(x)*y_omega;
      #ifndef TEST1
        d2_lambda_array[i+1] = ddfdll(x)+ddfdldx(x)*y_lambda*r2+ddfdxx(x)*y_lambda*y_lambda+dfdx(x)*d2_lambda_array[i];
        d2_omega_array[i+1] = ddfdww(x)+ddfdwdx(x)*y_omega*r2+ddfdxx(x)*y_omega*y_omega+dfdx(x)*d2_omega_array[i];
        d2_omega_lambda_array[i+1] = ddfdwdl(x)+ddfdldx(x)*y_omega+(ddfdwdx(x)+ddfdxx(x)*y_omega)*y_lambda+dfdx(x)*d2_omega_lambda_array[i];
      #endif
    }
    #ifdef VERBOSE
      for(i=0; i<N_der; i++){
        printf("\norbit_d[%ld] = ", i);
        PRINT(orbit_der[i]);
        printf("\n");
        printf("d1_lambda_array[%ld] = ", i);
        PRINT(d1_lambda_array[i]);
        printf("\n");
        printf("d1_omega_array[%ld] = ", i);
        PRINT(d1_omega_array[i]);
        printf("\n");
      }
    #endif
    //for(i=0; i<N_der; i++)
    i=N_der-1;
    #ifdef VERBOSE
      for(int p=0; p<N_der; p++){
        printf("w_d[%ld][%d] = ", i, p);
        PRINT(weightarray_der[i][p]);
        printf("\n");
      }
    #endif
    d1_lambda = GOODPROD(d1_lambda_array, weightarray_der[i], N_der)/tot_weight_der[i];
    d1_omega = GOODPROD(d1_omega_array, weightarray_der[i], N_der)/tot_weight_der[i];
      
    #ifdef VERBOSE
      PRINTTOCHECK("d1_l = ", d1_lambda);
      PRINTTOCHECK("d1_w = ", d1_omega);
    #endif
    
    #ifndef TEST1
      d2_lambda = GOODPROD(d2_lambda_array, weightarray_der[i], i+1)/tot_weight_der[i];
      d2_omega = GOODPROD(d2_omega_array, weightarray_der[i], i+1)/tot_weight_der[i];
      d_omega_lambda = GOODPROD(d2_omega_lambda_array, weightarray_der[i], i+1)/tot_weight_der[i];
      dw = -r1*d1_lambda/d1_omega;
      ddw = (d2_lambda-r2*d_omega_lambda*dw-d2_omega*dw)/d1_omega;
    #endif
			
    #ifndef VERBOSE
    #ifndef TEST1
      PRINT(omega); printf("%s ", buf);
      PRINT(lambda); printf("%s ", buf);
      PRINT(lambda_critical-lambda); printf("%s ", buf);
      PRINT(d1_lambda); printf("%s ", buf);
      PRINT(d1_omega); printf("%s ", buf);
      PRINT(d2_lambda); printf("%s ", buf);
      PRINT(d2_omega); printf("%s ", buf);
      PRINT(d_omega_lambda); printf("%s ", buf);
      PRINT(dw); printf("%s ", buf);
      PRINT(ddw); printf("%s \n", buf);
    #endif
    #endif
  }
  //}
  #ifdef TEST1
    printf("\n# column1 = step\n");
    printf("# column2 = error d_lambda\n");
    printf("# column3 = error d_omega\n\n");
    REAL step, currx, pertx, error_lambda, error_omega, currstep;
    int p=1; int ncomp = 7;
    //for(p=1; p<7; p++);
    for(step=(REAL)0.1; step>=0.0000001; step*=(REAL)0.1){
      currx = initx; 
      pertx = initx+pow(0.9,p)*(step); //0.9^p
      for(int t=1; t<=ncomp; t++){
        currx = map(currx);
        pertx = map(pertx);
      }
      currstep = pow(0.9,p)*step; 
      PRINT(currstep); printf("%s ", buf);
      error_lambda = pertx-currx-d1_lambda*pow(0.9,p)*step;
      error_omega = pertx-currx-d1_omega*pow(0.9,p)*step;
      PRINT(error_lambda); printf("%s ", buf);
      PRINT(error_omega); printf("%s \n", buf);
      p++;
    }
  #endif
}