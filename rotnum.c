#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#ifdef QUAD
#include <quadmath.h>
#endif
#include "define.h"
#include "goodsum.c"

long int N, n_t;
REAL rN, r0, r1, r2, r3, r5, r10, r2pi, rhalf;
REAL lambda, omega, target, initx;
REAL *orbit, *orbit_t, *tot_weight, *tot_weight_t;
REAL **weightarray, **weightarray_t;
char buf[60];

void PRINT(REAL x){
  #ifdef FLOAT
    printf(" %f ", x);
  #endif
  #ifdef DOUBLE
    printf(" %20.15f ", (double) x);
  #endif
  #ifdef LONG
    printf(" Lf ", (long double) x);
  #endif
  #ifdef QUAD
    quadmath_snprintf(buf, sizeof buf, "%+-#*.33Qe", 33, x);
  #endif
}

REAL weight(long int i){
  REAL frac, dummy;
  frac = (REAL)i/rN;
  if((0<i)&&(i<N)){
    dummy = EXP(-r1/(frac*(r1-frac)));
    return dummy;
  }
  return r0;
}

REAL step(REAL x){
  REAL dummy;
  dummy = lambda*SIN(r2pi*x)+omega;
  return dummy;
}

REAL rotnum(REAL omega_loc){
  REAL x, stepx, dummy;
  x = initx;
  omega = omega_loc;
  for(long int j=0; j<N; j++){
    stepx = step(x);
    orbit_t[j] = stepx;
    x = x+stepx;
    x = x-(REAL)FLOOR(x);
  }
  dummy = (GOODPROD(orbit_t, weightarray_t[n_t], N)/tot_weight_t[n_t])-target;
  return dummy;
}

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
  if (fb*fa>r0) nerror("Root must be bracketed in ZBRENT");
  fc=fb;
  for(iter=1;iter<=ITMAX;iter++){
    if(fb*fc>r0){c=a;fc=fa;e=d=b-a;}
    if(ABS(fc)<ABS(fb)){a=b;b=c;c=a;fa=fb;fb=fc;fc=fa;}
    tol1=r2*EPS*ABS(b)+(tol*rhalf);xm=(c-b)*rhalf;
    if(ABS(xm)<=tol1||fb==r0) return b;
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
    
  long int NMIN, NMAX, NMIN_t, NMAX_t;
  long int i, j, k;
  REAL x, stepx;
  REAL rotnum_cur, dist;
    
  if(argc != 7){
    printf("\nusage: %s NMIN, NMAX, NMIN_t, NMAX_t, initx, lambda\n", argv[0]);
    printf("\n(the golden mean = 0.61803398874989484820458683436563811772...) \n");
    printf("(the critical lambda = 0.15915494309189533576888376337251436203445964574...) \n\n");
    abort();
  }

  sscanf(argv[1], "%ld", &NMIN);
  sscanf(argv[2], "%ld", &NMAX);
  sscanf(argv[3], "%ld", &NMIN_t);
  sscanf(argv[4], "%ld", &NMAX_t);

  #ifdef FLOAT
  sscanf(argv[5], "%f", &initx);
  sscanf(argv[6], "%f", &lambda);
  #endif
  #ifdef DOUBLE
  sscanf(argv[5], "%lf", &initx);
  sscanf(argv[6], "%lf", &lambda);
  #endif
  #ifdef LONG
  sscanf(argv[5], "%Lf", &initx);
  sscanf(argv[6], "%Lf", &lambda);
  #endif
  #ifdef QUAD
  initx = strtoflt128(argv[5], NULL);
  lambda = strtoflt128(argv[6], NULL);
  #endif

  NMIN_t = NMAX_t; // no need to loop through for now
  printf("# NMIN %ld\n", NMIN);
  printf("# NMAX %ld\n", NMAX);
  printf("# NMIN_t %ld\n", NMIN_t);
  printf("# NMAX_t %ld\n", NMAX_t);
  printf("# initx "); PRINT(initx); printf("%s \n", buf);
  printf("# lambda "); PRINT(lambda); printf("%s \n", buf);

  r0 = (REAL)0;
  r1 = (REAL)1;
  r2 = (REAL)2;
  r3 = (REAL)3;
  r5 = (REAL)5;
  r10 = (REAL)10;
  r2pi = (REAL) (r2*PI);
  rhalf = (REAL) (r1/r2);
  target = ((REAL)(SQRT(r5)-r1)/r2);
  printf("# target "); PRINT(target); printf("%s \n", buf);

  // tongue

  orbit_t = (REAL *) calloc(sizeof(REAL), NMAX_t);
  tot_weight_t = (REAL *) calloc(sizeof(REAL), NMAX_t-NMIN_t+1);
  weightarray_t = (REAL **) calloc(sizeof(REAL *), NMAX_t-NMIN_t+1);

  for(n_t=0; n_t<NMAX_t-NMIN_t+1; n_t++){
    N = n_t+NMIN_t;
    rN = (REAL)N;
    weightarray_t[n_t] = (REAL *) calloc(sizeof(REAL), N);
    for(k=0; k<N; k++){
      weightarray_t[n_t][k] = weight(k);
    }
    tot_weight_t[n_t] = (GOODSUM(weightarray_t[n_t], N));
    omega = zbrent(rotnum, r0+r10*EPS, r1-r10*EPS, EPS*r10);
    free(weightarray_t[n_t]);
  }

  printf("# omega "); PRINT(omega); printf("%s \n", buf);

  free(tot_weight_t);
  free(weightarray_t);

  // rotnum

  orbit = (REAL *) calloc(sizeof(REAL), NMAX);
  tot_weight = (REAL *) calloc(sizeof(REAL), NMAX-NMIN+1);
  weightarray = (REAL **) calloc(sizeof(REAL *), NMAX-NMIN+1);

  for(j=0; j<NMAX; j++){
    stepx = step(x);
    orbit[j] = stepx;
    x = x+stepx;
    x = x-(REAL)FLOOR(x);
    #ifdef VERBOSE
      printf("\n orbit[%ld] = ", j);
      PRINT(orbit[j]);
    #endif
  }

  for(j=0; j<NMAX-NMIN+1; j++){
    N = j+NMIN;
    rN = (REAL)N;
    weightarray[j] = (REAL *) calloc(sizeof(REAL), N);
    for(k=0; k<j+NMIN; k++){
      weightarray[j][k] = (REAL) weight(k);
      #ifdef VERBOSE
        printf("\n w_array[%ld][%ld] = ",j,k);
        PRINT(weightarray[j][k]);
      #endif
    }
    tot_weight[j] = GOODSUM(weightarray[j], N);
    #ifdef VERBOSE
      printf("\n tot_w[%ld] = ", j);
      PRINT(tot_weight[j]);
    #endif
    rotnum_cur = GOODPROD(orbit, weightarray[j], N)/tot_weight[j];
    #ifdef VERBOSE
      printf("\n rotnum_cur = ");
      PRINT(rotnum_cur);
    #endif
    free(weightarray[j]);
    dist = ABS(rotnum_cur-target);
    #ifndef VERBOSE
      printf("%ld ", N);
      //PRINT(rotnum_cur);
      //printf("%s ", buf);
      PRINT(dist);
      printf("%s \n", buf);
    #endif
  }

  free(weightarray);
  free(orbit);
  free(tot_weight);	
}
