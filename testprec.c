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

int main(int argc, char** argv){
  long int i, N;
  REAL x;

  if(argc!=3){
    printf("usage: %s N, x \n", argv[0]);
    abort();
  }

  sscanf(argv[1], "%ld", &N);
	
  #ifdef FLOAT
    sscanf(argv[2], "%f", &x);
  #endif
  #ifdef DOUBLE
    sscanf(argv[2], "%lf", &x);
  #endif
  #ifdef LONG
    sscanf(argv[2], "%Lf", &x);
  #endif
  #ifdef QUAD
    x = strtoflt128(argv[2], NULL);
  #endif

  for(i=0; i<N; i++){
      x += x;
      x /=1.761348765013846510984561;
  }
      // print the result
      //printf("%ld   ", N);
      PRINT(x); printf("%s \n", buf);
    
  
}
