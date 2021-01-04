#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "define.h"
#include "goodsum.c"

#ifdef QUAD
#include <quadmath.h>
#endif
#include <assert.h>

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

int main(int argc, char** argv){
    long int i;
    REAL x;
    float d_value, p_value;
    for (int i = 0; i < 1000000; ++i){
        d_value = (REAL)10.03568375298674502645 / (REAL)(i);
        p_value = (REAL)0.010384015630856204856245 * (REAL)(i) + (REAL)100.0019348759081346502498;
    }
    PRINT(d_value); printf("%s ", buf);
    PRINT(p_value); printf("%s \n", buf);
}
