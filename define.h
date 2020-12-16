#ifdef QUAD
#define REAL __float128
#define PI ((REAL) M_PIq)
#define ABS(x) fabsq(x)
#define SIN(x) sinq(x)
#define COS(x) cosq(x)
#define EXP(x) expq(x)
#define ATAN(x) atanq(x)
#define SQRT(x) sqrtq(x)
#define FLOOR(x) floorq(x)
#define GOODSUM(x,n) goodsumq(x,n)
#define GOODPROD(x,z,n) goodprodq(x,z,n)
#endif

#ifdef LONG
#define REAL long double
#define PI ((long double) M_PI)
#define ABS(x) fabsl(x)
#define SIN(x) sinl(x)
#define COS(x) cosl(x)
#define EXP(x) expl(x)
#define ATAN(x) atanl(x)
#define SQRT(x) sqrtl(x)
#define FLOOR(x) floorl(x)
#define GOODSUM(x,n) goodsuml(x,n)
#define GOODPROD(x,z,n) goodprodl(x,z,n)
#endif


#ifdef DOUBLE
#define REAL double
#ifndef M_PI
#define M_PI 3.141592653589793238462643
#endif
#define PI ((double) M_PI)
#define ABS(x) fabs(x)
#define SIN(x) sin(x)
#define COS(x) cos(x)
#define EXP(x) exp(x)
#define LOG(x) log(x)
#define ATAN(x) atan(x)
#define SQRT(x) sqrt(x)
#define FLOOR(x) floor(x)
#define GOODSUM(x,n) goodsum(x,n)
#define GOODPROD(x,z,n) goodprod(x,z,n)
#endif

#ifdef FLOAT
#define REAL float
#define (PI (float) M_PI)
#define ABS(x) fabsf(x)
#define SIN(x) sinf(x)
#define COS(x) cosf(x)
#define EXP(x) expf(x)
#define ATAN(x) atanf(x)
#define FLOOR(x) floorf(x)
#define SQRT(x) sqrtf(x)
#define GOODSUM(x,n) goodsumf(x,n)
#define GOODPROD(x,z,n) goodprodf(x,z,n)
#endif
