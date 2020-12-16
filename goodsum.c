/* We provide implementations for several precisions of 
the corrected sums of an array according to 
the Kahn algorithm documented in Knuth. 
*/


float goodsumf( float *x,long  int n){
  float s_ant, c, y,s;
  long int i;

  s_ant = c = 0.0;
  for (i= 0; i < n ; i++){
    y  = x[i] -c;
    s = s_ant + y;
    c = (s - s_ant) - y;
    s_ant = s;
  }
  return(s);
}


double goodsum( double *x, long int n){
  double s_ant, c, y,s;
  long  int i;

  s_ant = c = 0.0;
  for (i = 0; i < n ; i++){
    y  = x[i] -c;
    s = s_ant + y;
    c = (s - s_ant) - y;
    s_ant = s;
  }
  return(s);
}

long double goodsuml( long double *x, long int n){
  long double s_ant, c, y,s;
  long int i;

  s_ant = c = 0.0;
  for (i= 0; i < n ; i++){
    y  = x[i] -c;
    s = s_ant + y;
    c = (s - s_ant) - y;
    s_ant = s;
  }
  return(s);
}

#ifdef QUAD
__float128 goodsumq( __float128 *x, long int n){
  __float128 s_ant, c, y,s;
  long int i;
  s_ant = c = 0.0;
  for (i= 0; i < n ; i++){
    y  = x[i] -c;
    s = s_ant + y;
    c = (s - s_ant) - y;
    s_ant = s;
  }
  return(s);
}
#endif


float goodprodf( float *x, float *z,long  int n){
  float s_ant, c, y,s;
  long int i;

  s_ant = c = 0.0;
  for (i= 0; i < n ; i++){
    y  = x[i]*z[i] -c;
    s = s_ant + y;
    c = (s - s_ant) - y;
    s_ant = s;
  }
  return(s);
}


double goodprod( double *x, double *z,long int n){
  double s_ant, c, y,s;
  long int i;

  s_ant = c = 0.0;
  for (i = 0; i < n ; i++){
    y  = x[i]*z[i] -c;
    s = s_ant + y;
    c = (s - s_ant) - y;
    s_ant = s;
  }
  return(s);
}

long double goodprodl( long double *x, long double *z, long int n){
  long double s_ant, c, y,s;
  long  int i;

  s_ant = c = 0.0;
  for (i= 0; i < n ; i++){
    y  = x[i]*z[i] -c;
    s = s_ant + y;
    c = (s - s_ant) - y;
    s_ant = s;
  }
  return(s);
}

#ifdef QUAD
__float128 goodprodq( __float128 *x, __float128 *z, long int n){
  __float128 s_ant, c, y,s;
  long  int i;
  s_ant = c = 0.0;
  for (i= 0; i < n ; i++){
    y  = x[i]*z[i] -c;
    s = s_ant + y;
    c = (s - s_ant) - y;
    s_ant = s;
  }
  return(s);
}
#endif
