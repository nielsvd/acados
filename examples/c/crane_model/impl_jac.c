/* This file was automatically generated by CasADi.
   The CasADi copyright holders make no ownership claim of its contents. */
#ifdef __cplusplus
extern "C" {
#endif

/* How to prefix internal symbols */
#ifdef CODEGEN_PREFIX
  #define NAMESPACE_CONCAT(NS, ID) _NAMESPACE_CONCAT(NS, ID)
  #define _NAMESPACE_CONCAT(NS, ID) NS ## ID
  #define CASADI_PREFIX(ID) NAMESPACE_CONCAT(CODEGEN_PREFIX, ID)
#else
  #define CASADI_PREFIX(ID) impl_jac_ ## ID
#endif

#include <math.h>

#ifndef casadi_real
#define casadi_real double
#endif

#define to_double(x) (double) x
#define to_int(x) (int) x
#define CASADI_CAST(x,y) (x) y

/* Pre-c99 compatibility */
#if __STDC_VERSION__ < 199901L
  #define fmin CASADI_PREFIX(fmin)
  casadi_real fmin(casadi_real x, casadi_real y) { return x<y ? x : y;}
  #define fmax CASADI_PREFIX(fmax)
  casadi_real fmax(casadi_real x, casadi_real y) { return x>y ? x : y;}
#endif

/* CasADi extensions */
#define sq CASADI_PREFIX(sq)
casadi_real sq(casadi_real x) { return x*x;}
#define sign CASADI_PREFIX(sign)
casadi_real CASADI_PREFIX(sign)(casadi_real x) { return x<0 ? -1 : x>0 ? 1 : x;}
#define twice CASADI_PREFIX(twice)
casadi_real twice(casadi_real x) { return x+x;}
#define if_else CASADI_PREFIX(if_else)
casadi_real if_else(casadi_real c, casadi_real x, casadi_real y) { return c!=0 ? x : y;}

/* Add prefix to internal symbols */
#define casadi_f0 CASADI_PREFIX(f0)
#define casadi_s0 CASADI_PREFIX(s0)
#define casadi_s1 CASADI_PREFIX(s1)
#define casadi_s2 CASADI_PREFIX(s2)
#define casadi_s3 CASADI_PREFIX(s3)
#define casadi_s4 CASADI_PREFIX(s4)

/* Printing routine */
#define PRINTF printf

/* Symbol visibility in DLLs */
#ifndef CASADI_SYMBOL_EXPORT
  #if defined(_WIN32) || defined(__WIN32__) || defined(__CYGWIN__)
    #if defined(STATIC_LINKED)
      #define CASADI_SYMBOL_EXPORT
    #else
      #define CASADI_SYMBOL_EXPORT __declspec(dllexport)
    #endif
  #elif defined(__GNUC__) && defined(GCC_HASCLASSVISIBILITY)
    #define CASADI_SYMBOL_EXPORT __attribute__ ((visibility ("default")))
  #else
    #define CASADI_SYMBOL_EXPORT
  #endif
#endif

static const int casadi_s0[8] = {4, 1, 0, 4, 0, 1, 2, 3};
static const int casadi_s1[5] = {1, 1, 0, 1, 0};
static const int casadi_s2[13] = {4, 4, 0, 0, 2, 3, 6, 2, 3, 0, 1, 2, 3};
static const int casadi_s3[11] = {4, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3};
static const int casadi_s4[6] = {4, 1, 0, 2, 2, 3};

/* impl_jacFun:(i0[4],i1[4],i2,i3)->(o0[4x4,6nz],o1[4x4,4nz],o2[4x1,2nz]) */
static int casadi_f0(const casadi_real** arg, casadi_real** res, int* iw, casadi_real* w, void* mem) {
  casadi_real a0, a1, a10, a11, a12, a13, a14, a15, a2, a3, a4, a5, a6, a7, a8, a9;
  a0=arg[0] ? arg[0][3] : 0;
  a1=sq(a0);
  a2=-8.0000000000000016e-02;
  a3=arg[0] ? arg[0][1] : 0;
  a4=cos(a3);
  a4=(a2*a4);
  a4=(a1*a4);
  a5=9.8100000000000009e-01;
  a6=cos(a3);
  a6=(a5*a6);
  a7=cos(a3);
  a7=(a6*a7);
  a8=sin(a3);
  a9=sin(a3);
  a9=(a5*a9);
  a9=(a8*a9);
  a7=(a7-a9);
  a4=(a4+a7);
  a7=arg[3] ? arg[3][0] : 0;
  a9=1.0000000000000001e-01;
  a10=(a7+a9);
  a11=cos(a3);
  a12=sq(a11);
  a12=(a9*a12);
  a10=(a10-a12);
  a4=(a4/a10);
  a12=sin(a3);
  a12=(a2*a12);
  a1=(a12*a1);
  a13=arg[2] ? arg[2][0] : 0;
  a1=(a1+a13);
  a6=(a6*a8);
  a1=(a1+a6);
  a1=(a1/a10);
  a1=(a1/a10);
  a11=(a11+a11);
  a6=sin(a3);
  a11=(a11*a6);
  a11=(a9*a11);
  a1=(a1*a11);
  a4=(a4-a1);
  a4=(-a4);
  if (res[0]!=0) res[0][0]=a4;
  a4=sq(a0);
  a1=cos(a3);
  a1=(a2*a1);
  a11=cos(a3);
  a11=(a1*a11);
  a6=sin(a3);
  a8=sin(a3);
  a2=(a2*a8);
  a2=(a6*a2);
  a11=(a11-a2);
  a11=(a4*a11);
  a2=sin(a3);
  a2=(a13*a2);
  a11=(a11-a2);
  a2=cos(a3);
  a2=(a5*a2);
  a11=(a11+a2);
  a2=9.8100000000000005e+00;
  a2=(a2*a7);
  a8=cos(a3);
  a8=(a2*a8);
  a11=(a11+a8);
  a8=8.0000000000000004e-01;
  a7=(a7+a9);
  a14=cos(a3);
  a15=sq(a14);
  a15=(a9*a15);
  a7=(a7-a15);
  a7=(a8*a7);
  a11=(a11/a7);
  a1=(a1*a6);
  a4=(a1*a4);
  a6=cos(a3);
  a13=(a13*a6);
  a4=(a4+a13);
  a13=sin(a3);
  a5=(a5*a13);
  a4=(a4+a5);
  a5=sin(a3);
  a2=(a2*a5);
  a4=(a4+a2);
  a4=(a4/a7);
  a4=(a4/a7);
  a14=(a14+a14);
  a3=sin(a3);
  a14=(a14*a3);
  a9=(a9*a14);
  a8=(a8*a9);
  a4=(a4*a8);
  a11=(a11-a4);
  a11=(-a11);
  if (res[0]!=0) res[0][1]=a11;
  a11=-1.;
  if (res[0]!=0) res[0][2]=a11;
  if (res[0]!=0) res[0][3]=a11;
  a11=(a0+a0);
  a12=(a12*a11);
  a12=(a12/a10);
  a12=(-a12);
  if (res[0]!=0) res[0][4]=a12;
  a0=(a0+a0);
  a1=(a1*a0);
  a1=(a1/a7);
  a1=(-a1);
  if (res[0]!=0) res[0][5]=a1;
  a1=1.;
  if (res[1]!=0) res[1][0]=a1;
  if (res[1]!=0) res[1][1]=a1;
  if (res[1]!=0) res[1][2]=a1;
  if (res[1]!=0) res[1][3]=a1;
  a10=(1./a10);
  a10=(-a10);
  if (res[2]!=0) res[2][0]=a10;
  a6=(a6/a7);
  a6=(-a6);
  if (res[2]!=0) res[2][1]=a6;
  return 0;
}

CASADI_SYMBOL_EXPORT int impl_jacFun(const casadi_real** arg, casadi_real** res, int* iw, casadi_real* w, void* mem){
  return casadi_f0(arg, res, iw, w, mem);
}

CASADI_SYMBOL_EXPORT void impl_jacFun_incref(void) {
}

CASADI_SYMBOL_EXPORT void impl_jacFun_decref(void) {
}

CASADI_SYMBOL_EXPORT int impl_jacFun_n_in(void) { return 4;}

CASADI_SYMBOL_EXPORT int impl_jacFun_n_out(void) { return 3;}

CASADI_SYMBOL_EXPORT const char* impl_jacFun_name_in(int i){
  switch (i) {
    case 0: return "i0";
    case 1: return "i1";
    case 2: return "i2";
    case 3: return "i3";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* impl_jacFun_name_out(int i){
  switch (i) {
    case 0: return "o0";
    case 1: return "o1";
    case 2: return "o2";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const int* impl_jacFun_sparsity_in(int i) {
  switch (i) {
    case 0: return casadi_s0;
    case 1: return casadi_s0;
    case 2: return casadi_s1;
    case 3: return casadi_s1;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const int* impl_jacFun_sparsity_out(int i) {
  switch (i) {
    case 0: return casadi_s2;
    case 1: return casadi_s3;
    case 2: return casadi_s4;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT int impl_jacFun_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w) {
  if (sz_arg) *sz_arg = 4;
  if (sz_res) *sz_res = 3;
  if (sz_iw) *sz_iw = 0;
  if (sz_w) *sz_w = 16;
  return 0;
}


#ifdef __cplusplus
} /* extern "C" */
#endif