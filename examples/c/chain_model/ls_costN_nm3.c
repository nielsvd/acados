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
  #define CASADI_PREFIX(ID) ls_costN_nm3_ ## ID
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

static const int casadi_s0[16] = {12, 1, 0, 12, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
static const int casadi_s1[7] = {3, 1, 0, 3, 0, 1, 2};
static const int casadi_s2[27] = {12, 12, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};

/* ls_costN_nm3:(i0[12],i1[3],i2[12])->(o0[12],o1[12x12,12nz]) */
static int casadi_f0(const casadi_real** arg, casadi_real** res, int* iw, casadi_real* w, void* mem) {
  casadi_real a0, a1;
  a0=arg[0] ? arg[0][0] : 0;
  a1=arg[2] ? arg[2][0] : 0;
  a0=(a0-a1);
  if (res[0]!=0) res[0][0]=a0;
  a0=arg[0] ? arg[0][1] : 0;
  a1=arg[2] ? arg[2][1] : 0;
  a0=(a0-a1);
  if (res[0]!=0) res[0][1]=a0;
  a0=arg[0] ? arg[0][2] : 0;
  a1=arg[2] ? arg[2][2] : 0;
  a0=(a0-a1);
  if (res[0]!=0) res[0][2]=a0;
  a0=arg[0] ? arg[0][3] : 0;
  a1=arg[2] ? arg[2][3] : 0;
  a0=(a0-a1);
  if (res[0]!=0) res[0][3]=a0;
  a0=arg[0] ? arg[0][4] : 0;
  a1=arg[2] ? arg[2][4] : 0;
  a0=(a0-a1);
  if (res[0]!=0) res[0][4]=a0;
  a0=arg[0] ? arg[0][5] : 0;
  a1=arg[2] ? arg[2][5] : 0;
  a0=(a0-a1);
  if (res[0]!=0) res[0][5]=a0;
  a0=arg[0] ? arg[0][6] : 0;
  a1=arg[2] ? arg[2][6] : 0;
  a0=(a0-a1);
  if (res[0]!=0) res[0][6]=a0;
  a0=arg[0] ? arg[0][7] : 0;
  a1=arg[2] ? arg[2][7] : 0;
  a0=(a0-a1);
  if (res[0]!=0) res[0][7]=a0;
  a0=arg[0] ? arg[0][8] : 0;
  a1=arg[2] ? arg[2][8] : 0;
  a0=(a0-a1);
  if (res[0]!=0) res[0][8]=a0;
  a0=arg[0] ? arg[0][9] : 0;
  a1=arg[2] ? arg[2][9] : 0;
  a0=(a0-a1);
  if (res[0]!=0) res[0][9]=a0;
  a0=arg[0] ? arg[0][10] : 0;
  a1=arg[2] ? arg[2][10] : 0;
  a0=(a0-a1);
  if (res[0]!=0) res[0][10]=a0;
  a0=arg[0] ? arg[0][11] : 0;
  a1=arg[2] ? arg[2][11] : 0;
  a0=(a0-a1);
  if (res[0]!=0) res[0][11]=a0;
  a0=1.;
  if (res[1]!=0) res[1][0]=a0;
  if (res[1]!=0) res[1][1]=a0;
  if (res[1]!=0) res[1][2]=a0;
  if (res[1]!=0) res[1][3]=a0;
  if (res[1]!=0) res[1][4]=a0;
  if (res[1]!=0) res[1][5]=a0;
  if (res[1]!=0) res[1][6]=a0;
  if (res[1]!=0) res[1][7]=a0;
  if (res[1]!=0) res[1][8]=a0;
  if (res[1]!=0) res[1][9]=a0;
  if (res[1]!=0) res[1][10]=a0;
  if (res[1]!=0) res[1][11]=a0;
  return 0;
}

CASADI_SYMBOL_EXPORT int ls_costN_nm3(const casadi_real** arg, casadi_real** res, int* iw, casadi_real* w, void* mem){
  return casadi_f0(arg, res, iw, w, mem);
}

CASADI_SYMBOL_EXPORT void ls_costN_nm3_incref(void) {
}

CASADI_SYMBOL_EXPORT void ls_costN_nm3_decref(void) {
}

CASADI_SYMBOL_EXPORT int ls_costN_nm3_n_in(void) { return 3;}

CASADI_SYMBOL_EXPORT int ls_costN_nm3_n_out(void) { return 2;}

CASADI_SYMBOL_EXPORT const char* ls_costN_nm3_name_in(int i){
  switch (i) {
    case 0: return "i0";
    case 1: return "i1";
    case 2: return "i2";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* ls_costN_nm3_name_out(int i){
  switch (i) {
    case 0: return "o0";
    case 1: return "o1";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const int* ls_costN_nm3_sparsity_in(int i) {
  switch (i) {
    case 0: return casadi_s0;
    case 1: return casadi_s1;
    case 2: return casadi_s0;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const int* ls_costN_nm3_sparsity_out(int i) {
  switch (i) {
    case 0: return casadi_s0;
    case 1: return casadi_s2;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT int ls_costN_nm3_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w) {
  if (sz_arg) *sz_arg = 3;
  if (sz_res) *sz_res = 2;
  if (sz_iw) *sz_iw = 0;
  if (sz_w) *sz_w = 2;
  return 0;
}


#ifdef __cplusplus
} /* extern "C" */
#endif
