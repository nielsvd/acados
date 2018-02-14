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
  #define CASADI_PREFIX(ID) vde_adj_chain_nm2_ ## ID
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

static const int casadi_s0[10] = {6, 1, 0, 6, 0, 1, 2, 3, 4, 5};
static const int casadi_s1[7] = {3, 1, 0, 3, 0, 1, 2};
static const int casadi_s2[13] = {9, 1, 0, 9, 0, 1, 2, 3, 4, 5, 6, 7, 8};

/* vde_adj_chain_nm2:(i0[6],i1[6],i2[3],i3[9])->(o0[9]) */
static int casadi_f0(const casadi_real** arg, casadi_real** res, int* iw, casadi_real* w, void* mem) {
  casadi_real a0;
  a0=0.;
  if (res[0]!=0) res[0][0]=a0;
  if (res[0]!=0) res[0][1]=a0;
  if (res[0]!=0) res[0][2]=a0;
  a0=arg[1] ? arg[1][0] : 0;
  if (res[0]!=0) res[0][3]=a0;
  a0=arg[1] ? arg[1][1] : 0;
  if (res[0]!=0) res[0][4]=a0;
  a0=arg[1] ? arg[1][2] : 0;
  if (res[0]!=0) res[0][5]=a0;
  a0=arg[1] ? arg[1][3] : 0;
  if (res[0]!=0) res[0][6]=a0;
  a0=arg[1] ? arg[1][4] : 0;
  if (res[0]!=0) res[0][7]=a0;
  a0=arg[1] ? arg[1][5] : 0;
  if (res[0]!=0) res[0][8]=a0;
  return 0;
}

CASADI_SYMBOL_EXPORT int vde_adj_chain_nm2(const casadi_real** arg, casadi_real** res, int* iw, casadi_real* w, void* mem){
  return casadi_f0(arg, res, iw, w, mem);
}

CASADI_SYMBOL_EXPORT void vde_adj_chain_nm2_incref(void) {
}

CASADI_SYMBOL_EXPORT void vde_adj_chain_nm2_decref(void) {
}

CASADI_SYMBOL_EXPORT int vde_adj_chain_nm2_n_in(void) { return 4;}

CASADI_SYMBOL_EXPORT int vde_adj_chain_nm2_n_out(void) { return 1;}

CASADI_SYMBOL_EXPORT const char* vde_adj_chain_nm2_name_in(int i){
  switch (i) {
    case 0: return "i0";
    case 1: return "i1";
    case 2: return "i2";
    case 3: return "i3";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* vde_adj_chain_nm2_name_out(int i){
  switch (i) {
    case 0: return "o0";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const int* vde_adj_chain_nm2_sparsity_in(int i) {
  switch (i) {
    case 0: return casadi_s0;
    case 1: return casadi_s0;
    case 2: return casadi_s1;
    case 3: return casadi_s2;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const int* vde_adj_chain_nm2_sparsity_out(int i) {
  switch (i) {
    case 0: return casadi_s2;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT int vde_adj_chain_nm2_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w) {
  if (sz_arg) *sz_arg = 4;
  if (sz_res) *sz_res = 1;
  if (sz_iw) *sz_iw = 0;
  if (sz_w) *sz_w = 1;
  return 0;
}


#ifdef __cplusplus
} /* extern "C" */
#endif
