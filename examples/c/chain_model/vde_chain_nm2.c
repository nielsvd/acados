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
  #define CASADI_PREFIX(ID) vde_chain_nm2_ ## ID
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

static const int casadi_s0[10] = {6, 1, 0, 6, 0, 1, 2, 3, 4, 5};
static const int casadi_s1[45] = {6, 6, 0, 6, 12, 18, 24, 30, 36, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5};
static const int casadi_s2[24] = {6, 3, 0, 6, 12, 18, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5};
static const int casadi_s3[7] = {3, 1, 0, 3, 0, 1, 2};
static const int casadi_s4[13] = {9, 1, 0, 9, 0, 1, 2, 3, 4, 5, 6, 7, 8};

/* vde_chain_nm2:(i0[6],i1[6x6],i2[6x3],i3[3],i4[9])->(o0[6],o1[6x6],o2[6x3]) */
static int casadi_f0(const casadi_real** arg, casadi_real** res, int* iw, casadi_real* w, void* mem) {
  casadi_real a0, a1, a2;
  a0=arg[0] ? arg[0][3] : 0;
  if (res[0]!=0) res[0][0]=a0;
  a0=arg[0] ? arg[0][4] : 0;
  if (res[0]!=0) res[0][1]=a0;
  a0=arg[0] ? arg[0][5] : 0;
  if (res[0]!=0) res[0][2]=a0;
  a0=arg[3] ? arg[3][0] : 0;
  if (res[0]!=0) res[0][3]=a0;
  a0=arg[3] ? arg[3][1] : 0;
  if (res[0]!=0) res[0][4]=a0;
  a0=arg[3] ? arg[3][2] : 0;
  if (res[0]!=0) res[0][5]=a0;
  a0=arg[1] ? arg[1][3] : 0;
  if (res[1]!=0) res[1][0]=a0;
  a0=arg[1] ? arg[1][4] : 0;
  if (res[1]!=0) res[1][1]=a0;
  a0=arg[1] ? arg[1][5] : 0;
  if (res[1]!=0) res[1][2]=a0;
  a0=0.;
  if (res[1]!=0) res[1][3]=a0;
  if (res[1]!=0) res[1][4]=a0;
  if (res[1]!=0) res[1][5]=a0;
  a1=arg[1] ? arg[1][9] : 0;
  if (res[1]!=0) res[1][6]=a1;
  a1=arg[1] ? arg[1][10] : 0;
  if (res[1]!=0) res[1][7]=a1;
  a1=arg[1] ? arg[1][11] : 0;
  if (res[1]!=0) res[1][8]=a1;
  if (res[1]!=0) res[1][9]=a0;
  if (res[1]!=0) res[1][10]=a0;
  if (res[1]!=0) res[1][11]=a0;
  a1=arg[1] ? arg[1][15] : 0;
  if (res[1]!=0) res[1][12]=a1;
  a1=arg[1] ? arg[1][16] : 0;
  if (res[1]!=0) res[1][13]=a1;
  a1=arg[1] ? arg[1][17] : 0;
  if (res[1]!=0) res[1][14]=a1;
  if (res[1]!=0) res[1][15]=a0;
  if (res[1]!=0) res[1][16]=a0;
  if (res[1]!=0) res[1][17]=a0;
  a1=arg[1] ? arg[1][21] : 0;
  if (res[1]!=0) res[1][18]=a1;
  a1=arg[1] ? arg[1][22] : 0;
  if (res[1]!=0) res[1][19]=a1;
  a1=arg[1] ? arg[1][23] : 0;
  if (res[1]!=0) res[1][20]=a1;
  if (res[1]!=0) res[1][21]=a0;
  if (res[1]!=0) res[1][22]=a0;
  if (res[1]!=0) res[1][23]=a0;
  a1=arg[1] ? arg[1][27] : 0;
  if (res[1]!=0) res[1][24]=a1;
  a1=arg[1] ? arg[1][28] : 0;
  if (res[1]!=0) res[1][25]=a1;
  a1=arg[1] ? arg[1][29] : 0;
  if (res[1]!=0) res[1][26]=a1;
  if (res[1]!=0) res[1][27]=a0;
  if (res[1]!=0) res[1][28]=a0;
  if (res[1]!=0) res[1][29]=a0;
  a1=arg[1] ? arg[1][33] : 0;
  if (res[1]!=0) res[1][30]=a1;
  a1=arg[1] ? arg[1][34] : 0;
  if (res[1]!=0) res[1][31]=a1;
  a1=arg[1] ? arg[1][35] : 0;
  if (res[1]!=0) res[1][32]=a1;
  if (res[1]!=0) res[1][33]=a0;
  if (res[1]!=0) res[1][34]=a0;
  if (res[1]!=0) res[1][35]=a0;
  a1=arg[2] ? arg[2][3] : 0;
  if (res[2]!=0) res[2][0]=a1;
  a1=arg[2] ? arg[2][4] : 0;
  if (res[2]!=0) res[2][1]=a1;
  a1=arg[2] ? arg[2][5] : 0;
  if (res[2]!=0) res[2][2]=a1;
  a1=1.;
  if (res[2]!=0) res[2][3]=a1;
  if (res[2]!=0) res[2][4]=a0;
  if (res[2]!=0) res[2][5]=a0;
  a2=arg[2] ? arg[2][9] : 0;
  if (res[2]!=0) res[2][6]=a2;
  a2=arg[2] ? arg[2][10] : 0;
  if (res[2]!=0) res[2][7]=a2;
  a2=arg[2] ? arg[2][11] : 0;
  if (res[2]!=0) res[2][8]=a2;
  if (res[2]!=0) res[2][9]=a0;
  if (res[2]!=0) res[2][10]=a1;
  if (res[2]!=0) res[2][11]=a0;
  a2=arg[2] ? arg[2][15] : 0;
  if (res[2]!=0) res[2][12]=a2;
  a2=arg[2] ? arg[2][16] : 0;
  if (res[2]!=0) res[2][13]=a2;
  a2=arg[2] ? arg[2][17] : 0;
  if (res[2]!=0) res[2][14]=a2;
  if (res[2]!=0) res[2][15]=a0;
  if (res[2]!=0) res[2][16]=a0;
  if (res[2]!=0) res[2][17]=a1;
  return 0;
}

CASADI_SYMBOL_EXPORT int vde_chain_nm2(const casadi_real** arg, casadi_real** res, int* iw, casadi_real* w, void* mem){
  return casadi_f0(arg, res, iw, w, mem);
}

CASADI_SYMBOL_EXPORT void vde_chain_nm2_incref(void) {
}

CASADI_SYMBOL_EXPORT void vde_chain_nm2_decref(void) {
}

CASADI_SYMBOL_EXPORT int vde_chain_nm2_n_in(void) { return 5;}

CASADI_SYMBOL_EXPORT int vde_chain_nm2_n_out(void) { return 3;}

CASADI_SYMBOL_EXPORT const char* vde_chain_nm2_name_in(int i){
  switch (i) {
    case 0: return "i0";
    case 1: return "i1";
    case 2: return "i2";
    case 3: return "i3";
    case 4: return "i4";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* vde_chain_nm2_name_out(int i){
  switch (i) {
    case 0: return "o0";
    case 1: return "o1";
    case 2: return "o2";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const int* vde_chain_nm2_sparsity_in(int i) {
  switch (i) {
    case 0: return casadi_s0;
    case 1: return casadi_s1;
    case 2: return casadi_s2;
    case 3: return casadi_s3;
    case 4: return casadi_s4;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const int* vde_chain_nm2_sparsity_out(int i) {
  switch (i) {
    case 0: return casadi_s0;
    case 1: return casadi_s1;
    case 2: return casadi_s2;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT int vde_chain_nm2_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w) {
  if (sz_arg) *sz_arg = 5;
  if (sz_res) *sz_res = 3;
  if (sz_iw) *sz_iw = 0;
  if (sz_w) *sz_w = 3;
  return 0;
}


#ifdef __cplusplus
} /* extern "C" */
#endif
