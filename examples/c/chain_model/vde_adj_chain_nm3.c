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
  #define CASADI_PREFIX(ID) vde_adj_chain_nm3_ ## ID
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
static const int casadi_s2[19] = {15, 1, 0, 15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14};

/* vde_adj_chain_nm3:(i0[12],i1[12],i2[3],i3[15])->(o0[15]) */
static int casadi_f0(const casadi_real** arg, casadi_real** res, int* iw, casadi_real* w, void* mem) {
  casadi_real a0, a1, a10, a11, a12, a13, a14, a15, a16, a17, a2, a3, a4, a5, a6, a7, a8, a9;
  a0=1.;
  a1=3.3000000000000002e-02;
  a2=arg[0] ? arg[0][0] : 0;
  a3=sq(a2);
  a4=arg[0] ? arg[0][1] : 0;
  a5=sq(a4);
  a3=(a3+a5);
  a5=arg[0] ? arg[0][2] : 0;
  a6=sq(a5);
  a3=(a3+a6);
  a3=sqrt(a3);
  a6=(a1/a3);
  a7=(a0-a6);
  a8=3.3333333333333336e+01;
  a9=arg[1] ? arg[1][3] : 0;
  a9=(a8*a9);
  a10=(a7*a9);
  a11=(a2+a2);
  a6=(a6/a3);
  a12=arg[1] ? arg[1][5] : 0;
  a12=(a8*a12);
  a13=(a5*a12);
  a14=arg[1] ? arg[1][4] : 0;
  a8=(a8*a14);
  a14=(a4*a8);
  a13=(a13+a14);
  a14=(a2*a9);
  a13=(a13+a14);
  a6=(a6*a13);
  a3=(a3+a3);
  a6=(a6/a3);
  a11=(a11*a6);
  a10=(a10+a11);
  a11=arg[0] ? arg[0][6] : 0;
  a11=(a11-a2);
  a2=sq(a11);
  a3=arg[0] ? arg[0][7] : 0;
  a3=(a3-a4);
  a13=sq(a3);
  a2=(a2+a13);
  a13=arg[0] ? arg[0][8] : 0;
  a13=(a13-a5);
  a14=sq(a13);
  a2=(a2+a14);
  a2=sqrt(a2);
  a1=(a1/a2);
  a0=(a0-a1);
  a14=(a0*a9);
  a15=(a11+a11);
  a1=(a1/a2);
  a16=(a13*a12);
  a17=(a3*a8);
  a16=(a16+a17);
  a11=(a11*a9);
  a16=(a16+a11);
  a1=(a1*a16);
  a2=(a2+a2);
  a1=(a1/a2);
  a15=(a15*a1);
  a14=(a14+a15);
  a10=(a10+a14);
  a10=(-a10);
  if (res[0]!=0) res[0][0]=a10;
  a10=(a7*a8);
  a4=(a4+a4);
  a4=(a4*a6);
  a10=(a10+a4);
  a8=(a0*a8);
  a3=(a3+a3);
  a3=(a3*a1);
  a8=(a8+a3);
  a10=(a10+a8);
  a10=(-a10);
  if (res[0]!=0) res[0][1]=a10;
  a7=(a7*a12);
  a5=(a5+a5);
  a5=(a5*a6);
  a7=(a7+a5);
  a0=(a0*a12);
  a13=(a13+a13);
  a13=(a13*a1);
  a0=(a0+a13);
  a7=(a7+a0);
  a7=(-a7);
  if (res[0]!=0) res[0][2]=a7;
  a7=arg[1] ? arg[1][0] : 0;
  if (res[0]!=0) res[0][3]=a7;
  a7=arg[1] ? arg[1][1] : 0;
  if (res[0]!=0) res[0][4]=a7;
  a7=arg[1] ? arg[1][2] : 0;
  if (res[0]!=0) res[0][5]=a7;
  if (res[0]!=0) res[0][6]=a14;
  if (res[0]!=0) res[0][7]=a8;
  if (res[0]!=0) res[0][8]=a0;
  a0=arg[1] ? arg[1][6] : 0;
  if (res[0]!=0) res[0][9]=a0;
  a0=arg[1] ? arg[1][7] : 0;
  if (res[0]!=0) res[0][10]=a0;
  a0=arg[1] ? arg[1][8] : 0;
  if (res[0]!=0) res[0][11]=a0;
  a0=arg[1] ? arg[1][9] : 0;
  if (res[0]!=0) res[0][12]=a0;
  a0=arg[1] ? arg[1][10] : 0;
  if (res[0]!=0) res[0][13]=a0;
  a0=arg[1] ? arg[1][11] : 0;
  if (res[0]!=0) res[0][14]=a0;
  return 0;
}

CASADI_SYMBOL_EXPORT int vde_adj_chain_nm3(const casadi_real** arg, casadi_real** res, int* iw, casadi_real* w, void* mem){
  return casadi_f0(arg, res, iw, w, mem);
}

CASADI_SYMBOL_EXPORT void vde_adj_chain_nm3_incref(void) {
}

CASADI_SYMBOL_EXPORT void vde_adj_chain_nm3_decref(void) {
}

CASADI_SYMBOL_EXPORT int vde_adj_chain_nm3_n_in(void) { return 4;}

CASADI_SYMBOL_EXPORT int vde_adj_chain_nm3_n_out(void) { return 1;}

CASADI_SYMBOL_EXPORT const char* vde_adj_chain_nm3_name_in(int i){
  switch (i) {
    case 0: return "i0";
    case 1: return "i1";
    case 2: return "i2";
    case 3: return "i3";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* vde_adj_chain_nm3_name_out(int i){
  switch (i) {
    case 0: return "o0";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const int* vde_adj_chain_nm3_sparsity_in(int i) {
  switch (i) {
    case 0: return casadi_s0;
    case 1: return casadi_s0;
    case 2: return casadi_s1;
    case 3: return casadi_s2;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const int* vde_adj_chain_nm3_sparsity_out(int i) {
  switch (i) {
    case 0: return casadi_s2;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT int vde_adj_chain_nm3_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w) {
  if (sz_arg) *sz_arg = 4;
  if (sz_res) *sz_res = 1;
  if (sz_iw) *sz_iw = 0;
  if (sz_w) *sz_w = 18;
  return 0;
}


#ifdef __cplusplus
} /* extern "C" */
#endif
