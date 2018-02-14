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
  #define CASADI_PREFIX(ID) vde_adj_chain_nm4_ ## ID
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

static const int casadi_s0[22] = {18, 1, 0, 18, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17};
static const int casadi_s1[7] = {3, 1, 0, 3, 0, 1, 2};
static const int casadi_s2[25] = {21, 1, 0, 21, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20};

/* vde_adj_chain_nm4:(i0[18],i1[18],i2[3],i3[21])->(o0[21]) */
static int casadi_f0(const casadi_real** arg, casadi_real** res, int* iw, casadi_real* w, void* mem) {
  casadi_real a0, a1, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a2, a20, a21, a22, a23, a24, a25, a26, a27, a3, a4, a5, a6, a7, a8, a9;
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
  a14=(a8*a14);
  a15=(a4*a14);
  a13=(a13+a15);
  a15=(a2*a9);
  a13=(a13+a15);
  a6=(a6*a13);
  a3=(a3+a3);
  a6=(a6/a3);
  a11=(a11*a6);
  a10=(a10+a11);
  a11=arg[0] ? arg[0][6] : 0;
  a2=(a11-a2);
  a3=sq(a2);
  a13=arg[0] ? arg[0][7] : 0;
  a15=(a13-a4);
  a16=sq(a15);
  a3=(a3+a16);
  a16=arg[0] ? arg[0][8] : 0;
  a17=(a16-a5);
  a18=sq(a17);
  a3=(a3+a18);
  a3=sqrt(a3);
  a18=(a1/a3);
  a19=(a0-a18);
  a20=arg[1] ? arg[1][9] : 0;
  a20=(a8*a20);
  a9=(a9-a20);
  a21=(a19*a9);
  a22=(a2+a2);
  a18=(a18/a3);
  a23=arg[1] ? arg[1][11] : 0;
  a23=(a8*a23);
  a24=(a12-a23);
  a25=(a17*a24);
  a26=arg[1] ? arg[1][10] : 0;
  a8=(a8*a26);
  a26=(a14-a8);
  a27=(a15*a26);
  a25=(a25+a27);
  a2=(a2*a9);
  a25=(a25+a2);
  a18=(a18*a25);
  a3=(a3+a3);
  a18=(a18/a3);
  a22=(a22*a18);
  a21=(a21+a22);
  a10=(a10+a21);
  a10=(-a10);
  if (res[0]!=0) res[0][0]=a10;
  a14=(a7*a14);
  a4=(a4+a4);
  a4=(a4*a6);
  a14=(a14+a4);
  a26=(a19*a26);
  a15=(a15+a15);
  a15=(a15*a18);
  a26=(a26+a15);
  a14=(a14+a26);
  a14=(-a14);
  if (res[0]!=0) res[0][1]=a14;
  a7=(a7*a12);
  a5=(a5+a5);
  a5=(a5*a6);
  a7=(a7+a5);
  a19=(a19*a24);
  a17=(a17+a17);
  a17=(a17*a18);
  a19=(a19+a17);
  a7=(a7+a19);
  a7=(-a7);
  if (res[0]!=0) res[0][2]=a7;
  a7=arg[1] ? arg[1][0] : 0;
  if (res[0]!=0) res[0][3]=a7;
  a7=arg[1] ? arg[1][1] : 0;
  if (res[0]!=0) res[0][4]=a7;
  a7=arg[1] ? arg[1][2] : 0;
  if (res[0]!=0) res[0][5]=a7;
  a7=arg[0] ? arg[0][12] : 0;
  a7=(a7-a11);
  a11=sq(a7);
  a17=arg[0] ? arg[0][13] : 0;
  a17=(a17-a13);
  a13=sq(a17);
  a11=(a11+a13);
  a13=arg[0] ? arg[0][14] : 0;
  a13=(a13-a16);
  a16=sq(a13);
  a11=(a11+a16);
  a11=sqrt(a11);
  a1=(a1/a11);
  a0=(a0-a1);
  a16=(a0*a20);
  a18=(a7+a7);
  a1=(a1/a11);
  a24=(a13*a23);
  a5=(a17*a8);
  a24=(a24+a5);
  a7=(a7*a20);
  a24=(a24+a7);
  a1=(a1*a24);
  a11=(a11+a11);
  a1=(a1/a11);
  a18=(a18*a1);
  a16=(a16+a18);
  a21=(a21-a16);
  if (res[0]!=0) res[0][6]=a21;
  a8=(a0*a8);
  a17=(a17+a17);
  a17=(a17*a1);
  a8=(a8+a17);
  a26=(a26-a8);
  if (res[0]!=0) res[0][7]=a26;
  a0=(a0*a23);
  a13=(a13+a13);
  a13=(a13*a1);
  a0=(a0+a13);
  a19=(a19-a0);
  if (res[0]!=0) res[0][8]=a19;
  a19=arg[1] ? arg[1][6] : 0;
  if (res[0]!=0) res[0][9]=a19;
  a19=arg[1] ? arg[1][7] : 0;
  if (res[0]!=0) res[0][10]=a19;
  a19=arg[1] ? arg[1][8] : 0;
  if (res[0]!=0) res[0][11]=a19;
  if (res[0]!=0) res[0][12]=a16;
  if (res[0]!=0) res[0][13]=a8;
  if (res[0]!=0) res[0][14]=a0;
  a0=arg[1] ? arg[1][12] : 0;
  if (res[0]!=0) res[0][15]=a0;
  a0=arg[1] ? arg[1][13] : 0;
  if (res[0]!=0) res[0][16]=a0;
  a0=arg[1] ? arg[1][14] : 0;
  if (res[0]!=0) res[0][17]=a0;
  a0=arg[1] ? arg[1][15] : 0;
  if (res[0]!=0) res[0][18]=a0;
  a0=arg[1] ? arg[1][16] : 0;
  if (res[0]!=0) res[0][19]=a0;
  a0=arg[1] ? arg[1][17] : 0;
  if (res[0]!=0) res[0][20]=a0;
  return 0;
}

CASADI_SYMBOL_EXPORT int vde_adj_chain_nm4(const casadi_real** arg, casadi_real** res, int* iw, casadi_real* w, void* mem){
  return casadi_f0(arg, res, iw, w, mem);
}

CASADI_SYMBOL_EXPORT void vde_adj_chain_nm4_incref(void) {
}

CASADI_SYMBOL_EXPORT void vde_adj_chain_nm4_decref(void) {
}

CASADI_SYMBOL_EXPORT int vde_adj_chain_nm4_n_in(void) { return 4;}

CASADI_SYMBOL_EXPORT int vde_adj_chain_nm4_n_out(void) { return 1;}

CASADI_SYMBOL_EXPORT const char* vde_adj_chain_nm4_name_in(int i){
  switch (i) {
    case 0: return "i0";
    case 1: return "i1";
    case 2: return "i2";
    case 3: return "i3";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* vde_adj_chain_nm4_name_out(int i){
  switch (i) {
    case 0: return "o0";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const int* vde_adj_chain_nm4_sparsity_in(int i) {
  switch (i) {
    case 0: return casadi_s0;
    case 1: return casadi_s0;
    case 2: return casadi_s1;
    case 3: return casadi_s2;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const int* vde_adj_chain_nm4_sparsity_out(int i) {
  switch (i) {
    case 0: return casadi_s2;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT int vde_adj_chain_nm4_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w) {
  if (sz_arg) *sz_arg = 4;
  if (sz_res) *sz_res = 1;
  if (sz_iw) *sz_iw = 0;
  if (sz_w) *sz_w = 28;
  return 0;
}


#ifdef __cplusplus
} /* extern "C" */
#endif
