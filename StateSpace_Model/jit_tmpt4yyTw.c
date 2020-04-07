/* This file was automatically generated by CasADi.
   The CasADi copyright holders make no ownership claim of its contents. */
#ifdef __cplusplus
extern "C" {
#endif

/* How to prefix internal symbols */
#ifdef CASADI_CODEGEN_PREFIX
  #define CASADI_NAMESPACE_CONCAT(NS, ID) _CASADI_NAMESPACE_CONCAT(NS, ID)
  #define _CASADI_NAMESPACE_CONCAT(NS, ID) NS ## ID
  #define CASADI_PREFIX(ID) CASADI_NAMESPACE_CONCAT(CODEGEN_PREFIX, ID)
#else
  #define CASADI_PREFIX(ID) jit_tmpt4yyTw_ ## ID
#endif

#include <math.h>

#ifndef casadi_real
#define casadi_real double
#endif

#ifndef casadi_int
#define casadi_int long long int
#endif

/* Add prefix to internal symbols */
#define casadi_f0 CASADI_PREFIX(f0)
#define casadi_f1 CASADI_PREFIX(f1)
#define casadi_f2 CASADI_PREFIX(f2)
#define casadi_s0 CASADI_PREFIX(s0)
#define casadi_s1 CASADI_PREFIX(s1)
#define casadi_s2 CASADI_PREFIX(s2)
#define casadi_s3 CASADI_PREFIX(s3)
#define casadi_sq CASADI_PREFIX(sq)

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

casadi_real casadi_sq(casadi_real x) { return x*x;}

static const casadi_int casadi_s0[3] = {0, 0, 0};
static const casadi_int casadi_s1[11] = {7, 1, 0, 7, 0, 1, 2, 3, 4, 5, 6};
static const casadi_int casadi_s2[4] = {0, 1, 0, 0};
static const casadi_int casadi_s3[15] = {11, 1, 0, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

/* dae:(t[],x[7],z[0],p[11],rx[0],rz[0],rp[0])->(ode[7],alg[0],quad[0],rode[0],ralg[0],rquad[0]) */
static int casadi_f0(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem) {
  casadi_real a0, a1, a10, a2, a3, a4, a5, a6, a7, a8, a9;
  a0=arg[1]? arg[1][1] : 0;
  if (res[0]!=0) res[0][0]=a0;
  a0=arg[3]? arg[3][5] : 0;
  a1=1.2705996954518719e+00;
  a1=(a0/a1);
  a2=arg[3]? arg[3][6] : 0;
  a1=(a1/a2);
  a3=arg[1]? arg[1][3] : 0;
  a1=(a1*a3);
  a4=7.8702993836651547e-01;
  a5=1.0986982031656709e+01;
  a6=(a0/a2);
  a6=(a5-a6);
  a4=(a4*a6);
  a6=arg[1]? arg[1][4] : 0;
  a4=(a4*a6);
  a1=(a1+a4);
  a4=8.6470837912087912e+00;
  a7=arg[1]? arg[1][0] : 0;
  a8=arg[3]? arg[3][9] : 0;
  a9=(a7+a8);
  a4=(a4*a9);
  a1=(a1+a4);
  if (res[0]!=0) res[0][1]=a1;
  a1=arg[3]? arg[3][7] : 0;
  a4=arg[3]? arg[3][8] : 0;
  a9=(a1/a4);
  a10=(a0/a2);
  a5=(a5-a10);
  a9=(a9*a5);
  a9=(a9*a6);
  a1=(a1*a0);
  a1=(a1/a2);
  a1=(a1/a4);
  a1=(a1*a3);
  a9=(a9+a1);
  a1=arg[1]? arg[1][2] : 0;
  a4=(a1/a4);
  a9=(a9-a4);
  if (res[0]!=0) res[0][2]=a9;
  a9=arg[3]? arg[3][4] : 0;
  a4=arg[3]? arg[3][0] : 0;
  a0=(a9*a4);
  a5=arg[3]? arg[3][1] : 0;
  a5=(a5+a9);
  a7=(a7+a8);
  a5=(a5*a7);
  a7=arg[3]? arg[3][2] : 0;
  a8=arg[1]? arg[1][5] : 0;
  a7=(a7*a8);
  a5=(a5+a7);
  a0=(a0-a5);
  a0=(a0+a1);
  a0=(a0/a2);
  a1=(a3/a2);
  a0=(a0-a1);
  if (res[0]!=0) res[0][3]=a0;
  a3=(a3-a6);
  a3=(a3/a2);
  if (res[0]!=0) res[0][4]=a3;
  a3=5.0000000000000000e-01;
  a2=arg[1]? arg[1][6] : 0;
  a6=arg[3]? arg[3][10] : 0;
  a0=(a2+a6);
  a0=(a4-a0);
  a1=2.5000000000000000e-01;
  a0=(a0/a1);
  a5=arg[3]? arg[3][3] : 0;
  a0=(a0-a5);
  a0=casadi_sq(a0);
  a7=5.0000000000000003e-02;
  a8=casadi_sq(a5);
  a8=(a7*a8);
  a0=(a0+a8);
  a0=sqrt(a0);
  a0=(a3*a0);
  a8=(a2+a6);
  a8=(a4-a8);
  a8=(a8/a1);
  a8=(a8+a5);
  a8=casadi_sq(a8);
  a5=casadi_sq(a5);
  a7=(a7*a5);
  a8=(a8+a7);
  a8=sqrt(a8);
  a3=(a3*a8);
  a0=(a0-a3);
  a3=(a2+a6);
  a3=(a4-a3);
  a3=(a3/a1);
  a0=(a0+a3);
  if (res[0]!=0) res[0][5]=a0;
  a2=(a2+a6);
  a4=(a4-a2);
  a4=(a4/a1);
  if (res[0]!=0) res[0][6]=a4;
  return 0;
}

CASADI_SYMBOL_EXPORT int dae(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem){
  return casadi_f0(arg, res, iw, w, mem);
}

CASADI_SYMBOL_EXPORT int dae_alloc_mem(void) {
  return 0;
}

CASADI_SYMBOL_EXPORT int dae_init_mem(int mem) {
  return 0;
}

CASADI_SYMBOL_EXPORT void dae_free_mem(int mem) {
}

CASADI_SYMBOL_EXPORT int dae_checkout(void) {
  return 0;
}

CASADI_SYMBOL_EXPORT void dae_release(int mem) {
}

CASADI_SYMBOL_EXPORT void dae_incref(void) {
}

CASADI_SYMBOL_EXPORT void dae_decref(void) {
}

CASADI_SYMBOL_EXPORT casadi_int dae_n_in(void) { return 7;}

CASADI_SYMBOL_EXPORT casadi_int dae_n_out(void) { return 6;}

CASADI_SYMBOL_EXPORT casadi_real dae_default_in(casadi_int i){
  switch (i) {
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* dae_name_in(casadi_int i){
  switch (i) {
    case 0: return "t";
    case 1: return "x";
    case 2: return "z";
    case 3: return "p";
    case 4: return "rx";
    case 5: return "rz";
    case 6: return "rp";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* dae_name_out(casadi_int i){
  switch (i) {
    case 0: return "ode";
    case 1: return "alg";
    case 2: return "quad";
    case 3: return "rode";
    case 4: return "ralg";
    case 5: return "rquad";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* dae_sparsity_in(casadi_int i) {
  switch (i) {
    case 0: return casadi_s0;
    case 1: return casadi_s1;
    case 2: return casadi_s2;
    case 3: return casadi_s3;
    case 4: return casadi_s2;
    case 5: return casadi_s2;
    case 6: return casadi_s2;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* dae_sparsity_out(casadi_int i) {
  switch (i) {
    case 0: return casadi_s1;
    case 1: return casadi_s2;
    case 2: return casadi_s2;
    case 3: return casadi_s2;
    case 4: return casadi_s2;
    case 5: return casadi_s2;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT int dae_work(casadi_int *sz_arg, casadi_int* sz_res, casadi_int *sz_iw, casadi_int *sz_w) {
  if (sz_arg) *sz_arg = 7;
  if (sz_res) *sz_res = 6;
  if (sz_iw) *sz_iw = 0;
  if (sz_w) *sz_w = 0;
  return 0;
}

/* f:(x[7],z[0],p[11],t[])->(ode[7],alg[0],quad[0]) */
static int casadi_f1(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem) {
  casadi_real a0, a1, a10, a2, a3, a4, a5, a6, a7, a8, a9;
  a0=arg[0]? arg[0][1] : 0;
  if (res[0]!=0) res[0][0]=a0;
  a0=arg[2]? arg[2][5] : 0;
  a1=1.2705996954518719e+00;
  a1=(a0/a1);
  a2=arg[2]? arg[2][6] : 0;
  a1=(a1/a2);
  a3=arg[0]? arg[0][3] : 0;
  a1=(a1*a3);
  a4=7.8702993836651547e-01;
  a5=1.0986982031656709e+01;
  a6=(a0/a2);
  a6=(a5-a6);
  a4=(a4*a6);
  a6=arg[0]? arg[0][4] : 0;
  a4=(a4*a6);
  a1=(a1+a4);
  a4=8.6470837912087912e+00;
  a7=arg[0]? arg[0][0] : 0;
  a8=arg[2]? arg[2][9] : 0;
  a9=(a7+a8);
  a4=(a4*a9);
  a1=(a1+a4);
  if (res[0]!=0) res[0][1]=a1;
  a1=arg[2]? arg[2][7] : 0;
  a4=arg[2]? arg[2][8] : 0;
  a9=(a1/a4);
  a10=(a0/a2);
  a5=(a5-a10);
  a9=(a9*a5);
  a9=(a9*a6);
  a1=(a1*a0);
  a1=(a1/a2);
  a1=(a1/a4);
  a1=(a1*a3);
  a9=(a9+a1);
  a1=arg[0]? arg[0][2] : 0;
  a4=(a1/a4);
  a9=(a9-a4);
  if (res[0]!=0) res[0][2]=a9;
  a9=arg[2]? arg[2][4] : 0;
  a4=arg[2]? arg[2][0] : 0;
  a0=(a9*a4);
  a5=arg[2]? arg[2][1] : 0;
  a5=(a5+a9);
  a7=(a7+a8);
  a5=(a5*a7);
  a7=arg[2]? arg[2][2] : 0;
  a8=arg[0]? arg[0][5] : 0;
  a7=(a7*a8);
  a5=(a5+a7);
  a0=(a0-a5);
  a0=(a0+a1);
  a0=(a0/a2);
  a1=(a3/a2);
  a0=(a0-a1);
  if (res[0]!=0) res[0][3]=a0;
  a3=(a3-a6);
  a3=(a3/a2);
  if (res[0]!=0) res[0][4]=a3;
  a3=5.0000000000000000e-01;
  a2=arg[0]? arg[0][6] : 0;
  a6=arg[2]? arg[2][10] : 0;
  a0=(a2+a6);
  a0=(a4-a0);
  a1=2.5000000000000000e-01;
  a0=(a0/a1);
  a5=arg[2]? arg[2][3] : 0;
  a0=(a0-a5);
  a0=casadi_sq(a0);
  a7=5.0000000000000003e-02;
  a8=casadi_sq(a5);
  a8=(a7*a8);
  a0=(a0+a8);
  a0=sqrt(a0);
  a0=(a3*a0);
  a8=(a2+a6);
  a8=(a4-a8);
  a8=(a8/a1);
  a8=(a8+a5);
  a8=casadi_sq(a8);
  a5=casadi_sq(a5);
  a7=(a7*a5);
  a8=(a8+a7);
  a8=sqrt(a8);
  a3=(a3*a8);
  a0=(a0-a3);
  a3=(a2+a6);
  a3=(a4-a3);
  a3=(a3/a1);
  a0=(a0+a3);
  if (res[0]!=0) res[0][5]=a0;
  a2=(a2+a6);
  a4=(a4-a2);
  a4=(a4/a1);
  if (res[0]!=0) res[0][6]=a4;
  return 0;
}

CASADI_SYMBOL_EXPORT int f(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem){
  return casadi_f1(arg, res, iw, w, mem);
}

CASADI_SYMBOL_EXPORT int f_alloc_mem(void) {
  return 0;
}

CASADI_SYMBOL_EXPORT int f_init_mem(int mem) {
  return 0;
}

CASADI_SYMBOL_EXPORT void f_free_mem(int mem) {
}

CASADI_SYMBOL_EXPORT int f_checkout(void) {
  return 0;
}

CASADI_SYMBOL_EXPORT void f_release(int mem) {
}

CASADI_SYMBOL_EXPORT void f_incref(void) {
}

CASADI_SYMBOL_EXPORT void f_decref(void) {
}

CASADI_SYMBOL_EXPORT casadi_int f_n_in(void) { return 4;}

CASADI_SYMBOL_EXPORT casadi_int f_n_out(void) { return 3;}

CASADI_SYMBOL_EXPORT casadi_real f_default_in(casadi_int i){
  switch (i) {
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* f_name_in(casadi_int i){
  switch (i) {
    case 0: return "x";
    case 1: return "z";
    case 2: return "p";
    case 3: return "t";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* f_name_out(casadi_int i){
  switch (i) {
    case 0: return "ode";
    case 1: return "alg";
    case 2: return "quad";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* f_sparsity_in(casadi_int i) {
  switch (i) {
    case 0: return casadi_s1;
    case 1: return casadi_s2;
    case 2: return casadi_s3;
    case 3: return casadi_s0;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* f_sparsity_out(casadi_int i) {
  switch (i) {
    case 0: return casadi_s1;
    case 1: return casadi_s2;
    case 2: return casadi_s2;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT int f_work(casadi_int *sz_arg, casadi_int* sz_res, casadi_int *sz_iw, casadi_int *sz_w) {
  if (sz_arg) *sz_arg = 4;
  if (sz_res) *sz_res = 3;
  if (sz_iw) *sz_iw = 0;
  if (sz_w) *sz_w = 0;
  return 0;
}

/* g:(rx[0],rz[0],rp[0],x[7],z[0],p[11],t[])->(rode[0],ralg[0],rquad[0]) */
static int casadi_f2(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem) {
  return 0;
}

CASADI_SYMBOL_EXPORT int g(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem){
  return casadi_f2(arg, res, iw, w, mem);
}

CASADI_SYMBOL_EXPORT int g_alloc_mem(void) {
  return 0;
}

CASADI_SYMBOL_EXPORT int g_init_mem(int mem) {
  return 0;
}

CASADI_SYMBOL_EXPORT void g_free_mem(int mem) {
}

CASADI_SYMBOL_EXPORT int g_checkout(void) {
  return 0;
}

CASADI_SYMBOL_EXPORT void g_release(int mem) {
}

CASADI_SYMBOL_EXPORT void g_incref(void) {
}

CASADI_SYMBOL_EXPORT void g_decref(void) {
}

CASADI_SYMBOL_EXPORT casadi_int g_n_in(void) { return 7;}

CASADI_SYMBOL_EXPORT casadi_int g_n_out(void) { return 3;}

CASADI_SYMBOL_EXPORT casadi_real g_default_in(casadi_int i){
  switch (i) {
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* g_name_in(casadi_int i){
  switch (i) {
    case 0: return "rx";
    case 1: return "rz";
    case 2: return "rp";
    case 3: return "x";
    case 4: return "z";
    case 5: return "p";
    case 6: return "t";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* g_name_out(casadi_int i){
  switch (i) {
    case 0: return "rode";
    case 1: return "ralg";
    case 2: return "rquad";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* g_sparsity_in(casadi_int i) {
  switch (i) {
    case 0: return casadi_s2;
    case 1: return casadi_s2;
    case 2: return casadi_s2;
    case 3: return casadi_s1;
    case 4: return casadi_s2;
    case 5: return casadi_s3;
    case 6: return casadi_s0;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* g_sparsity_out(casadi_int i) {
  switch (i) {
    case 0: return casadi_s2;
    case 1: return casadi_s2;
    case 2: return casadi_s2;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT int g_work(casadi_int *sz_arg, casadi_int* sz_res, casadi_int *sz_iw, casadi_int *sz_w) {
  if (sz_arg) *sz_arg = 7;
  if (sz_res) *sz_res = 3;
  if (sz_iw) *sz_iw = 0;
  if (sz_w) *sz_w = 0;
  return 0;
}


#ifdef __cplusplus
} /* extern "C" */
#endif
