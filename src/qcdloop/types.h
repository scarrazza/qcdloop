//
// QCDLoop 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch
//          Keith Ellis: keith.ellis@durham.ac.uk
//          Giulia Zanderighi: giulia.zanderighi@cern.ch

#pragma once

#ifdef __x86_64__
extern "C" { // for gcc4.7 compatibility
#include <quadmath.h>
}
#else
//#endif
//#ifdef __aarch64__
#include <stdlib.h>
#include <complex.h>
#include <tgmath.h>
#include <float.h>
typedef long double __float128;
using __complex128 = long double _Complex  ;
#define acosq acosl
#define acoshq acoshl 
#define asinq asinl
#define asinhq asinhl 
#define atanq atanl
#define atanhq atanhl 
#define cbrtq cbrtl
#define ceilq ceill
#define coshq coshl
#define cosq cosl
#define erfq erfl
#define erfcq erfcl
#define expq expl
#define expm1q expm1l
#define fabsq fabsl
#define floorq floorl
#define lgammaq lgammal
#define logq logl
#define log10q log10l
#define log2q log2l
#define log1pq log1pl
#define nearbyintq nearbyintl
#define rintq rintl
#define roundq roundl
#define sinhq sinhl
#define sinq sinl
#define sqrtq sqrtl
#define tanq tanl
#define tanhq tanhl
#define tgammaq tgammal
#define truncq truncl

#define atan2q atan2l
#define copysignq copysignl
#define hypotq hypotl
#define ldexpq ldexpl
#define nextafterq nextafterl
#define powq powl
#define remainderq remainderl


#define cabsq cabsl
#define crealq creall
#define cimagq cimagl
#define cargq cargl
#define conjq conjl
#define cacosq cacosl
#define cacoshq cacoshl
#define casinq casinl
#define casinhq casinhl
#define catanq catanl
#define catanhq catanhl
#define ccosq ccosl
#define ccoshq ccoshl
#define cexpq cexpl
#define clogq clogl
#define clog10q clog10l
#define csinq csinl
#define csinhq csinhl
#define csqrtq csqrtl
#define ctanq ctanl
#define ctanhq ctanhl
#define cpowq cpowl
#endif
#include <complex>

#define UNUSED(expr) (void)(expr)

/*!
 * Defines the standard types for:
 * double, complex, qdouble and qcomplex
 */
namespace ql
{
  enum { ep0, ep1, ep2 };
  typedef __float128 qdouble;
  typedef __complex128 qcomplex;
  typedef std::complex<double> complex;  
  enum Code { red = 31, green = 32, yellow = 33, blue = 34, def = 39};
}

namespace std
{
  //! implementation of operator<< for qdouble
  ostream& operator<<(std::ostream& out, ql::qdouble f);

  //! implementation of operator<< for qcomplex
  ostream& operator<<(std::ostream& out, ql::qcomplex f);

  //! implementation of operator<< for Code
  std::ostream& operator<<(std::ostream& os, ql::Code code);
}
