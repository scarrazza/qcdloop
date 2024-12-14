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
#define acos std::acos
#define acosh std::acoshl 
#define asin std::asin
#define asinh std::asinhl 
#define atan std::atan
#define atanh std::atanhl 
#define cbrt std::cbrt
#define ceil std::ceil
#define cosh std::cosh
#define cos std::cos
#define erf std::erf
#define erfc std::erfc
#define exp std::exp
#define expm1 std::expm1
#define fabs std::fabs
#define floor std::floor
#define lgamma std::lgamma
#define log std::log
#define log10 std::log10
#define log2 std::log2
#define log1p std::log1p
#define nearbyint std::nearbyint
#define rint std::rint
#define round std::round
#define sinh std::sinh
#define sin std::sin
#define sqrt std::sqrt
#define tan std::tan
#define tanh std::tanh
#define tgamma std::tgamma
#define trunc std::trunc

#define atan2 std::atan2
#define copysign std::copysign
#define hypot std::hypot
#define ldexp std::ldexp
#define nextafter std::nextafter
#define pow std::pow
#define remainder std::remainder


#define cabsq std::abs
#define crealq std::real
#define cimagq std::imag
#define cargq std::arg
#define conjq std::conj
#define cacosq std::acos
#define cacoshq std::acosh
#define casinq std::asin
#define casinhq std::asinh
#define catanq std::atan
#define catanhq std::atanh
#define ccosq std::cos
#define ccoshq std::cosh
#define cexpq std::exp
#define clogq std::log
#define clog10q std::log10
#define csinq std::sin
#define csinhq std::sinh
#define csqrtq std::sqrt
#define ctanq std::tan
#define ctanhq std::tanh
#define cpowq std::powl
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
