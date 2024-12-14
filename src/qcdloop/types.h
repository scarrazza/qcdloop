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
#endif
#ifdef __aarch64__
#include <stdlib.h>
#include <complex.h>
#include <tgmath.h>
#include <float.h>
typedef long double __float128;
using __complex128 = long double _Complex  ;
extern "C" {
__complex128 conjl(__complex128);
__float128 cimagl(__complex128);
__float128 creall(__complex128);
__complex128 csqrtl(__complex128);
__complex128 clogl(__complex128);
__float128 cabsl(__complex128);
__complex128 cpowl(__complex128,__complex128);
}
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

#ifdef __x86_64__
  //! implementation of operator<< for qdouble
  ostream& operator<<(std::ostream& out, ql::qdouble f);
#endif

  //! implementation of operator<< for qcomplex
  ostream& operator<<(std::ostream& out, ql::qcomplex f);

  //! implementation of operator<< for Code
  std::ostream& operator<<(std::ostream& os, ql::Code code);
}
