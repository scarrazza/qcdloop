//
// QCDLoop 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch
//          Keith Ellis: keith.ellis@durham.ac.uk
//          Giulia Zanderighi: giulia.zanderighi@cern.ch

#pragma once

#ifndef NOQUADMATH
extern "C" { // for gcc4.7 compatibility
#include <quadmath.h>
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
#ifdef NOQUADMATH
  typedef long double qdouble;
  typedef std::complex<long double> qcomplex;
#else
  typedef __float128 qdouble;
  typedef __complex128 qcomplex;
#endif
  typedef std::complex<double> complex;  
  enum Code { red = 31, green = 32, yellow = 33, blue = 34, def = 39};
}

namespace std
{
#ifndef NOQUADMATH
  //! implementation of operator<< for qdouble
  ostream& operator<<(std::ostream& out, ql::qdouble f);

  //! implementation of operator<< for qcomplex
  ostream& operator<<(std::ostream& out, ql::qcomplex f);
#endif

  //! implementation of operator<< for Code
  std::ostream& operator<<(std::ostream& os, ql::Code code);
}
