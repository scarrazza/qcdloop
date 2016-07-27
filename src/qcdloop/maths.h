//
// QCDLoop 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch
//          Keith Ellis: keith.ellis@durham.ac.uk
//          Giulia Zanderighi: giulia.zanderighi@cern.ch

#pragma once

#include "types.h"

/*!
 * Some basic math functions with inline (performance)
 * for double and quadruple precision.
 */
namespace ql
{
  // Logarithms
  inline double   Log(double const& x)   { return std::log(x); }
  inline qdouble  Log(qdouble const& x)  { return logq(x); }
  inline complex  Log(complex const& x)  { return std::log(x); }
  inline qcomplex Log(qcomplex const& x) { return clogq(x); }

  // Power
  inline double   Pow(double const& x, int const& a)   { return std::pow(x, a); }
  inline qdouble  Pow(qdouble const& x, int const& a)  { return powq(x,a); }
  inline complex  Pow(complex const& x, int const& a) { return std::pow(x,a); }
  inline qcomplex Pow(qcomplex const& x, int const& a){ return cpowq(x,a); }

  // Root
  inline double   Sqrt(double const& x)  { return std::sqrt(x); }
  inline qdouble  Sqrt(qdouble const& x) { return sqrtq(x); }
  inline complex  Sqrt(complex const& x) { return std::sqrt(x); }
  inline qcomplex Sqrt(qcomplex const& x){ return csqrtq(x); }

  // Absolute value
  inline double  Abs(double const& x)   { return std::fabs(x); }
  inline qdouble Abs(qdouble const& x)  { return fabsq(x);}
  inline double  Abs(complex const& x)  { return std::fabs(x);}
  inline qdouble Abs(qcomplex const& x) { return cabsq(x); }

  // Complex tools, imag, real and conj.
  inline double  Imag(double const& x)  { UNUSED(x); return 0; }
  inline qdouble Imag(qdouble const& x) { UNUSED(x); return qdouble(0); }
  inline double  Imag(complex const& x) { return x.imag(); }
  inline qdouble Imag(qcomplex const& x){ return cimagq(x);}

  inline double  Real(double const& x) { return x; }
  inline qdouble Real(qdouble const& x) { return x; }
  inline double  Real(complex const& x) { return x.real(); }
  inline qdouble Real(qcomplex const& x) { return crealq(x); }

  inline complex  Conjg(complex const& x) { return std::conj(x); }
  inline qcomplex Conjg(qcomplex const& x){ return conjq(x); }

  // Comparison and sign operations
  inline int Sign(double const& x) { return (double(0) < x) - (x < double(0)); }
  inline int Sign(qdouble const& x){ return (qdouble(0) < x) -(x < qdouble(0));}

  inline double Max(double const& a, double const& b)
  {
    if (Abs(a) > Abs(b)) return a;
    else return b;
  }
  inline qdouble Max(qdouble const& a, qdouble const& b)
  {
    if (Abs(a) > Abs(b)) return a;
    else return b;
  }
  inline complex  Max(complex const& a, complex const& b)
  {
    if (Abs(a) > Abs(b)) return a;
    else return b;
  }
  inline qcomplex Max(qcomplex const& a, qcomplex const& b)
  {
    if (Abs(a) > Abs(b)) return a;
    else return b;
  }

  inline double Min(double const& a, double const& b)
  {
    if (Abs(a) > Abs(b)) return b;
    else return a;
  }
  inline qdouble Min(qdouble const& a, qdouble const& b)
  {
    if (Abs(a) > Abs(b)) return b;
    else return a;
  }
  inline complex Min(complex const& a, complex const& b)
  {
    if (Abs(a) > Abs(b)) return b;
    else return a;
  }
  inline qcomplex Min(qcomplex const& a, qcomplex const& b)
  {
    if (Abs(a) > Abs(b)) return b;
    else return a;
  }

  // Heaviside function
  inline double Htheta(double const& x) { return 0.5*(1+Sign(x)); }
  inline qdouble Htheta(qdouble const& x){ return qdouble(0.5)*(qdouble(1)+Sign(x));}
}
