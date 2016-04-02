//
// QCDLoop 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch
//          Keith Ellis: keith.ellis@durham.ac.uk
//          Giulia Zanderighi: giulia.zanderighi@cern.ch

/*!
 * C/Fortran wrappers
 */

#pragma once

#include "qcdloop/types.h"
#include "qcdloop/config.h"
using namespace ql;

#ifdef __cplusplus
extern"C" {
#endif

  void qlcachesize_(int const& size);

  void qltadpole_(complex (&out)[3], double const& mu2, double const& m1);
  void qltadpolec_(complex (&out)[3], double const& mu2, complex const& m1);
  void qltadpoleq_(qcomplex (&out)[3], qdouble const& mu2, qdouble const& m1);
  void qltadpolecq_(qcomplex (&out)[3], qdouble const& mu2, qcomplex const& m1);

  void qlbubble_(complex (&out)[3], double const& mu2, double const& m1, double const& m2, double const& p1);
  void qlbubblec_(complex (&out)[3], double const& mu2, complex const& m1, complex const& m2, double const& p1);
  void qlbubbleq_(qcomplex (&out)[3], qdouble const& mu2, qdouble const& m1, qdouble const& m2, qdouble const& p1);
  void qlbubblecq_(qcomplex (&out)[3], qdouble const& mu2, qcomplex const& m1, qcomplex const& m2, qdouble const& p1);

  void qltriangle_(complex (&out)[3], double const& mu2, double const& m1, double const& m2, double const& m3, double const& p1, double const& p2, double const& p3);
  void qltrianglec_(complex (&out)[3], double const& mu2, complex const& m1, complex const& m2, complex const& m3, double const& p1, double const& p2, double const& p3);
  void qltriangleq_(qcomplex (&out)[3], qdouble const& mu2, qdouble const& m1, qdouble const& m2, qdouble const& m3, qdouble const& p1, qdouble const& p2, qdouble const& p3);
  void qltrianglecq_(qcomplex (&out)[3], qdouble const& mu2, qcomplex const& m1, qcomplex const& m2, qcomplex const& m3, qdouble const& p1, qdouble const& p2, qdouble const& p3);

  void qlbox_(complex (&out)[3], double const& mu2, double const& m1, double const& m2, double const& m3, double const& m4, double const& p1, double const& p2, double const& p3, double const& p4, double const& s12, double const& s23);
  void qlboxc_(complex (&out)[3], double const& mu2, complex const& m1, complex const& m2, complex const& m3, complex const& m4, double const& p1, double const& p2, double const& p3, double const& p4, double const& s12, double const& s23);
  void qlboxq_(qcomplex (&out)[3], qdouble const& mu2, qdouble const& m1, qdouble const& m2, qdouble const& m3, qdouble const& m4, qdouble const& p1, qdouble const& p2, qdouble const& p3, qdouble const& p4, qdouble const& s12, qdouble const& s23);
  void qlboxcq_(qcomplex (&out)[3], qdouble const& mu2, qcomplex const& m1, qcomplex const& m2, qcomplex const& m3, qcomplex const& m4, qdouble const& p1, qdouble const& p2, qdouble const& p3, qdouble const& p4, qdouble const& s12, qdouble const& s23);

#ifdef QL_NAMES

  void qlinit_();

  complex  qli1_(double const& m1, double const& mu2, int const& ep);
  complex  qli1c_(complex const& m1, double const& mu2, int const& ep);
  qcomplex qli1q_(qdouble const& m1, qdouble const& mu2, int const& ep);
  qcomplex qli1qc_(qcomplex const& m1, qdouble const& mu2, int const& ep);

  complex  qli2_(double const& p1, double const& m1, double const& m2, double const& mu2, int const& ep);
  complex  qli2c_(double const& p1, complex const& m1, complex const& m2, double const& mu2, int const& ep);
  qcomplex qli2q_(qdouble const& p1, qdouble const& m1, qdouble const& m2, qdouble const& mu2, int const& ep);
  qcomplex qli2qc_(qdouble const& p1, qcomplex const& m1, qcomplex const& m2, qdouble const& mu2, int const& ep);

  complex  qli3_(double const& p1, double const& p2, double const& p3, double const& m1, double const& m2, double const& m3, double const& mu2, int const& ep);
  complex  qli3c_(double const& p1, double const& p2, double const& p3, complex const& m1, complex const& m2, complex const& m3, double const& mu2, int const& ep);
  qcomplex qli3q_(qdouble const& p1, qdouble const& p2, qdouble const& p3, qdouble const& m1, qdouble const& m2, qdouble const& m3, qdouble const& mu2, int const& ep);
  qcomplex qli3qc_(qdouble const& p1, qdouble const& p2, qdouble const& p3, qcomplex const& m1, qcomplex const& m2, qcomplex const& m3, qdouble const& mu2, int const& ep);

  complex  qli4_(double const& p1, double const& p2, double const& p3, double const& p4, double const& s12, double const& s23, double const& m1, double const& m2, double const& m3, double const& m4, double const& mu2, int const& ep);
  complex  qli4c_(double const& p1, double const& p2, double const& p3, double const& p4, double const& s12, double const& s23, complex const& m1, complex const& m2, complex const& m3, complex const& m4, double const& mu2, int const& ep);
  qcomplex qli4q_(qdouble const& p1, qdouble const& p2, qdouble const& p3, qdouble const& p4, qdouble const& s12, qdouble const& s23, qdouble const& m1, qdouble const& m2, qdouble const& m3, qdouble const& m4, qdouble const& mu2, int const& ep);
  qcomplex qli4qc_(qdouble const& p1, qdouble const& p2, qdouble const& p3, qdouble const& p4, qdouble const& s12, qdouble const& s23, qcomplex const& m1, qcomplex const& m2, qcomplex const& m3, qcomplex const& m4, qdouble const& mu2, int const& ep);

  // some extra functions for MCMF
  bool qlzero_(double const& x);
  bool qlzeroq_(qdouble const& x);
  bool qlnonzero_(double const& x);
  bool qlnonzeroq_(qdouble const& x);
  complex cln_(complex const& x, double const& isig);
  qcomplex clnq_(qcomplex const& x, qdouble const& isig);

#endif

#ifdef __cplusplus
}
#endif

