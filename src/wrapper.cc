//
// QCDLoop 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch
//          Keith Ellis: keith.ellis@durham.ac.uk
//          Giulia Zanderighi: giulia.zanderighi@cern.ch

#include "qcdloop/wrapper.h"
#include "qcdloop/qcdloop.h"
#include <stdexcept>
#include <iostream>

#ifdef __cplusplus
extern"C" {
#endif

// result container
vector<complex> r(3);
vector<qcomplex> rq(3);

// topologies
TadPole<complex,double,double> td;
TadPole<complex,complex,double> tdc;
TadPole<qcomplex,qdouble,qdouble> tdq;
TadPole<qcomplex,qcomplex,qdouble> tdcq;

Bubble<complex,double,double> bb;
Bubble<complex,complex,double> bbc;
Bubble<qcomplex,qdouble,qdouble> bbq;
Bubble<qcomplex,qcomplex,qdouble> bbcq;

Triangle<complex,double,double> tr;
Triangle<complex,complex,double> trc;
Triangle<qcomplex,qdouble,qdouble> trq;
Triangle<qcomplex,qcomplex,qdouble> trcq;

Box<complex,double,double> bo;
Box<complex,complex,double> boc;
Box<qcomplex,qdouble,qdouble> boq;
Box<qcomplex,qcomplex,qdouble> bocq;

// global vectors for improving parsing speed
vector<double>   mI1(1);
vector<complex>  mI1c(1);
vector<qdouble>  mI1q(1);
vector<qcomplex> mI1cq(1);
vector<double>   mI2(2);
vector<complex>  mI2c(2);
vector<qdouble>  mI2q(2);
vector<qcomplex> mI2cq(2);
vector<double>   mI3(3);
vector<complex>  mI3c(3);
vector<qdouble>  mI3q(3);
vector<qcomplex> mI3cq(3);
vector<double>   mI4(4);
vector<complex>  mI4c(4);
vector<qdouble>  mI4q(4);
vector<qcomplex> mI4cq(4);
vector<double>   pI2(1);
vector<qdouble>  pI2q(1);
vector<double>   pI3(3);
vector<qdouble>  pI3q(3);
vector<double>   pI4(6);
vector<qdouble>  pI4q(6);

void qlcachesize_(const int &size)
{
  td.setCacheSize(size);
  tdc.setCacheSize(size);
  tdq.setCacheSize(size);
  tdcq.setCacheSize(size);

  bb.setCacheSize(size);
  bbc.setCacheSize(size);
  bbq.setCacheSize(size);
  bbcq.setCacheSize(size);

  tr.setCacheSize(size);
  trc.setCacheSize(size);
  trq.setCacheSize(size);
  trcq.setCacheSize(size);

  bo.setCacheSize(size);
  boc.setCacheSize(size);
  boq.setCacheSize(size);
  bocq.setCacheSize(size);
}

bool qlzero_(double const& x)
{
  return td.iszero(x);
}

bool qlzeroq_(qdouble const& x)
{
  return tdq.iszero(x);
}

bool qlnonzero_(double const& x)
{
  return !td.iszero(x);
}

bool qlnonzeroq_(qdouble const& x)
{
  return !tdq.iszero(x);
}

complex cln_(complex const& x, double const& isig)
{
  return td.cLn(x, isig);
}

qcomplex clnq_(qcomplex const& x, qdouble const& isig)
{
  return tdq.cLn(x, isig);
}

void qltadpole_(complex (&out)[3], double const& mu2, double const& m1)
{
  try
  {
    mI1[0] = m1;
    td.integral(r, mu2, mI1);
    out[0] = r[0];
    out[1] = r[1];
    out[2] = r[2];
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    exit(-1);
  }
}

void qltadpolec_(complex (&out)[3], double const& mu2, complex const& m1)
{
  try
  {
    mI1c[0] = m1;
    tdc.integral(r, mu2, mI1c);
    out[0] = r[0];
    out[1] = r[1];
    out[2] = r[2];
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    exit(-1);
  }
}

void qltadpoleq_(qcomplex (&out)[3], qdouble const& mu2, qdouble const& m1)
{
  try
  {
    mI1q[0] = m1;
    tdq.integral(rq, mu2, mI1q);
    out[0] = rq[0];
    out[1] = rq[1];
    out[2] = rq[2];
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    exit(-1);
  }
}

void qltadpolecq_(qcomplex (&out)[3], qdouble const& mu2, qcomplex const& m1)
{
  try
  {
    mI1cq[0] = m1;
    tdcq.integral(rq, mu2, mI1cq);
    out[0] = rq[0];
    out[1] = rq[1];
    out[2] = rq[2];
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    exit(-1);
  }
}

void qlbubble_(complex (&out)[3], double const& mu2, double const& m1, double const& m2, double const& p1)
{
  try
  {
    mI2[0] = m1;
    mI2[1] = m2;
    pI2[0] = p1;
    bb.integral(r, mu2, mI2, pI2);
    out[0] = r[0];
    out[1] = r[1];
    out[2] = r[2];
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    exit(-1);
  }
}

void qlbubblec_(complex (&out)[3], double const& mu2, complex const& m1, complex const& m2, double const& p1)
{
  try
  {
    mI2c[0] = m1;
    mI2c[1] = m2;
    pI2[0]  = p1;
    bbc.integral(r, mu2, mI2c, pI2);
    out[0] = r[0];
    out[1] = r[1];
    out[2] = r[2];
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    exit(-1);
  }
}

void qlbubbleq_(qcomplex (&out)[3], qdouble const& mu2, qdouble const& m1, qdouble const& m2, qdouble const& p1)
{
  try
  {
    mI2q[0] = m1;
    mI2q[1] = m2;
    pI2q[0]  = p1;
    bbq.integral(rq, mu2, mI2q, pI2q);
    out[0] = rq[0];
    out[1] = rq[1];
    out[2] = rq[2];
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    exit(-1);
  }
}

void qlbubblecq_(qcomplex (&out)[3], qdouble const& mu2, qcomplex const& m1, qcomplex const& m2, qdouble const& p1)
{
  try
  {
    mI2cq[0] = m1;
    mI2cq[1] = m2;
    pI2q[0]  = p1;
    bbcq.integral(rq, mu2, mI2cq, pI2q);
    out[0] = rq[0];
    out[1] = rq[1];
    out[2] = rq[2];
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    exit(-1);
  }
}

void qltriangle_(complex (&out)[3], double const& mu2, double const& m1, double const& m2, double const& m3, double const& p1, double const& p2, double const& p3)
{
  try
  {
    mI3[0] = m1;
    mI3[1] = m2;
    mI3[2] = m3;
    pI3[0] = p1;
    pI3[1] = p2;
    pI3[2] = p3;
    tr.integral(r, mu2, mI3, pI3);
    out[0] = r[0];
    out[1] = r[1];
    out[2] = r[2];
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    exit(-1);
  }
}

void qltrianglec_(complex (&out)[3], double const& mu2, complex const& m1, complex const& m2, complex const& m3, double const& p1, double const& p2, double const& p3)
{
  try
  {
    mI3c[0] = m1;
    mI3c[1] = m2;
    mI3c[2] = m3;
    pI3[0] = p1;
    pI3[1] = p2;
    pI3[2] = p3;
    trc.integral(r, mu2, mI3c, pI3);
    out[0] = r[0];
    out[1] = r[1];
    out[2] = r[2];
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    exit(-1);
  }
}

void qltriangleq_(qcomplex (&out)[3], qdouble const& mu2, qdouble const& m1, qdouble const& m2, qdouble const& m3, qdouble const& p1, qdouble const& p2, qdouble const& p3)
{
  try
  {
    mI3q[0] = m1;
    mI3q[1] = m2;
    mI3q[2] = m3;
    pI3q[0] = p1;
    pI3q[1] = p2;
    pI3q[2] = p3;
    trq.integral(rq, mu2, mI3q, pI3q);
    out[0] = rq[0];
    out[1] = rq[1];
    out[2] = rq[2];
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    exit(-1);
  }
}
void qltrianglecq_(qcomplex (&out)[3], qdouble const& mu2, qcomplex const& m1, qcomplex const& m2, qcomplex const& m3, qdouble const& p1, qdouble const& p2, qdouble const& p3)
{
  try
  {
    mI3cq[0] = m1;
    mI3cq[1] = m2;
    mI3cq[2] = m3;
    pI3q[0] = p1;
    pI3q[1] = p2;
    pI3q[2] = p3;
    trcq.integral(rq, mu2, mI3cq, pI3q);
    out[0] = rq[0];
    out[1] = rq[1];
    out[2] = rq[2];
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    exit(-1);
  }
}

void qlbox_(complex (&out)[3], double const& mu2, double const& m1, double const& m2, double const& m3, double const& m4, double const& p1, double const& p2, double const& p3, double const& p4, double const& s12, double const& s23)
{
  try
  {
    mI4[0] = m1;
    mI4[1] = m2;
    mI4[2] = m3;
    mI4[3] = m4;
    pI4[0] = p1;
    pI4[1] = p2;
    pI4[2] = p3;
    pI4[3] = p4;
    pI4[4] = s12;
    pI4[5] = s23;
    bo.integral(r, mu2, mI4, pI4);
    out[0] = r[0];
    out[1] = r[1];
    out[2] = r[2];
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    exit(-1);
  }
}

void qlboxc_(complex (&out)[3], double const& mu2, complex const& m1, complex const& m2, complex const& m3, complex const& m4, double const& p1, double const& p2, double const& p3, double const& p4, double const& s12, double const& s23)
{
  try
  {
    mI4c[0] = m1;
    mI4c[1] = m2;
    mI4c[2] = m3;
    mI4c[3] = m4;
    pI4[0] = p1;
    pI4[1] = p2;
    pI4[2] = p3;
    pI4[3] = p4;
    pI4[4] = s12;
    pI4[5] = s23;
    boc.integral(r, mu2, mI4c, pI4);
    out[0] = r[0];
    out[1] = r[1];
    out[2] = r[2];
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    exit(-1);
  }
}

void qlboxq_(qcomplex (&out)[3], qdouble const& mu2, qdouble const& m1, qdouble const& m2, qdouble const& m3, qdouble const& m4, qdouble const& p1, qdouble const& p2, qdouble const& p3, qdouble const& p4, qdouble const& s12, qdouble const& s23)
{
  try
  {
    mI4q[0] = m1;
    mI4q[1] = m2;
    mI4q[2] = m3;
    mI4q[3] = m4;
    pI4q[0] = p1;
    pI4q[1] = p2;
    pI4q[2] = p3;
    pI4q[3] = p4;
    pI4q[4] = s12;
    pI4q[5] = s23;
    boq.integral(rq, mu2, mI4q, pI4q);
    out[0] = rq[0];
    out[1] = rq[1];
    out[2] = rq[2];
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    exit(-1);
  }
}

void qlboxcq_(qcomplex (&out)[3], qdouble const& mu2, qcomplex const& m1, qcomplex const& m2, qcomplex const& m3, qcomplex const& m4, qdouble const& p1, qdouble const& p2, qdouble const& p3, qdouble const& p4, qdouble const& s12, qdouble const& s23)
{
  try
  {
    mI4cq[0] = m1;
    mI4cq[1] = m2;
    mI4cq[2] = m3;
    mI4cq[3] = m4;
    pI4q[0] = p1;
    pI4q[1] = p2;
    pI4q[2] = p3;
    pI4q[3] = p4;
    pI4q[4] = s12;
    pI4q[5] = s23;
    bocq.integral(rq, mu2, mI4cq, pI4q);
    out[0] = rq[0];
    out[1] = rq[1];
    out[2] = rq[2];
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    exit(-1);
  }
}

#ifdef QL_NAMES

  // for backward compatibility.
  void qlinit_()
  {
    std::cout << ql::yellow << "[QCDLoop warning]: this wrapper is not thread-safe." << std::endl;
    std::cout << "[QCDLoop suggestion]: consider developing object-oriented code." << ql::def << std::endl;
  }

  complex qli1_(double const& m1, double const& mu2, int const& ep)
  {
    try
    {
      mI1[0] = m1;
      td.integral(r, mu2, mI1);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    return r[abs(ep)];
  }

  complex qli1c_(complex const& m1, double const& mu2, int const& ep)
  {
    try
    {
      mI1c[0] = m1;
      tdc.integral(r, mu2, mI1c);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    return r[abs(ep)];
  }

  qcomplex qli1q_(qdouble const& m1, qdouble const& mu2, int const& ep)
  {
    try
    {
      mI1q[0] = m1;
      tdq.integral(rq, mu2, mI1q);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    return rq[abs(ep)];
  }

  qcomplex qli1qc_(qcomplex const& m1, qdouble const& mu2, int const& ep)
  {
    try
    {
      mI1cq[0] = m1;
      tdcq.integral(rq, mu2, mI1cq);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    return rq[abs(ep)];
  }

  complex qli2_(double const& p1, double const& m1, double const& m2, double const& mu2, int const& ep)
  {
    try
    {
      mI2[0] = m1;
      mI2[1] = m2;
      pI2[0] = p1;
      bb.integral(r, mu2, mI2, pI2);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }    
    return r[abs(ep)];
  }

  complex qli2c_(double const& p1, complex const& m1, complex const& m2, double const& mu2, int const& ep)
  {
    try
    {
      mI2c[0] = m1;
      mI2c[1] = m2;
      pI2[0] = p1;
      bbc.integral(r, mu2, mI2c, pI2);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    return r[abs(ep)];
  }

  qcomplex qli2q_(qdouble const& p1, qdouble const& m1, qdouble const& m2, qdouble const& mu2, int const& ep)
  {
    try
    {
      mI2q[0] = m1;
      mI2q[1] = m2;
      pI2q[0] = p1;
      bbq.integral(rq, mu2, mI2q, pI2q);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    return rq[abs(ep)];
  }

  qcomplex qli2qc_(qdouble const& p1, qcomplex const& m1, qcomplex const& m2, qdouble const& mu2, int const& ep)
  {
    try
    {
      mI2cq[0] = m1;
      mI2cq[1] = m2;
      pI2q[0] = p1;
      bbcq.integral(rq, mu2, mI2cq, pI2q);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    return rq[abs(ep)];
  }

  complex qli3_(double const& p1, double const& p2, double const& p3, double const& m1, double const& m2, double const& m3, double const& mu2, int const& ep)
  {
    try
    {
      mI3[0] = m1;
      mI3[1] = m2;
      mI3[2] = m3;
      pI3[0] = p1;
      pI3[1] = p2;
      pI3[2] = p3;
      tr.integral(r, mu2, mI3, pI3);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    return r[abs(ep)];
  }

  complex qli3c_(double const& p1, double const& p2, double const& p3, complex const& m1, complex const& m2, complex const& m3, double const& mu2, int const& ep)
  {
    try
    {
      mI3c[0] = m1;
      mI3c[1] = m2;
      mI3c[2] = m3;
      pI3[0] = p1;
      pI3[1] = p2;
      pI3[2] = p3;
      trc.integral(r, mu2, mI3c, pI3);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    return r[abs(ep)];
  }

  qcomplex qli3q_(qdouble const& p1, qdouble const& p2, qdouble const& p3, qdouble const& m1, qdouble const& m2, qdouble const& m3, qdouble const& mu2, int const& ep)
  {
    try
    {
      mI3q[0] = m1;
      mI3q[1] = m2;
      mI3q[2] = m3;
      pI3q[0] = p1;
      pI3q[1] = p2;
      pI3q[2] = p3;
      trq.integral(rq, mu2, mI3q, pI3q);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    return rq[abs(ep)];
  }

  qcomplex qli3qc_(qdouble const& p1, qdouble const& p2, qdouble const& p3, qcomplex const& m1, qcomplex const& m2, qcomplex const& m3, qdouble const& mu2, int const& ep)
  {
    try
    {
      mI3cq[0] = m1;
      mI3cq[1] = m2;
      mI3cq[2] = m3;
      pI3q[0] = p1;
      pI3q[1] = p2;
      pI3q[2] = p3;
      trcq.integral(rq, mu2, mI3cq, pI3q);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    return rq[abs(ep)];
  }

  complex qli4_(double const& p1, double const& p2, double const& p3, double const& p4, double const& s12, double const& s23, double const& m1, double const& m2, double const& m3, double const& m4, double const& mu2, int const& ep)
  {
    try
    {
      mI4[0] = m1;
      mI4[1] = m2;
      mI4[2] = m3;
      mI4[3] = m4;
      pI4[0] = p1;
      pI4[1] = p2;
      pI4[2] = p3;
      pI4[3] = p4;
      pI4[4] = s12;
      pI4[5] = s23;
      bo.integral(r, mu2, mI4, pI4);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }

    return r[abs(ep)];
  }

  complex qli4c_(double const& p1, double const& p2, double const& p3, double const& p4, double const& s12, double const& s23, complex const& m1, complex const& m2, complex const& m3, complex const& m4, double const& mu2, int const& ep)
  {
    try
    {
      mI4c[0] = m1;
      mI4c[1] = m2;
      mI4c[2] = m3;
      mI4c[3] = m4;
      pI4[0] = p1;
      pI4[1] = p2;
      pI4[2] = p3;
      pI4[3] = p4;
      pI4[4] = s12;
      pI4[5] = s23;
      boc.integral(r, mu2, mI4c, pI4);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    return r[abs(ep)];
  }

  qcomplex qli4q_(qdouble const& p1, qdouble const& p2, qdouble const& p3, qdouble const& p4, qdouble const& s12, qdouble const& s23, qdouble const& m1, qdouble const& m2, qdouble const& m3, qdouble const& m4, qdouble const& mu2, int const& ep)
  {
    try
    {
      mI4q[0] = m1;
      mI4q[1] = m2;
      mI4q[2] = m3;
      mI4q[3] = m4;
      pI4q[0] = p1;
      pI4q[1] = p2;
      pI4q[2] = p3;
      pI4q[3] = p4;
      pI4q[4] = s12;
      pI4q[5] = s23;
      boq.integral(rq, mu2, mI4q, pI4q);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    return rq[abs(ep)];
  }

  qcomplex qli4qc_(qdouble const& p1, qdouble const& p2, qdouble const& p3, qdouble const& p4, qdouble const& s12, qdouble const& s23, qcomplex const& m1, qcomplex const& m2, qcomplex const& m3, qcomplex const& m4, qdouble const& mu2, int const& ep)
  {
    try
    {
      mI4cq[0] = m1;
      mI4cq[1] = m2;
      mI4cq[2] = m3;
      mI4cq[3] = m4;
      pI4q[0] = p1;
      pI4q[1] = p2;
      pI4q[2] = p3;
      pI4q[3] = p4;
      pI4q[4] = s12;
      pI4q[5] = s23;
      bocq.integral(rq, mu2, mI4cq, pI4q);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    return rq[abs(ep)];
  }

#endif

#ifdef __cplusplus
}
#endif
