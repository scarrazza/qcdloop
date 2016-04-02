//
// QCDLoop 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch
//          Keith Ellis: keith.ellis@durham.ac.uk
//          Giulia Zanderighi: giulia.zanderighi@cern.ch

#include <iostream>
#include <iomanip>
#include <qcdloop/qcdloop.h>
#include <qcdloop/cache.h>

using std::vector;
using std::cout;
using std::endl;
using std::setprecision;
using std::scientific;
using ql::complex;
using ql::qcomplex;
using ql::qdouble;

int main()
{
  const double mu2 = ql::Pow(1.7,2.0);
  vector<double> p   = {};
  vector<double>   m = {5.0};
  vector<complex> cm = {{5.0,0.0}};
  vector<complex> res(3);

  ql::Timer tt;
  cout << scientific << setprecision(32);

  /**
   * DOUBLE PRECISION
   */

  /**
   * Tadpole
   */
  ql::TadPole<complex,double> tp;
  cout << ql::red << "\n====== Direct Tadpole Integral - double masses ======" << ql::def << endl;
  tt.start();
  for (int i = 0; i < 1e7; i++) tp.integral(res, mu2, m, p);
  tt.printTime(tt.stop());

  for (size_t i = 0; i < res.size(); i++)
    cout << "eps" << i << "\t" << res[i] << endl;

  ql::TadPole<complex,complex> ctp;
  cout << ql::red << "\n====== Direct Tadpole Integral - complex masses ======" << ql::def << endl;
  tt.start();
  for (int i = 0; i < 1e7; i++) ctp.integral(res, mu2, cm, p);
  tt.printTime(tt.stop());

  for (size_t i = 0; i < res.size(); i++)
    cout << "eps" << i << "\t" << res[i] << endl;

  /**
   * Bubble
   */
  m = {5.0, 2.0};
  cm= {{5.0,0.0}, {2.0,0.0}};
  p = {1.0};
  ql::Bubble<complex,double,double> bb;
  cout << ql::red << "\n====== Direct Bubble Integral - double masses ======" << ql::def << endl;
  tt.start();
  for (int i = 0; i < 1e7; i++) bb.integral(res, mu2, m, p);
  tt.printTime(tt.stop());

  for (size_t i = 0; i < res.size(); i++)
    cout << "eps" << i << "\t" << res[i] << endl;

  ql::Bubble<complex,complex,double> cbb;
  cout << ql::red << "\n====== Direct Bubble Integral - complex masses ======" << ql::def << endl;
  tt.start();
  for (int i = 0; i < 1e7; i++) cbb.integral(res, mu2, cm, p);
  tt.printTime(tt.stop());

  for (size_t i = 0; i < res.size(); i++)
    cout << "eps" << i << "\t" << res[i] << endl;

  /**
   * QUADRUPLE PRECISION
   */

  const qdouble mu2q = ql::Pow(1.7q,2.0q);
  vector<qdouble> pq   = {};
  vector<qdouble>   mq = {5.0q};
  vector<qcomplex> cmq = {{5.0q,0.0q}};
  vector<qcomplex> resq(3);

  /**
   * Tadpole
   */
  ql::TadPole<qcomplex,qdouble,qdouble> tpq;
  cout << ql::red << "\n====== Direct Tadpole Integral - qdouble masses ======" << ql::def << endl;
  tt.start();
  for (int i = 0; i < 1e7; i++) tpq.integral(resq, mu2q, mq, pq);
  tt.printTime(tt.stop());

  for (size_t i = 0; i < resq.size(); i++)
    cout << "eps" << i << "\t" << resq[i] << endl;

  ql::TadPole<qcomplex,qcomplex,qdouble> ctpq;
  cout << ql::red << "\n====== Direct Tadpole Integral - qcomplex masses ======" << ql::def << endl;
  tt.start();
  for (int i = 0; i < 1e7; i++) ctpq.integral(resq, mu2q, cmq, pq);
  tt.printTime(tt.stop());

  for (size_t i = 0; i < resq.size(); i++)
    cout << "eps" << i << "\t" << resq[i] << endl;

  /**
   * Bubble
   */
  mq = {5.0q, 2.0q};
  cmq= {{5.0q,0.0q},{2.0q,0.0q}};
  pq = {1.0q};
  ql::Bubble<qcomplex,qdouble,qdouble> bbq;
  cout << ql::red << "\n====== Direct Bubble Integral - qdouble masses ======" << ql::def << endl;
  tt.start();
  for (int i = 0; i < 1e7; i++) bbq.integral(resq, mu2q, mq, pq);
  tt.printTime(tt.stop());

  for (size_t i = 0; i < resq.size(); i++)
    cout << "eps" << i << "\t" << resq[i] << endl;

  ql::Bubble<qcomplex,qcomplex,qdouble> cbbq;
  cout << ql::red << "\n====== Direct Bubble Integral - qcomplex masses ======" << ql::def << endl;
  tt.start();
  for (int i = 0; i < 1e7; i++) cbbq.integral(resq, mu2q, cmq, pq);
  tt.printTime(tt.stop());

  for (size_t i = 0; i < resq.size(); i++)
    cout << "eps" << i << "\t" << resq[i] << endl;

  return 0;
}
