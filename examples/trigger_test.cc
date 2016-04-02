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
  double mu2 = ql::Pow(1.7,2.0);
  vector<double> p = {};
  vector<double> m = {5.0};
  vector<complex> res(3);
  ql::Timer tt;
  cout << scientific << setprecision(32);

  /**
   * DOUBLE PRECISION
   */

  /**
   * Tadpole
   */
  ql::TadPole<complex,double,double> tp;
  cout << ql::red << "\n====== Direct Tadpole Integral ======" << ql::def << endl;
  tt.start();
  for (int i = 0; i < 1e7; i++) tp.integral(res, mu2, m, p);
  tt.printTime(tt.stop());

  for (size_t i = 0; i < res.size(); i++)
    cout << "eps" << i << "\t" << res[i] << endl;

  // Second try with trigger
  ql::QCDLoop<complex,double,double> ql;
  cout << ql::blue << "\n====== Auto-Trigger Tadpole ======" << ql::def << endl;
  tt.start();
  for (int i = 0; i < 1e7; i++) ql.integral(res, mu2, m, p);
  tt.printTime(tt.stop());

  for (size_t i = 0; i < res.size(); i++)
    cout << "eps" << i << "\t" << res[i] << endl;    

  /**
   * Bubble
   */
  m = {5.0, 2.0};
  p = {1.0};
  ql::Bubble<complex,double> bb;
  cout << ql::red << "\n====== Direct Bubble Integral ======" << ql::def << endl;
  tt.start();
  for (int i = 0; i < 1e7; i++) bb.integral(res, mu2, m, p);
  tt.printTime(tt.stop());

  for (size_t i = 0; i < res.size(); i++)
    cout << "eps" << i << "\t" << res[i] << endl;

  // Second try with trigger
  cout << ql::blue << "\n====== Auto-Trigger Bubble ======" << ql::def << endl;
  tt.start();
  for (int i = 0; i < 1e7; i++) ql.integral(res, mu2, m, p);
  tt.printTime(tt.stop());

  for (size_t i = 0; i < res.size(); i++)
    cout << "eps" << i << "\t" << res[i] << endl;

  /**
   * Triangle
   */
  m = {5.0, 2.0, 3.0};
  p = {1.0, 2.0, 4.0};
  ql::Triangle<complex,double> tr;
  cout << ql::red << "\n====== Direct Triangle Integral ======" << ql::def << endl;
  tt.start();
  for (int i = 0; i < 1e7; i++) tr.integral(res, mu2, m, p);
  tt.printTime(tt.stop());

  for (size_t i = 0; i < res.size(); i++)
    cout << "eps" << i << "\t" << res[i] << endl;

  // Second try with trigger
  cout << ql::blue << "\n====== Auto-Trigger Triangle ======" << ql::def << endl;
  tt.start();
  for (int i = 0; i < 1e7; i++) ql.integral(res, mu2, m, p);
  tt.printTime(tt.stop());

  for (size_t i = 0; i < res.size(); i++)
    cout << "eps" << i << "\t" << res[i] << endl;

  /**
   * Box
   */
  m = {0.0, 0.0, 0.0, 0.0};
  p = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};

  ql::Box<complex,double> tb;
  cout << ql::red << "\n====== Direct Box Integral ======" << ql::def << endl;
  tt.start();
  for (int i = 0; i < 1e7; i++) tb.integral(res, mu2, m, p);
  tt.printTime(tt.stop());

  for (size_t i = 0; i < res.size(); i++)
    cout << "eps" << i << "\t" << res[i] << endl;

  // Second try with trigger
  cout << ql::blue << "\n====== Auto-Trigger Box ======" << ql::def << endl;
  tt.start();
  for (int i = 0; i < 1e7; i++) ql.integral(res, mu2, m, p);
  tt.printTime(tt.stop());

  for (size_t i = 0; i < res.size(); i++)
    cout << "eps" << i << "\t" << res[i] << endl;

  /**
    * QUADRUPLE PRECISION
    */
  cout << ql::green << "\n====== Quadruple precision ======" << ql::def << endl;
  qdouble mu2q = ql::Pow(1.7q,2.0q);
  vector<qdouble> pq = {};
  vector<qdouble> mq = {5.0q};
  vector<qcomplex> resq(3);

  ql::TadPole<qcomplex,qdouble,qdouble> tpq;
  cout << ql::red << "\n====== Direct Tadpole Integral ======" << ql::def << endl;
  tt.start();
  for (int i = 0; i < 1e7; i++) tpq.integral(resq, mu2q, mq, pq);
  tt.printTime(tt.stop());

  for (size_t i = 0; i < resq.size(); i++)
    cout << "eps" << i << "\t(" << crealq(resq[i]) << "," << cimagq(resq[i]) << endl;

  ql::QCDLoop<qcomplex,qdouble,qdouble> qlq;
  cout << ql::blue << "\n====== Auto-Trigger Tadpole ======" << ql::def << endl;
  tt.start();
  for (int i = 0; i < 1e7; i++) qlq.integral(resq, mu2q, mq, pq);
  tt.printTime(tt.stop());

  for (size_t i = 0; i < resq.size(); i++)
    cout << "eps" << i << "\t" << resq[i] << endl;

  /**
   * Bubble
   */
  mq = {5.0q, 2.0q};
  pq = {1.0q};
  ql::Bubble<qcomplex,qdouble,qdouble> bbq;
  cout << ql::red << "\n====== Direct Bubble Integral ======" << ql::def << endl;
  tt.start();
  for (int i = 0; i < 1e7; i++) bbq.integral(resq, mu2q, mq, pq);
  tt.printTime(tt.stop());

  for (size_t i = 0; i < res.size(); i++)
    cout << "eps" << i << "\t" << resq[i] << endl;

  // Second try with trigger
  cout << ql::blue << "\n====== Auto-Trigger Bubble ======" << ql::def << endl;
  tt.start();
  for (int i = 0; i < 1e7; i++) qlq.integral(resq, mu2q, mq, pq);
  tt.printTime(tt.stop());

  for (size_t i = 0; i < res.size(); i++)
    cout << "eps" << i << "\t" << resq[i] << endl;

  /**
   * Triangle
   */
  mq = {5.0q, 2.0q, 3.0q};
  pq = {1.0q, 2.0q, 4.0q};
  ql::Triangle<qcomplex,qdouble,qdouble> trq;
  cout << ql::red << "\n====== Direct Triangle Integral ======" << ql::def << endl;
  tt.start();
  for (int i = 0; i < 1e7; i++) trq.integral(resq, mu2q, mq, pq);
  tt.printTime(tt.stop());

  for (size_t i = 0; i < res.size(); i++)
    cout << "eps" << i << "\t" << resq[i] << endl;

  // Second try with trigger
  cout << ql::blue << "\n====== Auto-Trigger Triangle ======" << ql::def << endl;
  tt.start();
  for (int i = 0; i < 1e7; i++) qlq.integral(resq, mu2q, mq, pq);
  tt.printTime(tt.stop());

  for (size_t i = 0; i < res.size(); i++)
    cout << "eps" << i << "\t" << resq[i] << endl;

  /**
   * Box
   */
  mq = {0.0q, 0.0q, 0.0q, 0.0q};
  pq = {1.0q, 2.0q, 3.0q, 4.0q, 5.0q, 6.0q};

  ql::Box<qcomplex,qdouble,qdouble> tbq;
  cout << ql::red << "\n====== Direct Box Integral ======" << ql::def << endl;
  tt.start();
  for (int i = 0; i < 1e7; i++) tbq.integral(resq, mu2q, mq, pq);
  tt.printTime(tt.stop());

  for (size_t i = 0; i < resq.size(); i++)
    cout << "eps" << i << "\t" << resq[i] << endl;

  // Second try with trigger
  cout << ql::blue << "\n====== Auto-Trigger Box ======" << ql::def << endl;
  tt.start();
  for (int i = 0; i < 1e7; i++) qlq.integral(resq, mu2q, mq, pq);
  tt.printTime(tt.stop());

  for (size_t i = 0; i < resq.size(); i++)
    cout << "eps" << i << "\t" << resq[i] << endl;


  return 0;
}
