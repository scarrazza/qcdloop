
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
  const vector<double> p = {};
  const vector<double> m = {5.0};
  vector<complex> res(3);
  ql::Timer tt;
  cout << scientific << setprecision(32);

  /**
   * SAME INPUT
   */
  ql::TadPole<complex,double> tp;

  cout << ql::red << "\n====== Same input - Cache N=0 ======" << ql::def << endl;
  tp.setCacheSize(0);
  tt.start();
  for (int i = 0; i < 1e7; i++) tp.integral(res, mu2, m, p);
  tt.printTime(tt.stop());

  for (size_t i = 0; i < res.size(); i++)
    cout << "eps" << i << "\t" << res[i] << endl;

  cout << ql::red << "\n====== Same input - Cache N=1 ======" << ql::def << endl;
  tp.setCacheSize(1);
  tt.start();
  for (int i = 0; i < 1e7; i++) tp.integral(res, mu2, m, p);
  tt.printTime(tt.stop());

  for (size_t i = 0; i < res.size(); i++)
    cout << "eps" << i << "\t" << res[i] << endl;

  cout << ql::red << "\n====== Same input - Cache N=2 ======" << ql::def << endl;
  tp.setCacheSize(2);
  tt.start();
  for (int i = 0; i < 1e7; i++) tp.integral(res, mu2, m, p);
  tt.printTime(tt.stop());

  for (size_t i = 0; i < res.size(); i++)
    cout << "eps" << i << "\t" << res[i] << endl;

  /**
    * MULTIPLE INPUT
    */
  vector<double> mu2s = { ql::Pow(1.7,2.0), ql::Pow(1.5,2.0), ql::Pow(2.0,2.0), ql::Pow(1.2,2.0) };
  vector< vector<double> > ms = { {5.0}, {4.0}, {6.0}, {7.0}};
  int index = 0;

  cout << ql::blue << "\n====== 4 different input - Cache N=0 ======" << ql::def << endl;
  tp.setCacheSize(0); index = 0;
  tt.start();
  for (int i = 0; i < 1e7; i++)
    {
      tp.integral(res, mu2s[index], ms[index], p);
      index += 1; if (index == 4) index = 0;
    }
  tt.printTime(tt.stop());

  for (size_t i = 0; i < res.size(); i++)
    cout << "eps" << i << "\t" << res[i] << endl;

  cout << ql::blue << "\n====== 4 different input - Cache N=1 ======" << ql::def << endl;
  tp.setCacheSize(1); index = 0;
  tt.start();
  for (int i = 0; i < 1e7; i++)
    {
      tp.integral(res, mu2s[index], ms[index], p);
      index += 1; if (index == 4) index = 0;
    }
  tt.printTime(tt.stop());

  for (size_t i = 0; i < res.size(); i++)
    cout << "eps" << i << "\t" << res[i] << endl;

  cout << ql::blue << "\n====== 4 different input - Cache N=2 ======" << ql::def << endl;
  tp.setCacheSize(2); index = 0;
  tt.start();
  for (int i = 0; i < 1e7; i++)
    {
      tp.integral(res, mu2s[index], ms[index], p);
      index += 1; if (index == 4) index = 0;
    }
  tt.printTime(tt.stop());

  for (size_t i = 0; i < res.size(); i++)
    cout << "eps" << i << "\t" << res[i] << endl;

  cout << ql::blue << "\n====== 4 different input - Cache N=3 ======" << ql::def << endl;
  tp.setCacheSize(3); index = 0;
  tt.start();
  for (int i = 0; i < 1e7; i++)
    {
      tp.integral(res, mu2s[index], ms[index], p);
      index += 1; if (index == 4) index = 0;
    }
  tt.printTime(tt.stop());

  for (size_t i = 0; i < res.size(); i++)
    cout << "eps" << i << "\t" << res[i] << endl;

  cout << ql::blue << "\n====== 4 different input - Cache N=4 ======" << ql::def << endl;
  tp.setCacheSize(4); index = 0;
  tt.start();
  for (int i = 0; i < 1e7; i++)
    {
      tp.integral(res, mu2s[index], ms[index], p);
      index += 1; if (index == 4) index = 0;
    }
  tt.printTime(tt.stop());

  for (size_t i = 0; i < res.size(); i++)
    cout << "eps" << i << "\t" << res[i] << endl;

  cout << ql::blue << "\n====== 4 different input - Cache N=10 ======" << ql::def << endl;
  tp.setCacheSize(10); index = 0;
  tt.start();
  for (int i = 0; i < 1e7; i++)
    {
      tp.integral(res, mu2s[index], ms[index], p);
      index += 1; if (index == 4) index = 0;
    }
  tt.printTime(tt.stop());

  for (size_t i = 0; i < res.size(); i++)
    cout << "eps" << i << "\t" << res[i] << endl;

  return 0;
}
