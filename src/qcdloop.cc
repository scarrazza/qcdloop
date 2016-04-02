//
// QCDLoop 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch
//          Keith Ellis: keith.ellis@durham.ac.uk
//          Giulia Zanderighi: giulia.zanderighi@cern.ch

#include "qcdloop/qcdloop.h"
#include "qcdloop/tools.h"
#include "qcdloop/exceptions.h"
#include <iostream>

namespace ql
{
  template<typename TOutput, typename TMass, typename TScale>
  QCDLoop<TOutput,TMass,TScale>::QCDLoop():
    _tp(new TadPole<TOutput,TMass,TScale>()),
    _bb(new Bubble<TOutput,TMass,TScale>()),
    _tr(new Triangle<TOutput,TMass,TScale>()),
    _tb(new Box<TOutput,TMass,TScale>())
  {
  }

  template<typename TOutput, typename TMass, typename TScale>
  QCDLoop<TOutput,TMass,TScale>::~QCDLoop()
  {
    delete _tp;
    delete _bb;
    delete _tr;
    delete _tb;
  }

  /*!
   * \brief QCDLoop<TOutput, TScale, TMass>::I
   * \param res output object res[0,1,2] the coefficients in the Laurent series
   * \param mu2 is the square of the scale mu
   * \param m are the squares of the masses of the internal lines
   * \param p are the four-momentum suqared of the external lines
   */
  template<typename TOutput, typename TMass, typename TScale>
  void QCDLoop<TOutput,TMass,TScale>::integral(vector<TOutput> &res,
                                               TScale const& mu2,
                                               vector<TMass> const& m,
                                               vector<TScale> const& p) const
  {
    Topology<TOutput,TMass,TScale> *topos = nullptr;
    if (m.size() == 1 && p.size() == 0) topos = _tp;
    else if (m.size() == 2 && p.size() == 1) topos = _bb;
    else if (m.size() == 3 && p.size() == 3) topos = _tr;
    else if (m.size() == 4 && p.size() == 6) topos = _tb;
    else
      throw RangeError("QCDLoop","error topology not recognised");

    topos->integral(res, mu2, m, p);

    return;
  }

  template<typename TOutput, typename TMass, typename TScale>
  void QCDLoop<TOutput,TMass,TScale>::setCacheSize(int const& size)
  {
    _tp->setCacheSize(size);
    _bb->setCacheSize(size);
    _tr->setCacheSize(size);
    _tb->setCacheSize(size);
  }

  // explicity template declaration
  template class QCDLoop<complex,double,double>;
  template class QCDLoop<complex,complex,double>;
  template class QCDLoop<qcomplex,qdouble,qdouble>;
  template class QCDLoop<qcomplex,qcomplex,qdouble>;
}
