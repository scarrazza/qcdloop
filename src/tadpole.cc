//
// QCDLoop 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch
//          Keith Ellis: keith.ellis@durham.ac.uk
//          Giulia Zanderighi: giulia.zanderighi@cern.ch

#include "qcdloop/tadpole.h"
#include "qcdloop/tools.h"
#include "qcdloop/maths.h"
#include "qcdloop/exceptions.h"

namespace ql
{
  template<typename TOutput, typename TMass, typename TScale>
  TadPole<TOutput,TMass,TScale>::TadPole():
    Topology<TOutput,TMass,TScale>("TadPole")
  {
    this->_m.resize(1);
    this->_p.resize(0);
  }

  template<typename TOutput, typename TMass, typename TScale>
  TadPole<TOutput,TMass,TScale>::~TadPole()
  {
  }

  /*!
   * \brief Computes the TadPole integral defined as:
   * \f[
   * I_{1}^{D=4-2 \epsilon}(m^2)= m^2 \left( \frac{\mu^2}{m^2-i \epsilon}\right) \left[ \frac{1}{\epsilon} +1 \right] + O(\epsilon)
   *   \f]
   *
   * Implementation of the formulae of Denner and Dittmaier \cite Denner:2005nn.
   *
   * \param res output object res[0,1,2] the coefficients in the Laurent series
   * \param mu2 is the squre of the scale mu
   * \param m are the squares of the masses of the internal lines
   * \param p are the four-momentum squared of the external lines
   */
  template<typename TOutput, typename TMass, typename TScale>
  void TadPole<TOutput,TMass,TScale>::integral(vector<TOutput> &res,
                                               const TScale& mu2,
                                               vector<TMass> const& m,
                                               vector<TScale> const& p)
  {
    if (!this->checkCache(mu2,m,p))
      {
        if (mu2 < 0) throw RangeError("TadPole::integral","mu2 is negative!");

        std::fill(this->_val.begin(), this->_val.end(), this->_czero);
        if (!this->iszero(m[0]))
          {
            this->_val[1] = TOutput(m[0]);
            this->_val[0] = this->_val[1]*TOutput(Log(mu2/m[0])+this->_cone);
          }          
        this->storeCache(mu2,m,p);
      }

    if (res.size() != 3) { res.reserve(3); }
    std::copy(this->_val.begin(), this->_val.end(), res.begin());

    return;
  }

  // explicity template declaration
  template class TadPole<complex,double,double>;
  template class TadPole<complex,complex,double>;
  template class TadPole<qcomplex,qdouble,qdouble>;
  template class TadPole<qcomplex,qcomplex,qdouble>;
}
