//
// QCDLoop 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch
//          Keith Ellis: keith.ellis@durham.ac.uk
//          Giulia Zanderighi: giulia.zanderighi@cern.ch

#include "qcdloop/bubble.h"
#include "qcdloop/tools.h"
#include "qcdloop/maths.h"
#include "qcdloop/exceptions.h"
#include <iostream>
using std::cout;
using std::endl;

namespace ql
{
  template<typename TOutput, typename TMass, typename TScale>
  Bubble<TOutput,TMass,TScale>::Bubble():
    Topology<TOutput,TMass,TScale>("Bubble")
  {
    this->_m.resize(2);
    this->_p.resize(1);
  }

  template<typename TOutput, typename TMass, typename TScale>
  Bubble<TOutput,TMass,TScale>::~Bubble()
  {
  }

  /*!
   * The integral is defined as:
   * \f[
   * I_{2}^{D=4-2 \epsilon}(p^2; m_1^2, m_2^2)= \mu^{2 \epsilon} \left[ \frac{1}{\epsilon} - \int_{0}^{1} da \ln (-a (1-a) p^2 + am_{2}^2 + (1-a) m_{1}^{2} - i \epsilon ) \right]  + O(\epsilon)
   *   \f]
   * Implementation of the formulae of Denner and Dittmaier \cite Denner:2005nn.
   *
   * \param res output object res[0,1,2] the coefficients in the Laurent series
   * \param mu2 is the squre of the scale mu
   * \param m are the squares of the masses of the internal lines
   * \param p are the four-momentum squared of the external lines
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Bubble<TOutput,TMass,TScale>::integral(vector<TOutput> &res,
                                              const TScale& mu2,
                                              vector<TMass> const& m,
                                              vector<TScale> const& p)
  {
    if (!this->checkCache(mu2,m,p))
      {        
        if (mu2 < 0) throw RangeError("Bubble::integral","mu2 is negative!");

        // Normalization
        const TScale scalefac = Max(Max(Max(Abs(p[0]),mu2), Abs(m[0])), Abs(m[1]));
        const TMass m0 = Min(m[0],m[1])/scalefac;
        const TMass m1 = Max(m[0],m[1])/scalefac;
        const TScale p0 = p[0]/scalefac;
        const TScale musq = mu2/scalefac;

        if (this->iszero(p0) && this->iszero(m0) && this->iszero(m1))
          std::fill(this->_val.begin(), this->_val.end(), this->_czero);
        else if (this->iszero(Abs(p0/musq)) && this->iszero(Abs(m0/musq)) && this->iszero((m1/musq)))
          {
            cout << ql::yellow << "Bubble::integral : settings s=m1=m2=0 self-energy to zero\n";
            cout << "s,m0,m1 = "   << p0 << ", " << m0 << ", " << m1 << ql::def << endl;
            this->_val[0] = this->_val[2] = this->_czero;
            this->_val[1] = this->_cone;
          }        
        else if (this->iszero(m0/musq))
          {
            if (this->iszero(Abs((m1-p0)/musq))) BB1(this->_val,musq,m1);    // I(s;0,s) s = m1, DD(4.13)            
            else if (this->iszero(Abs(p0/musq))) BB2(this->_val,musq,m1);    // I(0;0,m2)
            else if (this->iszero(Abs(m1/musq))) BB3(this->_val,musq,m1-TMass(p0)); // I(s;0,0)
            else                                 BB4(this->_val,musq,m1,p0); // I(s;0,m2)
          }
        else if (this->iszero(Abs(p0/musq))) // deal with special case, s = 0
          BB5(this->_val, musq, m0, m1);
        else
          BB0(this->_val, musq, m0, m1, p0);          

        this->storeCache(mu2,m,p);
      }

    if (res.size() != 3) { res.reserve(3); }
    std::copy(this->_val.begin(), this->_val.end(), res.begin());
    return;
  }

  /*!
   * The integral is defined as in the general case but with the first term:
   * \f[
   * I_{2}(s; m_0^2, m_1^2) = \Delta + 2 - \ln \left( \frac{m_0 m_1}{m_0^2} \right) + \frac{m_0^2-m_1^2}{s} \ln \left(\frac{m_1}{m_0} \right) - \frac{m_0 m_1}{s} \left(\frac{1}{r}-r \right) \ln r
   *   \f]
   * with
   * \f[
   * x^2 + \frac{m_0^2 + m_1^2 - s - i \epsilon}{m_0 m_1} +1 = (x+r)\left(x+\frac{1}{r} \right)
   * \f]
   * Implementation of the formulae of Denner and Dittmaier \cite Denner:1991kt.
   *
   * \return the output object res[0,1,2] the coefficients in the Laurent series
   * \param mu2 is the squre of the scale mu
   * \param m0 is the square of the first mass of the internal line
   * \param m1 is the square of the second mass of the internal line
   * \param s is the square of the four-momentum external line
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Bubble<TOutput,TMass,TScale>::BB0(vector<TOutput> &res, TScale const& mu2, TMass const& m0, TMass const& m1, TScale const& s) const
  {
    // general case from Denner 0709.1075, Eq. (4.23)
    const TMass sqm0 = Sqrt(m0);
    const TMass sqm1 = Sqrt(m1);
    const TOutput bb = TOutput(m0+m1-s);
    const TOutput rtt= Sqrt(bb*bb - this->_cfour*TOutput(m1*m0));
    const TOutput x1 = this->_chalf*(bb+rtt)/(sqm0*sqm1);
    const TOutput x2 = this->_cone/x1;
    res[0] = this->_ctwo - Log(sqm0*sqm1/mu2) + (m0-m1)/s*Log(sqm1/sqm0) - sqm0*sqm1/s*(x2-x1)*this->cLn(x1, Sign(Real(x1-x2)));
    res[1] = this->_cone;
    res[2] = this->_czero;
  }

  /*!
   * The integral is defined as:
   * \f[
   * I_{2}^{D=4-2 \epsilon}(m^2; 0, m^2)=\left( \frac{\mu^{2}}{m^{2}} \right)^{\epsilon} \left[ \frac{1}{\epsilon} + 2 \right]  + O(\epsilon)
   *   \f]
   * Implementation of the formulae of Denner and Dittmaier \cite Denner:2005nn.
   *
   * \return the output object res[0,1,2] the coefficients in the Laurent series
   * \param mu2 is the squre of the scale mu
   * \param m is the squares of the mass of the internal line
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Bubble<TOutput,TMass,TScale>::BB1(vector<TOutput> &res, TScale const& mu2, TMass const& m) const
  {
    res[0] = TOutput(Log(mu2/m)) + this->_ctwo;
    res[1] = this->_cone;
    res[2] = this->_czero;
  }

  /*!
   * The integral is defined as:
   * \f[
   * I_{2}^{D=4-2 \epsilon}(0; 0, m^2)= \left( \frac{\mu^{2}}{m^{2}} \right)^{\epsilon} \left[ \frac{1}{\epsilon} + 1 \right]  + O(\epsilon)
   *   \f]
   * Implementation of the formulae of Denner and Dittmaier \cite Denner:2005nn.
   *
   * \return the output object res[0,1,2] the coefficients in the Laurent series
   * \param mu2 is the squre of the scale mu
   * \param m is the squares of the mass of the internal line
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Bubble<TOutput,TMass,TScale>::BB2(vector<TOutput> &res, TScale const& mu2, TMass const& m) const
  {
    res[0] = TOutput(Log(mu2/m)) + this->_cone;
    res[1] = this->_cone;
    res[2] = this->_czero;
  }

  /*!
   * The integral is defined as:
   * \f[
   * I_{2}^{D=4-2 \epsilon}(s; 0, 0)= \left( \frac{\mu^{2}}{-s-i \epsilon} \right)^{\epsilon} \left[ \frac{1}{\epsilon} + 2 \right]  + O(\epsilon)
   *   \f]
   * Implementation of the formulae of Denner and Dittmaier \cite Denner:2005nn.
   *
   * \return the output object res[0,1,2] the coefficients in the Laurent series
   * \param mu2 is the squre of the scale mu
   * \param s is the four-momentum squared of the external line
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Bubble<TOutput,TMass,TScale>::BB3(vector<TOutput> &res, TScale const& mu2, TMass const& s) const
  {
    res[0] = -this->cLn(s/mu2, -1) + this->_ctwo;
    res[1] = this->_cone;
    res[2] = this->_czero;
  }

  /*!
   * The integral is defined as:
   * \f[
   * I_{2}^{D=4-2 \epsilon}(s; 0, m^2)= \left( \frac{\mu^{2}}{m^2} \right)^{\epsilon} \left[ \frac{1}{\epsilon} + 2 + \frac{m^2-s}{s} \ln \left( \frac{m^2-s-i \epsilon}{m^2} \right) \right]  + O(\epsilon)
   *   \f]
   * Implementation of the formulae of Denner and Dittmaier \cite Denner:2005nn.
   *
   * \return the output object res[0,1,2] the coefficients in the Laurent series
   * \param mu2 is the squre of the scale mu
   * \param m is the squares of the mass of the internal line
   * \param s is the four-momentum squared of the external line
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Bubble<TOutput,TMass,TScale>::BB4(vector<TOutput> &res, TScale const& mu2, TMass const& m, TScale const& s) const
  {
    res[0] = -this->cLn((m-s)/mu2, -1) + this->_cone - this->fndd(0, TOutput(this->_one-m/s), 1);
    res[1] = this->_cone;
    res[2] = this->_czero;
  }

  /*!
   * The integral is defined as in the general case but with the first term:
   * \f[
   * I_{2}^{\epsilon=0}(0; m_0^2, m_1^2) = \ln \left( \frac{\mu}{m_0^2} \right) - f_0 \left( \frac{m0}{m0-m1} \right)
   *   \f]
   * Implementation of the formulae of Denner and Dittmaier \cite Denner:2005nn.
   *
   * \return the output object res[0,1,2] the coefficients in the Laurent series
   * \param mu2 is the squre of the scale mu
   * \param m0 is the square of the first mass of the internal line
   * \param m1 is the square of the second mass of the internal line
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Bubble<TOutput,TMass,TScale>::BB5(vector<TOutput> &res, TScale const& mu2, TMass const& m0, TMass const& m1) const
  {
    res[0] = TOutput(Log(mu2/m0)); // m1 = m0
    if (!this->iszero(Abs((m1-m0)/mu2)))
      res[0] -= this->fndd(0,TOutput(m0/(m0-m1)),1); // other root is formally infinite
    res[1] = this->_cone;
    res[2] = this->_czero;
  }

  // explicity tyoename declaration
  template class Bubble<complex,double,double>;
  template class Bubble<complex,complex,double>;
  template class Bubble<qcomplex,qdouble,qdouble>;
  template class Bubble<qcomplex,qcomplex,qdouble>;
}
