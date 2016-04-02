//
// QCDLoop 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch
//          Keith Ellis: keith.ellis@durham.ac.uk
//          Giulia Zanderighi: giulia.zanderighi@cern.ch

#include <iostream>
#include <cmath>
#include "qcdloop/tools.h"
#include "qcdloop/config.h"
#include "qcdloop/maths.h"
#include "qcdloop/exceptions.h"
using std::cout;
using std::endl;
using std::is_same;

namespace ql {    

  template<typename TOutput, typename TMass, typename TScale>
  Tools<TOutput,TMass,TScale>::Tools():
    _qlonshellcutoff(is_same<TScale,double>::value ? 1e-10 : 1e-20q),
    _pi(is_same<TScale,double>::value ? M_PI : M_PIq),
    _pi2   (_pi*_pi),
    _pio3  (_pi/TScale(3)),
    _pio6  (_pi/TScale(6)),
    _pi2o3 (_pio3*_pi),
    _pi2o6 (_pio6*_pi),
    _pi2o12(_pi2/TScale(12)),
    _zero (is_same<TScale,double>::value ? 0.0 : 0.0q),
    _half (is_same<TScale,double>::value ? 0.5 : 0.5q),
    _one  (is_same<TScale,double>::value ? 1.0 : 1.0q),
    _two  (is_same<TScale,double>::value ? 2.0 : 2.0q),
    _three(is_same<TScale,double>::value ? 3.0 : 3.0q),
    _four (is_same<TScale,double>::value ? 4.0 : 4.0q),
    _five (is_same<TScale,double>::value ? 5.0 : 5.0q),
    _six  (is_same<TScale,double>::value ? 6.0 : 6.0q),
    _ten  (is_same<TScale,double>::value ? 10.0 : 10.0q),
    _eps  (is_same<TScale,double>::value ? 1e-6 : 1e-12q),
    _eps4 (is_same<TScale,double>::value ? 1e-4 : 1e-4q),
    _eps7 (is_same<TScale,double>::value ? 1e-7 : 1e-7q),
    _eps10(is_same<TScale,double>::value ? 1e-10 : 1e-10q),
    _eps14(is_same<TScale,double>::value ? 1e-14 : 1e-14q),
    _eps15(is_same<TScale,double>::value ? 1e-15 : 1e-15q),
    _xloss(is_same<TScale,double>::value ? 0.125 : 0.125q),
    _neglig(is_same<TScale,double>::value ? 1e-14 : 1e-28),
    _reps  (is_same<TScale,double>::value ? 1e-16 : 1e-32q),
    _2ipi(TOutput{_zero,_two*_pi}),
    _ipio2(TOutput{_zero, _pi*_half}),
    _ipi(TOutput{_zero, _pi}),
    _czero (is_same<TOutput,complex>::value ? TOutput(0.0) : TOutput(0.0q)),
    _chalf (is_same<TOutput,complex>::value ? TOutput(0.5) : TOutput(0.5q)),
    _cone  (is_same<TOutput,complex>::value ? TOutput(1.0) : TOutput(1.0q)),
    _ctwo  (is_same<TOutput,complex>::value ? TOutput(2.0) : TOutput(2.0q)),
    _cthree(is_same<TOutput,complex>::value ? TOutput(3.0) : TOutput(3.0q)),
    _cfour (is_same<TOutput,complex>::value ? TOutput(4.0) : TOutput(4.0q)),
    _ieps (TOutput{_zero, _reps}),
    _ieps2(TOutput{_zero, _reps*_reps}),
    _ieps50(is_same<TOutput,complex>::value ? TOutput{_zero, 1e-50} : TOutput{_zero, 1e-50q})
  {
    // Show splash boot logo
    Splash::Show();

    // Allocate bernoulli coefficients
    if (is_same<TScale,double>::value)
      {
        _C = {
           0.4299669356081370,
           0.4097598753307711,
          -0.0185884366501460,
           0.0014575108406227,
          -0.0001430418444234,
           0.0000158841554188,
          -0.0000019078495939,
           0.0000002419518085,
          -0.0000000319334127,
           0.0000000043454506,
          -0.0000000006057848,
           0.0000000000861210,
          -0.0000000000124433,
           0.0000000000018226,
          -0.0000000000002701,
           0.0000000000000404,
          -0.0000000000000061,
           0.0000000000000009,
          -0.0000000000000001
        };

        _B = {
           0.02777777777777777777777777777777777777777778774E0,
          -0.000277777777777777777777777777777777777777777778E0,
           4.72411186696900982615268329554043839758125472E-6,
          -9.18577307466196355085243974132863021751910641E-8,
           1.89788699889709990720091730192740293750394761E-9,
          -4.06476164514422552680590938629196667454705711E-11,
           8.92169102045645255521798731675274885151428361E-13,
          -1.993929586072107568723644347793789705630694749E-14,
           4.51898002961991819165047655285559322839681901E-16,
          -1.035651761218124701448341154221865666596091238E-17,
           2.39521862102618674574028374300098038167894899E-19,
          -5.58178587432500933628307450562541990556705462E-21,
           1.309150755418321285812307399186592301749849833E-22,
          -3.087419802426740293242279764866462431595565203E-24,
           7.31597565270220342035790560925214859103339899E-26,
          -1.740845657234000740989055147759702545340841422E-27,
           4.15763564461389971961789962077522667348825413E-29,
          -9.96214848828462210319400670245583884985485196E-31,
           2.394034424896165300521167987893749562934279156E-32,
          -5.76834735536739008429179316187765424407233225E-34,
           1.393179479647007977827886603911548331732410612E-35,
          -3.372121965485089470468473635254930958979742891E-37,
           8.17820877756210262176477721487283426787618937E-39,
          -1.987010831152385925564820669234786567541858996E-40,
           4.83577851804055089628705937311537820769430091E-42
        };

      }
    else
      {
        _C = {
           0.4299669356081370q,
           0.4097598753307711q,
          -0.0185884366501460q,
           0.0014575108406227q,
          -0.0001430418444234q,
           0.0000158841554188q,
          -0.0000019078495939q,
           0.0000002419518085q,
          -0.0000000319334127q,
           0.0000000043454506q,
          -0.0000000006057848q,
           0.0000000000861210q,
          -0.0000000000124433q,
           0.0000000000018226q,
          -0.0000000000002701q,
           0.0000000000000404q,
          -0.0000000000000061q,
           0.0000000000000009q,
          -0.0000000000000001q
        };

        _B = {
          0.02777777777777777777777777777777777777777778774E0q,
         -0.000277777777777777777777777777777777777777777778E0q,
          4.72411186696900982615268329554043839758125472E-6q,
         -9.18577307466196355085243974132863021751910641E-8q,
          1.89788699889709990720091730192740293750394761E-9q,
         -4.06476164514422552680590938629196667454705711E-11q,
          8.92169102045645255521798731675274885151428361E-13q,
         -1.993929586072107568723644347793789705630694749E-14q,
          4.51898002961991819165047655285559322839681901E-16q,
         -1.035651761218124701448341154221865666596091238E-17q,
          2.39521862102618674574028374300098038167894899E-19q,
         -5.58178587432500933628307450562541990556705462E-21q,
          1.309150755418321285812307399186592301749849833E-22q,
         -3.087419802426740293242279764866462431595565203E-24q,
          7.31597565270220342035790560925214859103339899E-26q,
         -1.740845657234000740989055147759702545340841422E-27q,
          4.15763564461389971961789962077522667348825413E-29q,
         -9.96214848828462210319400670245583884985485196E-31q,
          2.394034424896165300521167987893749562934279156E-32q,
         -5.76834735536739008429179316187765424407233225E-34q,
          1.393179479647007977827886603911548331732410612E-35q,
         -3.372121965485089470468473635254930958979742891E-37q,
          8.17820877756210262176477721487283426787618937E-39q,
         -1.987010831152385925564820669234786567541858996E-40q,
          4.83577851804055089628705937311537820769430091E-42q
        };
      }

  }

  template<typename TOutput, typename TMass, typename TScale>
  Tools<TOutput,TMass,TScale>::~Tools()
  {
    _C.clear();
    _B.clear();
  }

  /*!
   * Compares the Abs(psq) to the cutoff
   * \param msq the input to be tested
   * \return true if |psq| < onshell cut-off
   */
  template<typename TOutput, typename TMass, typename TScale>
  bool Tools<TOutput,TMass,TScale>::iszero(TMass const& psq) const
  {
    return (Abs(psq) < _qlonshellcutoff) ? true : false;
  }

  /*!
   * Computes the log of a complex number z.
   * If the imag(z)=0 and real(z)<0 and extra ipi term is included.
   * \param z the complex argument of the logarithm
   * \param isig the sign of the imaginary part
   * \return the complex log(z)
   */
  template<typename TOutput, typename TMass, typename TScale>
  TOutput Tools<TOutput,TMass,TScale>::cLn(TOutput const& z, TScale const& isig) const
  {
    TOutput cln;
    if (Imag(z) == _zero && Real(z) <= _zero)
      cln = Log(-z) + TOutput{_zero, _pi*Sign(isig)};
    else
      cln = Log(z);
    return cln;
  }

  /*!
   * Computes the log of a real number x.
   * If the x<0 and extra ipi term is included.
   * \param x the real argument of the logarithm
   * \param isig the sign of the imaginary part
   * \return the complex log(x)
   */
  template<typename TOutput, typename TMass, typename TScale>
  TOutput Tools<TOutput,TMass,TScale>::cLn(TScale const& x, TScale const& isig) const
  {
    TOutput ln;
    if (x > 0)
      ln = TOutput(Log(x));
    else
      ln = TOutput(Log(-x)) + TOutput{_zero, _pi*Sign(isig)};
    return ln;
  }

  /*!
   * Implementation of the formulae of Denner and Dittmaier \cite Denner:2005nn.
   * \f[
   * f_{n}(x) = \ln \left( 1 - \frac{1}{x} \right) + \sum_{l=n+1}^{\infty} \frac{x^{n-l}}{l+1}
   *   \f]
   *
   * \param n the lower index of
   * \param x the argument of the function
   * \param iep the epsilon value
   * \return function DD from Eq. 4.11 of \cite Denner:2005nn.
   */
  template<typename TOutput, typename TMass, typename TScale>
  TOutput Tools<TOutput,TMass,TScale>::fndd(int const& n, TOutput const& x, TScale const& iep) const
  {
    const int infty = 16;    
    TOutput res = _czero;    
    if (Abs(x) < _ten)
      {
        if (!iszero(Abs(x-_cone)))
          res = (_cone-Pow(x,n+1))*(cLn(x-_cone, iep) - cLn(x, iep));

        for (int j = 0; j <= n; j++)
          res -= Pow(x, n-j)/(j+_one);
      }    
    else
      {        
        res = cLn(_cone-_cone/x, iep);
        for (int j = n+1; j <= n+infty; j++)
          res += Pow(x, n-j)/(j+_one);          
      }
    return res;
  }

  /*!
   * Computes \f${\rm Lnrat}(x,y) = \log(x-i \epsilon)-\log(y-i \epsilon)\f$
   * \param x TMass object for the numerator
   * \param y TMass object for the denumerator
   * \return returns the ratio of logs
   */
  template<typename TOutput, typename TMass, typename TScale>
  TOutput Tools<TOutput,TMass,TScale>::Lnrat(TOutput const& x, TOutput const& y) const
  {
    const TOutput r = x/y;
    if (iszero(Imag(r)))
      return TOutput(Log(Abs(r))) - _ipio2*TOutput(Sign(-Real(x))-Sign(-Real(y)));
    else
      return Log(r);
  }

  template<typename TOutput, typename TMass, typename TScale>
  TOutput Tools<TOutput,TMass,TScale>::Lnrat(TScale const& x, TScale const& y) const
  {
    return TOutput(Log(Abs(x/y))) - _ipio2*TOutput(Sign(-x)-Sign(-y));
  }

  /*!
   * The dilog function for real argument
   * \param x the argument of the dilog
   * \return the dilog
   */
  template<typename TOutput, typename TMass, typename TScale>
  TMass Tools<TOutput,TMass,TScale>::ddilog(TMass const& x) const
  {
    if ( x == _one )
      return _pi2o6;
    else if (x == -_one)
      return -_half*_pi2o6;

    const TMass T = -x;
    TMass Y, S, A;

    if (Real(T) <= -_two)
      {
        Y = -_one / (_one+T);
        S = _one;
        A = TMass(-_pi2o3 + _half*Real(Pow(Log(-T),2) - Pow(Log(_one+_one/T),2)));
      }
    else if (Real(T) < -_one)
      {
        Y = -_one - T;
        S = -_one;
        A = Log(-T);
        A = TMass(-_pi2o6 + A*(A+Log(_one+_one/T)));
      }
    else if (Real(T) <= -_half)
      {
        Y = (-_one-T)/T;
        S = _one;
        A = Log(-T);
        A = TMass(-_pi2o6 + A*(-_half*A+Log(_one+T)));
      }
    else if (Real(T) < _zero)
      {
        Y = -T/(_one+T);
        S = -_one;
        A = TMass(_half*Real(Pow(Log(_one+T),2)));
      }
    else if (Real(T) <= _one)
      {
        Y = T;
        S = _one;
        A = TMass(_zero);
      }
    else
      {
        Y = _one/T;
        S = -_one;
        A = TMass(_pi2o6 + _half*Real(Pow(Log(T),2)));
      }

    const TMass H = Y+Y-_one;
    const TMass ALFA = H+H;
    TMass B1 = _zero, B2 = _zero, B0;
    for (int i = 18; i >= 0; i--)
      {
        B0 = _C[i]+ALFA*B1-B2;
        B2 = B1;
        B1 = B0;
      }

    return -(S*(B0-H*B2)+A);
  }

  /*!
   * \param z the argument
   * \param isig the sign of z
   */
  template<typename TOutput, typename TMass, typename TScale>
  TOutput Tools<TOutput,TMass,TScale>::denspence(TOutput const& z, TScale const& isig) const
  {
    const TOutput z1 = _cone - z;
    const TScale az1 = Abs(z1);

    if (isig == _zero && Imag(z) == _zero && Abs(Real(z1)) < _qlonshellcutoff)
      cout << "denspence: argument on cut" << endl;

    if (az1 < _eps15)
      return TOutput{_pi2o6, _zero};
    else if (Real(z) < _half)
      {
        if (Abs(z) < _one)
          return li2series(z, isig);
        else
          return -_pi2o6 - _half*Pow(cLn(-z, -isig),2) - li2series(_one/z, -isig);
      }
    else
      {
        if (az1 < _one)
          return _pi2o6 - cLn(z, isig)*cLn(z1, -isig) - li2series(z1, -isig);
        else
          return _ctwo*_pi2o6 + _half*Pow(cLn(-z1, -isig),2) - cLn(z, isig)*cLn(z1, -isig) + li2series(_one/z1, isig);
      }
  }

  /*!
   * \param zrat1 first argument
   * \param zrat2 second argument
   * \param ieps1 sign for zrat1
   * \param ieps2 sign for zrat2
   * \return
   */
  template<typename TOutput, typename TMass, typename TScale>
  TOutput Tools<TOutput,TMass,TScale>::spencer(TOutput const& zrat1, TOutput const& zrat2, TScale const&ieps1, TScale const&ieps2) const
  {
    const TScale x1 = Real(zrat1);
    const TScale x2 = Real(zrat2);
    const TScale y1 = Imag(zrat1);
    const TScale y2 = Imag(zrat2);

    TOutput res, prod;
    if (iszero(y1) && iszero(y2))
      res = Li2omx(x1,x2,ieps1,ieps2);
    else
      {
        TOutput arg = zrat1*zrat2;
        const TScale ieps = _zero;
        if (Abs(arg) <= _one)
          {
            if (arg == _zero || arg == _one)
              prod = _czero;
            else
              {
                const TOutput lnarg = cLn(zrat1, ieps1) + cLn(zrat2, ieps2);
                const TOutput lnomarg = Log(_cone-arg);
                prod = lnarg*lnomarg;
              }
            res = TOutput(_pi2o6) - denspence(arg, ieps) - prod;
          }
        else if (Abs(arg) > _one)
          {
            arg = _cone/(zrat1*zrat2);
            const TOutput lnarg = -cLn(zrat1, ieps1)-cLn(zrat2, ieps2);
            const TOutput lnomarg = Log(_cone-arg);
            res = -TOutput(_pi2o6)+denspence(arg, ieps)+lnarg*lnomarg-_chalf*Pow(lnarg,2);
          }
      }
    return res;
  }

  /*!
   * \param z2 input arguments
   * \param im2 input signs.
   * \return the difference of cspence functions
   */
  template<typename TOutput, typename TMass, typename TScale>
  TOutput Tools<TOutput,TMass,TScale>::xspence(TOutput const (&z1)[2], TScale const (&im1)[2], TOutput const& z2, TScale const& im2) const
  {
    return cspence(z1[1],im1[1],z2,im2)-cspence(z1[0],im1[0],z2,im2);
  }

  /*!
   * \param z1 input argument
   * \param im1 sign of z1
   * \param z2 input argument
   * \param im2 sign of z2
   * \return the complex Spence's function
   */
  template<typename TOutput, typename TMass, typename TScale>
  TOutput Tools<TOutput,TMass,TScale>::cspence(TOutput const& z1, TScale const& im1, TOutput const& z2, TScale const& im2) const
  {
    TOutput cspence = _czero;
    const TOutput z12 = z1*z2;
    const TScale im12 = im2*Sign(Real(z1));

    if (Real(z12) > _half)
      {
        cspence = ltspence(1, z12, _zero);       
        const int etas = eta(z1, im1, z2, im2, im12);
        if (etas != 0) cspence += TOutput(etas)*cLn(_cone-z12, -im12)*_2ipi;
      }
    else if (Abs(z12) < _eps4)
      {
        cspence = TOutput(_pi2o6);        
        if (Abs(z12) > _eps14)
          cspence += -ltspence(0, z12, _zero) + (cLn(z1,im1) + cLn(z2,im2))*z12*(_cone + z12*(_chalf + z12*(_cone/_cthree + z12/_cfour)));
      }
    else
      cspence = TOutput(_pi2o6) - ltspence(0, z12, _zero) - (cLn(z1, im1) + cLn(z2, im2))*cLn(_cone-z12,_zero);      

    return cspence;
  }

  /*!
   * \param x numerator
   * \param y denominator
   * \return the Li2 ratio
   */
  template<typename TOutput, typename TMass, typename TScale>
  TOutput Tools<TOutput,TMass,TScale>::Li2omrat(TScale const& x, TScale const& y) const
  {
    const TScale omarg = x/y;
    const TScale arg = _one-omarg;
    if (arg > _one)
      return TOutput(_pi2o6-ddilog(omarg)) - Log(arg)*Lnrat(x,y);
    else
      return TOutput(ddilog(arg));
  }

  template<typename TOutput, typename TMass, typename TScale>
  TOutput Tools<TOutput,TMass,TScale>::Li2omrat(TOutput const& x, TOutput const& y,
                                                TScale const& ieps1, TScale const& ieps2) const
  {
    const TOutput omarg = x/y;
    const TOutput arg = _cone-omarg;
    const TScale isarg = Sign(Real(x)*ieps2-Real(y)*ieps1);

    if (Abs(arg) > _one)
      {
        const TScale isomarg = -isarg;
        return TOutput(_pi2o6) - denspence(omarg, isomarg) - cLn(omarg, isomarg)*cLn(arg, isarg);
      }
    else
      return denspence(arg, isarg);
  }

  /*!
   * expression for dilog(1-(v-i*ep)*(w-i*ep)/(x-i*ep)/(y-i*ep)) for real v,w,x and y
   * \param v numerator
   * \param w numerator
   * \param x denominator
   * \param y denominator
   * \return the dilog ratio
   */
  template<typename TOutput, typename TMass, typename TScale>
  TOutput Tools<TOutput,TMass,TScale>::Li2omx2(TScale const& v, TScale const& w, TScale const& x, TScale const& y) const
  {
    const TScale arg = (v*w)/(x*y);
    const TScale omarg = _one-arg;
    TOutput prod, Li2omx2;

    if (Abs(arg) <= _one)
      {
        if (Abs(arg) == _zero || Abs(arg) == _one)
          prod = _zero;
        else
          prod = (Lnrat(v,x)+Lnrat(w,y))*TOutput(Log(omarg));
        Li2omx2 = TOutput(_pi2o6-ddilog(arg))-prod;
      }
    else if (Abs(arg) > _one)
      {
        const TScale arg2 = (x*y)/(v*w);
        const TOutput lnarg = TOutput(-Lnrat(v,x)-Lnrat(w,y));
        const TOutput lnomarg = TOutput(Log(_one-arg2));
        Li2omx2 = -TOutput(_pi2o6-ddilog(arg2))+lnarg*lnomarg-_chalf*lnarg*lnarg;
      }

    return Li2omx2;
  }

  /*!
   * Calculates Li[2](1-(z1+ieps1)*(z2+ieps2)) for complex z1,z2
   * Using +Li2(1-z1*z2)                           for z1*z2<1
   * and   -Li2(1-1/(z1*z2))-1/2*(ln(z1)+ln(z2))^2 for z1*z2>1
   * \param z1 input argument
   * \param z2 input argument
   * \param ieps1 sign of z1
   * \param ieps2 sign of z2
   * \return Li2 of the 1-product
   */
  template<typename TOutput, typename TMass, typename TScale>
  TOutput Tools<TOutput,TMass,TScale>::Li2omx2(TOutput const& v, TOutput const& w, TOutput const& x, TOutput const& y, TScale const& ieps1, TScale const& ieps2) const
  {
    return cLi2omx2(v/x, w/y, ieps1, ieps2);
  }

  /*!
   * Calculates Li[2](1-(z1+ieps1)*(z2+ieps2)) for complex z1,z2
   * Using +Li2(1-z1*z2)                           for z1*z2<1
   * and   -Li2(1-1/(z1*z2))-1/2*(ln(z1)+ln(z2))^2 for z1*z2>1
   * \param z1 input argument
   * \param z2 input argument
   * \param ieps1 sign of z1
   * \param ieps2 sign of z2
   * \return Li2 of the 1-product
   */
  template<typename TOutput, typename TMass, typename TScale>
  TOutput Tools<TOutput,TMass,TScale>::cLi2omx2(TOutput const& z1, TOutput const& z2, TScale const& ieps1, TScale const& ieps2) const
  {
    const TOutput arg = z1*z2;
    const TScale ieps = Sign(Real(z2)*ieps1+Real(z1)*ieps2);

    TOutput prod, res;
    if (Abs(arg) <= _one)
      {
        if (arg == _czero || arg == _cone)
          prod = _czero;
        else
          prod = (cLn(z1, ieps1)+cLn(z2, ieps2)) * cLn(_cone-arg, -ieps);
        res = TOutput(_pi2o6) - denspence(arg, ieps) - prod;
      }
    else if (Abs(arg) > _one)
      {
        const TOutput arg2 = _cone/(z1*z2);
        const TOutput lnomarg = cLn(_cone-arg2, -ieps);
        const TOutput lnarg = -cLn(z1, ieps1) - cLn(z2, ieps2);
        res = TOutput(-_pi2o6) + denspence(arg2, ieps) + lnarg*lnomarg - _chalf*lnarg*lnarg;
      }
    return res;
  }

  /*!
   * Calculate Li[2](1-(x1+ieps1)*(x2+ieps2)) for real x1,x2
   * Using +Li2(1-x1*x2)                           for x1*x2<1
   * and   -Li2(1-1/(x1*x2))-1/2*(ln(x1)+ln(x2))^2 for x1*x2>1
   * \param x1 numerator
   * \param x2 denominator
   * \param ieps1 sign of x1
   * \param ieps2 sign of x2
   * \return the ratio Li2
   */
  template<typename TOutput, typename TMass, typename TScale>
  TOutput Tools<TOutput,TMass,TScale>::Li2omx(TMass const& x1, TMass const& x2, TScale const& ieps1, TScale const& ieps2) const
  {
    TOutput prod, Li2omx;

    TMass arg = x1*x2;
    const TScale ieps = Sign(Real(x2)*ieps1+Real(x1)*ieps2);
    if (Real(arg) <= _one)
      {
        if (Real(arg) == _one || Real(arg) == _zero)
          prod = _czero;
        else
          {
            const TOutput lnarg = cLn(Real(x1), ieps) + cLn(Real(x2), ieps2);
            const TOutput lnomarg = TOutput(Log(_one-arg));
            prod = lnarg*lnomarg;
          }
        Li2omx = TOutput(_pi2o6)-denspence(TOutput(arg), ieps)-prod;
      }
    else if (Real(arg) > _one)
      {
        arg = _one/(x1*x2);
        const TOutput lnarg = -cLn(Real(x1), ieps1)-cLn(Real(x2), ieps2);
        const TOutput lnomarg = TOutput(Log(_one-arg));
        Li2omx = -TOutput(_pi2o6)+denspence(TOutput(arg), ieps) + lnarg*lnomarg - _chalf*lnarg*lnarg;
      }
    return Li2omx;
  }

  /*!
   * Generalization of cLi2omx2 for 3 arguments.
   */
  template<typename TOutput, typename TMass, typename TScale>
  TOutput Tools<TOutput,TMass,TScale>::cLi2omx3(TOutput const& z1, TOutput const& z2, TOutput const& z3, TScale const& ieps1, TScale const& ieps2, TScale const& ieps3) const
  {
    const TOutput arg = z1*z2*z3;

    TScale ieps;
    if (iszero(Imag(arg)))
      ieps = Sign(Real(z2*z3)*ieps1+Real(z1*z3)*ieps2+Real(z1*z2)*ieps3);

    TOutput res = _czero;
    if (Abs(arg) <= _one)
      {
        TOutput prod;
        if (arg == _czero || arg == _cone)
          prod = _zero;
        else
          {
            const TOutput lnarg = cLn(z1, ieps1) + cLn(z2, ieps2) + cLn(z3, ieps3);
            const TOutput lnomarg = cLn(_cone-arg, _zero);
            prod = lnarg*lnomarg;
          }
        res = TOutput(_pi2o6)-denspence(arg, ieps)-prod;
      }
    else
      {
        const TOutput arg2 = _cone/(z1*z2*z3);
        const TOutput lnarg = -cLn(z1, ieps1) - cLn(z2,ieps2) - cLn(z3,ieps3);
        const TOutput lnomarg = cLn(_cone-arg2, _zero);
        res = - TOutput(_pi2o6) + denspence(arg2, ieps) + lnarg*lnomarg - _chalf*lnarg*lnarg;
      }

    return res;
  }

  /*!
   * \param x input argument
   * \param y input argument
   * \return
   */
  template<typename TOutput, typename TMass, typename TScale>
  TOutput Tools<TOutput,TMass,TScale>::L0(TMass const& x, TMass const& y) const
  {
    TOutput L0;
    const TMass denom = _one - x/y;
    if (Abs(denom) < _eps7)
      L0 = -_cone - TOutput(denom*(_half+denom/_three));
    else
      L0 = Lnrat(x,y)/TOutput(denom);

    return L0;
  }

  /*!
   * \param x input argument
   * \param y input argument
   * \return
   */
  template<typename TOutput, typename TMass, typename TScale>
  TOutput Tools<TOutput,TMass,TScale>::L1(TMass const& x, TMass const& y) const
  {
    TOutput L1;
    const TMass denom = _one - x/y;
    if (Abs(denom) < _eps7)
      L1 = -_cone*_chalf - TOutput(denom/_three*(_one+_three*denom/_four));
    else
      L1 = (L0(x,y)+_cone)/TOutput(denom);

    return L1;
  }

  /*!
   * Finite Triangle Li2 sum.
   */
  template<typename TOutput, typename TMass, typename TScale>
  TOutput Tools<TOutput,TMass,TScale>::R3int(TOutput const& a, TOutput const& s1, TOutput const& s2, TOutput const& t1, TOutput const& t2, TOutput const& t3, TOutput const& t4) const
  {
    const TOutput b = (s1+s2)*(s1-s2)-a;
    const TOutput c = s2*s2;
    const TOutput d = Sqrt((a-(s1+s2)*(s1+s2))*(a-(s1-s2)*(s1-s2)));

    TOutput y[2], s[2], res;
    solveabcd(a,b,c,d,y);
    solveabcd(a,t2,t3,t4,s);

    const TOutput y0 = -(t1+b*s[0])/t4;
    const TOutput dq0 = (y0-y[0]);
    const TOutput dq1 = (y0-y[1]);

    const TOutput OneOdq0 = _cone/dq0;
    const TOutput OneMy0 = _cone-y[0];
    const TScale SignImagOneOdq0 = Sign(Imag(OneOdq0));

    const TOutput OneOdq1 = _cone/dq1;
    const TOutput OneMy1 = _cone-y[1];
    const TScale SignImagOneOdq1 = Sign(Imag(OneOdq1));

    res =  cspence(-y[0], Sign(Imag(-y[0])), OneOdq0, SignImagOneOdq0)
          -cspence(OneMy0, Sign(Imag(OneMy0)), OneOdq0, SignImagOneOdq0)
          +cspence(-y[1], Sign(Imag(-y[1])), OneOdq1, SignImagOneOdq1)
          -cspence(OneMy1, Sign(Imag(OneMy1)), OneOdq1, SignImagOneOdq1);

    TOutput zz = y0*(a*y0+b);
    if (Abs(Real(zz))*_reps*_reps <= Abs(Imag(zz))*_neglig && Abs(Imag(zz)) <= Abs(Real(zz))*_neglig)
      zz = (TOutput(Real(zz))+c)/a;
    else
      zz = (zz+c)/a;

    // ajust complex logs
    TOutput extra = eta3(-y[0],-y[1],c/a) - eta3(dq0,dq1,zz);
    if (Real(a) < _zero && Imag(zz) < _zero) extra -= _2ipi;
    if (extra != _czero)
      {
        const TOutput arg4 = (y0-_cone)/y0;
        res += extra*cLn(arg4, Sign(Imag(arg4)));
      }

    return res;
  }

  /*!
   * Finite Triangle Li2 sum.
   */
  template<typename TOutput, typename TMass, typename TScale>
  TOutput Tools<TOutput,TMass,TScale>::R3int(TOutput const& a, TOutput const& s1, TOutput const& s2, TOutput const& t1) const
  {
    const TOutput b = (s1+s2)*(s1-s2)-a;
    const TOutput c = s2*s2;
    const TOutput d = Sqrt((a-(s1+s2)*(s1+s2))*(a-(s1-s2)*(s1-s2)));

    TOutput y[2], res;
    solveabcd(a,b,c,d,y);

    const TOutput y0 = t1;
    const TOutput dq0 = (y0-y[0]);
    const TOutput dq1 = (y0-y[1]);

    const TOutput OneOdq0 = _cone/dq0;
    const TOutput OneMy0 = _cone-y[0];
    const TScale SignImagOneOdq0 = Sign(Imag(OneOdq0));

    const TOutput OneOdq1 = _cone/dq1;
    const TOutput OneMy1 = _cone-y[1];
    const TScale SignImagOneOdq1 = Sign(Imag(OneOdq1));

    res =  cspence(-y[0], Sign(Imag(-y[0])), OneOdq0, SignImagOneOdq0)
          -cspence(OneMy0, Sign(Imag(OneMy0)), OneOdq0, SignImagOneOdq0)
          +cspence(-y[1], Sign(Imag(-y[1])), OneOdq1, SignImagOneOdq1)
          -cspence(OneMy1, Sign(Imag(OneMy1)), OneOdq1, SignImagOneOdq1);

    TOutput zz = y0*(a*y0+b);
    if (Abs(Real(zz))*_reps*_reps <= Abs(Imag(zz))*_neglig && Abs(Imag(zz)) <= Abs(Real(zz))*_neglig)
      zz = (TOutput(Real(zz))+c)/a;
    else
      zz = (zz+c)/a;

    // ajust complex logs
    TOutput extra = eta3(-y[0],-y[1],c/a) - eta3(dq0,dq1,zz);
    if (Real(a) < _zero && Imag(zz) < _zero) extra -= _2ipi;
    if (extra != _czero)
      {
        const TOutput arg4 = (y0-_cone)/y0;
        res += extra*cLn(arg4, Sign(Imag(arg4)));
      }      

    return res;
  }

  /*!
   * Finite Triangle Li2 sum.
   */
  template<typename TOutput, typename TMass, typename TScale>
  TOutput Tools<TOutput,TMass,TScale>::R2int(TOutput const& a, TOutput const& b, TOutput const& y0) const
  {
    const TOutput y1 = -b/a;
    const TOutput dq0 = y0-y1;

    const TOutput oneOdq0 = _cone/dq0;
    const TScale  SignImagOneOdq0 = Sign(Imag(oneOdq0));
    const TOutput oneMy1 = _cone-y1;

    TOutput res = cspence(-y1, Sign(Imag(-y1)), oneOdq0, SignImagOneOdq0)
                 -cspence(oneMy1, Sign(Imag(oneMy1)), oneOdq0, SignImagOneOdq0);

    TOutput extra = eta5(a,-y1,b,dq0,a*dq0);
    if (extra != _czero)
      {
        const TOutput arg4 = (y0-_cone)/y0;
        res += extra*cLn(arg4, Sign(Imag(arg4)));
      }
    return res;
  }

  /*!
   * \param y
   * \param z
   * \param ieps
   * \return
   */
  template<typename TOutput, typename TMass, typename TScale>
  TOutput Tools<TOutput,TMass,TScale>::Rint(TOutput const& y, TOutput const& z, TScale const& ieps) const
  {
    const TOutput omz = _cone-z;
    const TOutput ymone = y-_cone;
    const TOutput invymz = _cone/(y-z);
    const TOutput c[2] = { y*invymz, ymone*invymz };

    if (Imag(z) != _zero)
      {
        const TOutput c2ipi = TOutput{_zero, _two*_pi};
        const TOutput a = -z;
        const TOutput b = invymz;
        const TOutput ab= -z*invymz;
        const TOutput eta1 = c2ipi*TOutput( Htheta(-Imag(a))*Htheta(-Imag(b))*Htheta(Imag(ab))
                                     -Htheta(Imag(a))*Htheta(Imag(b))*Htheta(-Imag(ab)) );
        const TOutput a2 = omz;
        const TOutput ab2= omz*invymz;
        const TOutput eta2 = c2ipi*TOutput( Htheta(-Imag(a2))*Htheta(-Imag(b))*Htheta(Imag(ab2))
                                     -Htheta(Imag(a2))*Htheta(Imag(b))*Htheta(-Imag(ab2)) );

        TOutput logc1 = _czero, logc2 = _czero;
        if (eta1 != _czero) logc1 = TOutput(Log(c[0]));
        if (eta2 != _czero) logc2 = TOutput(Log(c[1]));

        return denspence(c[0],_zero) - denspence(c[1],_zero) + eta1*logc1 - eta2*logc2;
      }
    else
      {
        const TScale ieps1 = -ieps*Sign(Real(y));
        const TScale ieps2 = -ieps*Sign(Real(ymone));
        return denspence(c[0], ieps1) - denspence(c[1], ieps2);
      }
  }

  /*!
   * \brief Tools<TOutput, TMass, TScale>::R
   * \param q
   * \return
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Tools<TOutput,TMass,TScale>::R(TOutput &r, TOutput &d, TOutput const& q) const
  {
    d = Sqrt(q*q-_four);
    r  = q+d;
    TOutput r2 = q-d;
    TScale a = Abs(r), b = Abs(r2);
    if (b > a) { r = r2; d = -d; }
    a = Imag(q);
    b = Imag(r);
    if (a == _zero)
      {
        if (b <= _zero)
          r /= _ctwo;
        else
          {
            r = _ctwo/r;
            d = -d;
          }
      }
    else
      {
        const TScale ik = Sign(Real(a));
        const TScale ir = Sign(Real(b));
        if (ir == ik)
          r /= _ctwo;
        else
          {
            r = _ctwo/r;
            d = -d;
          }
      }
    return;
  }

  /*!
   * \brief Tools<TOutput, TMass, TScale>::qlZlogint
   * \param z
   * \param ieps
   * \return
   */
  template<typename TOutput, typename TMass, typename TScale>
  TOutput Tools<TOutput,TMass,TScale>::Zlogint(TOutput const& z, TScale const& ieps) const
  {
    const TOutput omz = _cone-z;
    return omz*(cLn(omz, ieps)-_cone)-(-z)*(cLn(-z, ieps)-_cone);
  }

  /*!
   * \param i_in
   * \param z_in
   * \param s
   * \return
   */
  template<typename TOutput, typename TMass, typename TScale>
  TOutput Tools<TOutput,TMass,TScale>::ltspence(int const& i_in, TOutput const& z_in, TScale const& s) const
  {
    TOutput z[2];
    z[i_in] = z_in;
    z[1-i_in] = _cone-z_in;

    TOutput ltspence;
    if (Real(z[0]) < _half)
      {
        if (Abs(z[0]) < _one)
          ltspence = ltli2series(z[1],s);
        else
          {
            const TOutput clnz = cLn(-z[0], -s);
            ltspence = TOutput(-_pi2o6) - _chalf*clnz*clnz - ltli2series(-z[1]/z[0], -s);
          }
      }
    else
      {
        const TScale az1 = Abs(z[1]);
        if (az1 < _eps15)
          ltspence = TOutput(_pi2o6);
        else if (az1 < _one)
          ltspence = TOutput(_pi2o6) - cLn(z[0],s)*cLn(z[1],-s) - ltli2series(z[0],-s);
        else
          {
            const TOutput clnz = cLn(-z[0], -s);
            ltspence = TOutput(_two*_pi2o6) + _chalf*clnz*clnz - cLn(z[0],s)*cLn(z[1],-s) + ltli2series(-z[0]/z[1],s);
          }
      }

    return ltspence;
  }

  /*!
   * \param z
   * \param isig
   * \return
   */
  template<typename TOutput, typename TMass, typename TScale>
  TOutput Tools<TOutput,TMass,TScale>::li2series(TOutput const& z, TScale const& isig) const
  {
    TOutput xm = -cLn(_cone-z, -isig);
    const TOutput x2 = xm*xm;
    TOutput res = xm - x2/_cfour;

    for (int j = 0; j < 25; j++)
      {
        xm *= x2;
        const TOutput n = res + xm*_B[j];
        if (n == res) return res;
        else res = n;
      }
    cout << "Tools::li2series: bad convergence" << endl;
    return _czero;
  }

  /*!
   * \param z1
   * \param s
   * \return
   */
  template<typename TOutput, typename TMass, typename TScale>
  TOutput Tools<TOutput,TMass,TScale>::ltli2series(TOutput const& z1, TScale const& s) const
  {
    TOutput xm = -cLn(z1, -s);
    const TOutput x2 = xm*xm;
    TOutput res = xm - x2/_cfour;

    for (int i = 0; i < 25; i++)
      {
        xm *= x2;
        const TOutput n = res + xm*_B[i];
        if (n == res) return res;
        else res = n;
      }
    cout << "Tools::ltli2series: bad convergence" << endl;
    return _czero;
  }

  /*!
   * \param a
   * \param b
   * \return
   */
  template<typename TOutput, typename TMass, typename TScale>
  TOutput Tools<TOutput,TMass,TScale>::eta2(TOutput const& a,TOutput const& b) const
  {
    const TScale ima = Imag(a);
    const TScale imb = Imag(b);
    const TScale imab = Imag(a*b);
    return _2ipi*(Htheta(-ima)*Htheta(-imb)*Htheta(imab)-Htheta(ima)*Htheta(imb)*Htheta(-imab));
  }

  /*!
   * \param a
   * \param b
   * \param c
   * \return
   */
  template<typename TOutput, typename TMass, typename TScale>
  TOutput Tools<TOutput,TMass,TScale>::eta3(TOutput const& a,TOutput const& b, TOutput const& c) const
  {
    TOutput res = _czero;
    const TScale ima = Sign(Imag(a));
    const TScale imb = Sign(Imag(b));
    const TScale imc = Sign(Imag(c));
    if (ima == imb && ima != imc) res = _2ipi*TOutput(imc);
    return res;
  }

  template<typename TOutput, typename TMass, typename TScale>
  TOutput Tools<TOutput,TMass,TScale>::eta5(TOutput const& a,TOutput const& b, TOutput const& c,TOutput const& d, TOutput const& e) const
  {
    TOutput res = _czero;
    const TScale ima = Sign(Imag(a));
    const TScale imb = Sign(Imag(b));
    const TScale imc = Sign(Imag(c));
    const TScale imd = Sign(Imag(d));
    const TScale ime = Sign(Imag(e));

    if (ima == imb)
      {
        if (ima == imd)
          {
            if (imc == ime)      res = _czero;
            else if (ima != imc) res = _2ipi*TOutput(imc);
            else                 res = _2ipi*TOutput(-ime);
          }
        else if (ima != imc) res = _2ipi*TOutput(imc);
        else res = _czero;
      }
    else if (ima == imd && ima != ime) res =  _2ipi*TOutput(-ime);
    else res = _czero;

    return res;
  }

  /*!
   * \param z2
   * \param im2
   * \return
   */
  template<typename TOutput, typename TMass, typename TScale>
  TOutput Tools<TOutput,TMass,TScale>::xetatilde(TOutput const (&z1)[2], TScale const (&im1)[2], TOutput const& z2, TScale const& im2, TOutput const (&l1)[2]) const
  {
    return l1[1]*TOutput(etatilde(z1[1], im1[1], z2, im2))
          -l1[0]*TOutput(etatilde(z1[0], im1[0], z2, im2));
  }

  /*!
   * \param z2
   * \param im2
   * \param im12
   * \return
   */
  template<typename TOutput, typename TMass, typename TScale>
  TOutput Tools<TOutput,TMass,TScale>::xeta(TOutput const (&z1)[2], TScale const (&im1)[2], TOutput const& z2, TScale const& im2, TScale const& im12, TOutput const (&l1)[2]) const
  {
    return l1[1]*TOutput(eta(z1[1], im1[1], z2, im2, im12))
          -l1[0]*TOutput(eta(z1[0], im1[0], z2, im2, im12));
  }

  /*!
   * \param z1
   * \param s1
   * \param z2
   * \param s2
   * \param s12
   * \return
   */
  template<typename TOutput, typename TMass, typename TScale>
  int Tools<TOutput,TMass,TScale>::eta(TOutput const& z1, TScale const& s1, TOutput const& z2, TScale const& s2, TScale const& s12) const
  {
    TScale im1 = Imag(z1), im2 = Imag(z2), im12 = Imag(z1*z2);
    if (im1 == _zero) im1 = s1;
    if (im2 == _zero) im2 = s2;
    if (im12 == _zero) im12 = s12;

    int eta;
    if (im1 < _zero && im2 < _zero && im12 > _zero)
      eta = 1;
    else if (im1 > _zero && im2 > _zero && im12 < _zero)
      eta = -1;
    else
      eta = 0;

    return eta;
  }

  /*!
   * \param c1
   * \param im1x
   * \param c2
   * \param im2x
   * \return
   */
  template<typename TOutput, typename TMass, typename TScale>
  int Tools<TOutput,TMass,TScale>::etatilde(TOutput const& c1, TScale const& im1x, TOutput const& c2, TScale const& im2x) const
  {
    int etatilde;
    TScale im1 = Imag(c1);
    TScale im2 = Imag(c2);

    if (im1 == _zero)
      im1 = im1x;

    if (im2 != _zero)
      etatilde = eta(c1, im1x, c2, _zero, _zero);
    else if ( Real(c2) > _zero)
      etatilde = 0;
    else if (im1 > _zero && im2x > _zero)
      etatilde = -1;
    else if (im1 < _zero && im2x < _zero)
      etatilde = 1;
    else
      etatilde = 0;

    return etatilde;
  }

  /*!
   * Calculate the K-function give in Eq. 2.7 of \cite Beenakker:1988jr
   * \f[
   *  K(p^2,m,m_p) = \frac{1-\sqrt{1-4m m_p / (z-(m-m_p)^2)}}{1+\sqrt{1-4m m_p / (z-(m-m_p)^2)}}
   * \f]
   * and fill x[0] = -K, x[1] = 1+K, x[2] = 1-K, the roots are allowed to be imaginary.
   * ieps gives the sign of the imaginary part of -K: 1 -> +i*eps
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Tools<TOutput,TMass,TScale>::kfn(TOutput (&res)[3], TScale& ieps, TMass const& xpi, TMass const& xm, TMass const& xmp) const
  {
    if (xm == TMass(_zero) || xmp == TMass(_zero))
      throw RangeError("Tools::ql","Error in kfn,xm,xmp");

    const TOutput xx1 = TOutput(xpi - (xm-xmp)*(xm-xmp));
    const TOutput rat = TOutput(xx1/(_four*xm*xmp));
    if (iszero(Real(rat)))
      {
        res[1] = -_ctwo*Sqrt(rat)*TOutput{_zero, _one} + _ctwo*rat;
        res[0] = _cone-res[1];
        res[2] = _ctwo-res[1];
      }
    else
      {
        const TOutput root = Sqrt((rat-_cone)/rat);
        const TOutput invopr = _cone/(_cone+root);
        res[0] = -invopr*invopr/rat;
        res[1] = _ctwo*invopr;
        res[2] = _ctwo*root*invopr;
      }
    ieps = _one;
  }

  /*!
   * Solution of a quadratic equation a*z^2+b*z+c=0.
   * \param a coefficient
   * \param b coefficient
   * \param c coefficient
   * \param l result [-im,+im]
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Tools<TOutput,TMass,TScale>::solveabc(TMass const& a, TMass const&b, TMass const& c, TOutput (&z)[2]) const
  {
    const TMass discr = b*b-TMass(_four)*a*c;

    if (iszero(a)) throw LogicException("Tools::solveabc","equation is not quadratic");

    if (iszero(Imag(discr)))
      {
        const TMass sgnb = Sign(Real(b));
        if (Real(discr) > 0)
            {
              const TMass q = -_half*(b+sgnb*Sqrt(discr));
              if (Real(b) > 0)
                {
                  z[0] = TOutput(c/q);
                  z[1] = TOutput(q/a);
                }
              else
                {
                  z[0] = TOutput(q/a);
                  z[1] = TOutput(c/q);
                }
            }
          else
            {
              z[1] = -(TOutput(b)+sgnb*Sqrt(TOutput(discr)))/(_ctwo*a);
              z[0] = Conjg(z[1]);
              if (Real(b) < 0)
                {
                  z[0] = z[1];
                  z[1] = Conjg(z[0]);
                }
            }
      }
    else
      {
        TOutput qq = -TOutput(b)+Sqrt(TOutput(discr));
        TOutput hh = -TOutput(b)-Sqrt(TOutput(discr));

        z[0] = qq*_chalf/TOutput(a);
        z[1] = (_ctwo*TOutput(c))/qq;

        if (Imag(z[0]) > _zero)
          {
            z[0] = hh*_chalf/TOutput(a);
            z[1] = (_ctwo*TOutput(c))/hh;
          }
      }
  }

  /*!
   * Solution of a quadratic equation a*z^2+b*z+c=0, with a give input discriminant.
   * \param a coefficient
   * \param b coefficient
   * \param c coefficient
   * \param d discriminant
   * \param l result [-im,+im]
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Tools<TOutput,TMass,TScale>::solveabcd(TOutput const& a, TOutput const&b, TOutput const& c, TOutput const& d, TOutput (&z)[2]) const
  {
    if (a == _czero)
      {
        if (b == _czero) throw LogicException("Tools::solveabcd","no possible solution");
        z[0] = -c/b; z[1] = z[0];
      }
    else if (c == _czero)
      {
        z[0] = d/a; z[1] = _czero;
      }
    else
      {
        const TOutput up = -b+d;
        const TOutput dn = -b-d;
        if (Abs(up) >= Abs(dn))
          {
            z[0] = _chalf*up/a;
            z[1] = _ctwo*c/up;
          }
        else
          {
            z[1] = _chalf*dn/a;
            z[0] = _ctwo*c/dn;
          }
      }
  }

  /*!
   * Solution of a quadratic equation a*z^2+b*z+c=0, with a give input discriminant.
   * \param a coefficient
   * \param b coefficient
   * \param c coefficient
   * \param d discriminant
   * \param l result [-im,+im]
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Tools<TOutput,TMass,TScale>::solveabcd(TOutput const& a, TOutput const&b, TOutput const& c, TOutput (&z)[2]) const
  {
    if (a == _czero)
      {
        if (b == _czero) throw LogicException("Tools::solveabcd","no possible solution");
        z[0] = -c/b; z[1] = z[0];
      }
    else if (c == _czero)
      {
        z[0] = _czero; z[1] = _czero;
      }
    else
      {
        const TOutput d = Sqrt(b*b - _four*a*c);
        const TOutput up = -b+d;
        const TOutput dn = -b-d;
        if (Abs(up) >= Abs(dn))
          {
            z[0] = _chalf*up/a;
            z[1] = _ctwo*c/up;
          }
        else
          {
            z[1] = _chalf*dn/a;
            z[0] = _ctwo*c/dn;
          }
      }
  }


  /*!
   * Calculate the function
   * \f[
   * R(p^2,m,m_p) = \frac{p_3^2+m_4^2-m_3^2+\sqrt{(p_3^2+m_4^2-m_3^2)^2-4 p_3^2 m_4^2}}{-p_3^2+m_4^2-m_3^2+\sqrt{(p_3^2+m_4^2-m_3^2)^2-4 p_3^2 m_4^2}}
   * \f]
   * where the roots are allowed to be imaginary.
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Tools<TOutput,TMass,TScale>::ratgam(TOutput &ratp, TOutput &ratm, TScale &ieps, TMass const& p3sq, TMass const& m3sq, TMass const& m4sq) const
  {
    const TOutput root = Sqrt(TOutput(Pow(p3sq-m3sq+m4sq,2) - _four*m4sq*p3sq));
    ratp = (TOutput(p3sq+m4sq-m3sq)+root)/(TOutput(-p3sq+m4sq-m3sq)+root);
    ratm = (TOutput(p3sq+m4sq-m3sq)-root)/(TOutput(-p3sq+m4sq-m3sq)-root);
    ieps = _zero;
  }

  /*!
   * Calculate the function
   * \f[
   * R = \frac{\sigma-i \epsilon}{\tau-i \epsilon}
   * \f]
   * where sigma and tau are real and ieps give the sign if i*pi.
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Tools<TOutput,TMass,TScale>::ratreal(TMass const& si, TMass const& ta, TMass &rat, TScale &ieps) const
  {
    rat = si/ta;
    if (Real(rat) > _zero)
      ieps = _zero;
    else if (Real(si) < _zero)
      ieps = -_one;
    else if (Real(ta) < _zero)
      ieps = _one;
    else if (Real(ta) == _zero)
      throw RangeError("Tools::ratreal", "error in ratreal");
  }

  template class Tools<complex,double,double>;
  template class Tools<complex,complex,double>;
  template class Tools<qcomplex,qdouble,qdouble>;
  template class Tools<qcomplex,qcomplex,qdouble>;

  Splash *Splash::_instance = nullptr;
  Splash::Splash()
  {
    cout << ql::blue << endl;
    cout << "      ____  __________  __                    "     << endl;
    cout << "     / __ \\/ ____/ __ \\/ /   ____  ____  ____ "   << endl;
    cout << "    / / / / /   / / / / /   / __ \\/ __ \\/ __ \\"  << endl;
    cout << "   / /_/ / /___/ /_/ / /___/ /_/ / /_/ / /_/ /"     << endl;
    cout << "   \\___\\_\\____/_____/_____/\\____/\\____/ .___/ "<< endl;
    cout << "                                     /_/      "       << endl;
    cout << "   ___git___: " << VERSION << " | __authors__: S.C., K.E., G.Z."<< ql::def << endl;
  }

}
