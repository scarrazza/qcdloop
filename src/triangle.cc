//
// QCDLoop 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch
//          Keith Ellis: keith.ellis@durham.ac.uk
//          Giulia Zanderighi: giulia.zanderighi@cern.ch

#include "qcdloop/triangle.h"
#include "qcdloop/tools.h"
#include "qcdloop/maths.h"
#include "qcdloop/exceptions.h"
#include <algorithm>
#include <iostream>
using std::sort;

namespace ql
{
  const int isort[6][6] = {
    {0,1,2,3,4,5},
    {1,2,0,4,5,3},
    {2,0,1,5,3,4},
    {1,0,2,3,5,4},
    {0,2,1,5,4,3},
    {2,1,0,4,3,5}
  };

  template<typename TOutput, typename TMass, typename TScale>
  Triangle<TOutput,TMass,TScale>::Triangle():
    Topology<TOutput,TMass,TScale>("Triangle")
  {
    this->_m.resize(3);
    this->_p.resize(3);
  }

  template<typename TOutput, typename TMass, typename TScale>
  Triangle<TOutput,TMass,TScale>::~Triangle()
  {
  }

  /*!
   * Sort arguments of triangle so that they are ordered in mass.
   * \param psq are the four-momentum squared of the external lines
   * \param msq are the squares of the masses of the internal lines
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Triangle<TOutput,TMass,TScale>::TriSort(TScale (&psq)[3], TMass (&msq)[3]) const
  {    
    const int x1[3] = {2,0,1};
    const int x2[3] = {1,2,0};
    TScale psqtmp[3];
    TMass  msqtmp[3];
    std::copy(std::begin(psq), std::end(psq), std::begin(psqtmp));
    std::copy(std::begin(msq), std::end(msq), std::begin(msqtmp));

    const TMass mmax = Max(msqtmp[0],Max(msqtmp[1],msqtmp[2]));
    if (mmax == msqtmp[0])
      for (int i = 0; i < 3; i++)
        {
          msq[x1[i]] = msqtmp[i];
          psq[x1[i]] = psqtmp[i];
        }
    else if (mmax == msqtmp[1])
      for (int i = 0; i < 3; i++)
        {
          msq[x2[i]] = msqtmp[i];
          psq[x2[i]] = psqtmp[i];
        }

    if (Abs(msq[0]) > Abs(msq[1]))
      {
        for (int i = 0; i < 2; i++)
          {
            msqtmp[i] = msq[i];
            psqtmp[i+1] = psq[i+1];
          }
        msq[0] = msqtmp[1];
        msq[1] = msqtmp[0];
        psq[1] = psqtmp[2];
        psq[2] = psqtmp[1];
      }      
  }

  /*!
   * Sort arguments of triangle so that |p3sq| > |p2sq| > |p1sq| and permute masses accordingly.
   * \param xpi i=0,2: mass^2, i=3,5: p^2
   * \param ypi abs(p3) < abs(p4) < abs(p5)
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Triangle<TOutput,TMass,TScale>::TriSort2(TMass const (&xpi)[6], TMass (&ypi)[6]) const
  {
    const TScale p1sq = Abs(xpi[3]);
    const TScale p2sq = Abs(xpi[4]);
    const TScale p3sq = Abs(xpi[5]);

    int j = 0;
    if      ( (p3sq >= p2sq) && (p2sq >= p1sq) ) j = 0;
    else if ( (p1sq >= p3sq) && (p3sq >= p2sq) ) j = 1;
    else if ( (p2sq >= p1sq) && (p1sq >= p3sq) ) j = 2;
    else if ( (p2sq >= p3sq) && (p3sq >= p1sq) ) j = 3;
    else if ( (p1sq >= p2sq) && (p2sq >= p3sq) ) j = 4;
    else if ( (p3sq >= p1sq) && (p1sq >= p2sq) ) j = 5;
    else j = 0;

    for (size_t k = 0; k < 6; k++)
      ypi[k] = xpi[isort[j][k]];
  }

  /*!
   * Sort an input vector of TScale objets based on its abs.
   * \param psq parameter to sort in ascending order
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Triangle<TOutput,TMass,TScale>::SnglSort(TScale (&psq)[3]) const
  {
    TScale absp[3] = { Abs(psq[0]), Abs(psq[1]), Abs(psq[2])};
    if (absp[0] > absp[1])
      {
        const TScale ptmp = psq[0], atmp = absp[0];
        psq[0] = psq[1]; absp[0] = absp[1];
        psq[1] = ptmp;   absp[1] = atmp;
      }

    if (absp[0] > absp[2])
      {
        const TScale ptmp = psq[0], atmp = absp[0];
        psq[0] = psq[2]; absp[0] = absp[2];
        psq[2] = ptmp;   absp[2] = atmp;
      }

    if (absp[1] > absp[2])
      {
        const TScale ptmp = psq[1];
        psq[1] = psq[2];
        psq[2] = ptmp;
      }
  }

  /*!
   * Computes the Triangle integral defined as:
   * \f[
   * I_{3}^{D}(p_1^2,p_2^2,p_3^2;m_1^2,m_2^2,m_3^2)= \frac{\mu^{4-D}}{i \pi^{D/2} r_{\Gamma}} \int d^Dl \frac{1}{(l^2-m_1^2+i \epsilon)((l+q_1)^2-m_2^2+i \epsilon)((l+q_2)^2-m_3^2+i\epsilon)}
   *   \f]
   *where \f$q_1=p_1,q_2=p_1+p_2\f$.
   *
   * Implementation of the formulae of Denner and Dittmaier \cite Denner:2005nn and
   * 't Hooft and Veltman \cite tHooft:1978xw.
   *
   * \param res output object res[0,1,2] the coefficients in the Laurent series
   * \param mu2 is the square of the scale mu
   * \param m are the squares of the masses of the internal lines
   * \param p are the four-momentum squared of the external lines
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Triangle<TOutput,TMass,TScale>::integral(vector<TOutput> &res,
                                                const TScale& mu2,
                                                vector<TMass> const& m,
                                                vector<TScale> const& p)
  {
    if (!this->checkCache(mu2,m,p))
      {
        if (mu2 < 0) throw RangeError("Triangle::integral","mu2 is negative!");

        // Normalization
        const TScale scalefac = Max(Abs(m[0]),Max(Abs(m[1]),Max(Abs(m[2]),Max(Abs(p[0]),Max(Abs(p[1]),Abs(p[2]))))));
        const TScale musq = mu2/scalefac;

        TMass msq[3];
        TScale psq[3];
        msq[0] = m[0]/scalefac;
        msq[1] = m[1]/scalefac;
        msq[2] = m[2]/scalefac;
        psq[0] = p[0]/scalefac;
        psq[1] = p[1]/scalefac;
        psq[2] = p[2]/scalefac;

        // Sort msq in ascending order
        TriSort(psq, msq);

        // if internal masses all 0, reorder abs(psq) in ascending order
        const bool iszeros[3] = {this->iszero(msq[0]),
                                 this->iszero(msq[1]),
                                 this->iszero(msq[2])};

        if (iszeros[0] && iszeros[1] && iszeros[2])
          SnglSort(psq);

        // calculate integral value
        const TMass Y01 = TMass(msq[0]+msq[1]-psq[0])/TMass(2);
        const TMass Y02 = TMass(msq[0]+msq[2]-psq[2])/TMass(2);
        const TMass Y12 = TMass(msq[1]+msq[2]-psq[1])/TMass(2);

        int massive = 0;
        for (size_t i = 0; i < 3; i++)
          if (!iszeros[i]) massive += 1;

        // building xpi
        const TMass xpi[6] = { msq[0], msq[1], msq[2], TMass(psq[0]), TMass(psq[1]), TMass(psq[2]) };

        if (massive == 3)      // three internal masses
          T0(this->_val, xpi, massive);
        else if (massive == 2) // two internal masses
          {
            if (this->iszero(Abs(Y01)) && this->iszero(Abs(Y02)))
              T6(this->_val, musq, msq[1], msq[2], psq[1]);
            else
              T0(this->_val, xpi, massive);
          }
        else if (massive == 1) // one internal masses
          {            
            if (!this->iszero(Abs(Y01)))
              T0(this->_val, xpi, massive);
            else if (this->iszero(Abs(Y02)) && this->iszero(Abs(Y12)))
              T5(this->_val, musq, msq[2]);
            else if (this->iszero(Abs(Y02)))
              T4(this->_val, musq, msq[2], psq[1]);
            else if (this->iszero(Abs(Y12)))
              T4(this->_val, musq, msq[2], psq[2]);
            else
              T3(this->_val, musq, msq[2], psq[1], psq[2]);
          }
        else // zero internal masses
          {          
            if (this->iszero(Abs(Y01)) && this->iszero(Abs(Y12)))
              T1(this->_val, musq, psq[2]);
            else if (this->iszero(Abs(Y01)))
              T2(this->_val, musq, psq[1], psq[2]);
            else
              T0(this->_val, xpi, massive);
          }

        this->_val[0] /= scalefac;
        this->_val[1] /= scalefac;
        this->_val[2] /= scalefac;

        this->storeCache(mu2,m,p);
      }

    if (res.size() != 3) { res.reserve(3); }
    std::copy(this->_val.begin(), this->_val.end(), res.begin());
    return;
  }

  /*!
   * Parses finite triangle integrals. Formulae from 't Hooft and Veltman \cite tHooft:1978xw
   * \param res output object res[0,1,2] the coefficients in the Laurent series
   * \param xpi an array with masses and momenta squared.
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Triangle<TOutput,TMass,TScale>::T0(vector<TOutput> &res,TMass const (&xpi)[6], int const& massive) const
  {    
    // Set poles to zero
    res[1] = res[2] = this->_czero;

    // Sort
    TMass ypi[6];
    TriSort2(xpi,ypi);

    const bool zypi3 = this->iszero(ypi[3]);
    const bool zypi4 = this->iszero(ypi[4]);

    // Trigger the finite topology
    if (zypi3 && zypi4 && this->iszero(ypi[5]))
      TIN0(res[0], ypi);
    else if (zypi3 && zypi4)
      TIN1(res[0], ypi, xpi, massive);
    else if (zypi3)
      TIN2(res[0], ypi, xpi, massive);
    else
      TIN3(res[0], ypi, xpi, massive);
  }

  /*!
   * Computes the finite triangle with all internal masses zero.
   * Formulae from 't Hooft and Veltman \cite tHooft:1978xw following the LoopTools implementation \cite Hahn:2006qw.
   * \param res the output object for the finite part.
   * \param xpi an array with masses and momenta squared
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Triangle<TOutput,TMass,TScale>::TIN0(TOutput &res, TMass const (&xpi)[6]) const
  {
    const TMass m1sq = xpi[0];
    const TMass m2sq = xpi[1];
    const TMass m3sq = xpi[2];

    if (this->iszero(m1sq-m2sq) && this->iszero(m2sq-m3sq))
      res = -this->_chalf/m1sq;
    else if (this->iszero(m1sq-m2sq))
      res = TOutput((m3sq*Log(m2sq/m3sq)+m3sq-m2sq)/Pow(m3sq-m2sq,2));
    else if (this->iszero(m2sq-m3sq))
      res = TOutput((m1sq*Log(m3sq/m1sq)+m1sq-m3sq)/Pow(m1sq-m3sq,2));
    else if (this->iszero(m3sq-m1sq))
      res = TOutput((m2sq*Log(m1sq/m2sq)+m2sq-m1sq)/Pow(m2sq-m1sq,2));
    else
      res = TOutput( m3sq*Log(m3sq/m1sq)/((m1sq-m3sq)*(m3sq-m2sq))
                    -m2sq*Log(m2sq/m1sq)/((m1sq-m2sq)*(m3sq-m2sq)));
  }

  /*!
   * Computes the finite triangle with all internal masses zero.
   * Formulae from 't Hooft and Veltman \cite tHooft:1978xw following the LoopTools and OneLoop implementation \cite Hahn:2006qw, \cite vanHameren:2010cp.
   * \param res the output object for the finite part.
   * \param xpi an array with masses and momenta squared
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Triangle<TOutput,TMass,TScale>::TIN1(TOutput &res, TMass const (&xpi)[6], TMass const (&sxpi)[6], int const& massive) const
  {    
    if (this->iszero(Imag(xpi[0])) && this->iszero(Imag(xpi[1])) && this->iszero(Imag(xpi[2])))
      {
        const TMass s1Ds1 = xpi[0];
        const TMass s2Ds2 = xpi[1];
        const TMass s3Ds3 = xpi[2];
        const TMass p3Dp3 = xpi[5];

        const TOutput x0 = TOutput((s1Ds1-s2Ds2)/p3Dp3);
        const TMass D23 = s2Ds2-s3Ds3;

        TOutput l[2];
        this->solveabc(p3Dp3,s3Ds3-s1Ds1-p3Dp3,s1Ds1,l);

        if (this->iszero(D23))
          res = TOutput(-(-this->Rint(x0,l[0],-1)-this->Rint(x0,l[1],1))/p3Dp3);
        else
          {
            const TOutput u0 = TOutput(s2Ds2/D23);
            const TScale ieps = Sign(Real(-D23));
            res = -(this->Rint(x0,u0,-ieps)-this->Rint(x0,l[0],-1)-this->Rint(x0,l[1],1))/p3Dp3;
          }
      }
    else
      {
        if (massive == 2)
          TINDNS2(res, sxpi);
        else if (massive == 1)
          TINDNS1(res, sxpi);
        else
          {
            if (Real(Kallen2(xpi[3],xpi[4],xpi[5])) < this->_zero) // never happens with real momenta (but just in case..)
              {
                const TOutput p2 = TOutput(xpi[5]);
                TOutput m[3] = {TOutput(xpi[0]),TOutput(xpi[1]),TOutput(xpi[2])};

                m[0] -= this->_ieps2*TOutput(Abs(Real(m[0])));
                m[1] -= this->_ieps2*TOutput(Abs(Real(m[1])));
                m[2] -= this->_ieps2*TOutput(Abs(Real(m[2])));

                const TOutput sm0 = Sqrt(m[0])-this->_ieps2;
                const TOutput sm2 = Sqrt(m[2])-this->_ieps2;

                const TOutput yy = -( (m[0]-m[1])-p2 )/p2;

                res = -(this->R3int(p2, sm0, sm2, yy) - this->R2int(m[1]-m[2], m[2], yy));
                res /= p2;
              }
            else
              TINDNS(res,xpi);
          }
      }
  }

  /*!
   * Computes the finite triangle with all internal masses zero.
   * Formulae from 't Hooft and Veltman \cite tHooft:1978xw following the LoopTools and OneLoop implementation \cite Hahn:2006qw, \cite vanHameren:2010cp.
   * \param res the output object for the finite part.
   * \param xpi an array with masses and momenta squared
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Triangle<TOutput,TMass,TScale>::TIN2(TOutput &res, TMass const (&xpi)[6], TMass const (&sxpi)[6], int const& massive) const
  {        
    if (this->iszero(Imag(xpi[0])) && this->iszero(Imag(xpi[1])) && this->iszero(Imag(xpi[2])))
      {
        const TMass s1Ds1 = xpi[0];
        const TMass s2Ds2 = xpi[1];
        const TMass s3Ds3 = xpi[2];
        const TMass p2Dp2 = xpi[4];
        const TMass p3Dp3 = xpi[5];

        const TOutput z0 = TOutput((s1Ds1-s2Ds2)/(p3Dp3-p2Dp2));
        TOutput zu[2], zl[2];
        this->solveabc(p2Dp2,s3Ds3-s2Ds2-p2Dp2,s2Ds2,zu);
        this->solveabc(p3Dp3,s3Ds3-s1Ds1-p3Dp3,s1Ds1,zl);

        if (this->iszero(p2Dp2-p3Dp3))
          {
            res = TOutput(-(this->Zlogint(zu[0],-1)+this->Zlogint(zu[1],1)
                           -this->Zlogint(zl[0],-1)+this->Zlogint(zl[1],1))/(s2Ds2-s1Ds1));
          }
        else
          {
            res = TOutput(-(this->Rint(z0,zu[0],-1)+this->Rint(z0,zu[1],1)
                           -this->Rint(z0,zl[0],-1)-this->Rint(z0,zl[1],1))/(p3Dp3-p2Dp2));
          }
      }
    else
      {
        if (massive == 2)
          TINDNS2(res, sxpi);
        else if (massive == 1)
          TINDNS1(res, sxpi);
        else
          {
            if (Real(Kallen2(xpi[3],xpi[4],xpi[5])) < this->_zero || xpi[4] != xpi[5])
              {
                const TOutput p[2] = {TOutput(xpi[4]),TOutput(xpi[5])};
                TOutput m[3] = {TOutput(xpi[0]), TOutput(xpi[1]), TOutput(xpi[2])};

                if (p[0] == p[1]) throw LogicException("Triangle::TIN2", "threshold singularity");

                m[0] -= this->_ieps2*TOutput(Abs(Real(m[0])));
                m[1] -= this->_ieps2*TOutput(Abs(Real(m[1])));
                m[2] -= this->_ieps2*TOutput(Abs(Real(m[2])));

                const TOutput sm0 = Sqrt(m[0])-this->_ieps2;
                const TOutput sm1 = Sqrt(m[1])-this->_ieps2;
                const TOutput sm2 = Sqrt(m[2])-this->_ieps2;

                const TOutput yy = ((m[0]-m[1])-p[1]+p[0])/(p[0]-p[1]);

                res = this->R3int(p[1],sm0,sm2,yy)-this->R3int(p[0],sm1,sm2,yy);
                res /= (p[0]-p[1]);
              }
            else
              TINDNS(res,xpi);
          }
      }
  }

  /*!
   * Computes the finite triangle with all internal masses zero.
   * Formulae from 't Hooft and Veltman \cite tHooft:1978xw following the LoopTools and OneLoop implementation \cite Hahn:2006qw, \cite vanHameren:2010cp.
   * \param res the output object for the finite part.
   * \param xpi an array with masses and momenta squared
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Triangle<TOutput,TMass,TScale>::TIN3(TOutput &res, TMass const (&xpi)[6], TMass const (&sxpi)[6], int const& massive) const
  {        
    if (this->iszero(Imag(xpi[0])) && this->iszero(Imag(xpi[1])) && this->iszero(Imag(xpi[2])))
      {
        TOutput Del2[3], y[3], z[2];
        TMass siDsi[3], piDpj[3][3], siDpj[3][3], kdel[3];
        int jp1[3] = {1,2,0}, jm1[3] = {2,0,1};

        for (int j = 0; j < 3; j++)
          {
            siDsi[j] = xpi[j];
            piDpj[j][j] = xpi[j+3];
            siDpj[j][j] = this->_half*(xpi[jp1[j]]-xpi[j]-piDpj[j][j]);
          }

        for (int j = 0; j < 3; j++)
          {
            piDpj[j][jp1[j]] = this->_half*(piDpj[jm1[j]][jm1[j]]-piDpj[j][j]-piDpj[jp1[j]][jp1[j]]);
            piDpj[jp1[j]][j] = piDpj[j][jp1[j]];
            siDpj[j][jm1[j]] = piDpj[jm1[j]][jm1[j]] + siDpj[jm1[j]][jm1[j]];
            siDpj[j][jp1[j]] =-piDpj[j][jp1[j]]+siDpj[jp1[j]][jp1[j]];
          }

        for (int j = 0; j < 3; j++)
          {
            Del2[j] = TOutput(piDpj[j][jp1[j]]*piDpj[j][jp1[j]] - piDpj[j][j]*piDpj[jp1[j]][jp1[j]]);
            Del2[j] = Sqrt(Del2[j]);
            kdel[j] = piDpj[j][j]*siDpj[jp1[j]][jp1[j]]-piDpj[j][jp1[j]]*siDpj[jp1[j]][j];
            y[j] = TOutput((siDpj[jp1[j]][j] + kdel[j]/Del2[j])/piDpj[j][j]);
          }

        res = this->_czero;
        for (int j = 0; j < 3; j++)
          {
            const TMass a = piDpj[j][j], b = -this->_two*siDpj[jp1[j]][j], c = siDsi[jp1[j]];
            this->solveabc(a,b,c,z);
            res += this->Rint(y[j],z[0],-1) + this->Rint(y[j],z[1],1);
          }

        res = -res/(this->_ctwo*Del2[0]);
      }
    else
      {
        if (massive == 2)
          TINDNS2(res, sxpi);
        else if (massive == 1)
          TINDNS1(res, sxpi);
        else
          {
            const TOutput K2 = Kallen2(xpi[3],xpi[4],xpi[5]);
            if (Real(K2) < this->_zero)
              {
                const TOutput p[3] = {TOutput(xpi[3]),TOutput(xpi[4]),TOutput(xpi[5])};
                TOutput m[3] = {TOutput(xpi[0]), TOutput(xpi[1]),TOutput(xpi[2])};

                m[0] -= this->_ieps2*TOutput(Abs(Real(m[0])));
                m[1] -= this->_ieps2*TOutput(Abs(Real(m[1])));
                m[2] -= this->_ieps2*TOutput(Abs(Real(m[2])));

                const TOutput alpha = Sqrt(K2)+this->_ieps2;
                const TOutput sm0 = Sqrt(m[0])-this->_ieps2;
                const TOutput sm1 = Sqrt(m[1])-this->_ieps2;
                const TOutput sm2 = Sqrt(m[2])-this->_ieps2;

                res = -(this->R3int(p[0], sm0, sm1,   m[1]-m[2]+p[1], p[2]-p[0]-p[1], p[1], alpha)
                       -this->R3int(p[2], sm0, sm2, -(m[0]-m[1])+p[2]-p[1], p[1]-p[0]-p[2], p[0], alpha)
                       +this->R3int(p[1], sm1, sm2, -(m[0]-m[1])+p[2]-p[1], p[0]+p[1]-p[2], p[0], alpha));
                res /= alpha;
              }
            else
              TINDNS(res,xpi);
          }
      }
  }  

  /*!
   * Computes the finite triangle, when Kallen2 > 0 and 3 massive particles.
   * Formulae from Denner, Nierste and Scharf \cite Denner:1991qq following the OneLoop implementation \cite vanHameren:2010cp.
   * \param res the output object for the finite part.
   * \param xpi an array with masses and momenta squared
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Triangle<TOutput,TMass,TScale>::TINDNS(TOutput &res, TMass const (&xpi)[6]) const
  {
    TOutput m1 = xpi[0];
    TOutput m2 = xpi[1];
    TOutput m3 = xpi[2];
    TOutput p1 = TOutput(xpi[3]);
    TOutput p2 = TOutput(xpi[4]);
    TOutput p3 = TOutput(xpi[5]);

    const TOutput sm1 = Sqrt(m1);
    const TOutput sm2 = Sqrt(m2);
    const TOutput sm3 = Sqrt(m3);

    TOutput k12 = this->_czero, k13 = this->_czero, k23 = this->_czero;
    if (m1+m2 != p1) k12 = (m1+m2-p1-p1*this->_ieps2)/(sm1*sm2);
    if (m1+m3 != p3) k13 = (m1+m3-p3-p3*this->_ieps2)/(sm1*sm3);
    if (m2+m3 != p2) k23 = (m2+m3-p2-p2*this->_ieps2)/(sm2*sm3);

    TOutput r12, r13, r23, d12, d13, d23;
    this->R(r12, d12, k12);
    this->R(r13, d13, k13);
    this->R(r23, d23, k23);

    const TOutput a = sm2/sm3 - k23 + r13*(k12 - sm2/sm1);
    if (a == this->_czero)
      {
        std::cout << "Triangle::TINDNS: threshold singularity, return 0" << std::endl;
        res = this->_czero;
        return;
      }
    const TOutput b = d13/sm2 + k12/sm3 - k23/sm1;
    const TOutput c = (sm1/sm3 - this->_one/r13)/(sm1*sm2);

    TOutput x[2];
    this->solveabcd(a,b,c,x);
    x[0] = -x[0];
    x[1] = -x[1];

    const TOutput qz[2]  = { x[0]*sm2, x[1]*sm2 };
    const TOutput qz2[2] = { x[0]*r13*sm2, x[1]*r13*sm2};
    const TOutput qz3[2] = { x[0]*r13, x[1]*r13};
    const TOutput oneOr12 = this->_cone/r12, oneOr23 = this->_cone/r23;
    const TOutput dqz = qz[0]-qz[1], dqz2 = qz2[0]-qz2[1], dqz3 = qz3[0]-qz3[1];
    const TScale siqz[2] = { (TScale) Sign(Imag(qz[0])), (TScale) Sign(Imag(qz[1]))};
    const TScale siqz2[2]= { (TScale) Sign(Imag(qz2[0])), (TScale) Sign(Imag(qz2[1]))};
    const TScale siqz3[2]= { (TScale) Sign(Imag(qz3[0])), (TScale) Sign(Imag(qz3[1]))};
    const TScale sigx[2]= { (TScale) Sign(Imag(x[0])), (TScale) Sign(Imag(x[1]))};

    res = (-this->xspence(qz,siqz,r12,Sign(Imag(r12)))/dqz
           -this->xspence(qz,siqz,oneOr12,Sign(Imag(oneOr12)))/dqz)*sm2
          +(this->xspence(qz2,siqz2,r23,Sign(Imag(r23)))/dqz2
          + this->xspence(qz2,siqz2,oneOr23,Sign(Imag(oneOr23)))/dqz2)*r13*sm2
          - this->xspence(qz3,siqz3,sm3,Sign(Imag(sm3)))/dqz3*r13
          + this->xspence(x,sigx,sm1,Sign(Imag(sm1)))/(x[0]-x[1]);

    if (x[1] != this->_czero)
      {
        const TOutput arg1 = qz3[0]/qz3[1];
        const TOutput arg2 = qz3[0]*qz3[1]/(sm3*sm3);
        const TOutput arg3 = x[0]/x[1];
        const TOutput arg4 = x[0]*x[1]/(sm1*sm1);

        TOutput log1 = this->cLn(arg2,Sign(Imag(arg2)));
        TOutput log2 = this->cLn(arg4,Sign(Imag(arg4)));
        if (Real(arg2) < this->_zero && Imag(arg2) < this->_zero) log1 += this->_2ipi;
        if (Real(arg4) < this->_zero && Imag(arg4) < this->_zero) log2 += this->_2ipi;

        res += (this->cLn(arg1,Sign(Imag(arg1)))/(this->_cone-arg1)*log1
              - this->cLn(arg3,Sign(Imag(arg3)))/(this->_cone-arg3)*log2)/(this->_ctwo*x[1]);
      }
    res /= (a*sm1*sm2*sm3);
  }

  /*!
   * Computes the finite triangle, with 2 massive particles
   * Formulae from Denner, Nierste and Scharf \cite Denner:1991qq following the OneLoop implementation \cite vanHameren:2010cp.
   * \param res the output object for the finite part.
   * \param xpi an array with masses and momenta squared
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Triangle<TOutput,TMass,TScale>::TINDNS2(TOutput &res, TMass const (&xpi)[6]) const
  {
    TOutput m2 = xpi[1];
    TOutput m4 = xpi[2];
    TOutput p2 = TOutput(xpi[3]);
    TOutput p3 = TOutput(xpi[5]);
    TOutput p23 = TOutput(xpi[4]);

    const TOutput sm2 = Sqrt(m2);
    const TOutput sm3 = Abs(sm2);
    const TOutput sm4 = Sqrt(m4);

    TOutput r23 = this->_czero, k24 = this->_czero, r34 = this->_czero;
    r23 = (m2-p2-p2*this->_ieps2)/(sm2*sm3);
    k24 = (m2+m4-p23-p23*this->_ieps2)/(sm2*sm4);
    r34 = (m4-p3-p3*this->_ieps2)/(sm3*sm4);

    TOutput r24, d24;
    this->R(r24, d24, k24);

    const TOutput a = r34/r24 - r23;
    if (a == this->_czero)
      {
        std::cout << "Triangle::TINDNS2: threshold singularity, return 0" << std::endl;
        res = this->_czero;
        return;
      }

    const TOutput b = -d24/sm3 + r34/sm2 - r23/sm4;
    const TOutput c = (sm4/sm2 - r24)/(sm3*sm4);

    TOutput x[2];
    this->solveabcd(a,b,c,x);
    x[0] = -x[0];
    x[1] = -x[1];

    const TOutput qz[2]  = { x[0]/r24, x[1]/r24};
    const TScale siqz[2] = { (TScale) Sign(Imag(qz[0])), (TScale) Sign(Imag(qz[1]))};

    res = -this->xspence(qz,siqz,sm2,Sign(Imag(sm2)))/(qz[0]-qz[1])/r24;

    if (x[1] != this->_czero)
      {
        const TOutput arg1 = qz[0]/qz[1];
        const TOutput arg2 = qz[0]*qz[1]/(sm2*sm2);
        const TOutput arg3 = x[0]/x[1];
        const TOutput arg4 = x[0]*x[1]/(sm4*sm4);

        TOutput log1 = this->cLn(arg2,Sign(Imag(arg2)));
        TOutput log2 = this->cLn(arg4,Sign(Imag(arg4)));
        if (Real(arg2) < this->_zero && Imag(arg2) < this->_zero) log1 += this->_2ipi;
        if (Real(arg4) < this->_zero && Imag(arg4) < this->_zero) log2 += this->_2ipi;

        res += (this->cLn(arg1,Sign(Imag(arg1)))/(this->_cone-arg1)*log1
              - this->cLn(arg3,Sign(Imag(arg3)))/(this->_cone-arg3)*log2)/(this->_ctwo*x[1]);
      }

    const TScale siqx[2] = { (TScale) Sign(Imag(x[0])), (TScale) Sign(Imag(x[1])) };
    res += this->xspence(x,siqx,sm4,Sign(Imag(sm4)))/(x[0]-x[1]);

    if (!this->iszero(Abs(r23)))
      {
        const TOutput arg = r23*sm3/r24;
        res += this->xspence(x,siqx,arg, Sign(Imag(arg)))/(x[0]-x[1]);
      }

    if (!this->iszero(Abs(r34)))
      {
        const TOutput arg = r34*sm3;
        res -= this->xspence(x,siqx,arg, Sign(Imag(arg)))/(x[0]-x[1]);
      }

    res /= (a*sm2*sm3*sm4);
  }

  /*!
   * Computes the finite triangle, with 1 massive particles.
   * Formulae from Denner, Nierste and Scharf \cite Denner:1991qq following the OneLoop implementation \cite vanHameren:2010cp.
   * \param res the output object for the finite part.
   * \param xpi an array with masses and momenta squared
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Triangle<TOutput,TMass,TScale>::TINDNS1(TOutput &res, TMass const (&xpi)[6]) const
  {
    TOutput m4 = xpi[2];
    TOutput p2 = TOutput(xpi[3]);
    TOutput p3 = TOutput(xpi[4]);
    TOutput p4 = TOutput(xpi[5]);
    TOutput p23 = p4;

    const TOutput sm4 = Sqrt(m4);
    const TOutput sm3 = Abs(sm4);
    const TOutput sm2 = sm3;

    TOutput r23 = this->_czero, r24 = this->_czero, r34 = this->_czero;
    r23 = (-p2-p2*this->_ieps2)/(sm2*sm3);
    r24 = (m4-p23-p23*this->_ieps2)/(sm2*sm4);
    r34 = (m4-p3-p3*this->_ieps2)/(sm3*sm4);

    const TOutput a = r34*r24 - r23;
    if (a == this->_czero)
      {
        std::cout << "Triangle::TINDNS1: threshold singularity, return 0" << std::endl;
        res = this->_czero;
        return;
      }

    const TOutput b = r24/sm3 + r34/sm2 - r23/sm4;
    const TOutput c = this->_cone/(sm2*sm3);

    TOutput x[2];
    this->solveabcd(a,b,c,x);
    x[0] = -x[0];
    x[1] = -x[1];

    const TScale siqx[2] = { (TScale) Sign(Imag(x[0])), (TScale) Sign(Imag(x[1])) };

    const TOutput arg1 = x[0]/x[1];
    const TOutput arg2 = x[0]*x[1]/(sm4*sm4);
    const TOutput arg3 = r23/(sm3*sm3);

    const TOutput log0 = this->cLn(arg1,Sign(Imag(arg1)))/(this->_cone-arg1);

    TOutput log1 = this->cLn(arg2,Sign(Imag(arg2)));
    TOutput log2 = this->cLn(arg3,Sign(Imag(arg3)));
    if (Real(arg2) < this->_zero && Imag(arg2) < this->_zero) log1 += this->_2ipi;
    if (Real(arg3) < this->_zero && Imag(arg3) < this->_zero) log2 += this->_2ipi;

    res = this->xspence(x,siqx,sm4,Sign(Imag(sm4)))/(x[0]-x[1])
         -log0*log1/(this->_ctwo*x[1])
         -log0*log2/x[1];

    if (!this->iszero(Abs(r24)))
      {
        const TOutput arg = r24*sm3;
        res += this->xspence(x,siqx,arg, Sign(Imag(arg)))/(x[0]-x[1]);
      }

    if (!this->iszero(Abs(r34)))
      {
        const TOutput arg = r34*sm3;
        res -= this->xspence(x,siqx,arg, Sign(Imag(arg)))/(x[0]-x[1]);
      }

    res /= (a*sm2*sm3*sm4);
  }

  /*!
   * Computes the Kallen function defined as:
   * \f[
   *  K(p_1,p_2,p_3) = \sqrt{p_1^2+p_2^2+p_3^2-2 (p_1 \cdot p_2 + p_2 \cdot p_3 +p_3 \cdot p_1)}
   * \f]
   * \param p1 four-momentum squared
   * \param p2 four-momentum squared
   * \param p3 four-momentum squared
   */
  template<typename TOutput, typename TMass, typename TScale>
  TOutput Triangle<TOutput,TMass,TScale>::Kallen(TOutput const& p1, TOutput const& p2, TOutput const& p3) const
  {
    return Sqrt(Kallen2(p1,p2,p3));
  }

  /*!
   * Computes the Kallen function defined as:
   * \f[
   *  K(p_1,p_2,p_3) = p_1^2+p_2^2+p_3^2-2 (p_1 \cdot p_2 + p_2 \cdot p_3 +p_3 \cdot p_1)
   * \f]
   * \param p1 four-momentum squared
   * \param p2 four-momentum squared
   * \param p3 four-momentum squared
   */
  template<typename TOutput, typename TMass, typename TScale>
  TOutput Triangle<TOutput,TMass,TScale>::Kallen2(TOutput const& p1, TOutput const& p2, TOutput const& p3) const
  {
    return TOutput(p1*p1+p2*p2+p3*p3-this->_ctwo*(p1*p2+p2*p3+p3*p1));
  }

  /*!
   * The integral is defined as:
   * \f[
   * I_{3}^{D=4-2 \epsilon}(0,0,p^2;0,0,0)= \frac{1}{p^2} \left( \frac{1}{\epsilon^2} + \frac{1}{\epsilon} \ln \left( \frac{\mu^2}{-p^2-i \epsilon} \right) + \frac{1}{2} \ln^2 \left( \frac{\mu^2}{-p^2-i \epsilon} \right) \right) + O(\epsilon)
   *   \f]
   * Implementation of the formulae from Beenakker et al. \cite Beenakker:2002nc.
   *
   * \param res output object res[0,1,2] the coefficients in the Laurent series
   * \param mu2 is the squre of the scale mu
   * \param p is the four-momentum squared of the external line
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Triangle<TOutput,TMass,TScale>::T1(vector<TOutput> &res, TScale const& mu2, TScale const& p) const
  {
    const TOutput wlogm = this->Lnrat(mu2, -p);
    res[2] = this->_cone/TOutput(p);
    res[1] = res[2]*wlogm;
    res[0] = this->_chalf*res[2]*wlogm*wlogm;
  }

  /*!
   * The integral is defined as:
   * \f[
   * I_{3}^{D=4-2 \epsilon}(0,p_1^2,p_2^2;0,0,0)= \frac{1}{p_1^2-p_2^2} \left\{ \frac{1}{\epsilon} \left[ \ln \left( \frac{\mu^2}{-p_1^2-i \epsilon} \right) - \ln \left( \frac{\mu^2}{-p_2^2-i \epsilon} \right) \right] + \frac{1}{2} \left[ \ln^2 \left( \frac{\mu^2}{-p_1^2-i \epsilon} \right) - \ln^2 \left( \frac{\mu^2}{-p_2^2-i \epsilon} \right) \right] \right\} + O(\epsilon)
   *   \f]
   * Implementation of the formulae from Beenakker et al. \cite Beenakker:2002nc.
   *
   * \param res output object res[0,1,2] the coefficients in the Laurent series
   * \param mu2 is the squre of the scale mu
   * \param p1 is the four-momentum squared of the external line
   * \param p2 is the four-momentum squared of the external line
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Triangle<TOutput,TMass,TScale>::T2(vector<TOutput> &res, TScale const& mu2, TScale const& p1, TScale const& p2) const
  {
    const TOutput wlog1 = this->Lnrat(mu2, -p1);
    const TOutput wlog2 = this->Lnrat(mu2, -p2);
    const TScale r = (p2-p1)/p1;
    res[2] = this->_czero;
    if (Abs(r) < this->_eps)
      {
        const TOutput ro2 = r/this->_ctwo;
        res[1] = -this->_cone/p1*(this->_cone-ro2);
        res[0] = res[1]*wlog1 + ro2/p1;
      }
    else
      {
        res[1] = (wlog1-wlog2)/TOutput(p1-p2);
        res[0] = this->_chalf*res[1]*(wlog1+wlog2);
      }
  }

  /*!
   * The integral is defined as:
   * \f[
   * I_{3}^{D=4-2 \epsilon}(0,p_1^2,p_2^2;0,0,m^2)= \frac{1}{p_1^2-p_2^2} \left( \frac{\mu^2}{m^2} \right)^\epsilon \left\{ \frac{1}{\epsilon} \ln \left( \frac{m^2-p_2^2}{m^2-p_1^2} \right) + {\rm Li}_2 \left( \frac{p_1^2}{m^2} \right) - {\rm Li}_2 \left( \frac{p_2^2}{m^2} \right) + \ln^2 \left( \frac{m^2-p_1^2}{m^2} \right) - \ln^2 \left( \frac{m^2-p_2^2}{m^2} \right) \right\} + O(\epsilon)
   *   \f]
   * Implementation of the formulae from Beenakker et al. \cite Beenakker:2002nc.
   *
   * \param res output object res[0,1,2] the coefficients in the Laurent series
   * \param mu2 is the squre of the scale mu
   * \param m is the square of the mass of the internal line
   * \param p1 is the four-momentum squared of the external line
   * \param p2 is the four-momentum squared of the external line
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Triangle<TOutput,TMass,TScale>::T3(vector<TOutput> &res, const TScale &mu2, const TMass &m, const TScale &p2, const TScale &p3) const
  {
    const TMass m2sqb = m-TMass(p2);
    const TMass m3sqb = m-TMass(p3);
    const TOutput dilog2 = this->Li2omrat(m2sqb, m);
    const TOutput dilog3 = this->Li2omrat(m3sqb, m);

    const TOutput wlog2  = this->Lnrat(m2sqb,m);
    const TOutput wlog3  = this->Lnrat(m3sqb,m);
    const TOutput wlogm  = this->Lnrat(mu2,m);
    const TMass r = (m3sqb-m2sqb)/m2sqb;

    res[2] = this->_czero;    
    if (Abs(r) < this->_eps)
      {        
        res[1] = (this->_cone-this->_chalf*r)/m2sqb;
        res[0] = (wlogm - (m+p2)/p2*wlog2);
        res[0] += -this->_chalf*(r*((m*m - this->_ctwo*p2*m-p2*p2)*wlog2 + p2*(m+p2+p2*wlogm)) / (p2*p2));
        res[0] /= m2sqb;        
      }
    else
      {        
        const TOutput fac = this->_cone/(p2-p3);
        res[1] = fac*(wlog3-wlog2);
        res[0] = res[1]*wlogm + fac*(dilog2-dilog3+(wlog2*wlog2-wlog3*wlog3));
      }
  }

  /*!
   * The integral is defined as:
   * \f[
   * I_{3}^{D=4-2 \epsilon}(0,m^2,p_2^2;0,0,m^2)= \left( \frac{\mu^2}{m^2} \right)^\epsilon \frac{1}{p_2^2-m^2} \left[ \frac{1}{2 \epsilon^2} + \frac{1}{\epsilon} \ln \left( \frac{m^2}{m^2-p_2^2} \right) + \frac{\pi^2}{12} + \frac{1}{2} \ln^2 \left( \frac{m^2}{m^2-p_2^2} \right) - {\rm Li}_2 \left( \frac{-p_2^2}{m^2-p_2^2} \right) \right] + O(\epsilon)
   *   \f]
   * Implementation of the formulae from Beenakker et al. \cite Beenakker:2002nc.
   *
   * \param res output object res[0,1,2] the coefficients in the Laurent series
   * \param mu2 is the squre of the scale mu
   * \param m is the square of the mass of the internal line
   * \param p is the four-momentum squared of the external line
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Triangle<TOutput,TMass,TScale>::T4(vector<TOutput> &res, const TScale &mu2, const TMass &m, const TScale &p2) const
  {
    const TOutput wlog  = this->Lnrat(m,m-p2);
    const TOutput wlogm = this->Lnrat(mu2, m);
    const TOutput fac   = this->_chalf/(p2-m);
    const TMass arg2    = -p2/(m-p2);
    const TMass omarg2  = TMass(this->_one)-arg2;
    const TOutput ct = TOutput(this->_pi2o6);

    TOutput dilog2;
    if (Abs(omarg2) < this->_zero)
      dilog2 = ct-TOutput(this->ddilog(omarg2))-Log(arg2)*wlog;
    else
      dilog2 = TOutput(this->ddilog(arg2));

    res[2] = fac;
    res[1] = res[2]*wlogm+fac*this->_ctwo*wlog;
    res[0] = -res[2]*this->_chalf*wlogm*wlogm+res[1]*wlogm + fac*(wlog*wlog+ct-this->_ctwo*dilog2);
  }

  /*!
   * The integral is defined as:
   * \f[
   * I_{3}^{D=4-2 \epsilon}(0,m^2,m^2;0,0,m^2)= \left( \frac{\mu^2}{m^2} \right)^\epsilon \frac{1}{m^2} \left( -\frac{1}{2 \epsilon} + 1 \right) + O(\epsilon)
   *   \f]
   * Implementation of the formulae from Beenakker et al. \cite Beenakker:2002nc.
   *
   * \param res output object res[0,1,2] the coefficients in the Laurent series
   * \param mu2 is the squre of the scale mu
   * \param m is the square of the mass of the internal line
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Triangle<TOutput,TMass,TScale>::T5(vector<TOutput> &res, const TScale &mu2, const TMass &m) const
  {
    const TOutput fac = this->_cone/m;
    const TOutput wlogm = this->Lnrat(mu2, m);
    res[2] = this->_czero;
    res[1] = -this->_chalf*fac;
    res[0] = res[1]*wlogm+fac;
  }

  /*!
   * The integral is defined as:
   * \f[
   * I_{3}^{D=4-2 \epsilon}(m_2^2,s,m_3^2;0,m_2^2,m_3^2)= \frac{\Gamma(1+\epsilon)\mu^\epsilon}{2\epsilon r_\Gamma} \int_0^1 d\gamma \frac{1}{\left[ \gamma m_2^2 + (1-\gamma) m_3^2 - \gamma (1-\gamma)s - i\epsilon \right]^{1+\epsilon}} + O(\epsilon)
   *   \f]
   * Implementation of the formulae from Beenakker et al. \cite Beenakker:2002nc.
   *
   * \param res output object res[0,1,2] the coefficients in the Laurent series
   * \param mu2 is the squre of the scale mu
   * \param m2 is the square of the mass of the internal line
   * \param m3 is the square of the mass of the internal line
   * \param p2 is the four-momentum squared of the external line
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Triangle<TOutput,TMass,TScale>::T6(vector<TOutput> &res, TScale const& mu2, TMass const& m2sq, TMass const& m3sq, TScale const& p2) const
  {
    const TMass m2 = Sqrt(m2sq);
    const TMass m3 = Sqrt(m3sq);

    TScale iepsd = 0;
    TOutput cxs[3];
    this->kfn(cxs, iepsd, p2, m2, m3);

    const TOutput xlog= this->cLn(cxs[0], iepsd);
    const TScale resx = Real(cxs[0]);
    const TScale imxs = Imag(cxs[0]);

    if (this->iszero(resx-this->_one) && this->iszero(imxs))
      {
        const TMass arg = mu2/(m2*m3);
        const TOutput fac = this->_chalf/(m2*m3);

        res[1] = fac;
        if (this->iszero(m2-m3))
          res[0] = fac*Log(arg);
        else
          res[0] = fac*(Log(arg)-this->_ctwo-(m3+m2)/(m3-m2)*Log(m2/m3));
      }
    else
      {
        const TMass arg = m2/m3;
        const TMass arg2= m2*m3;
        const TOutput logarg = TOutput(Log(arg));
        const TOutput fac = this->_cone/arg2*cxs[0]/(cxs[1]*cxs[2]);
        res[1] = -fac*xlog;
        res[0] = fac*(xlog*(-this->_chalf*xlog + Log(arg2/mu2))
                      - this->cLi2omx2(cxs[0], cxs[0], iepsd, iepsd)
                      + this->_chalf*logarg*logarg
                      + this->cLi2omx2(cxs[0], arg, iepsd, this->_zero)
                      + this->cLi2omx2(cxs[0], this->_cone/arg, iepsd, this->_zero));
      }
    res[2] = this->_czero;
  }

  // explicity tyoename declaration
  template class Triangle<complex,double,double>;
  template class Triangle<complex,complex,double>;
  template class Triangle<qcomplex,qdouble,qdouble>;
  template class Triangle<qcomplex,qcomplex,qdouble>;
}
