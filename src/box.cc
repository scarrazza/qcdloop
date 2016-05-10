//
// QCDLoop 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch
//          Keith Ellis: keith.ellis@durham.ac.uk
//          Giulia Zanderighi: giulia.zanderighi@cern.ch

#include "qcdloop/box.h"
#include "qcdloop/tools.h"
#include "qcdloop/maths.h"
#include "qcdloop/exceptions.h"
#include <iostream>

namespace ql
{
  // Constants for sort/swap
  const int swap_b0m[13][4] = { {3, 2, 1, 0},
                            {0, 3, 2, 1},
                            {1, 0, 3, 2},
                            {2, 1, 0, 3},
                            {7, 6, 5, 4},
                            {4, 7, 6, 5},
                            {5, 4, 7, 6},
                            {6, 5, 4, 7},
                            {9, 8, 9, 8},
                            {8, 9, 8, 9},
                            {10, 10, 10, 10},
                            {12, 11, 12, 11},
                            {11, 12, 11, 12},
                          };

  const int jsort_b0m[4] = {3,0,1,2};

  const int swap_b1m[13][5] = { {3, 2, 1, 0, 2},
                            {0, 3, 2, 1, 1},
                            {1, 0, 3, 2, 0},
                            {2, 1, 0, 3, 3},
                            {7, 6, 5, 4, 5},
                            {4, 7, 6, 5, 4},
                            {5, 4, 7, 6, 7},
                            {6, 5, 4, 7, 6},
                            {9, 8, 9, 8, 8},
                            {8, 9, 8, 9, 9},
                            {10, 10, 10, 10, 10},
                            {12, 11, 12, 11, 11},
                            {11, 12, 11, 12, 12},
                          };

  const int swap_b2m[13][5] = { {3, 2, 1, 0, 2},
                            {0, 3, 2, 1, 1},
                            {1, 0, 3, 2, 0},
                            {2, 1, 0, 3, 3},
                            {7, 6, 5, 4, 5},
                            {4, 7, 6, 5, 4},
                            {5, 4, 7, 6, 7},
                            {6, 5, 4, 7, 6},
                            {9, 8, 9, 8, 8},
                            {8, 9, 8, 9, 9},
                            {10, 10, 10, 10, 10},
                            {11, 11, 11, 11, 11},
                            {12, 12, 12, 12, 12},
                          };

  const int swap_b3m[13][4] = { {0, 3, 2, 1},
                            {1, 0, 3, 2},
                            {2, 1, 0, 3},
                            {3, 2, 1, 0},
                            {4, 7, 6, 5},
                            {5, 4, 7, 6},
                            {6, 5, 4, 7},
                            {7, 6, 5, 4},
                            {8, 9, 8, 9},
                            {9, 8, 9, 8},
                            {10, 10, 10, 10},
                            {11, 12, 11, 12},
                            {12, 11, 12, 11},
                          };


  template<typename TOutput, typename TMass, typename TScale>
  Box<TOutput,TMass,TScale>::Box():
    Topology<TOutput,TMass,TScale>("Box")
  {
    this->_m.resize(4);
    this->_p.resize(6);
  }

  template<typename TOutput, typename TMass, typename TScale>
  Box<TOutput,TMass,TScale>::~Box()
  {
  }

  template<typename TOutput, typename TMass, typename TScale>
  void Box<TOutput,TMass,TScale>::Ycalc(TMass (&Y)[4][4],TMass (&Yalt)[4][4], int const& massive, bool const& opposite) const
  {
    if (massive == 1)
      {
        //C---exchange (1<-->3)
        Yalt[0][0] = Y[2][2];
        Yalt[1][1] = Y[1][1];
        Yalt[2][2] = Y[0][0];
        Yalt[3][3] = Y[3][3];
        Yalt[0][1] = Yalt[1][0] = Y[1][2];
        Yalt[0][2] = Yalt[2][0] = Y[0][2];
        Yalt[0][3] = Yalt[3][0] = Y[2][3];
        Yalt[1][2] = Yalt[2][1] = Y[0][1];
        Yalt[1][3] = Yalt[3][1] = Y[1][3];
        Yalt[2][3] = Yalt[3][2] = Y[0][3];
      }
    else if (massive == 2 && opposite)
      {
        //C---exchange (2<-->4) .and (1<-->3)
        Yalt[0][0] = Y[2][2];
        Yalt[1][1] = Y[3][3];
        Yalt[2][2] = Y[0][0];
        Yalt[3][3] = Y[1][1];
        Yalt[0][1] = Yalt[1][0] = Y[2][3];
        Yalt[0][2] = Yalt[2][0] = Y[0][2];
        Yalt[0][3] = Yalt[3][0] = Y[1][2];
        Yalt[1][2] = Yalt[2][1] = Y[0][3];
        Yalt[1][3] = Yalt[3][1] = Y[1][3];
        Yalt[2][3] = Yalt[3][2] = Y[0][1];
      }
    else if (massive == 2 && !opposite)
      {
        //C---exchange (1<-->2)and(3<-->4)
        Yalt[0][0] = Y[1][1];
        Yalt[1][1] = Y[0][0];
        Yalt[2][2] = Y[3][3];
        Yalt[3][3] = Y[2][2];
        Yalt[0][1] = Yalt[1][0] = Y[0][1];
        Yalt[2][3] = Yalt[3][2] = Y[2][3];
        Yalt[0][2] = Yalt[2][0] = Y[1][3];
        Yalt[0][3] = Yalt[3][0] = Y[1][2];
        Yalt[1][2] = Yalt[2][1] = Y[0][3];
        Yalt[1][3] = Yalt[3][1] = Y[0][2];
      }
    else
      throw RangeError("Box::Ycalc","massive value not implemented");
  }

  /*!
   * Computes the Box integral defined as:
   * \f[
   * I_{4}^{D}(p_1^2,p_2^2,p_3^2,p_4^2;s_{12},s_{23};m_1^2,m_2^2,m_3^2,m_4^2)= \frac{\mu^{4-D}}{i \pi^{D/2} r_{\Gamma}} \int d^Dl \frac{1}{(l^2-m_1^2+i \epsilon)((l+q_1)^2-m_2^2+i \epsilon)((l+q_2)^2-m_3^2+i\epsilon)((l+q_4)^2-m_4^2+i\epsilon)}
   *   \f]
   * where \f$ q_1=p_1, q_2=p_1+p_2, q_3=p_1+p_2+p_3\f$ and \f$q_0=q_4=0\f$.
   *
   * Implementation of the formulae of Denner et al. \cite Denner:1991qq,
   * 't Hooft and Veltman \cite tHooft:1978xw, Bern et al. \cite Bern:1993kr.
   *
   * \param res output object res[0,1,2] the coefficients in the Laurent series
   * \param mu2 is the square of the scale mu
   * \param m are the squares of the masses of the internal lines
   * \param p are the four-momentum squared of the external lines
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Box<TOutput,TMass,TScale>::integral(vector<TOutput> &res,
                                           const TScale& mu2,
                                           vector<TMass> const& m,
                                           vector<TScale> const& p)
  {
    if (!this->checkCache(mu2,m,p))
      {
        if (mu2 < 0) throw RangeError("Box::integral","mu2 is negative!");

        // Normalization
        const TScale scalefac = Max(Abs(p[4]),Max(Abs(p[5]),Max(Abs(p[0]),Max(Abs(p[1]),Max(Abs(p[2]),Abs(p[3]))))));

        TMass xpi[13];
        xpi[0] = m[0]/scalefac;
        xpi[1] = m[1]/scalefac;
        xpi[2] = m[2]/scalefac;
        xpi[3] = m[3]/scalefac;
        xpi[4] = TMass(p[0]/scalefac);
        xpi[5] = TMass(p[1]/scalefac);
        xpi[6] = TMass(p[2]/scalefac);
        xpi[7] = TMass(p[3]/scalefac);
        xpi[8] = TMass(p[4]/scalefac);
        xpi[9] = TMass(p[5]/scalefac);
        xpi[10] = xpi[4]+xpi[5]+xpi[6]+xpi[7]-xpi[8]-xpi[9];
        xpi[11] =-xpi[4]+xpi[5]-xpi[6]+xpi[7]+xpi[8]+xpi[9];
        xpi[12] = xpi[4]-xpi[5]+xpi[6]-xpi[7]+xpi[8]+xpi[9];
        const TScale musq = mu2/scalefac;

        // Count number of internal masses
        int massive = 0;
        for (size_t i = 0; i < 4; i++)
          if (!this->iszero(Abs(xpi[i]))) massive += 1;

        // check cayley elements
        const TMass y13 = xpi[0] + xpi[2] - xpi[8];
        const TMass y24 = xpi[1] + xpi[3] - xpi[9];
        if (this->iszero(y13) || this->iszero(y24))
          {
            std::cout << "Box::integral: Modified Cayley elements y13 or y24=0" << std::endl;
            res[0] = res[1] = res[2] = this->_czero;
            return;
          }

        if (massive == 0)
          B0m(this->_val, xpi, musq);
        else if (massive == 1)
          B1m(this->_val, xpi, musq);        
        else if (massive == 2)
          B2m(this->_val, xpi, musq);
        else if (massive == 3)
          B3m(this->_val, xpi, musq);
        else if (massive == 4)
          B4m(this->_val, xpi);          

        this->_val[0] /= (scalefac*scalefac);
        this->_val[1] /= (scalefac*scalefac);
        this->_val[2] /= (scalefac*scalefac);

        this->storeCache(mu2,m,p);
      }

    if (res.size() != 3) { res.reserve(3); }
    res = this->_val;
    return;
  }

  /*!
   * This function trigger the topologies with 4-offshell external lines.
   * \param res output object res[0,1,2] the coefficients in the Laurent series
   * \param xpi array with masses and momenta squared
   * \param mu2 is the square of the scale mu
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Box<TOutput,TMass,TScale>::B0m(vector<TOutput> &res, const TMass (&xpi)[13], const TScale &mu2) const
  {    
    int offshell = 0, jsort0 = 0, jsort1 = 0, jsort2 = 0;
    for (int j = 0; j < 4; j++)
      {
        if (!this->iszero(xpi[j+4]))
          {
            offshell += 1;
            if (jsort1 == 0) jsort1 = j+1;
            else jsort2 = j+1;
          }
        else
          jsort0 = j;
      }    

    TMass xpiout[13];
    bool swapped = true;
    const int jdiff = jsort2 - jsort1;
    if ( offshell == 1 )
      {
        for (size_t j = 0; j < 13; j++)
          xpiout[swap_b0m[j][jsort1-1]] = xpi[j];
      }
    else if (offshell == 2 && (jdiff == 2 || jdiff == 1))
      {
        for (size_t j = 0; j < 13; j++)
          xpiout[swap_b0m[j][jsort2-1]] = xpi[j];
      }
    else if (offshell == 2 && jdiff == 3)
      {
        for (size_t j = 0; j < 13; j++)
          xpiout[swap_b0m[j][0]] = xpi[j];
      }
    else if (offshell == 3)
    {
      for (size_t j = 0; j < 13; j++)
        xpiout[swap_b0m[j][jsort_b0m[jsort0]]] = xpi[j];
    }
    else
      swapped = false;

    TMass Y[4][4];
    TMass xpi_in[13];
    if (swapped)
      std::copy(std::begin(xpiout), std::end(xpiout), std::begin(xpi_in));
    else
      std::copy(std::begin(xpi), std::end(xpi), std::begin(xpi_in));

    Y[0][0] = xpi_in[0];
    Y[1][1] = xpi_in[1];
    Y[2][2] = xpi_in[2];
    Y[3][3] = xpi_in[3];
    Y[0][1] = Y[1][0] = this->_half*(xpi_in[0]+xpi_in[1]-xpi_in[4]);
    Y[0][2] = Y[2][0] = this->_half*(xpi_in[0]+xpi_in[2]-xpi_in[8]);
    Y[0][3] = Y[3][0] = this->_half*(xpi_in[0]+xpi_in[3]-xpi_in[7]);
    Y[1][2] = Y[2][1] = this->_half*(xpi_in[1]+xpi_in[2]-xpi_in[5]);
    Y[1][3] = Y[3][1] = this->_half*(xpi_in[1]+xpi_in[3]-xpi_in[9]);
    Y[2][3] = Y[3][2] = this->_half*(xpi_in[2]+xpi_in[3]-xpi_in[6]);

    if (offshell == 4)
      BIN0(res, Y);
    else if (offshell == 3)
      B5(res, Y, mu2);
    else if (offshell == 2)
      {
        if (!this->iszero(xpi_in[5]))
          B3(res, Y, mu2);
        else if (!this->iszero(xpi_in[6]) && !this->iszero(xpi_in[7]))
          B4(res, Y, mu2);
      }
    else if (offshell == 1)
      B2(res, Y, mu2);
    else if (offshell == 0)
      B1(res, Y, mu2);
  }

  /*!
   * This function trigger the topologies with 3-offshell external lines.
   * \param res output object res[0,1,2] the coefficients in the Laurent series
   * \param xpi array with masses and momenta squared
   * \param mu2 is the square of the scale mu
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Box<TOutput,TMass,TScale>::B1m(vector<TOutput> &res, const TMass (&xpi)[13], const TScale &mu2) const
  {    
    int jsort = 0;
    for (int i = 0; i < 4; i++)
      if (!this->iszero(xpi[i]))
        jsort = i;

    TMass xpi_in[13];
    for (size_t i = 0; i < 13; i++)
      xpi_in[swap_b1m[i][jsort]] = xpi[i];

    TMass Y[4][4];
    Y[0][0] = xpi_in[0];
    Y[1][1] = xpi_in[1];
    Y[2][2] = xpi_in[2];
    Y[3][3] = xpi_in[3];
    Y[0][1] = Y[1][0] = this->_half*(xpi_in[0]+xpi_in[1]-xpi_in[4]);
    Y[0][2] = Y[2][0] = this->_half*(xpi_in[0]+xpi_in[2]-xpi_in[8]);
    Y[0][3] = Y[3][0] = this->_half*(xpi_in[0]+xpi_in[3]-xpi_in[7]);
    Y[1][2] = Y[2][1] = this->_half*(xpi_in[1]+xpi_in[2]-xpi_in[5]);
    Y[1][3] = Y[3][1] = this->_half*(xpi_in[1]+xpi_in[3]-xpi_in[9]);
    Y[2][3] = Y[3][2] = this->_half*(xpi_in[2]+xpi_in[3]-xpi_in[6]);

    if ( !this->iszero(Y[0][0]) || !this->iszero(Y[1][1]) || !this->iszero(Y[2][2]) )
      throw LogicException("Box::B1m","Wrong ordering.");

    const bool zY01 = this->iszero(Y[0][1]);
    const bool zY12 = this->iszero(Y[1][2]);
    const bool zY23 = this->iszero(Y[2][3]);
    const bool zY30 = this->iszero(Y[3][0]);

    if (zY01 && zY12 && zY23 && zY30)
      B6(res, Y, mu2);
    else if (zY01 && zY12 && zY23)
      B7(res, Y, mu2);
    else if (zY01 && zY12 && zY30)
      {
        TMass Yalt[4][4];
        Ycalc(Y,Yalt,1);
        B7(res, Yalt, mu2);
      }
    else if (zY01 && zY12)
      B8(res, Y, mu2);
    else if (zY01 && zY30)
      B9(res, Y, mu2);
    else if (zY12 && zY23)
      {
        TMass Yalt[4][4];
        Ycalc(Y,Yalt,1);
        B9(res, Yalt, mu2);
      }
    else if (zY01)
      B10(res, Y, mu2);
    else if (zY12)
      {
        TMass Yalt[4][4];
        Ycalc(Y,Yalt,1);
        B10(res, Yalt, mu2);
      }
    else
      BIN1(res, Y);      
  }

  /*!
   * This function trigger the topologies with 2-offshell external lines.
   * \param res output object res[0,1,2] the coefficients in the Laurent series
   * \param xpi array with masses and momenta squared
   * \param mu2 is the square of the scale mu
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Box<TOutput,TMass,TScale>::B2m(vector<TOutput> &res, const TMass (&xpi)[13], const TScale &mu2) const
  {    
    int jsort1 = -1, jsort2 = -1;
    for (int i = 0; i < 4; i++)
      if (!this->iszero(xpi[i]))
        {
          if (jsort1 == -1)
            jsort1 = i;
          else
            jsort2 = i;
        }

    int jdiff = jsort2-jsort1;
    TMass xpiout[13];

    if (jdiff == 1 || jdiff == 2)
      {
        for (size_t j = 0; j < 13; j++)
          xpiout[swap_b2m[j][jsort2]] = xpi[j];
      }
    else if (jdiff == 3)
      {
        for (size_t j = 0; j < 13; j++)
          xpiout[swap_b2m[j][0]] = xpi[j];
      }

    if (jdiff == 2)
      B2mo(res, xpiout, mu2);
    else
      B2ma(res, xpiout, mu2);      
  }

  /*!
   * This function trigger the topologies with 2-offshell external lines.
   * \param res output object res[0,1,2] the coefficients in the Laurent series
   * \param xpi array with masses and momenta squared
   * \param mu2 is the square of the scale mu
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Box<TOutput,TMass,TScale>::B2ma(vector<TOutput> &res, const TMass (&xpi)[13], const TScale &mu2) const
  {    
    TMass Y[4][4], Yalt[4][4];
    Y[0][0] = xpi[0];
    Y[1][1] = xpi[1];
    Y[2][2] = xpi[2];
    Y[3][3] = xpi[3];
    Y[0][1] = Y[1][0] = this->_half*(xpi[0]+xpi[1]-xpi[4]);
    Y[0][2] = Y[2][0] = this->_half*(xpi[0]+xpi[2]-xpi[8]);
    Y[0][3] = Y[3][0] = this->_half*(xpi[0]+xpi[3]-xpi[7]);
    Y[1][2] = Y[2][1] = this->_half*(xpi[1]+xpi[2]-xpi[5]);
    Y[1][3] = Y[3][1] = this->_half*(xpi[1]+xpi[3]-xpi[9]);
    Y[2][3] = Y[3][2] = this->_half*(xpi[2]+xpi[3]-xpi[6]);

    Ycalc(Y,Yalt,2,this->iszero(xpi[2]));

    const bool zY01 = this->iszero(Y[0][1]);
    const bool zY12 = this->iszero(Y[1][2]);
    const bool zY03 = this->iszero(Y[0][3]);

    if (zY01 && zY12 && zY03)
      B11(res, Y, mu2);
    else if (zY01 && zY12 && !zY03)
      B12(res, Y, mu2);
    else if (this->iszero(Yalt[0][1]) && this->iszero(Yalt[1][2]) && !this->iszero(Yalt[0][3]))
      B12(res, Yalt, mu2);
    else if (zY01 && !zY12 && !zY03)
      B13(res, Y, mu2);
    else
      BIN2(res, Y);      
  }

  /*!
   * This function trigger the topologies with 2-offshell external lines.
   * \param res output object res[0,1,2] the coefficients in the Laurent series
   * \param xpi array with masses and momenta squared
   * \param mu2 is the square of the scale mu
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Box<TOutput,TMass,TScale>::B2mo(vector<TOutput> &res, const TMass (&xpi)[13], const TScale &mu2) const
  {    
    //if (!this->iszero(xpi[0]) || !this->iszero(xpi[2]))
    //  throw RangeError("Box::B2mo","Error non zero m1sq/m3sq");

    TMass Y[4][4];
    Y[0][0] = xpi[0];
    Y[1][1] = xpi[1];
    Y[2][2] = xpi[2];
    Y[3][3] = xpi[3];
    Y[0][1] = Y[1][0] = this->_half*(xpi[0]+xpi[1]-xpi[4]);
    Y[0][2] = Y[2][0] = this->_half*(xpi[0]+xpi[2]-xpi[8]);
    Y[0][3] = Y[3][0] = this->_half*(xpi[0]+xpi[3]-xpi[7]);
    Y[1][2] = Y[2][1] = this->_half*(xpi[1]+xpi[2]-xpi[5]);
    Y[1][3] = Y[3][1] = this->_half*(xpi[1]+xpi[3]-xpi[9]);
    Y[2][3] = Y[3][2] = this->_half*(xpi[2]+xpi[3]-xpi[6]);

    const bool zY00 = this->iszero(Y[0][0]);
    const bool zY22 = this->iszero(Y[2][2]);
    const bool zY01 = this->iszero(Y[0][1]);
    const bool zY03 = this->iszero(Y[0][3]);
    const bool zY12 = this->iszero(Y[1][2]);
    const bool zY23 = this->iszero(Y[1][2]);

    if (zY00 && zY22 && zY01 && zY03 && zY12 && zY23)
      B14(res, Y, mu2);
    else if (zY00 && zY22 && zY01 && zY03)
      B15(res, Y, mu2);
    else if (zY00 && zY22 && zY12 && zY23)
      {
        TMass Yalt[4][4];
        Ycalc(Y,Yalt,2,this->iszero(xpi[2]));
        B15(res, Yalt, mu2);
      }
    else
      {
        TMass Yadj[4][4];
        const int swap23[4] = {0,2,1,3};
        for (int i = 0; i < 4; i++)
          for (int j = 0; j < 4; j++)
            Yadj[i][j] = Y[swap23[i]][swap23[j]];
        BIN2(res, Yadj);
      }      
  }

  /*!
   * This function trigger the topologies with 1-offshell external lines.
   * \param res output object res[0,1,2] the coefficients in the Laurent series
   * \param xpi array with masses and momenta squared
   * \param mu2 is the square of the scale mu
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Box<TOutput,TMass,TScale>::B3m(vector<TOutput> &res, const TMass (&xpi)[13], const TScale &mu2) const
  {    
    int jsort = 0;
    for (int i = 0; i < 4; i++)
      if (this->iszero(xpi[i])) jsort = i;

    TMass xpo[13];
    for (size_t i = 0; i < 13; i++)
      xpo[swap_b3m[i][jsort]] = xpi[i];

    TMass Y[4][4];
    Y[0][0] = xpo[0];
    Y[1][1] = xpo[1];
    Y[2][2] = xpo[2];
    Y[3][3] = xpo[3];
    Y[0][1] = Y[1][0] = this->_half*(xpo[0]+xpo[1]-xpo[4]);
    Y[0][2] = Y[2][0] = this->_half*(xpo[0]+xpo[2]-xpo[8]);
    Y[0][3] = Y[3][0] = this->_half*(xpo[0]+xpo[3]-xpo[7]);
    Y[1][2] = Y[2][1] = this->_half*(xpo[1]+xpo[2]-xpo[5]);
    Y[1][3] = Y[3][1] = this->_half*(xpo[1]+xpo[3]-xpo[9]);
    Y[2][3] = Y[3][2] = this->_half*(xpo[2]+xpo[3]-xpo[6]);

    if (this->iszero(Y[0][0]) && this->iszero(Y[0][1]) && this->iszero(Y[0][3]))
      B16(res, Y, mu2);
    else
      BIN3(res, Y);      
  }

  /*!
   * This function trigger the topologies with 0-offshell external lines.
   * \param res output object res[0,1,2] the coefficients in the Laurent series
   * \param xpi array with masses and momenta squared
   * \param mu2 is the square of the scale mu
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Box<TOutput,TMass,TScale>::B4m(vector<TOutput> &res, const TMass (&xpi)[13]) const
  {
    TMass Y[4][4];
    Y[0][0] = xpi[0];
    Y[1][1] = xpi[1];
    Y[2][2] = xpi[2];
    Y[3][3] = xpi[3];
    Y[0][1] = Y[1][0] = this->_half*(xpi[0]+xpi[1]-xpi[4]);
    Y[0][2] = Y[2][0] = this->_half*(xpi[0]+xpi[2]-xpi[8]);
    Y[0][3] = Y[3][0] = this->_half*(xpi[0]+xpi[3]-xpi[7]);
    Y[1][2] = Y[2][1] = this->_half*(xpi[1]+xpi[2]-xpi[5]);
    Y[1][3] = Y[3][1] = this->_half*(xpi[1]+xpi[3]-xpi[9]);
    Y[2][3] = Y[3][2] = this->_half*(xpi[2]+xpi[3]-xpi[6]);

    BIN4(res, Y);    
  }

  /*!
   * Finite box with zero masses. Formulae from \cite Denner:1991qq.
   * \param res output object res[0,1,2] the coefficients in the Laurent series, following the LoopTools implementation \cite Hahn:2006qw.
   * \param Y the modified Cayley matrix.
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Box<TOutput,TMass,TScale>::BIN0(vector<TOutput> &res, TMass const (&Y)[4][4]) const
  {
    const TMass m2 = Abs(Y[1][3]);
    const TMass k12 = this->_two*Y[0][1]/m2;
    const TMass k13 = this->_two*Y[0][2]/m2;
    const TMass k14 = this->_two*Y[0][3]/m2;
    const TMass k23 = this->_two*Y[1][2]/m2;
    const TMass k24 = this->_two*Y[1][3]/m2;
    const TMass k34 = this->_two*Y[2][3]/m2;

    const TOutput k12c = TOutput(k12-Max(Abs(k12),this->_one)*this->_ieps50);
    const TOutput k13c = TOutput(k13-Max(Abs(k13),this->_one)*this->_ieps50);
    const TOutput k14c = TOutput(k14-Max(Abs(k14),this->_one)*this->_ieps50);
    const TOutput k23c = TOutput(k23-Max(Abs(k23),this->_one)*this->_ieps50);
    const TOutput k24c = TOutput(k24-Max(Abs(k24),this->_one)*this->_ieps50)/k12c;
    const TOutput k34c = TOutput(k34-Max(Abs(k34),this->_one)*this->_ieps50)/k13c;

    const TMass a = k34*k24;
    const TMass b = k13*k24 + k12*k34 - k14*k23;
    const TOutput c = TOutput(k13*k12) + TOutput(k23)*this->_ieps50;
    const TOutput disc = Sqrt(TOutput(b*b) - TOutput(this->_four*a)*c);

    TOutput x4[2];
    x4[0] = this->_chalf*(TOutput(b) - disc)/TOutput(a);
    x4[1] = this->_chalf*(TOutput(b) + disc)/TOutput(a);
    if (Abs(x4[0]) > Abs(x4[1]))
      x4[1] = c/(TOutput(a)*x4[0]);
    else
      x4[0] = c/(TOutput(a)*x4[1]);

    const TScale imzero[2] = {this->_zero, this->_zero};
    const TOutput log_x40 = Log(x4[0]);
    const TOutput log_x41 = Log(x4[1]);

    res[0] = ( (log_x41-log_x40)*(-this->_chalf*(log_x41+log_x40)
              + Log(k12c) + Log(k13c) - Log(k23c) - Log(k14c))
              -this->xspence(x4,imzero,k34c,0) - this->xspence(x4,imzero,k24c,0)
              )/(m2*m2*disc);
    res[1] = res[2] = this->_czero;
  }

  /*!
   * Finite box with 1 non-zero mass. Formulae from \cite Denner:1991qq.
   * \param res output object res[0,1,2] the coefficients in the Laurent series, following the LoopTools implementation \cite Hahn:2006qw.
   * \param Y the modified Cayley matrix.
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Box<TOutput,TMass,TScale>::BIN1(vector<TOutput> &res, TMass const (&Y)[4][4]) const
  {
    const TMass m4 = Y[3][3];
    const TMass k34 = this->_two*Y[0][1]/m4;
    const TMass k23 = this->_two*Y[0][2]/m4;
    const TMass k13 = this->_two*Y[0][3]/m4;
    const TMass k24 = this->_two*Y[1][2]/m4;
    const TMass k14 = this->_two*Y[1][3]/m4;
    const TMass k12 = this->_two*Y[2][3]/m4;

    const TOutput k12c = TOutput(k12-Max(Abs(k12),this->_one)*this->_ieps50);
    const TOutput k13c = TOutput(k13-Max(Abs(k13),this->_one)*this->_ieps50);
    const TOutput k14c = TOutput(k14-Max(Abs(k14),this->_one)*this->_ieps50);
    const TOutput k23c = TOutput(k23-Max(Abs(k23),this->_one)*this->_ieps50);
    const TOutput k24c = TOutput(k24-Max(Abs(k24),this->_one)*this->_ieps50)/k12c;
    const TOutput k34c = TOutput(k34-Max(Abs(k34),this->_one)*this->_ieps50)/k13c;

    const TMass a = k34*k24;
    const TMass b = k13*k24 + k12*k34 - k14*k23;
    const TOutput c = TOutput(k13*k12) - TOutput(k23)*(this->_cone-this->_ieps50);
    const TOutput disc = Sqrt(TOutput(b*b) - this->_cfour*TOutput(a)*c);

    TOutput x4[2];
    x4[0] = this->_chalf*(TOutput(b) - disc)/TOutput(a);
    x4[1] = this->_chalf*(TOutput(b) + disc)/TOutput(a);
    if (Abs(x4[0]) > Abs(x4[1]))
      x4[1] = c/(TOutput(a)*x4[0]);
    else
      x4[0] = c/(TOutput(a)*x4[1]);

    const TScale imzero[2] = {this->_zero, this->_zero};
    res[0] = ( this->xspence(x4, imzero, k14c, this->_zero)
              -this->xspence(x4, imzero, k34c, this->_zero)
              -this->xspence(x4, imzero, k24c, this->_zero)
              +(Log(x4[1])-Log(x4[0]))*(Log(k12c)+Log(k13c)-Log(k23c))
              )/(TOutput(Pow(m4,2))*disc);
    res[1] = res[2] = this->_czero;
  }

  /*!
   * Finite box with 2 non-zero masses. Formulae from \cite Denner:1991qq.
   * \param res output object res[0,1,2] the coefficients in the Laurent series, following the LoopTools implementation \cite Hahn:2006qw.
   * \param Y the modified Cayley matrix.
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Box<TOutput,TMass,TScale>::BIN2(vector<TOutput> &res, TMass const (&Y)[4][4]) const
  {
    const TMass m3 = Y[2][2];
    const TMass m4 = Y[3][3];
    const TMass m = Sqrt(m3*m4);

    const TMass k12 = this->_two*Y[1][2]/m3;
    const TMass k13 = this->_two*Y[0][2]/m3;
    const TMass k14 = this->_two*Y[2][3]/m;
    const TMass k23 = this->_two*Y[0][1]/m3;
    const TMass k24 = this->_two*Y[1][3]/m;
    const TMass k34 = this->_two*Y[0][3]/m;

    const TOutput k12c = TOutput(k12-Max(Abs(k12),this->_one)*this->_ieps50);
    const TOutput k13c = TOutput(k13-Max(Abs(k13),this->_one)*this->_ieps50);
    const TOutput k23c = TOutput(k23-Max(Abs(k23),this->_one)*this->_ieps50);
    const TOutput k24c = TOutput(k24-Max(Abs(k24),this->_one)*this->_ieps50)/k12c;
    const TOutput k34c = TOutput(k34-Max(Abs(k34),this->_one)*this->_ieps50)/k13c;

    TOutput r14 =  this->_chalf*(TOutput(k14)+TOutput(Sign(Real(k14)))*Sqrt(TOutput((k14-this->_two)*(k14+this->_two))));
    r14 *= (this->_cone + this->_ieps50*TOutput(Sign(Real(this->_cone/r14 - r14))));

    const TMass a = k34*k24 - k23;
    const TMass b = k13*k24 + k12*k34 - k14*k23;
    const TOutput c = TOutput(k13*k12) - TOutput(k23)*(this->_cone-this->_ieps50);
    const TOutput disc = Sqrt(TOutput(b*b) - TOutput(this->_four*a)*c);

    TOutput x4[2];
    x4[0] = this->_chalf*(TOutput(b) - disc)/TOutput(a);
    x4[1] = this->_chalf*(TOutput(b) + disc)/TOutput(a);
    if (Abs(x4[0]) > Abs(x4[1]))
      x4[1] = c/(TOutput(a)*x4[0]);
    else
      x4[0] = c/(TOutput(a)*x4[1]);

    const TScale imzero[2] = {this->_zero, this->_zero};
    res[0] = ( this->xspence(x4, imzero, r14, this->_zero)
              +this->xspence(x4, imzero, this->_cone/r14, this->_zero)
              -this->xspence(x4, imzero, k34c, this->_zero)
              -this->xspence(x4, imzero, k24c, this->_zero)
              +(Log(x4[1]) - Log(x4[0]))*(Log(k12c) + Log(k13c) - Log(k23c))
             )/(m3*m*disc);
    res[2] = res[1] = this->_czero;
  }

  /*!
   * Finite box with 3 non-zero masses. Formulae from \cite Denner:1991qq.
   * \param res output object res[0,1,2] the coefficients in the Laurent series, following the LoopTools implementation \cite Hahn:2006qw.
   * \param Y the modified Cayley matrix.
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Box<TOutput,TMass,TScale>::BIN3(vector<TOutput> &res, TMass const (&Y)[4][4]) const
  {
    const TMass m2 = Y[1][1];
    const TMass m3 = Y[2][2];
    const TMass m4 = Y[3][3];
    const TMass m_0 = Sqrt(m3*m4);
    const TMass m_1 = Sqrt(m2*m3);
    const TMass m_2 = Sqrt(m2*m4);

    const TMass k12 = this->_two*Y[2][3]/m_0;
    const TMass k13 = this->_two*Y[0][2]/m3;
    const TMass k14 = this->_two*Y[1][2]/m_1;
    const TMass k23 = this->_two*Y[0][3]/m_0;
    const TMass k24 = this->_two*Y[1][3]/m_2;
    const TMass k34 = this->_two*Y[0][1]/m_1;

    int ir12 = 0, ir14 = 0, ir24 = 0;
    const TOutput r12 =  this->_chalf*(TOutput(k12)+TOutput(Sign(Real(k12)))*Sqrt(TOutput((k12-this->_two)*(k12+this->_two))));
    const TOutput r14 =  this->_chalf*(TOutput(k14)+TOutput(Sign(Real(k14)))*Sqrt(TOutput((k14-this->_two)*(k14+this->_two))));
    const TOutput r24 =  this->_chalf*(TOutput(k24)+TOutput(Sign(Real(k24)))*Sqrt(TOutput((k24-this->_two)*(k24+this->_two))));
    if (Real(k12) < -this->_two) ir12 = this->_ten*Sign(this->_one-Abs(r12));
    if (Real(k14) < -this->_two) ir14 = this->_ten*Sign(this->_one-Abs(r14));
    if (Real(k24) < -this->_two) ir24 = this->_ten*Sign(this->_one-Abs(r24));

    const TOutput q24 = r24-this->_cone/r24;
    const TOutput q12 = TOutput(k12)-r24*TOutput(k14);

    const TOutput a = TOutput(k34)/r24 - TOutput(k23);
    const TOutput b = TOutput(k12*k34) - TOutput(k13)*q24 - TOutput(k14*k23);
    const TOutput c = TOutput(k13)*q12 + r24*TOutput(k34) - TOutput(k23);
    const TOutput d = TOutput( (k12*k34 - k13*k24 - k14*k23)*(k12*k34 - k13*k24 - k14*k23)) -
                               this->_cfour*TOutput(k13*(k13 - k23*(k12 - k14*k24)) +
                                            k23*(k23 - k24*k34) + k34*(k34 - k13*k14));
    const TOutput discr = Sqrt(d);

    TOutput x4[2], x1[2], l4[2];
    x4[0] = this->_chalf*(b - discr)/a;
    x4[1] = this->_chalf*(b + discr)/a;
    if (Abs(x4[0]) > Abs(x4[1]))
      x4[1] = c/(a*x4[0]);
    else
      x4[0] = c/(a*x4[1]);

    const TOutput dd = -TOutput(k34)*r24 + TOutput(k23);
    TScale ix4[2], ix1[2];
    ix4[0] = Sign(Real(dd));
    ix4[1] = -ix4[0];

    x1[0] = x4[0]/r24;
    x1[1] = x4[1]/r24;
    ix1[0] = Sign(ix4[0]*Real(r24));
    ix1[1] = -ix1[0];

    const TOutput cc = this->cLn(Real(k13),-1);
    l4[0] = cc + this->cLn((q12 + q24*x4[0])/dd, Real(q24*ix4[0]/dd));
    l4[1] = cc + this->cLn((q12 + q24*x4[1])/dd, Real(q24*ix4[1]/dd));

    res[2] = res[1] = this->_czero;
    res[0] = (
          this->xspence(x4, ix4, r14, ir14) +
          this->xspence(x4, ix4, this->_cone/r14, -ir14) -
          this->xspence(x4, ix4, TOutput(k34/k13), -Real(k13)) -
          this->xspence(x1, ix1, r12, ir12) -
          this->xspence(x1, ix1, this->_cone/r12, -ir12) +
          this->xspence(x1, ix1, TOutput(k23/k13), -Real(k13)) -
          TOutput{this->_zero, this->_two*this->_pi}*this->xetatilde(x4, ix4, this->_cone/r24, -ir24, l4)
          )/(TOutput(m3*m_2)*discr);
  }

  /*!
   * Finite box with 4 non-zero masses. Formulae from \cite Denner:1991qq.
   * \param res output object res[0,1,2] the coefficients in the Laurent series, following the LoopTools implementation \cite Hahn:2006qw.
   * \param Y the modified Cayley matrix.
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Box<TOutput,TMass,TScale>::BIN4(vector<TOutput> &res, TMass const (&Y)[4][4]) const
  {
    TMass tmp;
    TMass M[4];
    for (int j = 0; j < 4; j++)
      M[j] = Y[j][j];

    TMass k12 = this->_two*Y[0][1]/Sqrt(M[0]*M[1]);
    TMass k23 = this->_two*Y[1][2]/Sqrt(M[1]*M[2]);
    TMass k34 = this->_two*Y[2][3]/Sqrt(M[2]*M[3]);
    TMass k14 = this->_two*Y[0][3]/Sqrt(M[0]*M[3]);
    TMass k13 = this->_two*Y[0][2]/Sqrt(M[0]*M[2]);
    TMass k24 = this->_two*Y[1][3]/Sqrt(M[1]*M[3]);

    if (Abs(k13) >= this->_two) { /*do nothing*/ }
    else if (Abs(k12) >= this->_two)
      {
        // 2 <-> 3
        tmp = k12;
        k12 = k13;
        k13 = tmp;
        tmp = k24;
        k24 = k34;
        k34 = tmp;
      }
    else if (Abs(k14) >= this->_two)
      {
        tmp = k13;
        k13 = k14;
        k14 = tmp;
        tmp = k23;
        k23 = k24;
        k24 = tmp;
      }
    else if (Abs(k23) >= this->_two)
      {
        tmp = k13;
        k13 = k23;
        k23 = tmp;
        tmp = k14;
        k14 = k24;
        k24 = tmp;
      }
    else if (Abs(k24) >= this->_two)
      {
        tmp = k12;
        k12 = k23;
        k23 = k34;
        k34 = k14;
        k14 = tmp;
        tmp = k13;
        k13 = k24;
        k24 = tmp;
      }
    else if (Abs(k34) >= this->_two)
      {
        tmp = k12;
        k12 = k24;
        k24 = tmp;
        tmp = k13;
        k13 = k34;
        k34 = tmp;
      }
    // else nothing found, all r_ij on the complex unit circle

    TMass kij[6];
    kij[0] = k12;
    kij[1] = k23;
    kij[2] = k34;
    kij[3] = k14;
    kij[4] = k13;
    kij[5] = k24;

    TOutput rij[6];
    TOutput r12 = rij[0] = this->_chalf*(TOutput(k12)+TOutput(Sign(Real(k12)))*Sqrt(TOutput((k12-this->_two)*(k12+this->_two))));
    TOutput r23 = rij[1] = this->_chalf*(TOutput(k23)+TOutput(Sign(Real(k23)))*Sqrt(TOutput((k23-this->_two)*(k23+this->_two))));
    TOutput r34 = rij[2] = this->_chalf*(TOutput(k34)+TOutput(Sign(Real(k34)))*Sqrt(TOutput((k34-this->_two)*(k34+this->_two))));
    TOutput r14 = rij[3] = this->_chalf*(TOutput(k14)+TOutput(Sign(Real(k14)))*Sqrt(TOutput((k14-this->_two)*(k14+this->_two))));
    TOutput r13 = rij[4] = this->_cone/(this->_chalf*(TOutput(k13)+TOutput(Sign(Real(k13)))*Sqrt(TOutput((k13-this->_two)*(k13+this->_two)))));
    TOutput r24 = rij[5] = this->_cone/(this->_chalf*(TOutput(k24)+TOutput(Sign(Real(k24)))*Sqrt(TOutput((k24-this->_two)*(k24+this->_two)))));

    TScale irij[6];
    for (int j = 0; j < 6; j++)
      {
        if (Imag(rij[j]) == this->_zero)
          {
            const TOutput ki = TOutput(kij[j]) - this->_ieps50;
            const TOutput kk = this->_chalf*(ki+TOutput(Sign(Real(ki)))*Sqrt((ki-this->_ctwo)*(ki+this->_ctwo)));
            irij[j] = Sign(Abs(rij[j])-this->_one)*Imag(kk);
          }
        else
          irij[j] = this->_zero;
      }

    TScale ir13 = irij[4], ir24 = irij[5];
    const TScale ir1324 = Sign(Real(r24))*ir13 - Sign(Real(r13))*ir24;

    const TOutput a = TOutput(k34)/r24 - TOutput(k23) + (TOutput(k12) - TOutput(k14)/r24)*r13;
    const TOutput b = (this->_cone/r13 - r13)*(this->_cone/r24 - r24) + TOutput(k12*k34) - TOutput(k14*k23);
    const TOutput c = TOutput(k34)*r24 - TOutput(k23) + (TOutput(k12) - TOutput(k14)*r24)/r13;
    const TOutput d = TOutput(k23) + (r24*TOutput(k14) - TOutput(k12))*r13 - r24*TOutput(k34);
    TOutput disc = Sqrt(b*b - this->_cfour*a*(c+d*this->_ieps50));

    TScale ix[2][4];
    ix[0][3] = Imag(this->_chalf/a*(b-disc));
    ix[1][3] = Imag(this->_chalf/a*(b+disc));

    disc = Sqrt(b*b - this->_cfour*a*c);
    TOutput x[2][4];
    x[0][3] = this->_chalf/a*(b-disc);
    x[1][3] = this->_chalf/a*(b+disc);
    if (Abs(x[0][3]) > Abs(x[1][3]))
      x[1][3] = c/(a*x[0][3]);
    else
      x[0][3] = c/(a*x[1][3]);

    x[0][0] = x[0][3]/r24;
    x[1][0] = x[1][3]/r24;
    x[0][1] = x[0][3]*r13/r24;
    x[1][1] = x[1][3]*r13/r24;
    x[0][2] = x[0][3]*r13;
    x[1][2] = x[1][3]*r13;

    const TScale s1 = Sign(Real(x[0][3]));
    const TScale s2 = Sign(Real(x[1][3]));
    ix[0][0] = ix[0][3]*Real(x[0][0])*s1;
    ix[1][0] = ix[1][3]*Real(x[1][0])*s2;
    ix[0][1] = ix[0][3]*Real(x[0][1])*s1;
    ix[1][1] = ix[1][3]*Real(x[1][1])*s2;
    ix[0][2] = ix[0][3]*Real(x[0][2])*s1;
    ix[1][2] = ix[1][3]*Real(x[1][2])*s2;

    res[0] = this->_czero;
    for (int j = 0; j < 4; j++)
      {
        const TOutput x_in[2] = { x[0][j], x[1][j] };
        const TScale ix_in[2] = { ix[0][j], ix[1][j] };
        res[0] += Pow(-this->_cone,j+1)*(
              this->xspence(x_in, ix_in, rij[j], irij[j]) +
              this->xspence(x_in, ix_in, this->_cone/rij[j], -irij[j])
              );
      }

    const TScale gamma = Sign(Real(a*(x[1][3] - x[0][3]))+this->_reps);
    TOutput l[2][4];
    for (int i = 0; i < 2; i++)
      for (int j = 0; j < 4; j++)
        l[i][j] = this->_czero;
    l[0][3] = this->_2ipi*TOutput(this->eta(r13, ir13, this->_cone/r24, -ir24, ir1324));
    l[1][3] = l[0][3];

    TOutput etas = this->_czero;
    if (Imag(r13) == 0)
      {
        r12 = TOutput(k12) - r24*TOutput(k14);
        r23 = TOutput(k23) - r24*TOutput(k34);
        r34 = TOutput(k34) - r13*TOutput(k14);
        r14 = TOutput(k23) - r13*TOutput(k12);
        const TOutput q13 = TOutput(k13) - this->_ctwo*r13;
        const TOutput q24 = TOutput(k24) - this->_ctwo*r24;

        TScale cc = gamma*Sign(Imag(r24)+ir24);
        l[0][0] = this->cLn(-x[0][0], -ix[0][0]) +
            this->cLn(r14 - q13/x[0][0], -this->_one) +
            this->cLn((r12 - q24*x[0][3])/d, cc);
        l[1][0] = this->cLn(-x[1][0], -ix[1][0]) +
            this->cLn(r14 - q13/x[1][0], -this->_one) +
            this->cLn((r12 - q24*x[1][3])/d, -cc);

        cc = gamma*Sign(Real(r13)*(Imag(r24)+ir24));
        l[0][1] = this->cLn(-x[0][1], -ix[0][1]) +
            this->cLn(r14 - q13/x[0][0], -this->_one) +
            this->cLn((r23 - q24*x[0][2])/d, cc);
        l[1][1] = this->cLn(-x[1][1], -ix[1][1]) +
            this->cLn(r14 - q13/x[1][0], -this->_one) +
            this->cLn((r23 - q24*x[1][2])/d, -cc);
        l[0][2] = this->cLn(-x[0][2], -ix[0][2]) +
            this->cLn(r34 - q13/x[0][3], -this->_one) +
            this->cLn((r23 - q24*x[0][2])/d, cc);
        l[1][2] = this->cLn(-x[1][2], -ix[1][2]) +
            this->cLn(r34 - q13/x[1][3], -this->_one) +
            this->cLn((r23 - q24*x[1][2])/d, -cc);

        const TOutput x_in[2] = {x[0][3], x[1][3]};
        const TScale ix_in[2] = {ix[0][3], ix[1][3]};
        const TOutput l_in_a[2] = {l[0][2], l[1][2]};
        const TOutput l_in_b[2] = {l[0][0], l[1][0]};
        const TOutput l_in_c[2] = {l[0][1], l[1][1]};
        const TOutput l_in_d[2] = {l[0][3], l[1][3]};

        etas = this->xetatilde(x_in, ix_in, r13, ir13, l_in_a) +
            this->xetatilde(x_in, ix_in, this->_cone/r24, -ir24, l_in_b) -
            this->xetatilde(x_in, ix_in, r13/r24, ir1324, l_in_c) +
            this->xetatilde(x_in, ix_in, -r13/r24, -ir1324, l_in_d);
      }
    else
      {
        for (int j = 0; j < 3; j++)
          {
            l[0][j] = Log(-x[0][j]) + this->cLn(TOutput(kij[j])-this->_cone/x[0][j]-x[0][j], -Real(x[0][j]*b*TOutput(gamma)));
            l[1][j] = Log(-x[1][j]) + this->cLn(TOutput(kij[j])-this->_cone/x[1][j]-x[1][j], -Real(x[1][j]*b*TOutput(gamma)));
          }

        const TOutput x_in[2] = {x[0][3], x[1][3]};
        const TScale ix_in[2] = {ix[0][3], ix[1][3]};
        const TOutput l_in_a[2] = {l[0][2], l[1][2]};
        const TOutput l_in_b[2] = {l[0][0], l[1][0]};
        const TOutput l_in_c[2] = {l[0][1], l[1][1]};
        const TOutput l_in_d[2] = {l[0][3], l[1][3]};

        etas = this->xeta(x_in, ix_in, r13, ir13, ix[0][2], l_in_a) +
            this->xeta(x_in, ix_in, this->_cone/r24, -ir24, ix[0][0], l_in_b) -
            this->xeta(x_in, ix_in, r13/r24, ir1324, ix[0][1], l_in_c) +
            this->xeta(x_in, ix_in, -r13/r24, -ir1324, ix[0][3], l_in_d)*
            (this->_cone - TOutput(Sign(Real(b))*gamma));
      }

    res[0] = (res[0] - this->_2ipi*TOutput(etas) +
        (l[1][1]-l[0][1])*l[0][3])/(TOutput(Sqrt(M[0]*M[1]*M[2]*M[3]))*disc) ;

    res[2] = res[1] = this->_czero;
  }

  /*!
   * The integral is defined as:
   * \f[
   * I_4^{D=4-2\epsilon}(0,0,0,0;s_{12},s_{23};0,0,0,0) = \frac{\mu^{2\epsilon}}{s_{12}s_{23}}\left[ \frac{2}{\epsilon^2}\left( (-s_{12}-i\epsilon)^{-\epsilon} + (-s_{23}-i \epsilon)^{-\epsilon} \right) - \ln^2 \left( \frac{-s_{12}-i\epsilon}{-s_{23}-i\epsilon} \right) - \pi^2 \right] + \mathcal{O}(\epsilon)
   * \f]
   * Implementation of the formulae from Bern et al. \cite Bern:1993kr.
   *
   * \param res output object res[0,1,2] the coefficients in the Laurent series
   * \param mu2 the energy scale squared.
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Box<TOutput,TMass,TScale>::B1(vector<TOutput> &res, TMass const (&Y)[4][4], TScale const& mu2) const
  {
    const TMass si = this->_two*Y[0][2];
    const TMass ta = this->_two*Y[1][3];
    const TOutput fac = this->_cone/(si*ta);
    const TOutput lnrat_tamu2 = this->Lnrat(ta,mu2);
    const TOutput lnrat_simu2 = this->Lnrat(si, mu2);
    const TOutput lnrat_tasi = this->Lnrat(ta,si);

    res[2] = fac*this->_ctwo*this->_ctwo;
    res[1] = fac*this->_ctwo*(-lnrat_tamu2 - lnrat_simu2);
    res[0] = fac*(lnrat_tamu2*lnrat_tamu2+lnrat_simu2*lnrat_simu2-lnrat_tasi*lnrat_tasi-this->_pi2);
  }

  /*!
   * The integral is defined as:
   * \f[
   * I_4^{D=4-2\epsilon}(0,0,0,p_4^2;s_{12},s_{23};0,0,0,0) = \frac{\mu^{2\epsilon}}{s_{12}s_{23}}\left[ \frac{2}{\epsilon^2}\left( (-s_{12})^{-\epsilon} + (-s_{23})^{-\epsilon} -(-p_4^2)^{-\epsilon} \right) - 2 {\rm Li}_2 \left( 1-\frac{p_4^2}{s_{12}} \right) - 2 {\rm Li}_2 \left( 1-\frac{p_4^2}{s_{23}} \right) -  \ln^2 \left( \frac{-s_{12}}{-s_{23}} \right) - \frac{\pi^2}{3} \right] + \mathcal{O}(\epsilon)
   * \f]
   * Implementation of the formulae from Bern et al. \cite Bern:1993kr.
   *
   * \param res output object res[0,1,2] the coefficients in the Laurent series
   * \param mu2 the energy scale squared.
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Box<TOutput,TMass,TScale>::B2(vector<TOutput> &res, TMass const (&Y)[4][4], TScale const& mu2) const
  {
    const TMass si = this->_two*Y[0][2];
    const TMass ta = this->_two*Y[1][3];
    const TMass mp4sq = this->_two*Y[0][3];
    const TOutput fac = this->_cone/(si*ta);
    const TOutput ln_mp4_mu2 = this->Lnrat(mp4sq,mu2);
    const TOutput ln_ta_mu2  = this->Lnrat(ta,mu2);
    const TOutput ln_si_mu2  = this->Lnrat(si,mu2);
    const TOutput ln_mp4_ta  = this->Lnrat(mp4sq,ta);
    const TOutput ln_mp4_si  = this->Lnrat(mp4sq,si);
    const TOutput ln_ta_si   = this->Lnrat(ta,si);    

    res[2] = fac*this->_ctwo;
    res[1] = res[2]*(ln_mp4_mu2-ln_ta_mu2-ln_si_mu2);
    res[0] = fac*(-ln_mp4_mu2*ln_mp4_mu2+ln_ta_mu2*ln_ta_mu2+ln_si_mu2*ln_si_mu2
                  +this->_ctwo*(this->Li2omrat(ta,mp4sq)+this->Li2omrat(si,mp4sq)-this->_pi2o6)
                  + ln_mp4_ta*ln_mp4_ta + ln_mp4_si*ln_mp4_si - ln_ta_si*ln_ta_si);
  }

  /*!
   * The integral is defined as:
   * \f[
   * I_4^{D=4-2\epsilon}(0,p_2^2,0,p_4^2;s_{12},s_{23};0,0,0,0) = \frac{\mu^{2\epsilon}}{s_{12}s_{23}-p_4^2 p_2^2} \\
   *    \left[ \frac{2}{\epsilon^2}\left( (-s_{12})^{-\epsilon} + (-s_{23})^{-\epsilon}-(-p_2^2)^{-\epsilon}-(-p_4^2)^{-\epsilon} \right) \\
   * - 2 {\rm Li}_2 \left( 1-\frac{p_2^2}{s_{12}} \right) - 2 {\rm Li}_2 \left( 1-\frac{p_2^2}{s_{23}} \right) - 2 {\rm Li}_2 \left( 1-\frac{p_4^2}{s_{12}} \right) - 2 {\rm Li}_2 \left( 1-\frac{p_4^2}{s_{23}} \right) \\
   * + 2 {\rm Li}_2 \left( 1-\frac{p_4^2 p_2^2}{s_{23}s_{12}} \right) -  \ln^2 \left( \frac{-s_{12}}{-s_{23}} \right) \right] + \mathcal{O}(\epsilon)
   * \f]
   * Implementation of the formulae from Bern et al. \cite Bern:1993kr.
   *
   * \param res output object res[0,1,2] the coefficients in the Laurent series
   * \param mu2 the energy scale squared.
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Box<TOutput,TMass,TScale>::B3(vector<TOutput> &res, TMass const (&Y)[4][4], TScale const& mu2) const
  {
    const TMass si = this->_two*Y[0][2];
    const TMass ta = this->_two*Y[1][3];
    const TMass mp4sq = this->_two*Y[0][3];
    const TMass mp2sq = this->_two*Y[1][2];
    const TMass r = this->_one-mp2sq*mp4sq/(si*ta);
    const TScale SignRealsi = Sign(Real(si));
    const TScale SignRealmp2sq = Sign(Real(mp2sq));

    // use expansion only in cases where signs are not ++-- or --++
    const bool landau = ( (SignRealsi == Sign(Real(ta))) &&
                          (SignRealmp2sq == Sign(Real(mp4sq))) &&
                          (SignRealsi != SignRealmp2sq) );

    if ((Abs(r) < this->_eps) && (landau == false))
      {
        // expanded case
        const TOutput fac = this->_cone/(si*ta);
        const TOutput ln_si_mu2 = this->Lnrat(si,mu2);
        const TOutput ln_ta_mp4 = this->Lnrat(ta,mp4sq);
        const TOutput l0_mp4_ta = this->L0(mp4sq,ta);
        const TOutput l0_mp4_si = this->L0(mp4sq,si);
        const TOutput l1_mp4_ta = this->L1(mp4sq,ta);
        const TOutput l1_mp4_si = this->L1(mp4sq,si);

        res[2] = this->_czero;
        res[1] = -(this->_ctwo+r)*fac;
        res[0] = fac*( this->_ctwo-this->_chalf*r + (this->_ctwo+r)*(ln_si_mu2+ln_ta_mp4)
                      +this->_ctwo*(l0_mp4_ta+l0_mp4_si)+r*(l1_mp4_ta+l1_mp4_si));
      }
    else
      {
        const TOutput fac = this->_cone/(si*ta-mp2sq*mp4sq);
        const TOutput ln_mp2_mu2 = this->Lnrat(mp2sq,mu2);
        const TOutput ln_mp4_mu2 = this->Lnrat(mp4sq,mu2);
        const TOutput ln_si_mu2  = this->Lnrat(si,mu2);
        const TOutput ln_ta_mu2  = this->Lnrat(ta,mu2);
        const TOutput ln_s_t = this->Lnrat(si, ta);
        const TOutput li2_1 = this->Li2omrat(mp2sq, si);
        const TOutput li2_2 = this->Li2omrat(mp2sq, ta);
        const TOutput li2_3 = this->Li2omrat(mp4sq, si);
        const TOutput li2_4 = this->Li2omrat(mp4sq, ta);
        const TOutput li2_5 = this->Li2omx2(mp2sq,mp4sq,si,ta);

        res[2] = this->_czero;
        res[1] = fac*this->_ctwo*(this->Lnrat(mp2sq,si) + this->Lnrat(mp4sq,ta));
        res[0] = fac*(ln_si_mu2*ln_si_mu2 + ln_ta_mu2*ln_ta_mu2
                      -ln_mp2_mu2*ln_mp2_mu2 - ln_mp4_mu2*ln_mp4_mu2
                      +this->_ctwo*(li2_5-li2_1-li2_2-li2_3-li2_4-this->_chalf*ln_s_t*ln_s_t));
      }
  }

  /*!
   * The integral is defined as:
   * \f[
   * I_4^{D=4-2\epsilon}(0,0,p_3^2,p_4^2;s_{12},s_{23};0,0,0,0) = \frac{\mu^{2\epsilon}}{s_{12}s_{23}} \\
   *     \left[ \frac{2}{\epsilon^2}\left( (-s_{12})^{-\epsilon} + (-s_{23})^{-\epsilon}-(-p_3^2)^{-\epsilon}-(-p_4^2)^{-\epsilon} \right) + \frac{1}{\epsilon^2} \left( (-p_3^2)^{-\epsilon}(-p_4)^{-\epsilon} \right) / (-s_{12})^{-\epsilon} \\
   * - 2 {\rm Li}_2 \left( 1-\frac{p_3^2}{s_{23}} \right) - 2 {\rm Li}_2 \left( 1-\frac{p_4^2}{s_{23}} \right) -  \ln^2 \left( \frac{-s_{12}}{-s_{23}} \right) \right] + \mathcal{O}(\epsilon)
   * \f]
   * Implementation of the formulae from Bern et al. \cite Bern:1993kr.
   *
   * \param res output object res[0,1,2] the coefficients in the Laurent series
   * \param mu2 the energy scale squared.
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Box<TOutput,TMass,TScale>::B4(vector<TOutput> &res, TMass const (&Y)[4][4], TScale const& mu2) const
  {
    const TMass si = this->_two*Y[0][2];
    const TMass ta = this->_two*Y[1][3];
    const TMass mp4sq = this->_two*Y[0][3];
    const TMass mp3sq = this->_two*Y[2][3];
    const TOutput fac = this->_cone/(si*ta);
    const TOutput ln_ta_mu2 = this->Lnrat(ta,mu2);
    const TOutput ln_si_mu2 = this->Lnrat(si,mu2);
    const TOutput ln_mp3_mu2 = this->Lnrat(mp3sq, mu2);
    const TOutput ln_mp4_mu2 = this->Lnrat(mp4sq, mu2);
    const TOutput ln_si_ta = this->Lnrat(si,ta);
    const TOutput ln_si_mp3 = this->Lnrat(si,mp3sq);

    res[2] = fac;
    res[1] = -fac*(ln_si_mp3 + this->Lnrat(ta, mp4sq) + ln_ta_mu2);
    res[0] = fac*(ln_ta_mu2*ln_ta_mu2
                  + this->_chalf*ln_si_mu2*ln_si_mu2
                  - this->_chalf*ln_mp3_mu2*ln_mp3_mu2
                  - this->_chalf*ln_mp4_mu2*ln_mp4_mu2
                  + this->_ctwo*(-this->Li2omrat(mp3sq,ta)-this->Li2omrat(mp4sq,ta)+
                                 this->_chalf*(ln_si_mp3*this->Lnrat(si,mp4sq) - ln_si_ta*ln_si_ta) )
                  );
  }

  /*!
   * The integral is defined as:
   * \f[
   * I_4^{D=4-2\epsilon}(0,p_2^2,p_3^2,p_4^2;s_{12},s_{23};0,0,0,0) = \frac{\mu^{2\epsilon}}{s_{12}s_{23}-p_2^2 p_4^2} \\
   *     \left[ \frac{2}{\epsilon^2}\left( (-s_{12})^{-\epsilon} + (-s_{23})^{-\epsilon}-(-p_2^2)^{-\epsilon}-(-p_3^2)^{-\epsilon} -(-p_4^2)^{-\epsilon} \right) \\
   *  + \frac{1}{\epsilon^2} \left( (-p_2^2)^{-\epsilon}(-p_3)^{-\epsilon} \right) / (-s_{23})^{-\epsilon}+ \frac{1}{\epsilon^2} \left( (-p_3^2)^{-\epsilon}(-p_4)^{-\epsilon} \right) / (-s_{12})^{-\epsilon} \\
   * - 2 {\rm Li}_2 \left( 1-\frac{p_2^2}{s_{12}} \right) - 2 {\rm Li}_2 \left( 1-\frac{p_4^2}{s_{23}} \right) + 2 {\rm Li}_2 \left( 1-\frac{p_2^2 p_4^2}{s_{12} s_{23}} \right) -  \ln^2 \left( \frac{-s_{12}}{-s_{23}} \right) \right] + \mathcal{O}(\epsilon)
   * \f]
   * Implementation of the formulae from Bern et al. \cite Bern:1993kr.
   *
   * \param res output object res[0,1,2] the coefficients in the Laurent series
   * \param mu2 the energy scale squared.
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Box<TOutput,TMass,TScale>::B5(vector<TOutput> &res, TMass const (&Y)[4][4], TScale const& mu2) const
  {
    const TMass si = this->_two*Y[0][2];
    const TMass ta = this->_two*Y[1][3];
    const TMass mp2sq = this->_two*Y[1][2];
    const TMass mp3sq = this->_two*Y[2][3];
    const TMass mp4sq = this->_two*Y[0][3];
    const TMass r = this->_one - mp2sq*mp4sq/(si*ta);
    const TScale SignRealsi = Sign(Real(si));
    const TScale SignRealmp2sq = Sign(Real(mp2sq));

    // use expansion only in cases where signs are not ++-- or --++
    const bool landau = ( (SignRealsi == Sign(Real(ta))) &&
                          (SignRealmp2sq == Sign(Real(mp4sq))) &&
                          (SignRealsi != SignRealmp2sq));

    if ( Abs(r) < this->_eps && landau == false)
      {
        const TOutput l0_mp4_ta = this->L0(mp4sq,ta);
        const TOutput l1_mp4_ta = this->L1(mp4sq,ta);

        // expanded case
        res[2] = this->_czero;
        res[1] = -(this->_cone+this->_chalf*r)/(si*ta);
        res[0] = res[1]*(this->Lnrat(mu2,si)+this->Lnrat(mp3sq,ta)-this->_ctwo
                         -(this->_cone+mp4sq/ta)*l0_mp4_ta)
                         +(r/(si*ta))*(l1_mp4_ta-l0_mp4_ta-this->_cone);
      }
    else
      {
        // general case
        const TOutput fac = this->_cone/(si*ta-mp2sq*mp4sq);
        const TOutput li2_1 = this->Li2omrat(mp2sq,si);
        const TOutput li2_2 = this->Li2omrat(mp4sq,ta);
        const TOutput li2_3 = this->Li2omx2(mp2sq,mp4sq,si,ta);
        const TOutput ln_ta_mp2 = this->Lnrat(ta,mp2sq);
        const TOutput ln_si_mp4 = this->Lnrat(si,mp4sq);
        const TOutput ln_si_ta = this->Lnrat(si,ta);

        res[2] = this->_czero;
        res[1] = -ln_ta_mp2 - ln_si_mp4;
        res[0] = -this->_chalf*(ln_ta_mp2*ln_ta_mp2 + ln_si_mp4*ln_si_mp4)
                                -(this->Lnrat(mp3sq,ta)+this->Lnrat(mu2,ta))*ln_ta_mp2
                                -(this->Lnrat(mp3sq,si)+this->Lnrat(mu2,si))*ln_si_mp4
                                -this->_ctwo*(li2_1+li2_2-li2_3)-ln_si_ta*ln_si_ta;
        res[1] *= fac;
        res[0] *= fac;
      }
  }

  /*!
   * The integral is defined as:
   * \f[
   * I_4^{D=4-2\epsilon}(0,0,m^2,m^2;s_{12},s_{23};0,0,0,m^2) = -\frac{1}{s_{12} (m^2-s_{23})} \left( \frac{\mu^2}{m^2} \right)^\epsilon \\
   *     \left[ \frac{2}{\epsilon^2} - \frac{1}{\epsilon} \left( 2 \ln \left( \frac{m^2-s_{23}}{m^2} \right) + \ln \left( \frac{-s_{12}}{m^2} \right) \right) + 2 \ln \left( \frac{m^2-s_{23}}{m^2} \right) \ln \left( \frac{-s_{12}}{m^2} \right) - \frac{\pi^2}{2} \right] + \mathcal{O}(\epsilon)
   * \f]
   * Implementation of the formulae from Beenakker et al. \cite Beenakker:1988bq.
   *
   * \param res output object res[0,1,2] the coefficients in the Laurent series
   * \param mu2 the energy scale squared.
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Box<TOutput,TMass,TScale>::B6(vector<TOutput> &res, TMass const (&Y)[4][4], TScale const& mu2) const
  {
    const TMass si = this->_two*Y[0][2];
    const TMass tabar = this->_two*Y[1][3];
    const TMass msq = Y[3][3];
    const TOutput wlogs = this->Lnrat(si, msq);
    const TOutput wlogt = this->Lnrat(tabar, msq);
    const TOutput wlogm = this->Lnrat(mu2, msq);

    res[2] = this->_ctwo;
    res[1] = this->_ctwo*(wlogm-wlogt)-wlogs;
    res[0] = wlogm*wlogm-wlogm*(this->_ctwo*wlogt+wlogs) + this->_ctwo*wlogt*wlogs-this->_chalf*this->_pi2;

    const TOutput d = TOutput(si*tabar);
    for (size_t i = 0; i < 3; i++)
      res[i] /= d;
  }

  /*!
   * The integral is defined as:
   * \f[
   * I_4^{D=4-2\epsilon}(0,0,m^2,p_2^2;s_{12},s_{23};0,0,0,m^2) = \left( \frac{\mu^2}{m^2} \right)^\epsilon \frac{1}{s_{12} (s_{23}-m^2)} \\
   *     \left[ \frac{3}{2 \epsilon^2} - \frac{1}{\epsilon} \left\{ 2 \ln \left( 1-\frac{s_{23}}{m^2} \right) + \ln \left( \frac{-s_{12}}{m^2} \right) - \ln \left( 1-\frac{p_4^2}{m^2} \right) \right\} \\
   *   -2 {\rm Li}_2 \left( 1 - \frac{m^2-p_4^2}{m^2-s_{23}} \right) + 2 \ln \left( \frac{-s_{12}}{m^2} \right) \ln \left( 1-\frac{s_{23}}{m^2} \right) - \ln^2 \left( 1 - \frac{p_4^2}{m^2} \right) -\frac{5\pi^2}{12} \right] + \mathcal{O}(\epsilon)
   * \f]
   * Implementation of the formulae from Beenakker et al. \cite Beenakker:1988bq.
   *
   * \param res output object res[0,1,2] the coefficients in the Laurent series
   * \param mu2 the energy scale squared.
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Box<TOutput,TMass,TScale>::B7(vector<TOutput> &res, TMass const (&Y)[4][4], TScale const& mu2) const
  {
    const TMass tabar = this->_two*Y[1][3];
    const TMass p4sqbar = this->_two*Y[0][3];
    const TMass si = this->_two*Y[0][2];
    const TMass msq = Y[3][3];
    const TOutput wlogs = this->Lnrat(si, msq);
    const TOutput wlogt = this->Lnrat(tabar, msq);
    const TOutput wlogp = this->Lnrat(p4sqbar, msq);
    const TOutput wlogm = this->Lnrat(mu2, msq);

    res[2] = TOutput(this->_three/this->_two);
    res[1] = TOutput(this->_three/this->_two)*wlogm - this->_ctwo*wlogt - wlogs + wlogp;
    res[0] = this->_ctwo*wlogs*wlogt - wlogp*wlogp - TOutput(this->_five*this->_pi2o12)
             +this->_three/this->_four*wlogm*wlogm+wlogm*(-this->_two*wlogt-wlogs+wlogp)
             -this->_two*this->Li2omrat(p4sqbar, tabar);

    const TOutput d = TOutput(si*tabar);
    for (size_t i = 0; i < 3; i++)
      res[i] /= d;
  }

  /*!
   * The integral is defined as:
   * \f[
   * I_4^{D=4-2\epsilon}(0,0,p_3^2,p_4^2;s_{12},s_{23};0,0,0,m^2) = \frac{1}{s_{12} (s_{23}-m^2)} \left[ \frac{1}{\epsilon^2} - \frac{1}{\epsilon} \left[ \ln \frac{-s_{12}}{\mu^2} + \ln \frac{(m^2-s_{23}^2)}{(m^2-p_3^2)(m^2-p_4^2)} \right] \\
   *   - 2 {\rm Li}_2 \left( 1 - \frac{m^2-p_3^2}{m^2-s_{23}} \right) - 2 {\rm Li}_2 \left( 1-\frac{m^2-p_4^2}{m^2-s_{23}} \right) - {\rm Li}_2 \left( 1 + \frac{(m^2-p_3^2)(m^2-p_4^2)}{s_{12} m^2} \right) \\
   *   - \frac{\pi^2}{6} + \frac{1}{2} \ln^2 \left( \frac{-s_{12}}{\mu^2} \right) - \frac{1}{2} \ln^2 \left( \frac{-s_{12}}{m^2} \right) + 2 \ln \left( \frac{-s_{12}}{\mu^2} \right) \ln \left( \frac{m^2-s_{23}}{m^2} \right) \\
   *   - \ln \left( \frac{m^2-p_3^2}{\mu^2} \right) \ln \left( \frac{m^2-p_3^2}{m^2} \right) - \ln \left( \frac{m^2-p_4^2}{\mu^2} \right) \ln \left( \frac{m^2-p_4^2}{m^2} \right)  \right] + \mathcal{O}(\epsilon)
   * \f]
   * Implementation of the formulae from Beenakker et al. \cite Beenakker:1988bq.
   *
   * \param res output object res[0,1,2] the coefficients in the Laurent series
   * \param mu2 the energy scale squared.
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Box<TOutput,TMass,TScale>::B8(vector<TOutput> &res, TMass const (&Y)[4][4], TScale const& mu2) const
  {
    const TMass msq = Y[3][3];
    const TMass tabar = this->_two*Y[1][3];
    const TMass si    = this->_two*Y[0][2];
    const TMass p3sqbar = this->_two*Y[2][3];
    const TMass p4sqbar = this->_two*Y[0][3];
    const TOutput wlogs = this->Lnrat(si, mu2);
    const TOutput wlogp3 = this->Lnrat(p3sqbar, tabar);
    const TOutput wlogp4 = this->Lnrat(p4sqbar, tabar);

    const TOutput dilog3 = this->Li2omrat(p3sqbar, tabar);
    const TOutput dilog4 = this->Li2omrat(p4sqbar, tabar);
    const TOutput dilog34 = this->Li2omx2(p3sqbar,p4sqbar,si,msq);
    const TOutput ln_si_mu2 = this->Lnrat(si, mu2);

    res[2] = this->_cone;
    res[1] = wlogp3 + wlogp4 - wlogs;
    res[0] = -this->_ctwo*dilog3 - this->_ctwo*dilog4 - dilog34
             -TOutput(this->_pi2o6)+this->_chalf*(ln_si_mu2*ln_si_mu2-Pow(this->Lnrat(si, msq),2))
             +this->_ctwo*ln_si_mu2*this->Lnrat(tabar, msq)
             -this->Lnrat(p3sqbar,mu2)*this->Lnrat(p3sqbar,msq)
             -this->Lnrat(p4sqbar,mu2)*this->Lnrat(p4sqbar,msq);

    const TOutput d = TOutput(si*tabar);
    for (size_t i = 0; i < 3; i++)
      res[i] /= d;
  }

  /*!
   * The integral is defined as:
   * \f[
   * I_4^{D=4-2\epsilon}(0,p_2^2,p_3^2,m^2;s_{12},s_{23};0,0,0,m^2) = \frac{1}{s_{12} (s_{23}-m^2)} \left[ \frac{1}{2 \epsilon^2} - \frac{1}{\epsilon} \left( \frac{s_{12} (m^2-s_{23})}{p_2^2 \mu m} \right) \\
   *  + {\rm Li}_2 \left(1+\frac{(m^2-p_3^2)(m^2-s_{23})}{m^2 p_2^2} \right) + 2 {\rm Li}_2 \left( 1-\frac{s_{12}}{p_2^2} \right) + \frac{\pi^2}{12} + \ln^2 \left( \frac{s_{12}(m^2-s_{23})}{p_2^2 \mu m} \right) \right] + \mathcal{O}(\epsilon)
   * \f]
   * Implementation of the formulae from Ellis et al. \cite Ellis:2007qk.
   *
   * \param res output object res[0,1,2] the coefficients in the Laurent series
   * \param mu2 the energy scale squared.
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Box<TOutput,TMass,TScale>::B9(vector<TOutput> &res, TMass const (&Y)[4][4], TScale const& mu2) const
  {
    const TMass msq = Y[3][3];
    const TMass mean = Sqrt(TMass(mu2*msq));
    const TMass tabar = this->_two*Y[1][3];
    const TMass si = this->_two*Y[0][2];
    const TMass m3sqbar = this->_two*Y[2][3];
    const TMass mp2sq = this->_two*Y[1][2];
    const TOutput fac = TOutput(si*tabar);

    const TOutput wlogt = this->Lnrat(tabar,mean);
    const TOutput wlog2 = this->Lnrat(si,mp2sq);

    const TOutput dilog1 = this->Li2omx2(m3sqbar, tabar, mp2sq, msq);
    const TOutput dilog2 = this->Li2omrat(si,mp2sq);

    res[2] = this->_chalf;
    res[1] = -wlogt-wlog2;
    res[0] = dilog1+this->_two*dilog2 + TOutput(this->_pi2o12) + (wlogt+wlog2)*(wlogt+wlog2);

    for (size_t i = 0; i < 3; i++)
      res[i] /= fac;
  }

  /*!
   * The integral is defined as:
   * \f[
   * I_4^{D=4-2\epsilon}(0,p_2^2,p_3^2,p_4^2;s_{12},s_{23};0,0,0,m^2) = \frac{1}{(s_{12}s_{23}-m^2 s_{12} - p_2^2 p_4^2 + m^2 p_2^2)} \\
   *  \left[ \frac{1}{\epsilon} \ln \left( \frac{(m^2-p_4^2) p_2^2}{(m^2-s_{23})s_{12})} \right) + {\rm Li}_2 \left( 1 + \frac{(m^2-p_3^2)(m^2-s_{23})}{p_2^2 m^2} \right) - {\rm Li}_2 \left( 1 + \frac{(m^2-p_3^2)(m^2-p_4^2)}{s_{12} m^2} \right)  \\
   *   +2 {\rm Li}_2 \left( \right) - 2 {\rm Li}_2 \left( 1-\frac{p_2^2}{s_{12}} \right) + 2 {\rm Li}_2 \left( 1-\frac{p_2 (m^2-p_4^2)}{s_{12}(m^2-s_{23})} \right) \\
   *   +2 \ln \left( \frac{\mu m}{m^2-s_{23}} \right) \ln \left( \frac{(m^2-p_4^2) p_2^2}{(m^2-s_{23}) s_{12}} \right) \right] + \mathcal{O}(\epsilon)
   * \f]
   * Implementation of the formulae from Ellis et al. \cite Ellis:2007qk.
   *
   * \param res output object res[0,1,2] the coefficients in the Laurent series
   * \param mu2 the energy scale squared.
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Box<TOutput,TMass,TScale>::B10(vector<TOutput> &res, TMass const (&Y)[4][4], TScale const& mu2) const
  {
    const TMass msq = Y[3][3];
    const TMass si = this->_two*Y[0][2];
    const TMass tabar = this->_two*Y[1][3];
    const TMass m4sqbar = this->_two*Y[0][3];
    const TMass m3sqbar = this->_two*Y[2][3];
    const TMass mp2sq = this->_two*Y[1][2];
    const TMass mean = Sqrt(TMass(mu2*msq));

    const TOutput fac = TOutput(si*tabar - mp2sq*m4sqbar);
    const TOutput wlogsmu = this->Lnrat(si, mu2);
    const TOutput wlogtmu = this->Lnrat(tabar, mu2);
    const TOutput wlog2mu = this->Lnrat(mp2sq, mu2);
    const TOutput wlog4mu = this->Lnrat(m4sqbar,mu2);

    const TOutput dilog1 = this->Li2omrat(mp2sq,si);
    const TOutput dilog2 = this->Li2omrat(tabar, m4sqbar);
    const TOutput dilog3 = this->Li2omx2(mp2sq,m4sqbar,si,tabar);
    const TOutput dilog4 = this->Li2omx2(m3sqbar, tabar, mp2sq, msq);
    const TOutput dilog5 = this->Li2omx2(m3sqbar, m4sqbar, si, msq);

    res[2] = this->_czero;
    res[1] = wlog2mu + wlog4mu - wlogsmu - wlogtmu;
    res[0] = dilog4 - dilog5
        -this->_two*dilog1 + this->_two*dilog2 + this->_two*dilog3
        +this->_two*res[1]*this->Lnrat(mean, tabar);

    for (size_t i = 0; i < 3; i++)
      res[i] /= fac;
  }

  /*!
   * The integral is defined as:
   * \f[
   * I_4^{D=4-2\epsilon}(0,m_3^2,p_3^2,m_4^2;s_{12},s_{23};0,0,m_3^2,m_4^2) = \frac{1}{(m_3^2-s_{12})(m_4^2 - s_{23})} \\
   *  \left[ \frac{1}{\epsilon^2} - \frac{1}{\epsilon} \ln \left( \frac{(m^2-s_{23})(m_3^2-s_{12})}{m_3 m_4 \mu^2} \right) + 2 \ln \left( \frac{m_3^2-s_{12}}{m_3 \mu} \right) \ln \left( \frac{m_4^2-s_{23}}{m_4 \mu} \right)  \\
   *  - \frac{\pi^2}{2} + \ln^2 \frac{m_3}{m_4} - \frac{1}{2} \ln^2 \left( \frac{\gamma^+_{34}}{\gamma^+_{34} - 1} \right) - \frac{1}{2} \ln \left( \frac{\gamma^-_{34}}{\gamma^-_{34} - 1} \right) \right] + \mathcal{O}(\epsilon)
   * \f]
   * Implementation of the formulae from Ellis et al. \cite Ellis:2007qk.
   *
   * \param res output object res[0,1,2] the coefficients in the Laurent series
   * \param mu2 the energy scale squared.
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Box<TOutput,TMass,TScale>::B11(vector<TOutput> &res, TMass const (&Y)[4][4], TScale const& mu2) const
  {
    const TMass m3sq = Y[2][2];
    const TMass m4sq = Y[3][3];
    const TMass sibar = this->_two*Y[0][2];
    const TMass tabar = this->_two*Y[1][3];
    const TMass p3sq = -(this->_two*Y[2][3]-Y[2][2]-Y[3][3]);
    const TMass m3mu = Sqrt(m3sq*mu2);
    const TMass m4mu = Sqrt(m4sq*mu2);

    const TOutput wlogt = this->Lnrat(tabar, m4mu);
    const TOutput wlogs = this->Lnrat(sibar, m3mu);

    TOutput root;
    TMass x43p, x43pm1, x43m, x43mm1;
    if (this->iszero(p3sq))
      {
        root = this->_cone;
        x43p = -this->_one;
        x43pm1 = -this->_one;
        x43m = m3sq;
        x43mm1 = m4sq;
      }
    else
      {
        root = Sqrt(TOutput(Pow(p3sq+m3sq-m4sq,2)-this->_four*m3sq*p3sq));
        const TOutput ga43p   = TOutput(+p3sq+m3sq-m4sq) + root;
        const TOutput ga43pm1 = TOutput(-p3sq+m3sq-m4sq) + root;
        const TOutput ga43m   = TOutput(+p3sq+m3sq-m4sq) - root;
        const TOutput ga43mm1 = TOutput(-p3sq+m3sq-m4sq) - root;

        x43p = -Real(ga43p);
        x43pm1 = -Real(ga43pm1);
        x43m = Real(ga43m);
        x43mm1 = Real(ga43mm1);
      }

    // deal with real roots
    TScale ieps2;
    TOutput ln43p, ln43m, rat2p, rat2m;
    if (this->iszero(Imag(root)))
      {
        ln43p = this->Lnrat(x43p, x43pm1);
        ln43m = this->Lnrat(x43m, x43mm1);
      }
    else
      {
        this->ratgam(rat2p, rat2m, ieps2, p3sq, m4sq, m3sq);
        ln43p = this->cLn(rat2p, ieps2);
        ln43m = this->cLn(rat2m, ieps2);
      }

    TOutput intbit;
    if (this->iszero(p3sq))
      intbit = -this->_chalf*TOutput(Pow(Log(m3sq/m4sq),2));
    else
      intbit = -this->_chalf*(Pow(ln43p,2)+Pow(ln43m,2));

    res[2] = this->_cone;
    res[1] = -wlogt-wlogs;
    res[0] = intbit
        + this->_ctwo*wlogt*wlogs-TOutput(this->_half*this->_pi2)
        + TOutput(Pow(Log(m3sq/m4sq),2)/this->_four);

    const TOutput d = TOutput(sibar*tabar);
    for (size_t i = 0; i < 3; i++)
      res[i] /= d;
  }

  /*!
   * The integral is defined as:
   * \f[
   * I_4^{D=4-2\epsilon}(0,m_3^2,p_3^2,p_4^2;s_{12},s_{23};0,0,m_3^2,m_4^2) = \frac{1}{(s_{12}-m_3^2)(s_{23}-m_4^2)} \\
   *  \left[ \frac{1}{2 \epsilon^2} - \frac{1}{\epsilon} \ln \left( \frac{(m_4^2-s_{23})(m_3^2-s_{12})}{(m_4-p_4^2) m_3 \mu} \right) + 2 \ln \left( \frac{m_4^2-s_{23}}{m_3 \mu} \right) \ln \left( \frac{m_3^2-s_{12}}{m_3 \mu} \right)  \\
   *  - \ln^2 \left( \frac{m_4^2-p_4^2}{m_3 \mu}\right) -\frac{\pi^2}{12} + \ln \left( \frac{m_4^2-p_4^2}{m_3^2-s_{12}} \right) \ln \left( \frac{m_4^2}{m_3^2} \right) - \frac{1}{2} \ln^2 \left( \frac{\gamma^+_{34}}{ \gamma^+_{34}-1 }\right) - \frac{1}{2} \ln^2 \left( \frac{\gamma^-_{34}}{\gamma^-_{34}-1} \right) \\
   *  - 2 {\rm Li}_2 \left( 1 - \frac{(m_4^2-p_4^2)}{(m_4^2-s_{23})} \right) - {\rm Li}_2 \left( 1 - \frac{(m_4-p_4^2) \gamma^+_{43}}{(m_3^2-s_{12})(\gamma^+_{43}-1)} \right)- {\rm Li}_2 \left( 1 - \frac{(m_4-p_4^2) \gamma^-_{43}}{(m_3^2-s_{12})(\gamma^-_{43}-1)} \right) \right] + \mathcal{O}(\epsilon)
   * \f]
   * Implementation of the formulae from Ellis et al. \cite Ellis:2007qk.
   *
   * \param res output object res[0,1,2] the coefficients in the Laurent series
   * \param mu2 the energy scale squared.
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Box<TOutput,TMass,TScale>::B12(vector<TOutput> &res, TMass const (&Y)[4][4], TScale const& mu2) const
  {
    const TMass m3sq = Y[2][2];
    const TMass m4sq = Y[3][3];
    const TMass sibar = this->_two*Y[0][2];
    const TMass tabar = this->_two*Y[1][3];
    const TMass m4sqbar = this->_two*Y[0][3];
    const TMass p3sq = -(this->_two*Y[2][3]-Y[2][2]-Y[3][3]);

    const TMass mean = Sqrt(mu2*Real(m3sq));
    const TOutput fac = TOutput(sibar*tabar);

    const TOutput wlogsmu = this->Lnrat(sibar, mean);
    const TOutput wlogtmu = this->Lnrat(tabar, mean);
    const TOutput wlog4mu = this->Lnrat(m4sqbar, mean);
    const TOutput wlog = wlogsmu + wlogtmu - wlog4mu;

    TOutput root;
    TMass x43p, x43pm1, x43m, x43mm1;
    if (this->iszero(p3sq))
      {
        root = this->_cone;
        x43p = -this->_one;
        x43pm1 = -this->_one;
        x43m = m3sq;
        x43mm1 = m4sq;
      }
    else
      {
        root = Sqrt(TOutput(Pow(p3sq+m3sq-m4sq,2)-this->_four*m3sq*p3sq));
        const TOutput ga43p   = TOutput(+p3sq+m3sq-m4sq) + root;
        const TOutput ga43pm1 = TOutput(-p3sq+m3sq-m4sq) + root;
        const TOutput ga43m   = TOutput(+p3sq+m3sq-m4sq) - root;
        const TOutput ga43mm1 = TOutput(-p3sq+m3sq-m4sq) - root;

        x43p = -Real(ga43p);
        x43pm1 = -Real(ga43pm1);
        x43m = Real(ga43m);
        x43mm1 = Real(ga43mm1);
      }

    const TOutput dilog1 = this->Li2omrat(m4sqbar, tabar);
    TOutput dilog2, dilog3;

    // deal with real roots
    TScale ieps2, ieps1;
    TMass rat1;
    TOutput ln43p, ln43m, rat2p, rat2m, zrat1;
    if (this->iszero(Imag(root)))
      {
        ln43p = this->Lnrat(x43p, x43pm1);
        ln43m = this->Lnrat(x43m, x43mm1);
        dilog2 = this->Li2omx2(m4sqbar, x43p, sibar, x43pm1);
        dilog3 = this->Li2omx2(m4sqbar, x43m, sibar, x43mm1);
      }
    else
      {
        this->ratreal(m4sqbar, sibar, rat1, ieps1);
        this->ratgam(rat2p, rat2m, ieps2, p3sq, m4sq, m3sq);
        zrat1 = TOutput(rat1);
        ln43p = this->cLn(rat2p, ieps2);
        ln43m = this->cLn(rat2m, ieps2);

        dilog2 = this->spencer(zrat1, rat2p, ieps1, ieps2);
        dilog3 = this->spencer(zrat1, rat2m, ieps1, ieps2);
      }

    res[2] = this->_chalf;
    res[1] = -wlog;
    res[0] = -TOutput(this->_pi2o12)
        +this->_ctwo*wlogsmu*wlogtmu-wlog4mu*wlog4mu
        +(wlog4mu-wlogsmu)*Log(m4sq/m3sq)-this->_half*(ln43p*ln43p + ln43m*ln43m)
        - this->_two*dilog1-dilog2-dilog3;

    for (size_t i = 0; i < 3; i++)
      res[i] /= fac;
  }

  /*!
   * The integral is defined as:
   * \f[
   * I_4^{D=4-2\epsilon}(0,p_2^2,p_3^2,p_4^2;s_{12},s_{23};0,0,m_3^2,m_4^2) = \frac{1}{\Delta} \left[ \frac{1}{\epsilon} \ln \left( \frac{(m_3^2-p_2^2)(m_4^2-p_4^2)}{(m_3^2-s_{12})(m_4^2-s_{23})} \right) \\
   * - 2 {\rm Li}_2 \left(1-\frac{(m_3^2-p_2^2)}{(m_3^2-s_{12})} \right) - {\rm Li}_2 \left( 1 - \frac{(m_3^2-p_2^2)\gamma^+_{34}}{(m_4^2-s_{23})(\gamma_{34}^+ - 1)} \right) - {\rm Li}_2 \left( 1 - \frac{(m_3^2-p_2^2)\gamma^-_{34}}{(m_4^2-s_{23})(\gamma_{34}^- - 1)} \right) \\
   * - 2 {\rm Li}_2 \left(1-\frac{(m_4^2-p_4^2)}{(m_4^2-s_{23})} \right) - {\rm Li}_2 \left( 1 - \frac{(m_4^2-p_4^2)\gamma^+_{43}}{(m_3^2-s_{12})(\gamma_{43}^+ - 1)} \right) - {\rm Li}_2 \left( 1 - \frac{(m_4^2-p_4^2)\gamma^-_{43}}{(m_2^2-s_{12})(\gamma_{43}^- - 1)} \right) \\
   * + 2 {\rm Li}_2 \left(1-\frac{(m_3^2-p_2^2)(m_4^2-p_4^2)}{(m_3^2-s_{12})(m_4^2-s_{23})} \right) + 2 \ln \left( \frac{m_3^2-s_{12}}{\mu^2} \right) \ln \left( \frac{m_4^2-s_{23}}{\mu^2} \right) \\
   * - \ln^2 \left( \frac{m_3^2-p_2^2}{\mu^2} \right) -ln^2 \left( \frac{m_4^2-p_4^2}{\mu^2} \right) + \ln \left( \frac{m_3^2-p_2^2}{m_4^2-s_{23}} \right) \ln \left( \frac{m_3^2}{\mu^2} \right) + \ln \left( \frac{m_4^2-p_4^2}{m_3^2 - s_{12}} \right) \ln \left( \frac{m_4^2}{\mu^2} \right) \\
   * -\frac{1}{2} \ln^2 \left( \frac{\gamma_{34}^+}{\gamma_{34}^+-1} \right) -\frac{1}{2} \ln^2 \left( \frac{\gamma_{34}^-}{\gamma_{34}^--1} \right) \right] + \mathcal{O}(\epsilon)
   * \f]
   * Implementation of the formulae from Ellis et al. \cite Ellis:2007qk.
   *
   * \param res output object res[0,1,2] the coefficients in the Laurent series
   * \param mu2 the energy scale squared.
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Box<TOutput,TMass,TScale>::B13(vector<TOutput> &res, TMass const (&Y)[4][4], TScale const& mu2) const
  {
    const TMass m3sq = Y[2][2];
    const TMass m4sq = Y[3][3];
    const TMass sibar = this->_two*Y[0][2];
    const TMass tabar = this->_two*Y[1][3];
    const TMass m4sqbar = this->_two*Y[0][3];
    const TMass m3sqbar = this->_two*Y[1][2];
    const TMass p3sq = -(this->_two*Y[2][3]-Y[2][2]-Y[3][3]);

    const TOutput fac = TOutput(sibar*tabar - m3sqbar*m4sqbar);
    const TOutput wlogsmu = this->Lnrat(sibar, TMass(mu2));
    const TOutput wlogtmu = this->Lnrat(tabar, TMass(mu2));
    const TOutput wlog3mu = this->Lnrat(m3sqbar, TMass(mu2));
    const TOutput wlog4mu = this->Lnrat(m4sqbar, TMass(mu2));
    const TOutput dilog1 = this->Li2omrat(m3sqbar, sibar);
    const TOutput dilog4 = this->Li2omrat(m4sqbar, tabar);
    const TOutput dilog7 = this->Li2omx2(m3sqbar, m4sqbar, sibar, tabar);

    TOutput root, ga34p, ga34pm1, ga34m, ga34mm1;
    TOutput ga43p, ga43pm1, ga43m, ga43mm1;
    TMass x34p, x34pm1, x34m, x34mm1;
    TMass x43p, x43pm1, x43m, x43mm1;
    if (this->iszero(p3sq))
      {
        root = this->_cone;
        x34p = -this->_one;
        x34pm1 = -this->_one;
        x34m = m4sq;
        x34mm1 = m3sq;

        x43p = m3sq;
        x43pm1 = m4sq;
        x43m = -this->_one;
        x43mm1 = -this->_one;
      }
    else
      {
        root = Sqrt(TOutput(Pow(p3sq-m3sq+m4sq,2)-this->_four*m4sq*p3sq));
        ga34p   = TOutput(+p3sq+m4sq-m3sq)+root;
        ga34pm1 = TOutput(-p3sq+m4sq-m3sq)+root;
        ga34m   = TOutput(+p3sq+m4sq-m3sq)-root;
        ga34mm1 = TOutput(-p3sq+m4sq-m3sq)-root;

        ga43p   = TOutput(+p3sq+m3sq-m4sq)+root;
        ga43pm1 = TOutput(-p3sq+m3sq-m4sq)+root;
        ga43m   = TOutput(+p3sq+m3sq-m4sq)-root;
        ga43mm1 = TOutput(-p3sq+m3sq-m4sq)-root;

        x34p = -Real(ga34p);
        x34pm1 = -Real(ga34pm1);
        x34m = Real(ga34m);
        x34mm1 = Real(ga34mm1);

        x43p = -Real(ga43p);
        x43pm1 = -Real(ga43pm1);
        x43m = Real(ga43m);
        x43mm1 = Real(ga43mm1);
      }

    TMass rat3t, rat4s;
    TOutput ln43p, ln43m, dilog2, dilog3, dilog5, dilog6;
    TOutput zrat3t, zrat4s, rat34p, rat34m, rat43p, rat43m;
    TScale ieps3t, ieps4s, ieps34, ieps43;
    if (this->iszero(Imag(root)))
      {
        ln43p = this->Lnrat(x43p, x43pm1);
        ln43m = this->Lnrat(x43m, x43mm1);

        dilog2 = this->Li2omx2(m3sqbar, x34p, tabar, x34pm1);
        dilog3 = this->Li2omx2(m3sqbar, x34m, tabar, x34mm1);
        dilog5 = this->Li2omx2(m4sqbar, x43p, sibar, x43pm1);
        dilog6 = this->Li2omx2(m4sqbar, x43m, sibar, x43mm1);
      }
    else
      {
        this->ratreal(m3sqbar, tabar, rat3t, ieps3t);
        this->ratreal(m4sqbar, sibar, rat4s, ieps4s);

        this->ratgam(rat34p, rat34m, ieps34, p3sq, m3sq, m4sq);
        this->ratgam(rat43p, rat43m, ieps43, p3sq, m4sq, m3sq);

        zrat3t = TOutput(rat3t);
        zrat4s = TOutput(rat4s);

        dilog2 = this->spencer(zrat3t, rat34p, ieps3t, ieps34);
        dilog3 = this->spencer(zrat3t, rat34m, ieps3t, ieps34);
        dilog5 = this->spencer(zrat4s, rat43p, ieps4s, ieps43);
        dilog6 = this->spencer(zrat4s, rat43m, ieps4s, ieps43);

        ln43p = this->cLn(rat43p, this->_zero);
        ln43m = this->cLn(rat43m, this->_zero);
      }


    res[2] = this->_czero;
    res[1] = wlog3mu+wlog4mu-wlogsmu-wlogtmu;
    res[0] =
        -this->_ctwo*dilog1 - dilog2 - dilog3
        -this->_ctwo*dilog4 - dilog5 - dilog6
        +this->_ctwo*dilog7
        +this->_ctwo*wlogsmu*wlogtmu-wlog3mu*wlog3mu-wlog4mu*wlog4mu
        +(wlog3mu-wlogtmu)*Log(m3sq/mu2)
        +(wlog4mu-wlogsmu)*Log(m4sq/mu2)
        -this->_chalf*(ln43p*ln43p+ln43m*ln43m);

    for (size_t i = 0; i < 3; i++)
      res[i] /= fac;
  }

  /*!
   * The integral is defined as:
   * \f[
   * I_4^{D=4-2\epsilon}(m_2^2,m_2^2,m_4^2,m_4^2;t,s;0,m_2^2,0,m_4^2) = \frac{-2}{m_2 m_4 t} \frac{x_s \ln x_s}{1-x_s^2} \left[ \frac{1}{\epsilon} + \ln \left( \frac{\mu^2}{-t} \right) \right],\, s-(m_2-m_4)^2 \neq 0 \\
   *  = \frac{1}{m_2 m_4 t} \left[ \frac{1}{\epsilon} + \ln \left( \frac{\mu^2}{-t} \right) \right], s-(m_2-m_4)^2 = 0.
   * \f]
   * Implementation of the formulae from Beenakker et al. \cite Beenakker:1988jr.
   *
   * \param res output object res[0,1,2] the coefficients in the Laurent series
   * \param mu2 the energy scale squared.
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Box<TOutput,TMass,TScale>::B14(vector<TOutput> &res, TMass const (&Y)[4][4], TScale const& mu2) const
  {
    const TMass m2sq = Y[1][1];
    const TMass m4sq = Y[3][3];
    const TMass ta = this->_two*Y[0][2];
    const TMass si = this->_two*Y[1][3]-Y[1][1]-Y[3][3];
    const TMass m2 = Sqrt(m2sq);
    const TMass m4 = Sqrt(m4sq);

    const TOutput wlogtmu = this->Lnrat(mu2, ta);

    TScale ieps = 0;
    TOutput cxs[3];
    this->kfn(cxs, ieps, -si, m2, m4);
    const TScale xs = Real(cxs[0]);
    const TScale imxs = Imag(cxs[0]);

    TOutput fac;
    if ( this->iszero(xs-this->_one) && this->iszero(imxs))
      fac = TOutput(-xs/(m2*m4*ta));
    else
      {
        const TOutput xlog = this->cLn(cxs[0], ieps);
        fac = TOutput(this->_two/(m2*m4*ta))*cxs[0]/(cxs[1]*cxs[2])*xlog;
      }
    res[2] = this->_czero;
    res[1] = fac;
    res[0] = fac*wlogtmu;
  }

  /*!
   * The integral is defined as:
   * \f[
   * I_4^{D=4-2\epsilon}(m_2^2,p_2^2,p_3^2,m_4^2;t,s;0,m_2^2,0,m_4^2) = \\
   *  \frac{x_s}{m_2 m_4 t (1-x_s^2)} \left\{ \ln x_s \left[ -\frac{1}{\epsilon} - \frac{1}{2} \ln x_s - \ln \left( \frac{\mu^2}{m_2 m_4} \right) - \ln \left( \frac{m_2^2-p_2^2}{-t} \right) - \ln \left( \frac{m_4^2-p_3^2}{-t} \right) \right] \\
   *  - {\rm Li}_2 (1-x_s^2) + \frac{1}{2} ln^2 y + \sum_{\rho=\pm1} {\rm Li}_2 (1-x_s y^\rho) \right\}
   * \f]
   * Implementation of the formulae from Beenakker et al. \cite Beenakker:1988jr.
   *
   * \param res output object res[0,1,2] the coefficients in the Laurent series
   * \param mu2 the energy scale squared.
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Box<TOutput,TMass,TScale>::B15(vector<TOutput> &res, TMass const (&Y)[4][4], TScale const& mu2) const
  {
    const TMass m2sq = Y[1][1];
    const TMass m4sq = Y[3][3];
    const TMass m2sqbar = this->_two*Y[1][2];
    const TMass m4sqbar = this->_two*Y[2][3];
    const TMass si = this->_two*Y[1][3]-Y[1][1]-Y[3][3];
    const TMass ta = this->_two*Y[0][2];
    const TMass m2 = Sqrt(m2sq);
    const TMass m4 = Sqrt(m4sq);

    TScale ieps = 0;
    TOutput cxs[3];
    this->kfn(cxs, ieps, -si, m2, m4);
    const TOutput xs = cxs[0];

    TOutput fac;
    if (this->iszero(m2sqbar) && !this->iszero(m4sqbar))
      {
        TMass yi;
        TScale iepyi;
        this->ratreal(m4*m2sqbar, m2*m4sqbar, yi, iepyi);
        TOutput cyi = TOutput(yi);
        fac = xs/(this->_cone-xs*xs)/TOutput(-m2*m4*ta);
        const TOutput xlog = this->cLn(xs, ieps);
        res[2] = this->_czero;
        res[1] = -xlog;
        res[0] = xlog*(-xlog-TOutput(Log(mu2/m4sq))
            -this->_ctwo*this->Lnrat(m4sqbar, ta))
            -this->cLi2omx2(xs,xs,ieps,ieps)
            +this->cLi2omx2(xs, cyi, ieps, iepyi)
            -this->cLi2omx2(this->_cone/xs, cyi, -ieps, iepyi);

        for (size_t i = 0; i < res.size(); i++)
          res[i] *= TOutput(fac);

        return;
      }
    else if (this->iszero(m4sqbar) && !this->iszero(m2sqbar))
      {
        TMass yy;
        TScale iepsyy;
        this->ratreal(m2*m4sqbar, m4*m2sqbar, yy, iepsyy);
        TOutput cyy = TOutput(yy);
        fac = xs/(this->_cone-xs*xs)/TOutput(-m2*m4*ta);
        const TOutput xlog = this->cLn(xs, ieps);
        res[2] = this->_czero;
        res[1] = -xlog;
        res[0] = xlog*(-xlog-TOutput(Log(mu2/m2sq))
                       -this->_ctwo*this->Lnrat(m2sqbar, ta))
            -this->cLi2omx2(xs,xs,ieps,ieps)
            +this->cLi2omx2(xs,cyy,ieps,iepsyy)
            -this->cLi2omx2(this->_cone/xs, cyy, -ieps, iepsyy);

        for (size_t i = 0; i < res.size(); i++)
          res[i] *= TOutput(fac);

        return;
      }
    else if (this->iszero(m4sqbar) && this->iszero(m2sqbar))
      throw RangeError("Box::B15","wrong kinematics, this is really B14.");

    TMass yy;
    TScale iepsyy;
    this->ratreal(m2*m4sqbar, m4*m2sqbar, yy, iepsyy);
    const TScale rexs = Real(xs);
    const TScale imxs = Imag(xs);

    if (this->iszero(rexs-this->_one) && this->iszero(imxs))
      {
        fac = TOutput(this->_chalf/(m2*m4*ta));
        res[2] = this->_czero;
        res[1] = this->_cone;
        res[0] = TOutput(Log(mu2/(m2*m4)))
            -this->Lnrat(m2sqbar, ta)-this->Lnrat(m4sqbar, ta)-this->_ctwo
            -TOutput((this->_one+yy)/(this->_one-yy))*this->Lnrat(m2*m4sqbar, m4*m2sqbar);
      }
    else
      {
        fac = xs/(this->_cone-xs*xs)/TOutput(-m2*m4*ta);
        const TOutput xlog = this->cLn(xs, ieps);
        res[2] = this->_czero;
        res[1] = -xlog;
        res[0] = xlog*(-this->_chalf*xlog-TOutput(Log(mu2/(m2*m4)))
                       -this->Lnrat(m2sqbar,ta)-this->Lnrat(m4sqbar, ta))
            -this->cLi2omx2(xs,xs,ieps,ieps)
            +this->_chalf*Pow(this->Lnrat(m2*m4sqbar, m4*m2sqbar),2)
            +this->cLi2omx2(xs,TOutput(yy), ieps, iepsyy)
            +this->cLi2omx2(xs,TOutput(this->_one/yy), ieps, -iepsyy);
      }

    for (size_t i = 0; i < res.size(); i++)
      res[i] *= TOutput(fac);
  }

  /*!
   * The integral is defined as:
   * \f[
   * I_4^{D=4-2\epsilon}(m_2^2,p_2^2,p_3^2,m_4^2;t,s;0,m_2^2,m_3^2,m_4^2) = frac{x_s}{m_2 m_4 (t-m_3^2)(t-x_s^2)} \\
   *  \left\{ - \frac{\ln x_s}{ \epsilon} - 2 \ln x_s \ln \left( \frac{m_3 \mu}{m_3^2-t} \right) + \ln^2 x_2 + \ln^2 x_3 - {\rm Li}_2 (1-x_s^2) \\
   *  {\rm Li}_2 (1-x_s x_2 x_3) + {\rm Li}_2 \left( 1- \frac{x_s}{x_2 x_3} \right) + {\rm Li}_2 \left( 1- \frac{x_s x_2}{x_3} \right) + {\rm Li}_2 \left( 1- \frac{x_s x_3}{x_2} \right) \right\}
   * \f]
   * Implementation of the formulae from Beenakker et al. \cite Beenakker:1988jr.
   *
   * \param res output object res[0,1,2] the coefficients in the Laurent series
   * \param mu2 the energy scale squared.
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Box<TOutput,TMass,TScale>::B16(vector<TOutput> &res, TMass const (&Y)[4][4], TScale const& mu2) const
  {
    const TMass m2sq = Y[1][1];
    const TMass m3sq = Y[2][2];
    const TMass m4sq = Y[3][3];
    const TMass tabar = this->_two*Y[0][2];
    const TMass si = this->_two*Y[1][3]-Y[1][1]-Y[3][3];
    const TMass mp2sq = this->_two*Y[1][2]-m3sq-m2sq;
    const TMass mp3sq = this->_two*Y[2][3]-m3sq-m4sq;
    const TMass m2 = Sqrt(m2sq);
    const TMass m3 = Sqrt(m3sq);
    const TMass m4 = Sqrt(m4sq);
    const TMass mean = Sqrt(Real(m3sq)*mu2);

    TScale ieps = 0, iep2 = 0, iep3 = 0;
    TOutput cxs[3], cx2[3], cx3[3];
    this->kfn(cxs, ieps, -si, m2, m4);
    this->kfn(cx2, iep2, -mp2sq, m2, m3);
    this->kfn(cx3, iep3, -mp3sq, m3, m4);

    const TOutput xs = cxs[0];
    const TScale rexs = Real(xs);
    const TScale imxs = Imag(xs);

    TOutput fac;
    if (this->iszero(rexs-this->_one) && this->iszero(imxs))
      {
        fac = TOutput(-this->_half/(m2*m4*tabar));
        res[2] = this->_czero;
        res[1] = this->_cone;
        res[0] = this->_two*this->Lnrat(mean, tabar)-this->_ctwo;

        if (this->iszero(Real(cx2[0]-cx3[0])) && this->iszero(Imag(cx2[0]-cx3[0]))
            && this->iszero(Real(cx2[0]-this->_one)) && this->iszero(Imag(cx2[0])))
          res[0] += this->_cfour;
        else if (this->iszero(Real(cx2[0]-cx3[0])) && this->iszero(Imag(cx2[0]-cx3[0])))
          res[0] += this->_ctwo + this->_ctwo*(cx2[0]*cx2[0] + this->_cone)*this->cLn(cx2[0], iep2)/(cx2[0]*cx2[0]-this->_cone);
        else
          {
            const TOutput cln_cx2_iep2 = this->cLn(cx2[0], iep2);
            const TOutput cln_cx3_iep3 = this->cLn(cx3[0], iep3);
            res[0] += -(this->_cone+cx2[0]*cx3[0])/(this->_cone-cx2[0]*cx3[0])
                *(cln_cx2_iep2 + cln_cx3_iep3)
                -(this->_cone+cx2[0]/cx3[0])/(this->_cone-cx2[0]/cx3[0])
                *(cln_cx2_iep2 - cln_cx3_iep3);
          }
      }
    else
      {
        fac = TOutput(-this->_one/(m2*m4*tabar))*cxs[0]/(this->_cone-cxs[0]*cxs[0]);
        const TOutput xlog = this->cLn(xs, ieps);
        const TOutput cln_cx2_iep2 = this->cLn(cx2[0], iep2);
        const TOutput cln_cx3_iep3 = this->cLn(cx3[0], iep3);
        res[2] = this->_czero;
        res[1] = -xlog;
        res[0] = -this->_ctwo*xlog*this->Lnrat(mean, tabar)
            +cln_cx2_iep2*cln_cx2_iep2 + cln_cx3_iep3*cln_cx3_iep3
            -this->cLi2omx2(xs,xs,ieps,ieps)
            +this->cLi2omx3(cxs[0],cx2[0],cx3[0],ieps,iep2,iep3)
            +this->cLi2omx3(cxs[0],this->_cone/cx2[0], this->_cone/cx3[0], ieps,-iep2,-iep3)
            +this->cLi2omx3(cxs[0],cx2[0],this->_cone/cx3[0], ieps, iep2,-iep3)
            +this->cLi2omx3(cxs[0],this->_cone/cx2[0], cx3[0], ieps, -iep2, iep3);            
      }

    for (size_t i = 0; i < 3; i++)
      res[i] *= TOutput(fac);
  }

  // explicity tyoename declaration
  template class Box<complex,double,double>;
  template class Box<complex,complex,double>;
  template class Box<qcomplex,qdouble,qdouble>;
  template class Box<qcomplex,qcomplex,qdouble>;
}
