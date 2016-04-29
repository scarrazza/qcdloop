//
// QCDLoop 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch
//          Keith Ellis: keith.ellis@durham.ac.uk
//          Giulia Zanderighi: giulia.zanderighi@cern.ch

#pragma once

#include "topology.h"

namespace ql
{
  /**
   * @brief The Triangle integral.
   *
   * Parses automatically the topology and computes the integral
   */
  template<typename TOutput = complex, typename TMass = double, typename TScale = double>
  class Triangle: public Topology<TOutput,TMass,TScale>
  {
  public:
    Triangle();  //!< The Constructor.
    ~Triangle(); //!< The Destructor.

    //! Computes the Triangle integral
    void integral(vector<TOutput> &res, TScale const& mu2, vector<TMass> const& m, vector<TScale> const& p);

    //! General case triangle integral I(p1,p2,p3;m1,m2,m3)
    void T0(vector<TOutput> &res, TMass const (&xpi)[6], int const& massive) const;

    //! Divergent triangle integral I(0,0,p2;0,0,0)
    void T1(vector<TOutput> &res, TScale const& mu2, TScale const& p) const;

    //! Divergent triangle integral I(0,p1,p2;0,0,0)
    void T2(vector<TOutput> &res, TScale const& mu2, TScale const& p2, TScale const& p3) const;

    //! Divergent triangle integral I(0,p1,p2;0,0,m2)
    void T3(vector<TOutput> &res, TScale const& mu2, TMass const& m, TScale const& p1, TScale const& p2) const;

    //! Divergent triangle integral I(0,m2,p2;0,0,m2)
    void T4(vector<TOutput> &res, TScale const& mu2, TMass const& m, TScale const& p2) const;

    //! Divergent triangle integral I(0,m2,m2;0,0,m2)
    void T5(vector<TOutput> &res, TScale const& mu2, TMass const& m) const;

    //! Divergent triangle integral I(m2,s,m3;0,m2,m3)
    void T6(vector<TOutput> &res, TScale const& mu2, TMass const& m2, TMass const& m3, TScale const& p2) const;

  private:
    //! Sort arguments of triangle, with msq in ascending order
    void TriSort(TScale (&psq)[3], TMass (&msq)[3]) const;

    //! Sort arguments of triangle, with msq in ascending order
    void TriSort2(TMass const (&xpi)[6], TMass (&ypi)[6]) const;

    //! Sort on the basis of abs value.
    void SnglSort(TScale (&psq)[3]) const;

    //! Calculation of triangle with p1sq=p2sq=p3sq=0
    void TIN0(TOutput &res, TMass const (&xpi)[6]) const;

    //! Calculation of triangle with p1sq=p2sq=0
    void TIN1(TOutput &res, TMass const (&xpi)[6], TMass const (&sxpi)[6], int const& massive) const;

    //! Calculation of triangle with p1sq=0
    void TIN2(TOutput &res, TMass const (&xpi)[6], TMass const (&sxpi)[6], int const& massive) const;

    //! Calculation of triangle with all p non-zero
    void TIN3(TOutput &res, TMass const (&xpi)[6], TMass const (&sxpi)[6], int const& massive) const;

    //! Calculation of triangle with all p and complex m non-zero
    void TINDNS(TOutput &res, TMass const (&xpi)[6]) const;

    //! Calculation of triangle with 1 massive particles
    void TINDNS1(TOutput &res, TMass const (&xpi)[6]) const;

    //! Calculation of triangle with 2 massive particles
    void TINDNS2(TOutput &res, TMass const (&xpi)[6]) const;

    //! Calculate the kallen formula
    TOutput Kallen(TOutput const& p1, TOutput const& p2, TOutput const& p3) const;

    //! Calculate the kallen^2 formula
    TOutput Kallen2(TOutput const& p1, TOutput const& p2, TOutput const& p3) const;

  };
}
