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
   * @brief The Box integral.
   *
   * Parses automatically the topology and computes the integral
   */
  template<typename TOutput = complex, typename TMass = double, typename TScale = double>
  class Box: public Topology<TOutput,TMass,TScale>
  {
  public:
    Box();  //!< The Constructor.
    ~Box(); //!< The Destructor.

    //! Computes the tadpole integral
    void integral(vector<TOutput> &res, TScale const& mu2, vector<TMass> const& m, vector<TScale> const& p);

    //! Divergent box I(0,0,0,0;s12,s23;0,0,0,0)
    void B1(vector<TOutput> &res, TMass const (&Y)[4][4], TScale const& mu2) const;

    //! Divergent box I(0,0,0,p2;s12,s23;0,0,0,0)
    void B2(vector<TOutput> &res, TMass const (&Y)[4][4], TScale const& mu2) const;

    //! Divergent box I(0,pd2,0,pq2;s12,s23;0,0,0,0)
    void B3(vector<TOutput> &res, TMass const (&Y)[4][4], TScale const& mu2) const;

    //! Divergent box I(0,0,pt2,pq2;s12,s23;0,0,0,0)
    void B4(vector<TOutput> &res, TMass const (&Y)[4][4], TScale const& mu2) const;

    //! Divergent box I(0,pd2,pt2,pq2;s12,s23;0,0,0,0)
    void B5(vector<TOutput> &res, TMass const (&Y)[4][4], TScale const& mu2) const;

    //! Divergent box I(0,0,m2,m2;s12,s23;0,0,0,m2)
    void B6(vector<TOutput> &res, TMass const (&Y)[4][4], TScale const& mu2) const;

    //! Divergent box I(0,0,m2,pq2;s12,s23;0,0,0,m2)
    void B7(vector<TOutput> &res, TMass const (&Y)[4][4], TScale const& mu2) const;

    //! Divergent box I(0,0,pt2,pq2;s12,s23;0,0,0,m2)
    void B8(vector<TOutput> &res, TMass const (&Y)[4][4], TScale const& mu2) const;

    //! Divergent box I(0,p2,p3,m2;s12,s23;0,0,0,m2)
    void B9(vector<TOutput> &res, TMass const (&Y)[4][4], TScale const& mu2) const;

    //! Divergent box I(0,p2,p3,p4;s12,s23;0,0,0,m2)
    void B10(vector<TOutput> &res, TMass const (&Y)[4][4], TScale const& mu2) const;

    //! Divergent box I(0,m3,pt2,m4;s12,s23;0,0,m3,m4)
    void B11(vector<TOutput> &res, TMass const (&Y)[4][4], TScale const& mu2) const;

    //! Divergent box I(0,m3,pt2,pq2;s12,s23;0,0,m3,m4)
    void B12(vector<TOutput> &res, TMass const (&Y)[4][4], TScale const& mu2) const;

    //! Divergent box I(0,p2,p3,p4;s12,s23;0,0,m3,m4)
    void B13(vector<TOutput> &res, TMass const (&Y)[4][4], TScale const& mu2) const;

    //! Divergent box I(m2,m2,m4,m4;s12,s23;0,m2,0,m4)
    void B14(vector<TOutput> &res, TMass const (&Y)[4][4], TScale const& mu2) const;

    //! Divergent box I(m2,p2,p3,m4;s12,s23;0,m2,0,m4)
    void B15(vector<TOutput> &res, TMass const (&Y)[4][4], TScale const& mu2) const;

    //! Divergent box I(m2,p2,p3,m4;s12,s23;0,m2,m3,m4)
    void B16(vector<TOutput> &res, TMass const (&Y)[4][4], TScale const& mu2) const;

  private:  
    //! Calculate Y from xpi
    void Ycalc(TMass (&Y)[4][4],TMass (&Yalt)[4][4], int const& massive, bool const& opposite = false) const;

    //! Trigger for 4-offshell external lines
    void B0m(vector<TOutput> &res, TMass const (&xpi)[13], TScale const& mu2) const;

    //! Trigger for 3-offshell external lines
    void B1m(vector<TOutput> &res, TMass const (&xpi)[13], TScale const& mu2) const;

    //! Trigger for 2-offshell external lines
    void B2m(vector<TOutput> &res, TMass const (&xpi)[13], TScale const& mu2) const;
    void B2mo(vector<TOutput> &res, TMass const (&xpi)[13], TScale const& mu2) const;
    void B2ma(vector<TOutput> &res, TMass const (&xpi)[13], TScale const& mu2) const;

    //! Trigger for 1-offshell external lines
    void B3m(vector<TOutput> &res, TMass const (&xpi)[13], TScale const& mu2) const;

    //! Trigger for 0-offshell external lines
    void B4m(vector<TOutput> &res, TMass const (&xpi)[13]) const;

    //! Finite box, 4-offshell external lines
    void BIN0(vector<TOutput> &res, TMass const (&Y)[4][4]) const;

    //! Finite box, 3-offshell external lines
    void BIN1(vector<TOutput> &res, TMass const (&Y)[4][4]) const;

    //! Finite box, 2-offshell external lines
    void BIN2(vector<TOutput> &res, TMass const (&Y)[4][4]) const;

    //! Finite box, 1-offshell external lines
    void BIN3(vector<TOutput> &res, TMass const (&Y)[4][4]) const;

    //! Finite box, 0-offshell external lines
    void BIN4(vector<TOutput> &res, TMass const (&Y)[4][4]) const;
  };
}
