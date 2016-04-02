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
   * @brief The Bubble integral.
   *
   * Parses automatically the topology and computes the integral
   */
  template<typename TOutput = complex, typename TMass = double, typename TScale = double>
  class Bubble: public Topology<TOutput,TMass,TScale>
  {
  public:
    Bubble();  //!< The Constructor.
    ~Bubble(); //!< The Destructor.

    //! Computes the Bubble integral automatically
    void integral(vector<TOutput> &res, TScale const& mu2, vector<TMass> const& m, vector<TScale> const& p);

    //! General configuration I(s;m0,m1)
    void BB0(vector<TOutput> &res, TScale const& mu2, TMass const& m0, TMass const& m1, TScale const& s) const;

    //! Special configuration I(m2;0,m2)
    void BB1(vector<TOutput> &res, TScale const& mu2, TMass const& m) const;

    //! Special configuration I(0;0,m2)
    void BB2(vector<TOutput> &res, TScale const& mu2, TMass const& m) const;

    //! Special configuration I(s;0,0)
    void BB3(vector<TOutput> &res, TScale const& mu2, TMass const& s) const;

    //! Special configuration I(s;0,m2)
    void BB4(vector<TOutput> &res, TScale const& mu2, TMass const& m, TScale const& s) const;

    //! Special configuration I(0;m0,m1)
    void BB5(vector<TOutput> &res, TScale const& mu2, TMass const& m0, TMass const& m1) const;
    };
}
