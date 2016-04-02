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
   * @brief The TadPole integral.
   *
   * Parses automatically the topology and computes the integral
   */
  template<typename TOutput = complex, typename TMass = double, typename TScale = double>
  class TadPole: public Topology<TOutput,TMass,TScale>
  {
  public:
    TadPole();  //!< The Constructor.
    ~TadPole(); //!< The Destructor.

    //! Computes the tadpole integral
    void integral(vector<TOutput> &res, TScale const& mu2, vector<TMass> const& m, vector<TScale> const& p = {});
  };
}
