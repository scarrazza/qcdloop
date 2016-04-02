//
// QCDLoop 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch
//          Keith Ellis: keith.ellis@durham.ac.uk
//          Giulia Zanderighi: giulia.zanderighi@cern.ch

#pragma once

#include <vector>
#include <memory>
#include "types.h"
#include "tadpole.h"
#include "bubble.h"
#include "triangle.h"
#include "box.h"
#include "maths.h"
#include "tools.h"
#include "timer.h"
using std::vector;

namespace ql
{
  /**
   * @brief The one-loop scalar integral trigger.
   *
   * Parses automatically the topology and computes the integral
   */
  template<typename TOutput = complex, typename TMass = double, typename TScale = double>
  class QCDLoop
  {
  public:
    QCDLoop();  //!< The Constructor.
    ~QCDLoop(); //!< The Destructor.

    //! Computes the one-loop scalar integral
    void integral(vector<TOutput> &res,
                  TScale const& mu2,
                  vector<TMass> const& m,
                  vector<TScale> const& p = {}) const;

    //! Set the Cache size
    void setCacheSize(int const& size);

  private:
    TadPole<TOutput,TMass,TScale>*  _tp;
    Bubble<TOutput,TMass,TScale>*   _bb;
    Triangle<TOutput,TMass,TScale>* _tr;
    Box<TOutput,TMass,TScale>*      _tb;
  };
}

