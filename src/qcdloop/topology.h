//
// QCDLoop 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch
//          Keith Ellis: keith.ellis@durham.ac.uk
//          Giulia Zanderighi: giulia.zanderighi@cern.ch

#pragma once

#include <vector>
#include <string>
#include "types.h"
#include "tools.h"
#include "cache.h"
using std::vector;
using std::string;

namespace ql
{  
  /**
   * @brief The Topology abstract class.
   *
   * The pure virtual class for topology extentions.
   */
  template<typename TOutput = complex, typename TMass = double, typename TScale = double>
  class Topology: public Tools<TOutput,TMass,TScale>, public LRUCache<size_t, vector<TOutput> >
  {
  public:
    Topology(string name = "None");
    Topology(const Topology &obj);
    virtual ~Topology();

    //! Pure virtual function for integral computation.
    virtual void integral(vector<TOutput>&, TScale const&,vector<TMass> const&,vector<TScale> const&) = 0;

    // Get methods
    string const& getName() const { return _name;} //!< Get topology name

  protected:
    //!< Check stored cached results
    bool checkCache(TScale const&,vector<TMass> const&,vector<TScale> const&);

    //!< Store results with LRU caching for size > 1
    void storeCache(TScale const&,vector<TMass> const&,vector<TScale> const&);

    string _name;
    size_t _key;
    TScale _mu2;
    vector<TMass> _m;
    vector<TScale> _p;
    vector<TOutput>  _val;
    ContainerHasher<TMass,TScale> *_ch;
  };

}
