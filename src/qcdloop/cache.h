//
// QCDLoop 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch
//          Keith Ellis: keith.ellis@durham.ac.uk
//          Giulia Zanderighi: giulia.zanderighi@cern.ch

#pragma once

#include <unordered_map>
#include <list>
#include <cstddef>
#include <stdexcept>
#include <vector>
#include "types.h"
using std::vector;

namespace ql {

  /**
   * @brief The Hasher algorithm container for key generation
   */
  template <typename TMass, typename TScale>
  class ContainerHasher
  {
  public:
    ContainerHasher();
    size_t genkey(TScale const& mu2, vector<TMass> const& m, vector<TScale> const& p) const;
  };

  /**
   * @brief The LRU Cache class.
   *
   * Implements the LRU algorithm for caching.
   */
  template<typename Tkey, typename Tvalue>
  class LRUCache
  {
  public:
    /*!
     * \brief LRUCache constructor
     * \param size the caching size
     */
    LRUCache(int const& size = 1);
    LRUCache(const LRUCache &obj);
    ~LRUCache();

    typedef typename std::pair<Tkey, Tvalue> key_value_pair_t;
    typedef typename std::list<key_value_pair_t>::iterator list_iterator_t;

    //! Set the Cache size
    void setCacheSize(int const& size);

    //! Get the Cache size
    int const& getCacheSize() const { return _size; }

    //! Store the cached data
    void store(Tkey const& key, Tvalue const& value);

    //! Get the cached data
    bool get(Tkey const& key, Tvalue & out);

  private:
    int _size;
    std::list<key_value_pair_t> _cache_list;
    std::unordered_map<Tkey, list_iterator_t> _cache_map;
  };
}
