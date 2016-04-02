//
// QCDLoop 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch
//          Keith Ellis: keith.ellis@durham.ac.uk
//          Giulia Zanderighi: giulia.zanderighi@cern.ch

#include "qcdloop/cache.h"


namespace std {

  template <typename T>
  inline void hash_combine(std::size_t & seed, const T & v)
  {
    hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  }

  template <>
  struct hash<ql::qdouble> : public __hash_base<size_t, ql::qdouble>
  {
    inline size_t operator()(const ql::qdouble & x) const noexcept
    {
      return x != 0.0q ? std::_Hash_impl::hash(x) : 0;
    }
  };

  template <> struct hash<ql::complex> : public __hash_base<size_t, ql::complex>
  {
    inline size_t operator()(const ql::complex & x) const noexcept
    {
      return x != ql::complex{0.0,0.0} ? std::_Hash_impl::hash(x) : 0;
    }
  };

  template <> struct hash<ql::qcomplex> : public __hash_base<size_t, ql::qcomplex>
  {
    inline size_t operator()(const ql::qcomplex & x) const noexcept
    {
      return x != ql::qcomplex{0.0q,0.0q} ? std::_Hash_impl::hash(x) : 0;
    }
  };
}


namespace ql
{
  template<typename TMass, typename TScale>
  ContainerHasher<TMass,TScale>::ContainerHasher()
  {
  }

  template<typename TMass, typename TScale>
  size_t ContainerHasher<TMass,TScale>::genkey(TScale const& mu2,
                                               vector<TMass> const& m,
                                               vector<TScale> const& p) const
  {
    size_t seed = 0;
    std::hash_combine(seed, mu2);
    for (typename std::vector<TMass>::const_iterator it = m.begin(), end = m.end(); it != end; ++it)
      std::hash_combine(seed, *it);
    for (typename std::vector<TScale>::const_iterator it = p.begin(), end = p.end(); it != end; ++it)
      std::hash_combine(seed, *it);
    return seed;
  }

  // explicity template declaration
  template class ContainerHasher<double,double>;
  template class ContainerHasher<complex,double>;
  template class ContainerHasher<qdouble,qdouble>;
  template class ContainerHasher<qcomplex,qdouble>;

  /*!
   * The LRU Cache constructor
   */
  template<typename Tkey, typename Tvalue>
  LRUCache<Tkey,Tvalue>::LRUCache(int const& size): _size(size)
  {
  }

  template<typename Tkey, typename Tvalue>
  LRUCache<Tkey,Tvalue>::LRUCache(const LRUCache &obj):
    _size(obj._size),
    _cache_list(obj._cache_list),
    _cache_map(obj._cache_map)
  {
  }

  template<typename Tkey, typename Tvalue>
  LRUCache<Tkey,Tvalue>::~LRUCache() {
    _cache_list.clear();
    _cache_map.clear();
  }

  template<typename Tkey, typename Tvalue>
  void LRUCache<Tkey,Tvalue>::setCacheSize(int const& size)
  {
    _size = size;
    _cache_list.clear();
    _cache_map.clear();
  }

  template<typename Tkey, typename Tvalue>
  void LRUCache<Tkey,Tvalue>::store(Tkey const& key, Tvalue const& value)
  {
    auto it = _cache_map.find(key);
    if (it != _cache_map.end()) {
            _cache_list.erase(it->second);
            _cache_map.erase(it);
    }

    _cache_list.push_front(key_value_pair_t(key, value));
    _cache_map[key] = _cache_list.begin();

    if ((int)_cache_map.size() > _size) {
            auto last = _cache_list.end();
            last--;
            _cache_map.erase(last->first);
            _cache_list.pop_back();
    }
  }

  template<typename Tkey, typename Tvalue>
  bool LRUCache<Tkey,Tvalue>::get(Tkey const& key, Tvalue & out)
  {
    auto it = _cache_map.find(key);
    if (it == _cache_map.end()) {
      return false;
    } else {
      _cache_list.splice(_cache_list.begin(), _cache_list, it->second);
      out = it->second->second;
      return true;
    }
  }

  // explicity template declaration
  template class LRUCache<size_t,vector<complex> >;
  template class LRUCache<size_t,vector<qcomplex> >;
}
