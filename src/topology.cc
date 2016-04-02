//
// QCDLoop 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch
//          Keith Ellis: keith.ellis@durham.ac.uk
//          Giulia Zanderighi: giulia.zanderighi@cern.ch

#include "qcdloop/topology.h"
#include "qcdloop/tools.h"
#include "qcdloop/exceptions.h"

namespace ql {

  /*!
   * \brief Topology<TOutput, TMass>::Topology
   * \param name
   */
  template<typename TOutput, typename TMass, typename TScale>
  Topology<TOutput,TMass,TScale>::Topology(string name):
    Tools<TOutput,TMass,TScale>(),
    LRUCache<size_t,vector<TOutput> >(1),
    _name(name),
    _key(),
    _mu2(this->_zero),
    _ch(new ContainerHasher<TMass, TScale>())
  {
    _val.resize(3);
  }

  /*!
   * \brief Topology<TOutput, TMass>::Topology
   * \param obj
   */
  template<typename TOutput, typename TMass, typename TScale>
  Topology<TOutput,TMass,TScale>::Topology(const Topology &obj):
    Tools<TOutput,TMass,TScale>(obj),
    LRUCache<size_t,vector<TOutput> >(obj),
    _name(obj._name),
    _key(obj._key),
    _mu2(obj._mu2),
    _m(obj._m),
    _p(obj._p),
    _val(obj._val),
    _ch(obj._ch)
  {
  }

  /*!
   * \brief Topology<TOutput, TMass>::~Topology
   */
  template<typename TOutput, typename TMass, typename TScale>
  Topology<TOutput,TMass,TScale>::~Topology()
  {
    _val.clear();
    _m.clear();
    _p.clear();
    delete _ch;
  }

  /*!
   * \brief Topology<TOutput, TMass>::checkCache
   * \param mu2
   * \param m
   * \param p
   * \return true if the cache already contains the cofiguration
   */
  template<typename TOutput, typename TMass, typename TScale>
  bool Topology<TOutput,TMass,TScale>::checkCache(const TScale &mu2, const vector<TMass> &m, const vector<TScale> &p)
  {
    const int cz = this->getCacheSize();
    if (cz == 1)
      {
        if (_mu2 == mu2 && _m == m && _p == p) return true;
      }
    else if (cz > 1)
      {
        _key = _ch->genkey(mu2, m, p);
        if (this->get(_key,_val)) return true;
      }

    return false;
  }

  /*!
   * \brief Topology<TOutput, TMass>::storeCache
   */
  template<typename TOutput, typename TMass, typename TScale>
  void Topology<TOutput,TMass,TScale>::storeCache(const TScale &mu2, const vector<TMass> &m, const vector<TScale> &p)
  {    
    const int cz = this->getCacheSize();
    if (cz == 1)
      {
        _mu2 = mu2;
        std::copy(m.begin(), m.end(), _m.begin());
        std::copy(p.begin(), p.end(), _p.begin());
      }
    else if (cz > 1) this->store(_key,_val);
  }

  // explicity template declaration
  template class Topology<complex,double,double>;
  template class Topology<complex,complex,double>;
  template class Topology<qcomplex,qdouble,qdouble>;
  template class Topology<qcomplex,qcomplex,qdouble>;
}
