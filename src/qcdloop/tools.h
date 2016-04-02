//
// QCDLoop 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch
//          Keith Ellis: keith.ellis@durham.ac.uk
//          Giulia Zanderighi: giulia.zanderighi@cern.ch

#pragma once

#include "types.h"
#include <vector>
using std::vector;

namespace ql
{
  /*!
   * \brief Standard Tools for the library.
   *
   * The container of common functions used by Topology.
   */
  template<typename TOutput = complex, typename TMass = double, typename TScale = double>
  class Tools
  {
  public:
    Tools();  //!< The Constructor.
    ~Tools(); //!< The Destructor.

    //! Check for zeros.
    bool iszero(TMass const& psq) const;

    /* Logarithm functions*/
    //! Log of complex argument with explicit sign for imag part.
    TOutput cLn(TOutput const& z, TScale const& isig) const;
    TOutput cLn(TScale  const& x, TScale const& isig) const;

    //! The fndd function.
    TOutput fndd(int const& n, TOutput const& x, TScale const& iep) const;

    //! Computes the ratio of logs.
    TOutput Lnrat(TOutput const& x, TOutput const& y) const;
    TOutput Lnrat(TScale const& x, TScale const& y) const;

    /* Polylog functions & variations */
    //! Computes the dilog function for real argument.
    TMass ddilog(TMass const& x) const;

    //! Spence's function.
    TOutput denspence(TOutput const& z, TScale const& isig) const;

    //! Spence's function for ratio arguments.
    TOutput spencer(TOutput const& zrat1, TOutput const& zrat2, TScale const&ieps1, TScale const&ieps2) const;

    //! Difference of complex Spence's function.
    TOutput xspence(TOutput const (&z1)[2], TScale const (&im1)[2], TOutput const& z2, TScale const& im2) const;

    //! Complex Spence's function.
    TOutput cspence(TOutput const& z1, TScale const& im1, TOutput const& z2, TScale const& im2) const;

    //! Expression for dilog(1-(x-i*ep)/(y-i*ep)).
    TOutput Li2omrat(TScale const& x, TScale const& y) const;
    TOutput Li2omrat(TOutput const& x, TOutput const& y, TScale const& ieps1 = -1, TScale const& ieps2 = -1) const;

    //! Expression for dilog(1-(v-i*ep)*(w-i*ep)/(x-i*ep)/(y-i*ep)).
    TOutput Li2omx2(TScale const& v, TScale const& w, TScale const& x, TScale const& y) const;
    TOutput Li2omx2(TOutput const& v, TOutput const& w, TOutput const& x, TOutput const& y, TScale const& ieps1 = -1, TScale const& ieps2 = -1) const;
    TOutput cLi2omx2(TOutput const& z1, TOutput const& z2, TScale const& ieps1 = -1, TScale const& ieps2 = -1) const;

    //! Calculate Li[2](1-(x1+ieps1)*(x2+ieps2))
    TOutput Li2omx(TMass const& x1, TMass const& x2, TScale const& ieps1, TScale const& ieps2) const;

    //! Calculate Li[2](1-(z1+ieps1)*(z2+ieps2)*(z3+ieps3))
    TOutput cLi2omx3(TOutput const& z1, TOutput const& z2, TOutput const& z3, TScale const& ieps1, TScale const& ieps2, TScale const& ieps3) const;

    /* Special Li2 functions */
    TOutput L0(TMass const& x, TMass const& y) const;
    TOutput L1(TMass const& x, TMass const& y) const;
    TOutput R3int(TOutput const& p, TOutput const& s1, TOutput const& s2, TOutput const& t1) const;
    TOutput R3int(TOutput const& p, TOutput const& s1, TOutput const& s2, TOutput const& t1, TOutput const& t2, TOutput const& t3, TOutput const& t4) const;
    TOutput R2int(TOutput const& a, TOutput const& b, TOutput const& y0) const;
    TOutput Rint(TOutput const& y, TOutput const& z, TScale const& ieps) const;
    void    R(TOutput &r, TOutput &d, TOutput const& q) const;
    TOutput Zlogint(TOutput const& z, TScale const& ieps) const;

    /* Auxiliary Li2 functions */
    TOutput ltspence(int const& i_in, TOutput const& z_in, TScale const& s) const;
    TOutput li2series(TOutput const& z, TScale const& isig) const;
    TOutput ltli2series(TOutput const& z1, TScale const& s) const;

    /* eta functions */
    TOutput eta2(TOutput const& a,TOutput const& b) const;
    TOutput eta3(TOutput const& a,TOutput const& b, TOutput const& c) const;
    TOutput eta5(TOutput const& a,TOutput const& b, TOutput const& c, TOutput const& d, TOutput const& e) const;
    TOutput xetatilde(TOutput const (&z1)[2], TScale const (&im1)[2], TOutput const& z2, TScale const& im2, TOutput const (&l1)[2]) const;
    TOutput xeta(TOutput const (&z1)[2], TScale const (&im1)[2], TOutput const& z2, TScale const& im2, TScale const& im12, TOutput const (&l1)[2]) const;
    int eta(TOutput const& z1, TScale const& s1, TOutput const& z2, TScale const& s2, TScale const& s12) const;
    int etatilde(TOutput const& c1, TScale const& im1x, TOutput const& c2, TScale const& im2x) const;

    /* Extra functions */
    //! The K-function
    void kfn(TOutput (&res)[3], TScale& ieps, TMass const& xpi, TMass const& xm, TMass const& xmp) const;

    //! Solution of the quadratic equation
    void solveabc(TMass const& a, TMass const&b, TMass const& c, TOutput (&z)[2]) const;

    //! Solution of the quadratic equation passing the discriminant.
    void solveabcd(TOutput const& a, TOutput const&b, TOutput const& c, TOutput const& d, TOutput (&z)[2]) const;
    void solveabcd(TOutput const& a, TOutput const&b, TOutput const& c, TOutput (&z)[2]) const;

    //! Ratio function.
    void ratgam(TOutput &ratp, TOutput &ratm, TScale &ieps, TMass const& p3sq, TMass const& m3sq, TMass const& m4sq) const;

    //! Ratio function.
    void ratreal(TMass const& si, TMass const& ta, TMass &rat, TScale &ieps) const;

  private:
    //! zero cutoff value
    TScale _qlonshellcutoff;

    //! ddilog coefficients
    vector<TScale> _C;

    //! even-n Bernoulli numbers.
    vector<TScale> _B;

  protected:
    // Pi constants
    TScale _pi;
    TScale _pi2;
    TScale _pio3;
    TScale _pio6;
    TScale _pi2o3;
    TScale _pi2o6;
    TScale _pi2o12;

    // Real constants
    TScale _zero;
    TScale _half;
    TScale _one;
    TScale _two;
    TScale _three;
    TScale _four;
    TScale _five;
    TScale _six;
    TScale _ten;

    // Real Epsilons
    TScale _eps;
    TScale _eps4;
    TScale _eps7;
    TScale _eps10;
    TScale _eps14;
    TScale _eps15;
    TScale _xloss;
    TScale _neglig;
    TScale _reps;

    // Pi complex contants
    TOutput _2ipi;
    TOutput _ipio2;
    TOutput _ipi;

    // Complex constants
    TOutput _czero;
    TOutput _chalf;
    TOutput _cone;
    TOutput _ctwo;
    TOutput _cthree;
    TOutput _cfour;

    // Complex epsilons
    TOutput _ieps;
    TOutput _ieps2;
    TOutput _ieps50;
  };

  /*!
   * \brief The Splash class.
   *
   * Splash allocates a singleton which prints to screen
   * the QCDLoop welcome message when it is constructed.
   */
  class Splash
  {
  public:
    //! Prints splash screen
    static void Show() { getInstance(); }

  private:
    Splash(); //!< The Constructor

    //! Allocates the Splash singleton
    static Splash* getInstance()
    {
      if (!_instance)
        _instance = new Splash();
      return _instance;
    }
    static Splash* _instance;
  };
}
