//
// QCDLoop 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch
//          Keith Ellis: keith.ellis@durham.ac.uk
//          Giulia Zanderighi: giulia.zanderighi@cern.ch

#include "qcdloop/types.h"

namespace std
{
#if defined(__x86_64__) || defined(__i386__) || defined(__PPC64__) || defined(__ppc64__) || defined(_ARCH_PPC64)
  ostream& operator<<(std::ostream& out, ql::qdouble f)
  {
     char buf[200];
     std::ostringstream format;
     format << "%." << (std::min)(190L, out.precision()) << "Qe";
     quadmath_snprintf(buf, 200, format.str().c_str(), f);
     out << buf;
     return out;
   }

  ostream& operator<<(std::ostream& out, ql::qcomplex f)
  {
     out << "(" << crealq(f) << "," << cimagq(f) << ")";
     return out;
  }
#endif
#if defined(__aarch64__)
  ostream& operator<<(std::ostream& out, ql::qcomplex f)
  {
     out << "(" << creall(f) << "," << cimagl(f) << ")";
     return out;
  }
#endif
 
  ostream& operator<<(std::ostream& os, ql::Code code) 
  {      
    return os << "\033[" << static_cast<int>(code) << "m";
  }
}
