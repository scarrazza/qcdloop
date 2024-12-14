//
// QCDLoop 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch
//          Keith Ellis: keith.ellis@durham.ac.uk
//          Giulia Zanderighi: giulia.zanderighi@cern.ch

#include "qcdloop/types.h"

namespace std
{
  ostream& operator<<(std::ostream& out, ql::qdouble f)
  {
     char buf[200];
     std::ostringstream format;
#ifdef __x86_64__
     format << "%." << (std::min)(190L, out.precision()) << "Qe";
     quadmath_snprintf(buf, 200, format.str().c_str(), f);
#endif
#ifdef __aarch64__
     format << "%." << (std::min)(190L, out.precision()) << "Le";
     snprintf(buf, 200, format.str().c_str(), f);
#endif
     out << buf;
     return out;
   }

  ostream& operator<<(std::ostream& out, ql::qcomplex f)
  {
     out << "(" << crealq(f) << "," << cimagq(f) << ")";
     return out;
  }
 
  ostream& operator<<(std::ostream& os, ql::Code code) 
  {      
    return os << "\033[" << static_cast<int>(code) << "m";
  }
}
