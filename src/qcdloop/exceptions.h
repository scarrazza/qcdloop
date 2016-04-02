//
// QCDLoop 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch
//          Keith Ellis: keith.ellis@durham.ac.uk
//          Giulia Zanderighi: giulia.zanderighi@cern.ch

#pragma once

#include <exception>
#include <stdexcept>

namespace ql
{
  /// Error to be thrown when out of the valid range.
  class LogicException: public std::logic_error
  {
  public:
    LogicException(const std::string& tag, const std::string& what) : std::logic_error(tag + ": " + what) {}
  };

  /// Error to be thrown when out of the valid range.
  class RangeError: public LogicException
  {
  public:
    RangeError(const std::string& tag, const std::string& what) : LogicException(tag,what) {}
  };

  /// Error to be thrown when out of the valid range.
  class LengthError: public LogicException
  {
  public:
    LengthError(const std::string& tag, const std::string& what): LogicException(tag,what) {}
  };

  /// Error to be thrown when out of the valid range.
  class ConvergenceError: public LogicException
  {
  public:
    ConvergenceError(const std::string& tag, const std::string& what) : LogicException(tag,what) {}
  };

}
