//
// QCDLoop 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch
//          Keith Ellis: keith.ellis@durham.ac.uk
//          Giulia Zanderighi: giulia.zanderighi@cern.ch

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

namespace ql
{
  /*!
   * \brief The Timer class.
   *
   * Computes the calculation time.
   */
  class Timer {

   public:
    //! Starts the timer.
    void start(){ gettimeofday(&startTime, NULL); }

    //! Stops the timer.
    double stop()
    {
      timeval endTime;
      long seconds, useconds;
      double duration;
      gettimeofday(&endTime, NULL);
      seconds  = endTime.tv_sec  - startTime.tv_sec;
      useconds = endTime.tv_usec - startTime.tv_usec;
      duration = seconds + useconds/1E6f;
      return duration;
    }

    /*!
     * \brief Prints enlapsed time
     * \param duration input time from Timer::stop()
     */
    static void printTime(double const& duration)
    {
      printf("Elapsed Time: %5.6f seconds\n", duration);
    }

  private:
   timeval startTime;
  };

}
