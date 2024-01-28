// Copyright (C) 2012, 2013 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__BASE_GSL_HPP
#define MC__BASE_GSL_HPP

#include <iostream>
#include <vector>
#include <sys/time.h>

#include "gsl/gsl_errno.h"
#include "gsl/gsl_odeiv2.h"

namespace mc
{
//! @brief C++ base class computing solutions of parametric ODEs
////////////////////////////////////////////////////////////////////////
//! mc::BASE_ODE is a C++ base class for computing solutions of
//! parametric ordinary differential equations (ODEs).
////////////////////////////////////////////////////////////////////////
class BASE_ODE
{
protected:
  //! @brief current time
  double _t;

  //! @brief current stage
  unsigned int _istg;

public:
  /** @defgroup ODEBND Enclosing the Reachable Set of Parametric Nonlinear ODEs
   *  @{
   */
  //! @brief Integrator status
  enum STATUS{
     NORMAL=0,	//!< Normal execution
     FAILURE,	//!< Integration breakdown (bounds explosion)
     FATAL	//!< Interruption due to errors in third-party libraries
  };

  //! @brief Default constructor
  BASE_ODE()
    {}

  //! @brief Default destructor
  virtual ~BASE_ODE()
    {}

  //! @brief Structure storing integration statistics
  struct Stats
  {
    //! @brief Constructor
    Stats():
      cputime(0.), numRHS(0), numJAC(0)
      {}
    //! @brief Constructor
    void reset()
      { cputime = 0.; numRHS = numJAC = 0; }

    //! @brief CPU time
    double cputime;
    //! @brief Number of right-hand side (RHS) evaluations
    unsigned int numRHS;
    //! @brief Number of Jacobian evaluations
    unsigned int numJAC;
  } stats;

  //! @brief Return last successful integration time
  double time() const
    { return _t; }
  /** @} */

protected:

  //! @brief Function to initialize GSL driver
  void _init_stats();

  //! @brief Function to finalize statistics for GSL
  void _final_stats();

  //! @brief Function to display intermediate results
  void _print_stats
    ( const Stats&stats, std::ostream&os=std::cout ) const;

  //! @brief Private methods to block default compiler methods
  BASE_ODE(const BASE_ODE&);
  BASE_ODE& operator=(const BASE_ODE&);
};

inline void
BASE_ODE::_init_stats()
{
  // Initialize statistics
  stats.reset();
  timeval time;
  gettimeofday(&time, 0) ;
  stats.cputime = - time.tv_sec - time.tv_usec*1e-6;
}

inline void
BASE_ODE::_final_stats()
{
  // Get final CPU time
  timeval time;
  gettimeofday(&time, 0);
  stats.cputime += time.tv_sec + time.tv_usec*1e-6;
}

inline void
BASE_ODE::_print_stats
( const Stats&stats, std::ostream&os ) const
{
  // Statistics
  os << " No EVALATIONS" << "   RHS: " << stats.numRHS
                         << "   JAC: " << stats.numJAC
     << std::endl
     << " CPU TIME (SEC)     " << std::fixed << std::left
                               << std::setprecision(5) << stats.cputime
     << std::endl;
  return;
}

//! @brief C++ base class computing solutions of parametric ODEs using GSL
////////////////////////////////////////////////////////////////////////
//! mc::BASE_GSL is a C++ base class for computing solutions of
//! parametric ordinary differential equations (ODEs) using GSL.
////////////////////////////////////////////////////////////////////////
class BASE_GSL: public BASE_ODE
{
protected:
  //! @brief full GSL state
  double *_vec_state;

public:
  /** @defgroup ODEBND Enclosing the Reachable Set of Parametric Nonlinear ODEs
   *  @{
   */
  //! @brief Integrator status
  enum STATUS{
     NORMAL=0,	//!< Normal execution
     FAILURE,	//!< Integration breakdown (bounds explosion)
     FATAL	//!< Interruption due to errors in third-party libraries
  };

  //! @brief Default constructor
  BASE_GSL
    (): BASE_ODE()
    {
      // Initalize GSL arrays
      _vec_state = 0;
    }

  //! @brief Default destructor
  virtual ~BASE_GSL()
    {     
      // Free GSL arrays
      delete[] _vec_state;
    }

  //! @brief GSL options
  struct Options
  {
    //! @brief Constructor
    Options():
      INTMETH(RKF45), H0(1e-2), HMIN(0e0), NMAX(0), RTOL(1e-6), ATOL(1e-6)
      {}
    //! @brief Assignment operator
    template <typename U> Options& operator=
      ( U&options ){
        INTMETH   = options.INTMETH;
        H0        = options.H0;
        HMIN      = options.HMIN;
        NMAX      = options.NMAX;
        RTOL      = options.RTOL;
        ATOL      = options.ATOL;
        return *this;
      }
    //! @brief Enumeration of numerical integration algorithms
    enum INTEGRATION_METHOD{
      RKF45=0,		//!< Explicit embedded Runge-Kutta-Fehlberg (4,5) method (non-stiff systems) [Default]
      RK8PD,		//!< Explicit embedded Runge-Kutta Prince-Dormand (8,9) method (non-stiff systems)
      MSADAMS,		//!< Variable-coefficient linear multistep Adams method in Nordsieck form (non-stiff systems)
      MSBDF		//!< Variable-coefficient linear multistep backward differentiation formula (BDF) method in Nordsieck form (stiff systems)
    } options;
    //! @brief Numerical integration method
    INTEGRATION_METHOD INTMETH;
    //! @brief Initial step-size (Default: 1e-2)
    double H0;
    //! @brief Minimum step-size (Default: 0e0)
    double HMIN;
    //! @brief Maximum number of steps in a time stage (Default: 0)
    unsigned int NMAX;
    //! @brief Relative integration tolerance (Default: 1e-6)
    double RTOL;
    //! @brief Absolute integration tolerance (Default: 1e-6)
    double ATOL;
  };
  /** @} */

protected:

  //! @brief Private methods to block default compiler methods
  BASE_GSL(const BASE_GSL&);
  BASE_GSL& operator=(const BASE_GSL&);
};

} // end namescape mc

#endif

