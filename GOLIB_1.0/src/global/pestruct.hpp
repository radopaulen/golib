// Copyright (C) 2012, 2014 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__PESTRUCT_HPP
#define MC__PESTRUCT_HPP

#include "odestruct.hpp"

namespace mc
{

//! @brief C++ base class for dynamic optimization problem formulation
class PESTRUCT: public ODESTRUCT
{
public:
  //! @brief Infinity
  static double INF;
  //! @brief Enumeration type for constraints
  enum t_CTR{
    EQ=0,	//!< Equality constraint
    LE,		//!< Inequality constraint
    GE		//!< Inequality constraint
  };

protected:
  //! @brief Number of output variables
  unsigned int _ny;
  //! @brief Number of constraints
  unsigned int _nc;
  
public:
  //! @brief Default constructor
  PESTRUCT
    ( const unsigned int np, const unsigned int nx, const unsigned int ny,
      const unsigned int nc=0, const unsigned int ni=0 )
    : ODESTRUCT(np,nx,ni), _ny(ny), _nc(nc)
    {}

  //! @brief Const member function retreiving the number of output variables
  unsigned int ny() const
    { return _ny; }
  //! @brief Const member function returning number of constraints
  unsigned int nc() const
    { return _nc; }

private:
  //! @brief Private methods to block default compiler methods
  PESTRUCT();
};

double PESTRUCT::INF = 2e19;

} // end namespace mc

#endif
