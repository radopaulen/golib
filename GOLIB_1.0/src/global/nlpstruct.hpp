// Copyright (C) 2012, 2014 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__NLPSTRUCT_HPP
#define MC__NLPSTRUCT_HPP

namespace mc
{

//! @brief C++ base class for nonlinear optimization problem formulation
class NLPSTRUCT
{
public:
  //! @brief Infinity
  static double INF;
  //! @brief Enumeration type for objective function
  enum t_OBJ{
    MIN=0,	//!< Minimization
    MAX		//!< Maximization
  };
  //! @brief Enumeration type for constraints
  enum t_CTR{
    EQ=0,	//!< Equality constraint
    LE,		//!< Inequality constraint
    GE		//!< Inequality constraint
  };

protected:
  //! @brief Number of decision variables
  unsigned int _np;
  //! @brief Number of constraints
  unsigned int _nc;
  
public:
  //! @brief Default constructor
  NLPSTRUCT
    ( const unsigned int np, const unsigned int nc )
    : _np(np), _nc(nc)
    {}

  //! @brief Const member function returning number of decision variables
  unsigned int np() const
    { return _np; }
  //! @brief Const member function returning number of constraints
  unsigned int nc() const
    { return _nc; }

private:
  //! @brief Private methods to block default compiler methods
  NLPSTRUCT();
};

double NLPSTRUCT::INF = 2e19;

} // end namespace mc

#endif
