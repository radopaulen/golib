// Copyright (C) 2009-2013 Benoit Chachuat, Imperial College London.
// All Rights Reserved.

/*!
\page page_FFUNC Tree Decomposition of Factorable Functions
\author Mario E. Villanueva & Benoit Chachuat
\version 0.0
\date 2013
\bug No known bugs.

TBC.

*/

#ifndef MC__FFUNC_HPP
#define MC__FFUNC_HPP

#include <iostream>
#include <string>
#include <sstream>
#include <stdarg.h>
#include <set>
#include <list>
#include <utility>
#include <algorithm>

#include "mcfunc.hpp"
#include "ffdep.hpp"

#define  MC__FFUNC_DEBUG

namespace mc
{
class FFOp;
class FFTree;

//! @brief Structure defining the numeric field of a factorable program variable
////////////////////////////////////////////////////////////////////////
//! mc::FFNum is a C++ structure defining the numeric field of a
//! variable in a factorable function, which can either a real scalar
//! (double), or an integer scalar (int)
////////////////////////////////////////////////////////////////////////
struct FFNum
{
  //! @brief Enumeration type for numeric variables in factorable function
  enum TYPE{
    INT=0,	//!< Integer value
    REAL	//!< Real value
  };
  //! @brief Variable type
  TYPE t;
  //! @brief Integer/real variable value
  union{
    int n;
    double x;
  }; 

  //! @brief Constructor for an integer variable
  FFNum( const int i=0 ):
    t(INT), n(i)
    {}
  //! @brief Constructor for a real variable
  FFNum( const double d ):
    t(REAL), x(d)
    {}

  //! @brief Constructor for an integer scalar
  FFNum& operator=
    ( const int i )
    { t = INT; n = i; return *this; }
  //! @brief Constructor for a real scalar
  FFNum& operator=
    ( const double d )
    { t = REAL; x = d; return *this; }
  //! @brief Copy constructor
  FFNum& operator=
    ( const FFNum&num )
    { t = num.t; t==REAL? x=num.x: n; return *this; }
};

//! @brief Structure comparing values of scalars in factorable functions for equality
////////////////////////////////////////////////////////////////////////
//! mc::eq_FFNum is a C++ structure comparing the numeric field of a
//! FFVar object in a factorable function for equality.
////////////////////////////////////////////////////////////////////////
struct eq_FFNum
////////////////////////////////////////////////////////////////////////
{
  bool operator()
    ( const FFNum*Num1, const FFNum*Num2 ) const
    {
      if( Num1->t != Num2->t ) return false;
      switch( Num1->t ){
        case FFNum::INT:  return Num1->n==Num2->n? true: false;
        case FFNum::REAL: return isequal( Num1->x, Num2->x );
      }
    }
};

//! @brief Structure comparing values of scalars in factorable functions for strict inequality
////////////////////////////////////////////////////////////////////////
//! mc::lt_FFNum is a C++ structure comparing the numeric field of a
//! FFVar object in a factorable function for strict inequality.
////////////////////////////////////////////////////////////////////////
struct lt_FFNum
////////////////////////////////////////////////////////////////////////
{
  bool operator()
    ( const FFNum*Num1, const FFNum*Num2 ) const
    {
      if( Num1->t < Num2->t ) return true;
      if( Num1->t > Num2->t ) return false;
      switch( Num1->t ){
        case FFNum::INT:  return Num1->n<Num2->n? true: false;
        case FFNum::REAL: return !isequal( Num1->x, Num2->x ) && Num1->x<Num2->x? true: false;
      }
      return false;
    }
};

//! @brief Class defining variables in a factorable function
////////////////////////////////////////////////////////////////////////
//! mc::FFVar is a C++ class defining variables in the factored form of
//! a factorable function.
////////////////////////////////////////////////////////////////////////
class FFVar
////////////////////////////////////////////////////////////////////////
{
  // friends of this class with other classes/structures/operators
  friend class FFTree;
  friend struct lt_FFVar;
  friend std::ostream& operator<< ( std::ostream&, const FFTree& );
  friend std::ostream& operator<< ( std::ostream&, const FFOp& );

  // friends of this class for operator and function overloading
  friend std::ostream& operator<< ( std::ostream&, const FFVar& );
  friend FFVar operator+ ( const FFVar& );
  friend FFVar operator+ ( const FFVar&, const FFVar& );
  template <typename V> friend FFVar operator+ ( const V&, const FFVar& );
  template <typename V> friend FFVar operator+ ( const FFVar&, const V& );
  friend FFVar operator- ( const FFVar& );
  friend FFVar operator- ( const FFVar&, const FFVar& );
  template <typename V> friend FFVar operator- ( const V&, const FFVar& );
  template <typename V> friend FFVar operator- ( const FFVar&, const V& );
  friend FFVar operator* ( const FFVar&, const FFVar& );
  template <typename V> friend FFVar operator* ( const V&, const FFVar& );
  template <typename V> friend FFVar operator* ( const FFVar&, const V& );
  friend FFVar operator/ ( const FFVar&, const FFVar& );
  template <typename V> friend FFVar operator/ ( const V&, const FFVar& );
  template <typename V> friend FFVar operator/ ( const FFVar&, const V& );
  friend FFVar max ( const FFVar&, const FFVar& );
  friend FFVar max ( const unsigned int, const FFVar* );
  template <typename V> friend FFVar max ( const V&, const FFVar& );
  template <typename V> friend FFVar max ( const FFVar&, const V& );
  friend FFVar min ( const FFVar&, const FFVar& );
  friend FFVar min ( const unsigned int, const FFVar* );
  template <typename V> friend FFVar min ( const V&, const FFVar& );
  template <typename V> friend FFVar min ( const FFVar&, const V& );
  friend FFVar inv   ( const FFVar& );
  friend FFVar sqr   ( const FFVar& );
  friend FFVar exp   ( const FFVar& );
  friend FFVar log   ( const FFVar& );
  friend FFVar sqrt  ( const FFVar& );
  friend FFVar fabs  ( const FFVar& );
  friend FFVar cos   ( const FFVar& );
  friend FFVar sin   ( const FFVar& );
  friend FFVar tan   ( const FFVar& );
  friend FFVar acos  ( const FFVar& );
  friend FFVar asin  ( const FFVar& );
  friend FFVar atan  ( const FFVar& );
  friend FFVar erf   ( const FFVar& );
  friend FFVar erfc  ( const FFVar& );
  friend FFVar fstep ( const FFVar& );
  friend FFVar bstep ( const FFVar& );
  friend FFVar pow   ( const FFVar&, const int );
  template <typename V> friend FFVar pow ( const FFVar&, const V& );
  friend FFVar pow   ( const double, const FFVar& );

public:

  // other operator overloadings
  FFVar& operator= ( const FFVar& );
  FFVar& operator= ( const int );
  FFVar& operator= ( const double );
  template <typename V> FFVar& operator+= ( const V& );
  template <typename V> FFVar& operator-= ( const V& );
  template <typename V> FFVar& operator*= ( const V& );
  template <typename V> FFVar& operator/= ( const V& );

  typedef std::list< FFOp* > t_Ops;
  //typedef typename t_Ops::iterator it_Ops;
  //typedef typename t_Ops::const_iterator cit_Ops;
  typedef typename std::pair< FFOp*, t_Ops > pt_Ops;

  /** @defgroup FFunc Tree Decomposition of Factorable Functions
   *  @{
   */
  //! @brief Index for 'free' variables in factorable function
  static const int NOREF = -33;
  //! @brief Enumeration type for variables in factorable function
  enum TYPE{
    VAR=0,	//!< Original variable
    AUX,	//!< Auxiliary variable
    CINT,	//!< Integer constant
    CREAL	//!< Real constant
  };
  //! @brief Typedef for variable identifier in factorable function
  typedef std::pair< TYPE, long > pt_idVar;
  /** @} */

private:
  //! @brief Identifier (type and index)
  pt_idVar _id;
  //! @brief Numeric field (integer or real)
  FFNum _num;
  //! @brief Dependence and linearity
  FFDep _dep;
  //! @brief Pointer to underlying factorable function - NULL for variable identifier NOREF
  FFTree *_tree;
  //! @brief Pointer to parent (_ops.first) and children (_ops.second) operations - _ops.first=NULL for unreferenced constants
  pt_Ops _ops;

public:

  /** @ingroup FFunc
   *  @{
   */
  //! @brief Constructor for variable <a>ix</a> in factorable function <a>*tree</a>
  FFVar
    ( FFTree*tree, const unsigned long ix );

  //! @brief Set the variable <a>ix</a> in the factorable function <a>*tree</a>.
  FFVar& set
    ( FFTree*tree, const unsigned long ix )
    { *this = FFVar( tree, ix ); return *this; }

  //! @brief Constructor for integer constant
  FFVar
    ( const int i=0 )
    : _id( CINT, NOREF ), _num( i ), _dep( i ), _tree( 0 )
    { _ops.first = 0; }

  //! @brief Constructor for real parameter
  FFVar
    ( const double d )
    : _id( CREAL, NOREF ), _num( d ), _dep( d ), _tree( 0 )
    { _ops.first = 0; }

  //! @brief Copy constructor
  FFVar
    ( const FFVar&Var )
    : _id( Var._id ), _num( Var._num ), _dep( Var._dep ), _tree( Var._tree ), _ops( Var._ops )
    {}
  /** @} */

private:

  //! @brief Constructor for auxiliary variable in factorable function <a>*tree</a> defined from operation <a>*Op</a>
  FFVar
    ( FFTree*tree, const FFDep&dep, FFOp*op=0 );

  //! @brief Constructor for a variable with identifier <a>id</a> in factorable function <a>*tree</a>
  FFVar
    ( FFTree*tree, const pt_idVar&id );
    
public:

  /** @ingroup FFunc
   *  @{
   */
  //! @brief Get variable identifier
  const std::pair<TYPE,long> id() const
    { return _id; }

  //! @brief Get reference to variable identifier
  std::pair<TYPE,long>& id()
    { return _id; }

  //! @brief Get const reference to variable numeric field
  const FFNum& num() const
    { return _num; }

  //! @brief Get const reference to variable dependencies
  const FFDep& dep() const
    { return _dep; }

  //! @brief Get reference to variable dependencies
  FFDep& dep()
    { return _dep; }

  //! @brief Get const pointer to defining operation
  const pt_Ops ops() const
    { return _ops; }

  //! @brief Get pointer to defining operation
  pt_Ops& ops()
    { return _ops; }

  //! @brief Get const pointer to factorable function
  const FFTree* tree() const
    { return _tree; }

  //! @brief Get pointer to factorable function
  FFTree*& tree()
    { return _tree; }

  //! @brief Get variable name
  std::string name() const
    { return _name(_id); }
  /** @} */

private:

  //! @brief Return string with variable name for identifier <a>id</a>
  static std::string _name
    ( const std::pair<TYPE,long> id )
    {
      std::ostringstream ovar;
      id.first==VAR? ovar<<"X": ovar<<"Z";
      ovar<<id.second;
      return ovar.str();
    }
};

//! @brief Structure comparing variable identifiers in a factorable function for ordering in set FFTree::_Vars
////////////////////////////////////////////////////////////////////////
//! mc::lt_FFVar is a C++ structure comparing variable identifiers in a
//! factorable function for ordering in set FFTree::_Vars.
////////////////////////////////////////////////////////////////////////
struct lt_FFVar
////////////////////////////////////////////////////////////////////////
{
  bool operator()
    ( const FFVar*Var1, const FFVar*Var2 ) const
    {
      // Order variables/constants w.r.t. their types first
      if( Var1->_id.first < Var2->_id.first ) return true;
      if( Var1->_id.first > Var2->_id.first ) return false;
      // If variables, order w.r.t. their index next
      switch( Var1->_id.first ){
        case FFVar::VAR:
        case FFVar::AUX:
          if( Var1->_id.second < Var2->_id.second ) return true;
          if( Var1->_id.second > Var2->_id.second ) return false;
          break;
        case FFVar::CINT: case FFVar::CREAL:
          lt_FFNum ltNum;
          return ltNum( &Var1->_num, &Var2->_num );
          break;
      }
      return false;
    }
};

//! @brief Class defining operations in a factorable function
////////////////////////////////////////////////////////////////////////
//! mc::FFOp is a C++ class defining operations in the factored form of
//! a factorable function.
////////////////////////////////////////////////////////////////////////
class FFOp
////////////////////////////////////////////////////////////////////////
{
public:

  /** @ingroup FP
   *  @{
   */
  //! @brief Enumeration type for unary and binary operations
  enum TYPE{
    CNST=0, VAR,
    PLUS, NEG, MINUS, TIMES, SCALE, DIV,
    EXP, LOG, SQRT, SQR, IPOW, POW, SIN, COS, TAN, ASIN, ACOS, ATAN,
    FABS, ERF, FSTEP, MINF, MAXF
  };

  //! @brief Constructor
  FFOp( TYPE top, FFVar*lop=0, FFVar*rop=0, FFVar*res=0 );

  //! @brief Destructor
  ~FFOp()
    {}

  //! @brief Type of operation
  TYPE type;
  //! @brief Pointer to operation result
  FFVar* pres;
  //! @brief Pointer to left operand
  FFVar* plop;
  //! @brief Pointer to right operand
  FFVar* prop;

  //! @brief Propagate subset of operations participating in subtree
  void propagate_subtree
    ( std::list<FFOp*>&Ops );
  //! @brief Propagate script for tree decomposition using DOT and display to <a>os</a>
  void generate_dot_script
    ( std::ostream&os ) const;
  //! @brief Append script for current operation using DOT to <a>os</a>
  void append_dot_script
    ( std::ostream&os ) const;
  //! @brief Append script for factor <a>fname</a> using DOT to <a>os</a>
  void append_dot_script_factor
    ( std::ostream&os, const std::string&fname, const bool unary,
      const unsigned int fontsize, const bool dotted=false ) const;
  //! @brief Append script for current variable using DOT to <a>os</a>
  void append_dot_script_variable
    ( std::ostream&os, const bool constant, const unsigned int fontsize ) const;

  //! @brief Flag operation as visited or not
  void flag
    ( const bool visited=true ) const
    { _visited = visited; }
  //! @brief Retreive operation status (visited or not)
  const bool stat() const
    { return _visited; }
  //! @brief Retreive/set operation status (visited or not)
  const bool& stat()
    { return _visited; }
  //! @brief Return whether or not operation is univariate
  bool is_univariate();
  
  /** @} */

private:
  //! @brief Whether a cut has been visited (e.g. in a tree navigation)
  mutable bool _visited;
};

//! @brief C++ structure for comparing operations in a factorable program for ordering in set FFTree::_Ops
////////////////////////////////////////////////////////////////////////
//! mc::lt_FFOp is a C++ structure for comparing operations in a
//! factorable program based on their types and operands for ordering
//! in set FFTree::_Ops.
////////////////////////////////////////////////////////////////////////
struct lt_FFOp
////////////////////////////////////////////////////////////////////////
{
  bool operator()
    ( const FFOp*Op1, const FFOp*Op2 ) const
    {
      // Sort by type of operation first
      if( Op1->type < Op2->type ) return true;
      if( Op1->type > Op2->type ) return false;

      // Sort by variable type next
      lt_FFVar ltVar;
      if( !Op1->plop ) return ltVar( Op1->pres, Op2->pres );
      if( ltVar( Op1->plop, Op2->plop ) ) return true;
      if( ltVar( Op2->plop, Op1->plop ) ) return false;
      if( Op1->prop ) return ltVar( Op1->prop, Op2->prop );
      return false;
    }
};

//! @brief C++ structure for comparing operations in a factorable program
////////////////////////////////////////////////////////////////////////
//! mc::lt_FFOp is a C++ structure for comparing operations in a
//! factorable program based on their types only.
////////////////////////////////////////////////////////////////////////
struct range_FFOp
////////////////////////////////////////////////////////////////////////
{
  bool operator()
    ( const FFOp*Op1, const FFOp*Op2 ) const
    { return ( Op1->type < Op2->type ); }
};

//! @brief C++ class representing the tree decomposition of factorable functions
////////////////////////////////////////////////////////////////////////
//! mc::FFTree is a C++ class representing the tree decomposition of
//! factorable functions and enabling basic manipulations on that tree.
////////////////////////////////////////////////////////////////////////
class FFTree
////////////////////////////////////////////////////////////////////////
{
  // friends of this class with other classes/structures/operators
  friend class FFVar;
  friend class FFOp;
  friend FFVar operator+ ( const FFVar& );
  friend FFVar operator+ ( const FFVar&, const FFVar& );
  template <typename V> friend FFVar operator+ ( const V&, const FFVar& );
  template <typename V> friend FFVar operator+ ( const FFVar&, const V& );
  friend FFVar operator- ( const FFVar& );
  friend FFVar operator- ( const FFVar&, const FFVar& );
  template <typename V> friend FFVar operator- ( const V&, const FFVar& );
  template <typename V> friend FFVar operator- ( const FFVar&, const V& );
  friend FFVar operator* ( const FFVar&, const FFVar& );
  template <typename V> friend FFVar operator* ( const V&, const FFVar& );
  template <typename V> friend FFVar operator* ( const FFVar&, const V& );
  friend FFVar operator/ ( const FFVar&, const FFVar& );
  template <typename V> friend FFVar operator/ ( const V&, const FFVar& );
  template <typename V> friend FFVar operator/ ( const FFVar&, const V& );
  friend FFVar max ( const FFVar&, const FFVar& );
  friend FFVar max ( const unsigned int, const FFVar* );
  template <typename V> friend FFVar max ( const V&, const FFVar& );
  template <typename V> friend FFVar max ( const FFVar&, const V& );
  friend FFVar min ( const FFVar&, const FFVar& );
  friend FFVar min ( const unsigned int, const FFVar* );
  template <typename V> friend FFVar min ( const V&, const FFVar& );
  template <typename V> friend FFVar min ( const FFVar&, const V& );
  friend FFVar inv   ( const FFVar& );
  friend FFVar sqr   ( const FFVar& );
  friend FFVar exp   ( const FFVar& );
  friend FFVar log   ( const FFVar& );
  friend FFVar sqrt  ( const FFVar& );
  friend FFVar fabs  ( const FFVar& );
  friend FFVar cos   ( const FFVar& );
  friend FFVar sin   ( const FFVar& );
  friend FFVar tan   ( const FFVar& );
  friend FFVar acos  ( const FFVar& );
  friend FFVar asin  ( const FFVar& );
  friend FFVar atan  ( const FFVar& );
  friend FFVar erf   ( const FFVar& );
  friend FFVar erfc  ( const FFVar& );
  friend FFVar fstep ( const FFVar& );
  friend FFVar bstep ( const FFVar& );
  friend FFVar pow   ( const FFVar&, const int );
  template <typename V> friend FFVar pow ( const FFVar&, const V& );
  friend FFVar pow   ( const double, const FFVar& );

  // friends of this class for operator and function overloading
  friend std::ostream& operator<< ( std::ostream&, const FFTree& );

public:
  /** @ingroup FFunc
   *  @{
   */
  typedef typename FFVar::pt_idVar pt_idVar;
  typedef std::set< FFVar*, lt_FFVar > t_Vars;
  typedef std::set< FFOp*,  lt_FFOp >  t_Ops;
  typedef typename t_Vars::iterator it_Vars;
  typedef typename t_Vars::const_iterator cit_Vars;
  typedef typename t_Ops::iterator  it_Ops;
  typedef typename t_Ops::const_iterator  cit_Ops;
  /** @} */

protected:
  //! @brief Number of original variables in factorable function
  unsigned long _nvar;
  //! @brief Number of auxiliary variables in factorable function
  unsigned long _naux;

  //! @brief Set of variables in factorable function
  t_Vars _Vars;
  //! @brief Set of operations in factorable function
  t_Ops  _Ops;

public:
  /** @ingroup FFunc
   *  @{
   */
  //! @brief Default Constructor
  FFTree():
    _nvar( 0 ), _naux( 0 ) 
    {}

  //! @brief Destructor
  virtual ~FFTree()
    { clear(); }

  //! @brief Exceptions of mc::FFTree
  class Exceptions
  {
  public:
    //! @brief Enumeration type for exception handling
    enum TYPE{
      INIT = 1,		//!< Error due to a invalid FFTree pointer in initialization of FFVar
      TREE,		//!< Error due to an operation between variables linked to different factorable functions
      INTERN = -1, 	//!< Internal error
      UNDEF = -33, 	//!< Error due to calling a function/feature not yet implemented in MC++
    };
    //! @brief Constructor for error <a>ierr</a>
    Exceptions( TYPE ierr ) : _ierr( ierr ){}
    //! @brief Inline function returning the error flag
    int ierr(){ return _ierr; }
    //! @brief Error description
    std::string what(){
      switch( _ierr ){
      case INIT:
        return "mc::FFTree\t Inavlid mc::FFTree passed in initilization of mc:FFVar object";
      case TREE:
        return "mc::FFTree\t Operation between mc::FFVar objects linked to different mc::FFTree objects";
      case UNDEF:
        return "mc::FFTree\t Feature not yet implemented in MC++";
      default:
        return "mc::FFDep\t Undocumented error";
      }
    }
  private:
    TYPE _ierr;
  };

/*
  //! @brief Options of mc::FFTree
  struct Options
  {
    //! @brief Constructor
    Options():
      SOLVER_OUTPUT_FILE("")
      {}
    //! @brief Name of output file for optimization model
    std::string SOLVER_OUTPUT_FILE;
  } options;
*/

  //! @brief Number of original variables in factorable function
  unsigned long nvar() const
    { return _nvar; }
  
  //! @brief Number of auxiliary variables in factorable function
  unsigned long naux() const
    { return _naux; }
  
  //! @brief Reference to set of (all) variables in factorable function
  const t_Vars& Vars() const
    { return _Vars; }

  //! @brief Clear factorable function (all variables and operations)
  void clear()
    { _clear_variables(); _clear_operations(); _naux=_nvar=0; }

  //! @brief Extract list of nodes corresponding to the factors <a>*F</a>
  std::list<FFOp*> subtree
    ( const unsigned int nVar, const FFVar*pVar ) const;

  //! @brief Output list of nodes in <a>Ops</a> to <a>os</a>
  static void output
    ( std::list<FFOp*>&Ops, std::ostream&os=std::cout );

  //! @brief Generate script for tree visualization of factors <a>*F</a> using DOT
  void dot_script
    ( const unsigned int nVar, const FFVar*pVar, std::ostream&os=std::cout ) const;

  //! @brief Propagate constraints in the factorable program to reduce variable range, for the constraint <a>constr</a> or all constraints in case <a>constr=0</a>; on return, <a>false</a> is an indication of an infeasible program
  //bool propagate_constraints
  //  ( FFOp*constr=0 );
  /** @} */
   
protected:
  //! @brief Erase operation <a>op</a> in set <a>_Ops</a>
  bool _remove_operation
    ( FFOp*op );
  //! @brief Erase all operations in set <a>_Ops</a>
  void _clear_operations()
    { it_Ops ito = _Ops.begin();
      for( ; ito != _Ops.end(); ++ito ) delete *ito;
      _Ops.clear(); }
  //! @brief Reset all operations in set <a>_Ops</a>
  void _reset_operations() const
    { it_Ops ito = _Ops.begin();
      for( ; ito != _Ops.end(); ++ito ) (*ito)->flag( false ); }
  //! @brief Looks for the operation of type <a>top</a> with left and right operands <a>lop</a>, <a>rop</a> in set <a>_Ops</a> and adds it if not found
  FFOp* _insert_operation
    ( const typename FFOp::TYPE top, FFVar*lop, FFVar*rop=0 );
  //! @brief Looks for the binary operation of type <a>top</a> with left and right operands <a>Var1</a>, <a>Var2</a> in set <a>_Ops</a> and adds it if not found; also adds new auxiliary variable in set <a>_Vars</a> and update list of dependencies in both operands in <a>_Vars</a>
  static FFVar& _insert_binary_operation
    ( const typename FFOp::TYPE top, const FFDep&dep, const FFVar&Var1, const FFVar&Var2 );
  //! @brief Looks for the binary operation of type <a>top</a> with left and right operands <a>Cst1</a>, <a>Var2</a> in set <a>_Ops</a> and adds it if not found; also adds new auxiliary variable as well as constant <a>Cst1</a> in set <a>_Vars</a> and update list of dependencies in both operands in <a>_Vars</a>
  template <typename U> static FFVar&
  _insert_binary_operation
    ( const typename FFOp::TYPE top, const FFDep&dep, const U&Cst1, const FFVar&Var2 );
  //! @brief Looks for the binary operation of type <a>top</a> with left and right operands <a>Var1</a>, <a>Cst2</a> in set <a>_Ops</a> and adds it if not found; also adds new auxiliary variable as well as constant <a>Cst2</a> in set <a>_Vars</a> and update list of dependencies in both operands in <a>_Vars</a>
  template <typename U> static FFVar&
  _insert_binary_operation
    ( const typename FFOp::TYPE top, const FFDep&dep, const FFVar&Var1, const U&Cst2 );
  //! @brief Looks for the unary operation of type <a>top</a> with operand <a>Var1</a>, <a>Var2</a> in set <a>_Ops</a> and adds it if not found; also adds new auxiliary variable in set <a>_Vars</a> and update list of dependencies in both operands in <a>_Vars</a>
  static FFVar& _insert_unary_operation
    ( const typename FFOp::TYPE top, const FFDep&dep, const FFVar&Var );

  //! @brief Adds the auxiliary variable with dependency <a>dep</a> from operation <a>op</a>
  FFVar* _add_auxiliary
    ( const FFDep&dep, FFOp*pOp );
  //! @brief Looks for the real constant <a>x</a> and adds it if not found
  FFVar* _add_constant
    ( const double x );
  //! @brief Looks for the integer constant <a>n</a> and adds it if not found
  FFVar* _add_constant
    ( const int n );

  //! @brief Erase all variables in _Vars
  void _clear_variables()
    { it_Vars itv = _Vars.begin();
      for( ; itv != _Vars.end(); ++itv ) delete *itv;
      _Vars.clear(); }
  //! @brief Appends the auxiliary variable <a>pAux</a> and define it in _Ops with type <a>tOp</a>
  void _append_aux
    ( FFVar*pAux, typename FFOp::TYPE tOp );
  //! @brief Appends new auxiliary variable
  virtual void _append_aux
    ( FFVar*pAux );
  //! @brief Appends new original variable
  virtual void _append_var
    ( FFVar*pVar );
  //! @brief Search for the variable with identify <a>id</a> in <a>_Vars</a>
  FFVar* _find_var
    ( const typename FFVar::pt_idVar&id );

private:
  //! @brief Private methods to block default compiler methods
  FFTree(const FFTree&);
  FFTree& operator=(const FFTree&);

};

////////////////////////////////// FFNum ///////////////////////////////////////

inline std::ostream&
operator <<
( std::ostream&out, const FFNum&Num )
{
  switch( Num.t ){
    case FFNum::INT:
      out << Num.n << " (I)"; break;
    case FFNum::REAL:
      out << Num.x << " (D)"; break;
  }
  return out;
}

////////////////////////////////// FFVar ///////////////////////////////////////

inline FFVar::FFVar
( FFTree*tree, const unsigned long ix )
: _id( VAR, ix ), _num( 0./0. ), _dep(), _tree( tree )
{ 
  if( !_tree ) throw typename FFTree::Exceptions( FFTree::Exceptions::INIT );
  // Keep track of variable count
  if( _tree->_nvar <= ix ) _tree->_nvar = ix+1;
  // Initialize dependence
  _dep.indep(ix);
  // Search for existing variable with same identifier in set FFTree::_Vars
  // If found, remove variable and corresponding operation in FFTree::_Ops
  FFVar* pVar = new FFVar( *this );
  typename FFTree::it_Vars iVar = _tree->_Vars.find( pVar );
  if( iVar!=_tree->_Vars.end() ){
    _tree->_Ops.erase( (*iVar)->_ops.first );
    delete (*iVar)->_ops.first;
    (*iVar)->_ops.second.clear();
    delete *iVar;
    _tree->_Vars.erase( iVar );
  }
  // Insert new variable in set FFTree::_Vars and corresponding operation in set FFTree::_Ops
  FFOp* pOp = new FFOp( FFOp::VAR, 0, 0, pVar );
  _tree->_Ops.insert( pOp );
  pVar->_ops.first = _ops.first = pOp;
  _tree->_append_var( pVar );
}

inline FFVar::FFVar
( FFTree*tree, const FFDep&dep, FFOp*op )
: _id( AUX, tree->_naux++ ), _num( 0./0. ), _dep(dep), _tree( tree )
{ _ops.first = op; }

inline FFVar::FFVar
( FFTree*tree, const pt_idVar&id )
: _id( id ), _num( 0 ), _dep(), _tree( tree )
{ _ops.first = 0; }

inline std::ostream&
operator <<
( std::ostream&out, const FFVar&Var)
{
  out << Var.name();
   // << " <= " << std::left << Var._num << "\t(" << Var._tree << ")";
  return out;
}

inline FFVar&
FFVar::operator=
( const FFVar&Var )
{
  if( this == &Var ) return *this;
  _id   = Var._id;
  _num  = Var._num;
  _dep  = Var._dep;
  _tree = Var._tree;
  _ops  = Var._ops;
  return *this;
}

inline FFVar&
FFVar::operator=
( const int n )
{
  _id   = std::make_pair(CINT,NOREF);
  _num  = n;
  _dep  = n;
  _tree = 0;
  _ops.first = 0;
  _ops.second.clear();
  return *this;
}

inline FFVar&
FFVar::operator=
( const double x )
{
  _id   = std::make_pair(CREAL,NOREF);
  _num  = x;
  _dep  = x;
  _tree = 0;
  _ops.first = 0;
  _ops.second.clear();
  return *this;
}

inline FFVar
operator+
( const FFVar&Var )
{
  return Var;
}

template <typename U> inline FFVar&
FFVar::operator+=
( const U&Var )
{
  FFVar VarNew = *this + Var;
  *this = VarNew;
  return *this;
}

inline FFVar
operator+
( const FFVar&Var1, const FFVar&Var2 )
{ 
  // Case either or both operands are (unreferenced) numeric constants
  if( Var1._id.second == FFVar::NOREF && Var2._id.second == FFVar::NOREF ){
    switch( Var1._num.t ){
      case FFNum::INT:
      switch( Var2._num.t ){
        case FFNum::INT:   return( Var1._num.n + Var2._num.n );
        case FFNum::REAL:  return( Var1._num.n + Var2._num.x );
      }
      case FFNum::REAL:
      switch( Var2._num.t ){
        case FFNum::INT:   return( Var1._num.x + Var2._num.n );
        case FFNum::REAL:  return( Var1._num.x + Var2._num.x );
      }
    }
  }
  if( Var1._id.second == FFVar::NOREF ){
    switch( Var1._num.t ){
      case FFNum::INT:   return( Var2 + Var1._num.n );
      case FFNum::REAL:  return( Var2 + Var1._num.x );
    }
  }
  if( Var2._id.second == FFVar::NOREF ){
    switch( Var2._num.t ){
      case FFNum::INT:   return( Var1 + Var2._num.n );
      case FFNum::REAL:  return( Var1 + Var2._num.x );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  return FFTree::_insert_binary_operation( FFOp::PLUS, Var1._dep+Var2._dep, Var1, Var2 );
}

template <typename U> inline FFVar
operator+
( const U&Cst1, const FFVar&Var2 )
{
  // Case right operand is (unreferenced) numeric constrants
  if( Var2._id.second == FFVar::NOREF ){
    switch( Var2._num.t ){
      case FFNum::INT:   return( Cst1 + Var2._num.n );
      case FFNum::REAL:  return( Cst1 + Var2._num.x );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  // Also append constant Cst1 if not defined
  return FFTree::_insert_binary_operation( FFOp::PLUS, Cst1+Var2._dep, Cst1, Var2 );
}

template <typename U> inline FFVar
operator+
( const FFVar&Var1, const U&Cst2 )
{
  return( Cst2 + Var1 );
}

inline FFVar
operator-
( const FFVar&Var )
{
  // Case right operand is a numeric constant
  if( Var._id.second == FFVar::NOREF ){
    switch( Var._num.t ){
      case FFNum::INT:   return( -Var._num.n );
      case FFNum::REAL:  return( -Var._num.x );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  return FFTree::_insert_unary_operation( FFOp::NEG, -Var._dep, Var );
}

template <typename U> inline FFVar&
FFVar::operator-=
( const U&Var )
{
  FFVar VarNew = *this - Var;
  *this = VarNew;
  return *this;
}

inline FFVar
operator-
( const FFVar&Var1, const FFVar&Var2 )
{
  if( &Var1 == &Var2 ) return 0.;

  // Case either or both operands are numeric constants
  if( Var1._id.second == FFVar::NOREF
   && Var2._id.second == FFVar::NOREF ){
    switch( Var1._num.t ){
      case FFNum::INT:
      switch( Var2._num.t ){
        case FFNum::INT:   return( Var1._num.n - Var2._num.n );
        case FFNum::REAL:  return( Var1._num.n - Var2._num.x );
      }
      case FFNum::REAL:
      switch( Var2._num.t ){
        case FFNum::INT:   return( Var1._num.x - Var2._num.n );
        case FFNum::REAL:  return( Var1._num.x - Var2._num.x );
      }
    }
  }

  if( Var1._id.second == FFVar::NOREF ){
    switch( Var1._num.t ){
      case FFNum::INT:   return( (double)Var1._num.n - Var2 );
      case FFNum::REAL:  return( Var1._num.x - Var2 );
    }
  }
  
  if( Var2._id.second == FFVar::NOREF ){
    switch( Var2._num.t ){
      case FFNum::INT:   return( Var1 - (double)Var2._num.n );
      case FFNum::REAL:  return( Var1 - Var2._num.x );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  return FFTree::_insert_binary_operation( FFOp::MINUS, Var1._dep-Var2._dep, Var1, Var2 );
}

template <typename U> inline FFVar
operator-
( const FFVar&Var1, const U&Cst2 )
{
  return( Var1 + (-Cst2) );
}

template <typename U> inline FFVar
operator-
( const U&Cst1, const FFVar&Var2 )
{
  // Case right operand is a numeric constant
  if( Var2._id.second == FFVar::NOREF ){
    switch( Var2._num.t ){
      case FFNum::INT:   return( Cst1 - Var2._num.n );
      case FFNum::REAL:  return( Cst1 - Var2._num.x );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  // Also append constant Cst1 if not defined
  return FFTree::_insert_binary_operation( FFOp::MINUS, Cst1-Var2._dep, Cst1, Var2 );
}

template <typename U> inline FFVar&
FFVar::operator*=
( const U&Var )
{
  FFVar VarNew = *this * Var;
  *this = VarNew;
  return *this;
}

inline FFVar
operator*
( const FFVar&Var1, const FFVar&Var2 )
{
  if( &Var1 == &Var2 ) return sqr(Var1);

  // Case either or both operands are numeric constants
  if( Var1._id.second == FFVar::NOREF && Var2._id.second == FFVar::NOREF ){
    switch( Var1._num.t ){
      case FFNum::INT:
      switch( Var2._num.t ){
        case FFNum::INT:   return( Var1._num.n * Var2._num.n );
        case FFNum::REAL:  return( Var1._num.n * Var2._num.x );
      }
      case FFNum::REAL:
      switch( Var2._num.t ){
        case FFNum::INT:   return( Var1._num.x * Var2._num.n );
        case FFNum::REAL:  return( Var1._num.x * Var2._num.x );
      }
    }
  }

  if( Var1._id.second == FFVar::NOREF ){
    switch( Var1._num.t ){
      case FFNum::INT:   return( (double)Var1._num.n * Var2 );
      case FFNum::REAL:  return( Var1._num.x * Var2 );
    }
  }
  
  if( Var2._id.second == FFVar::NOREF ){
    switch( Var2._num.t ){
      case FFNum::INT:   return( Var1 * (double)Var2._num.n );
      case FFNum::REAL:  return( Var1 * Var2._num.x );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  return FFTree::_insert_binary_operation( FFOp::TIMES, Var1._dep*Var2._dep, Var1, Var2 );
}

template <typename U> inline FFVar
operator*
( const FFVar&Var1, const U&Cst2 )
{
  return( Cst2 * Var1 );
}

template <typename U> inline FFVar
operator*
( const U&Cst1, const FFVar&Var2 )
{
  // Case right operand is a numeric constant
  if( Var2._id.second == FFVar::NOREF ){
    switch( Var2._num.t ){
      case FFNum::INT:   return( Cst1 * Var2._num.n );
      case FFNum::REAL:  return( Cst1 * Var2._num.x );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  // Also append constant Cst1 if not defined
  return FFTree::_insert_binary_operation( FFOp::SCALE, Cst1*Var2._dep, Cst1, Var2 );
}

template <typename U> inline FFVar&
FFVar::operator/=
( const U&Var )
{
  FFVar VarNew = *this / Var;
  *this = VarNew;
  return *this;
}

inline FFVar
operator/
( const FFVar&Var1, const FFVar&Var2 )
{
  if( &Var1 == &Var2 ) return 1.;

  // Case either or both operands are numeric constants
  if( Var1._id.second == FFVar::NOREF && Var2._id.second == FFVar::NOREF ){
    switch( Var1._num.t ){
      case FFNum::INT:
      switch( Var2._num.t ){
        case FFNum::INT:   return( Var1._num.n / Var2._num.n );
        case FFNum::REAL:  return( Var1._num.n / Var2._num.x );
      }
      case FFNum::REAL:
      switch( Var2._num.t ){
        case FFNum::INT:   return( Var1._num.x / Var2._num.n );
        case FFNum::REAL:  return( Var1._num.x / Var2._num.x );
      }
    }
  }

  if( Var1._id.second == FFVar::NOREF ){
    switch( Var1._num.t ){
      case FFNum::INT:   return( Var1._num.n / Var2 );
      case FFNum::REAL:  return( Var1._num.x / Var2 );
    }
  }
  
  if( Var2._id.second == FFVar::NOREF ){
    switch( Var2._num.t ){
      case FFNum::INT:   return( Var1 / Var2._num.n );
      case FFNum::REAL:  return( Var1 / Var2._num.x );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  return FFTree::_insert_binary_operation( FFOp::DIV, Var1._dep/Var2._dep, Var1, Var2 );
}

template <typename U> inline FFVar
operator/
( const FFVar&Var1, const U&Cst2 )
{
  return( ( 1 / Cst2 ) * Var1 );
}

template <typename U> inline FFVar
operator/
( const U&Cst1, const FFVar&Var2 )
{
  // Case right operand is a numeric constant
  if( Var2._id.second == FFVar::NOREF ){
    switch( Var2._num.t ){
      case FFNum::INT:   return( Cst1 / Var2._num.n );
      case FFNum::REAL:  return( Cst1 / Var2._num.x );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  // Also append constant Cst1 if not defined
  return FFTree::_insert_binary_operation( FFOp::DIV, Cst1/Var2._dep, Cst1, Var2 );
}

inline FFVar
inv
( const FFVar&Var )
{
  return( 1 / Var );
}

inline FFVar
max
( const FFVar&Var1, const FFVar&Var2 )
{
  if( &Var1 == &Var2 ) return Var1;

  // Case either or both operands are numeric constants
  if( Var1._id.second == FFVar::NOREF && Var2._id.second == FFVar::NOREF ){
    switch( Var1._num.t ){
      case FFNum::INT:
      switch( Var2._num.t ){
        case FFNum::INT:   return( Var1._num.n>Var2._num.n? Var1._num.n: Var2._num.n );
        case FFNum::REAL:  return( std::max((double)Var1._num.n,Var2._num.x) );
      }
      case FFNum::REAL:
      switch( Var2._num.t ){
        case FFNum::INT:   return( std::max(Var1._num.x,(double)Var2._num.n) );
        case FFNum::REAL:  return( std::max(Var1._num.x,Var2._num.x) );
      }
    }
  }

  if( Var1._id.second == FFVar::NOREF ){
    switch( Var1._num.t ){
      case FFNum::INT:   return( max((double)Var1._num.n,Var2) );
      case FFNum::REAL:  return( max(Var1._num.x,Var2) );
    }
  }
  
  if( Var2._id.second == FFVar::NOREF ){
    switch( Var2._num.t ){
      case FFNum::INT:   return( max((double)Var2._num.n,Var1) );
      case FFNum::REAL:  return( max(Var2._num.x,Var1) );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  return FFTree::_insert_binary_operation( FFOp::MAXF, max(Var1._dep,Var2._dep), Var1, Var2 );
}

template <typename U> inline FFVar
max
( const FFVar&Var1, const U&Cst2 )
{
  return( max( Cst2, Var1 ) );
}

template <typename U> inline FFVar
max
( const U&Cst1, const FFVar&Var2 )
{
  // Case right operand is a numeric constant
  if( Var2._id.second == FFVar::NOREF ){
    switch( Var2._num.t ){
      case FFNum::INT:   return( std::max(Cst1,(double)Var2._num.n) );
      case FFNum::REAL:  return( std::max(Cst1,Var2._num.x) );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  // Also append constant Cst1 if not defined
  return FFTree::_insert_binary_operation( FFOp::MAXF, max((double)Cst1,Var2._dep), Cst1, Var2 );
}

inline FFVar
max
( const unsigned nVar, const FFVar*pVar )
{
  if( nVar<=0 || !pVar ) return( 0 );
  
  FFVar VarR = pVar[0];
  for( unsigned int i=1; i<nVar; i++ ) VarR = max( VarR, pVar[i] );
  return( VarR );
}

inline FFVar
min
( const FFVar&Var1, const FFVar&Var2 )
{
  if( &Var1 == &Var2 ) return Var1;

  // Case either or both operands are numeric constants
  if( Var1._id.second == FFVar::NOREF && Var2._id.second == FFVar::NOREF ){
    switch( Var1._num.t ){
      case FFNum::INT:
      switch( Var2._num.t ){
        case FFNum::INT:   return( Var1._num.n<Var2._num.n? Var1._num.n: Var2._num.n );
        case FFNum::REAL:  return( std::min((double)Var1._num.n,Var2._num.x) );
      }
      case FFNum::REAL:
      switch( Var2._num.t ){
        case FFNum::INT:   return( std::min(Var1._num.x,(double)Var2._num.n) );
        case FFNum::REAL:  return( std::min(Var1._num.x,Var2._num.x) );
      }
    }
  }

  if( Var1._id.second == FFVar::NOREF ){
    switch( Var1._num.t ){
      case FFNum::INT:   return( min((double)Var1._num.n,Var2) );
      case FFNum::REAL:  return( min(Var1._num.x,Var2) );
    }
  }
  
  if( Var2._id.second == FFVar::NOREF ){
    switch( Var2._num.t ){
      case FFNum::INT:   return( min((double)Var2._num.n,Var1) );
      case FFNum::REAL:  return( min(Var2._num.x,Var1) );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  return FFTree::_insert_binary_operation( FFOp::MINF, max(Var1._dep,Var2._dep), Var1, Var2 );
}

template <typename U> inline FFVar
min
( const FFVar&Var1, const U&Cst2 )
{
  return( min( Cst2, Var1 ) );
}

template <typename U> inline FFVar
min
( const U&Cst1, const FFVar&Var2 )
{
  // Case right operand is a numeric constant
  if( Var2._id.second == FFVar::NOREF ){
    switch( Var2._num.t ){
      case FFNum::INT:   return( std::min(Cst1,(double)Var2._num.n) );
      case FFNum::REAL:  return( std::min(Cst1,Var2._num.x) );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  // Also append constant Cst1 if not defined
  return FFTree::_insert_binary_operation( FFOp::MINF, min((double)Cst1,Var2._dep), Cst1, Var2 );
}

inline FFVar
min
( const unsigned nVar, const FFVar*pVar )
{
  if( nVar<=0 || !pVar ) return( 0 );
  
  FFVar VarR = pVar[0];
  for( unsigned int i=1; i<nVar; i++ ) VarR = min( VarR, pVar[i] );
  return( VarR );
}

inline FFVar
exp
( const FFVar&Var )
{
  // Case operand is a numeric constant
  if( Var._id.second == FFVar::NOREF ){
    switch( Var._num.t ){
      case FFNum::INT:   return( std::exp( Var._num.n ) );
      case FFNum::REAL:  return( std::exp( Var._num.x ) );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  return FFTree::_insert_unary_operation( FFOp::EXP, exp(Var._dep), Var );
}

inline FFVar
log
( const FFVar&Var )
{
  // Case operand is a numeric constant
  if( Var._id.second == FFVar::NOREF ){
    switch( Var._num.t ){
      case FFNum::INT:   return( std::log( Var._num.n ) );
      case FFNum::REAL:  return( std::log( Var._num.x ) );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  return FFTree::_insert_unary_operation( FFOp::LOG, log(Var._dep), Var );
}

inline FFVar
sqr
( const FFVar&Var )
{
  // Case operand is a numeric constant
  if( Var._id.second == FFVar::NOREF ){
    switch( Var._num.t ){
      case FFNum::INT:   return( mc::sqr( Var._num.n ) );
      case FFNum::REAL:  return( mc::sqr( Var._num.x ) );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  return FFTree::_insert_unary_operation( FFOp::SQR, sqr(Var._dep), Var );
}

inline FFVar
sqrt
( const FFVar&Var )
{
  // Case operand is a numeric constant
  if( Var._id.second == FFVar::NOREF ){
    switch( Var._num.t ){
      case FFNum::INT:   return( std::sqrt( Var._num.n ) );
      case FFNum::REAL:  return( std::sqrt( Var._num.x ) );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  return FFTree::_insert_unary_operation( FFOp::SQRT, sqrt(Var._dep), Var );
}

inline FFVar
fabs
( const FFVar&Var )
{
  // Case operand is a numeric constant
  if( Var._id.second == FFVar::NOREF ){
    switch( Var._num.t ){
      case FFNum::INT:   return( std::fabs( Var._num.n ) );
      case FFNum::REAL:  return( std::fabs( Var._num.x ) );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  return FFTree::_insert_unary_operation( FFOp::FABS, fabs(Var._dep), Var );
}

inline FFVar
pow
( const FFVar&Var, const int iExp )
{
  // Case operand is a numeric constant
  if( Var._id.second == FFVar::NOREF ){
    switch( Var._num.t ){
      case FFNum::INT:   return( std::pow( Var._num.n, iExp ) );
      case FFNum::REAL:  return( std::pow( Var._num.x, iExp ) );
    }
  }

  // Case integer exponent is 0,1, or 2
  if( iExp == 0 ) return( 1. );
  if( iExp == 1 ) return Var;
  if( iExp == 2 ) return sqr(Var);

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  // Also append constant iExp if not defined
  return FFTree::_insert_binary_operation( FFOp::IPOW, pow(Var._dep,iExp), Var, iExp );
}

template <typename U> inline FFVar
pow
( const FFVar&Var1, const U&Var2 )
{
  return exp( Var2 * log( Var1 ) );
}

inline FFVar
pow
( const double Var1, const FFVar&Var2 )
{
  return exp( Var2 * std::log( Var1 ) );
}

inline FFVar
cos
( const FFVar&Var )
{
  // Case operand is a numeric constant
  if( Var._id.second == FFVar::NOREF ){
    switch( Var._num.t ){
      case FFNum::INT:   return( std::cos( Var._num.n ) );
      case FFNum::REAL:  return( std::cos( Var._num.x ) );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  return FFTree::_insert_unary_operation( FFOp::COS, cos(Var._dep), Var );
}

inline FFVar
sin
( const FFVar&Var )
{
  // Case operand is a numeric constant
  if( Var._id.second == FFVar::NOREF ){
    switch( Var._num.t ){
      case FFNum::INT:   return( std::sin( Var._num.n ) );
      case FFNum::REAL:  return( std::sin( Var._num.x ) );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  return FFTree::_insert_unary_operation( FFOp::SIN, sin(Var._dep), Var );
}

inline FFVar
tan
( const FFVar&Var )
{
  // Case operand is a numeric constant
  if( Var._id.second == FFVar::NOREF ){
    switch( Var._num.t ){
      case FFNum::INT:   return( std::tan( Var._num.n ) );
      case FFNum::REAL:  return( std::tan( Var._num.x ) );

    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  return FFTree::_insert_unary_operation( FFOp::TAN, tan(Var._dep), Var );
}

inline FFVar
asin
( const FFVar&Var )
{
  // Case operand is a numeric constant
  if( Var._id.second == FFVar::NOREF ){
    switch( Var._num.t ){
      case FFNum::INT:   return( std::asin( Var._num.n ) );
      case FFNum::REAL:  return( std::asin( Var._num.x ) );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  return FFTree::_insert_unary_operation( FFOp::ASIN, asin(Var._dep), Var );
}

inline FFVar
acos
( const FFVar&Var )
{
  // Case operand is a numeric constant
  if( Var._id.second == FFVar::NOREF ){
    switch( Var._num.t ){
      case FFNum::INT:   return( std::acos( Var._num.n ) );
      case FFNum::REAL:  return( std::acos( Var._num.x ) );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  return FFTree::_insert_unary_operation( FFOp::ACOS, acos(Var._dep), Var );
}

inline FFVar
atan
( const FFVar&Var )
{
  // Case operand is a numeric constant
  if( Var._id.second == FFVar::NOREF ){
    switch( Var._num.t ){
      case FFNum::INT:   return( std::atan( Var._num.n ) );
      case FFNum::REAL:  return( std::atan( Var._num.x ) );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  return FFTree::_insert_unary_operation( FFOp::ATAN, atan(Var._dep), Var );
}

inline FFVar
erf
( const FFVar&Var )
{
  // Case operand is a numeric constant
  if( Var._id.second == FFVar::NOREF ){
    switch( Var._num.t ){
      case FFNum::INT:   return( ::erf( Var._num.n ) );
      case FFNum::REAL:  return( ::erf( Var._num.x ) );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  return FFTree::_insert_unary_operation( FFOp::ERF, erf(Var._dep), Var );
}

inline FFVar
erfc
( const FFVar&Var )
{
  return ( 1. - erf( Var ) );
}

inline FFVar
fstep
( const FFVar&Var )
{
  // Case operand is a numeric constant
  if( Var._id.second == FFVar::NOREF ){
    switch( Var._num.t ){
      case FFNum::INT:   return( mc::fstep( Var._num.n ) );
      case FFNum::REAL:  return( mc::fstep( Var._num.x ) );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if operation does not exist already)
  return FFTree::_insert_unary_operation( FFOp::FSTEP, fstep(Var._dep), Var );
}

inline FFVar
bstep
( const FFVar&Var )
{
  return ( fstep( -Var ) );
}


////////////////////////////////// FFOp ////////////////////////////////////////

inline
FFOp::FFOp
( TYPE top, FFVar*lop, FFVar*rop, FFVar*res ):
  type( top ), pres( res ), _visited(false)
{
  // Order operands for commutative binary operators
  if( !lop
    || ( top != PLUS && top != TIMES && top != SCALE )
    || lt_FFVar()( lop, rop ) )
    { plop = lop; prop = rop; }
  else
    { plop = rop; prop = lop; }
}

inline std::ostream&
operator <<
( std::ostream&out, const FFOp&Op)
{
  switch( Op.type ){
    case FFOp::CNST:  out << Op.pres->num(); break;
    case FFOp::VAR:   out << "VARIABLE"; break;
    case FFOp::PLUS:  out << FFVar::_name( Op.plop->id() ) << " + " << FFVar::_name( Op.prop->id() ); break;
    case FFOp::NEG:   out << "- " << FFVar::_name( Op.plop->id() ); break;
    case FFOp::MINUS: out << FFVar::_name( Op.plop->id() ) << " - " << FFVar::_name( Op.prop->id() ); break;
    case FFOp::TIMES:
    case FFOp::SCALE: out << FFVar::_name( Op.plop->id() ) << " * " << FFVar::_name( Op.prop->id() ); break;
    case FFOp::DIV:   out << FFVar::_name( Op.plop->id() ) << " / " << FFVar::_name( Op.prop->id() ); break;
    case FFOp::MINF:  out << "MIN( " << FFVar::_name( Op.plop->id() ) << ", " << FFVar::_name( Op.prop->id() ) << " )"; break;
    case FFOp::MAXF:  out << "MAX( " << FFVar::_name( Op.plop->id() ) << ", " << FFVar::_name( Op.prop->id() ) << " )"; break;
    case FFOp::EXP:   out << "EXP( " << FFVar::_name( Op.plop->id() ) << " )"; break;
    case FFOp::LOG:   out << "LOG( " << FFVar::_name( Op.plop->id() ) << " )"; break;
    case FFOp::SQR:   out << "SQR( " << FFVar::_name( Op.plop->id() ) << " )"; break;
    case FFOp::SQRT:  out << "SQRT( " << FFVar::_name( Op.plop->id() ) << " )"; break;
    case FFOp::FABS:  out << "FABS( " << FFVar::_name( Op.plop->id() ) << " )"; break;
    case FFOp::IPOW:  out << "POW( " << FFVar::_name( Op.plop->id() ) << ", " << FFVar::_name( Op.prop->id() ) << " )"; break;
    case FFOp::COS:   out << "COS( " << FFVar::_name( Op.plop->id() ) << " )"; break;
    case FFOp::SIN:   out << "SIN( " << FFVar::_name( Op.plop->id() ) << " )"; break;
    case FFOp::TAN:   out << "TAN( " << FFVar::_name( Op.plop->id() ) << " )"; break;
    case FFOp::ASIN:  out << "ASIN( " << FFVar::_name( Op.plop->id() ) << " )"; break;
    case FFOp::ACOS:  out << "ACOS( " << FFVar::_name( Op.plop->id() ) << " )"; break;
    case FFOp::ATAN:  out << "ATAN( " << FFVar::_name( Op.plop->id() ) << " )"; break;
    case FFOp::ERF:   out << "ERF( " << FFVar::_name( Op.plop->id() ) << " )"; break;
    case FFOp::FSTEP: out << "FSTEP( " << FFVar::_name( Op.plop->id() ) << " )"; break;
    default:;
  } 
  return out;
}

inline void
FFOp::propagate_subtree
( std::list<FFOp*>&Ops )
{
  if( _visited ) return;
  _visited = true;

  //Ops.push_front( this );
  if( plop && plop->ops().first ) plop->ops().first->propagate_subtree( Ops );
  if( prop && prop->ops().first ) prop->ops().first->propagate_subtree( Ops );
  Ops.push_back( this );
}

inline void
FFOp::generate_dot_script
( std::ostream&os ) const
{
  if( _visited ) return;
  _visited = true;

  if( plop && plop->ops().first ) plop->ops().first->generate_dot_script( os );
  if( prop && prop->ops().first ) prop->ops().first->generate_dot_script( os );
  append_dot_script( os );
}

inline void
FFOp::append_dot_script
( std::ostream&os ) const
{
  std::ostringstream var_color; var_color << "red";
  std::ostringstream aux_color; aux_color << "blue";
  std::ostringstream op_color;  op_color  << "black";
  switch( type ){
  case FFOp::VAR:   return append_dot_script_variable( os, false, 14 );
  case FFOp::CNST:  return append_dot_script_variable( os, true,  14 );
  case FFOp::PLUS:  return append_dot_script_factor( os, " + ",   false, 18 );
  case FFOp::NEG:   return append_dot_script_factor( os, " - ",   true,  18 );
  case FFOp::MINUS: return append_dot_script_factor( os, " - ",   false, 18 );
  case FFOp::SCALE: return append_dot_script_factor( os, " x ",   false, 18 );
  case FFOp::TIMES: return append_dot_script_factor( os, " x ",   false, 18 );
  case FFOp::DIV:   return append_dot_script_factor( os, " / ",   false, 18 );
  case FFOp::MINF:  return append_dot_script_factor( os, "min",   false, 14 );
  case FFOp::MAXF:  return append_dot_script_factor( os, "max",   false, 14 );
  case FFOp::IPOW:  return append_dot_script_factor( os, "pow",   false, 14, true );
  case FFOp::EXP:   return append_dot_script_factor( os, "exp",   true,  14 );
  case FFOp::LOG:   return append_dot_script_factor( os, "log",   true,  14 );
  case FFOp::FABS:  return append_dot_script_factor( os, "fabs",  true,  14 );
  case FFOp::SQR:   return append_dot_script_factor( os, "sqr",   true,  14 );
  case FFOp::SQRT:  return append_dot_script_factor( os, "sqrt",  true,  14 );
  case FFOp::COS:   return append_dot_script_factor( os, "cos",   true,  14 );
  case FFOp::SIN:   return append_dot_script_factor( os, "sin",   true,  14 );
  case FFOp::TAN:   return append_dot_script_factor( os, "tan",   true,  14 );
  case FFOp::ACOS:  return append_dot_script_factor( os, "acos",  true,  14 );
  case FFOp::ASIN:  return append_dot_script_factor( os, "asin",  true,  14 );
  case FFOp::ATAN:  return append_dot_script_factor( os, "atan",  true,  14 );
  case FFOp::ERF:   return append_dot_script_factor( os, "erf",   true,  14 );
  case FFOp::FSTEP: return append_dot_script_factor( os, "fstep", true,  14 );
  default: os << "/* a factor was not displayed */\n";
  }
  return;
}

inline void
FFOp::append_dot_script_factor
( std::ostream&os, const std::string&fname, const bool unary,
  const unsigned int fontsize, const bool dotted ) const
{
  std::ostringstream aux_color; aux_color << "blue";
  std::ostringstream op_color;  op_color  << "black";

  os << "  " << pres->name() << " [shape=box,fontname=\"Arial\",color=" << aux_color.str().c_str()
     << ",peripheries=" << (pres->ops().second.empty()?2:1) << "];\n";
  std::ostringstream oop; oop << pres->name() << "_op";
  os << "  " << oop.str().c_str() << " [shape=ellipse,fontname=\"Arial\",fontsize="
     << fontsize << ",label=\"" << fname.c_str() << "\",color=" << op_color.str().c_str()
     << "];\n";
  os << "  " << oop.str().c_str() << " -> " << pres->name() << ";\n";
  os << "  " << plop->name() << " -> " << oop.str().c_str() << ";\n";
  if( unary ) return;
  os << "  " << prop->name() << " -> " << oop.str().c_str()
     << (dotted? "[style=dotted];\n": ";\n");
}

inline void
FFOp::append_dot_script_variable
( std::ostream&os, const bool constant, const unsigned int fontsize ) const
{
  std::ostringstream var_color; var_color << "red";
  std::ostringstream aux_color; aux_color << "blue";

  os << "  " << pres->name() << " [shape=box,fontname=\"Arial\",color="
     << (constant? aux_color.str().c_str(): var_color.str().c_str())
     << ",peripheries=" << (pres->ops().second.empty()?2:1) << "];\n";
  if( !constant) return;
  std::ostringstream ocst; ocst << pres->name() << "_0";
  os << "  " << ocst.str().c_str() << " [shape=ellipse,fontname=\"Arial\",label=\""
     << pres->num() << "\",color=" << aux_color.str().c_str() << "];\n";
  os << "  " << ocst.str().c_str() << " -> " << pres->name()
     << " [style=dotted,color=" << aux_color.str().c_str() << "];\n";
}

inline bool
FFOp::is_univariate
()
{
  switch( type ){
  case FFOp::PLUS: case FFOp::MINUS: case FFOp::TIMES: case FFOp::DIV:
  case FFOp::MINF: case FFOp::MAXF:
    return false;
  case FFOp::NEG:  case FFOp::SCALE: case FFOp::IPOW:  case FFOp::EXP:   
  case FFOp::LOG:  case FFOp::FABS:  case FFOp::SQR:   case FFOp::SQRT:  
  case FFOp::COS:  case FFOp::SIN:   case FFOp::TAN:   case FFOp::ACOS:  
  case FFOp::ASIN: case FFOp::ATAN:  case FFOp::ERF:   case FFOp::FSTEP:
  default:
    return true;
  }
}


///////////////////////////////// FFTree //////////////////////////////////////

inline std::ostream&
operator <<
( std::ostream&out, const FFTree&Tree)
{
  typename FFTree::t_Vars Vars = Tree._Vars;
  typename FFTree::it_Vars itv = Vars.begin();

  out << ( Tree._nvar? "\nVARIABLES:\n": "\nNO VARIABLES\n" );
  for( ; itv!=Vars.end() && (*itv)->_id.first<=FFVar::VAR; ++itv ){
    out << "  " << **itv;
    out << "\t => {";
    typename FFVar::t_Ops Ops = (*itv)->_ops.second;
    typename FFVar::t_Ops::iterator ito = Ops.begin();
    for( ; ito!=Ops.end(); ++ito ) out << " " << *(*ito)->pres;
    out << " }" << std::endl;
  }

  out << ( Tree._naux? "\nINTERMEDIATES:\n": "\nNO INTERMEDIATES\n" );
  for( ; itv!=Vars.end(); ++itv ){
    out << "  " << **itv;
    if( (*itv)->_ops.first ) out << "\t" << "<=  " << *((*itv)->_ops.first);
    out << "\t => {";
    typename FFVar::t_Ops Ops = (*itv)->_ops.second;
    typename FFVar::t_Ops::iterator ito = Ops.begin();
    for( ; ito!=Ops.end(); ++ito ) out << " " << *(*ito)->pres;
    out << " }" << std::endl;
  }

  return out;
}

inline FFVar&
FFTree::_insert_binary_operation
( const typename FFOp::TYPE top, const FFDep&dep, const FFVar&Var1, const FFVar&Var2 )
{
  if( Var1._tree != Var2._tree ) throw Exceptions( Exceptions::TREE );
  FFTree *pTree = Var1._tree;
  FFVar *pVar1 = Var1._ops.first->pres, *pVar2 = Var2._ops.first->pres;

  FFOp *pOp = pTree->_insert_operation( top, pVar1, pVar2 );
  if( pOp->pres ) return *pOp->pres;
  pVar1->_ops.second.push_back( pOp );
  pVar2->_ops.second.push_back( pOp );
  pOp->pres = pTree->_add_auxiliary( dep, pOp );
  return *pOp->pres;
}

template <typename U> inline FFVar&
FFTree::_insert_binary_operation
( const typename FFOp::TYPE top, const FFDep&dep, const U&Cst1, const FFVar&Var2 )
{
  FFTree *pTree = Var2._tree;
  FFVar *pVar2 = Var2._ops.first->pres;
  FFVar* pCst1 = pTree->_add_constant( Cst1 );
  FFVar *pVar1 = pCst1->_ops.first->pres;

  FFOp *pOp = pTree->_insert_operation( top, pVar1, pVar2 );
  if( pOp->pres ) return *pOp->pres;
  pVar1->_ops.second.push_back( pOp );
  pVar2->_ops.second.push_back( pOp );
  pOp->pres = pTree->_add_auxiliary( dep, pOp );
  return *pOp->pres;
}

template <typename U> inline FFVar&
FFTree::_insert_binary_operation
( const typename FFOp::TYPE top, const FFDep&dep, const FFVar&Var1, const U&Cst2 )
{
  FFTree *pTree = Var1._tree;
  FFVar *pVar1 = Var1._ops.first->pres;
  FFVar* pCst2 = pTree->_add_constant( Cst2 );
  FFVar *pVar2 = pCst2->_ops.first->pres;

  FFOp *pOp = pTree->_insert_operation( top, pVar1, pVar2 );
  if( pOp->pres ) return *pOp->pres;
  pVar1->_ops.second.push_back( pOp );
  pVar2->_ops.second.push_back( pOp );
  pOp->pres = pTree->_add_auxiliary( dep, pOp );
  return *pOp->pres;
}

inline FFVar&
FFTree::_insert_unary_operation
( const typename FFOp::TYPE top, const FFDep&dep, const FFVar&Var )
{
  FFTree* pTree = Var._tree;
  FFVar *pVar = Var._ops.first->pres;

  FFOp* pOp = pTree->_insert_operation( top, Var._ops.first->pres );
  if( pOp->pres ) return *pOp->pres;
  pVar->_ops.second.push_back( pOp );
  pOp->pres = pTree->_add_auxiliary( dep, pOp );
  return *pOp->pres;
}

inline FFOp*
FFTree::_insert_operation
( const typename FFOp::TYPE top, FFVar*lop, FFVar*rop )
{
  FFOp* op = new FFOp( top, lop, rop );
  typename FFTree::it_Ops itop = _Ops.find( op );
  if( itop!=_Ops.end() ){ delete op; return *itop; }
  _Ops.insert( op );
  return op;
}

inline bool
FFTree::_remove_operation
( FFOp* op )
{
  typename FFTree::it_Ops itop = _Ops.find( op );
  if( itop==_Ops.end() ) return false;
  delete op;
  _Ops.erase( itop );
  return true;
}
 
inline FFVar*
FFTree::_add_auxiliary
( const FFDep&dep, FFOp*pOp )
{
  pOp->pres = new FFVar( this, dep, pOp );
  _append_aux( pOp->pres );
  return pOp->pres;
}

inline FFVar*
FFTree::_add_constant
( const double x )
{
  // Check if real constant x already defined in _Vars
  FFVar* pAux = new FFVar( x );
  it_Vars iAux = _Vars.find( pAux );
  if( iAux!=_Vars.end() ){ delete pAux; return *iAux; }

  // Otherwise, append constant x
  _append_aux( pAux, FFOp::CNST );
  return pAux;
}

inline FFVar*
FFTree::_add_constant
( const int n )
{
  // Check if real constant x already defined in _Vars
  FFVar* pAux = new FFVar( n );
  it_Vars iAux = _Vars.find( pAux );
  if( iAux!=_Vars.end() ){ delete pAux; return *iAux; }

  // Otherwise, append constant n
  _append_aux( pAux, FFOp::CNST );
  return pAux;
}

inline void
FFTree::_append_aux
( FFVar*pAux, typename FFOp::TYPE tOp )
{
  FFOp*pOp = new FFOp( tOp, 0, 0, pAux );
  _Ops.insert( pOp );
  pAux->tree() = this;
  pAux->ops().first = pOp;
  pAux->id().second = _naux++;
  _append_aux( pAux );
}

inline void
FFTree::_append_aux
( FFVar*pAux )
{
  _Vars.insert( pAux );
}

inline void
FFTree::_append_var
( FFVar*pVar )
{
  _Vars.insert( pVar );
}   

inline FFVar*
FFTree::_find_var
( const typename FFVar::pt_idVar&id )
{
  FFVar* pVar = new FFVar( this, id );
  it_Vars iVar = _Vars.find( pVar );
  delete pVar;
  return( iVar==_Vars.end()? 0: *iVar );
}

inline std::list<FFOp*>
FFTree::subtree
( const unsigned int nVar, const FFVar*pVar ) const
{
  std::list<FFOp*> Ops;
  _reset_operations();
  for( unsigned int iVar=0; iVar<nVar; iVar++ ){
    pVar[iVar].ops().first->propagate_subtree( Ops );
  }
  return Ops;
}

inline void
FFTree::output
( std::list<FFOp*>&Ops, std::ostream&os )
{
  os << ( !Ops.empty()? "\nFACTORS:\n": "\nNO FACTORS\n" );
  typename FFVar::t_Ops::iterator ito = Ops.begin();
  for( ; ito!=Ops.end(); ++ito )
    os << "  " << *(*ito)->pres << "\t" << "<=  " << **ito << std::endl;
}

inline void
FFTree::dot_script
( const unsigned int nVar, const FFVar*pVar, std::ostream&os ) const
{
  _reset_operations();
  os << "\ndigraph G {\n";
  for( unsigned int iVar=0; iVar<nVar; iVar++ ){
    pVar[iVar].ops().first->generate_dot_script( os );
  }
  os << "}\n";
}

} // namespace mc

#endif
