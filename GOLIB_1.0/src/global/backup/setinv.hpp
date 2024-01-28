// Copyright (C) 2012 Benoit Chachuat (b.chachuat@imperial.ac.uk)
// All rights reserved.

// This code is provided "as is", without any warranty of any kind,
// either expressed or implied, including but not limited to, any implied
// warranty of merchantibility or fitness for any purpose. In no event
// will any party who distributed the code be liable for damages or for
// any claim(s) by any other party, including but not limited to, any
// lost profits, lost monies, lost data or data rendered inaccurate,
// losses sustained by third parties, or any other special, incidental or
// consequential damages arising out of the use or inability to use the
// program, even if the possibility of such damages has been advised
// against. The entire risk as to the quality, the performance, and the
// fitness of the program for any particular purpose lies with the party
// using the code.

// This code, and any derivative of this code, may not be used in a
// commercial package without the prior explicit written permission of
// the authors. Verbatim copies of this code may be made and distributed
// in any medium, provided that this copyright notice is not removed or
// altered in any way. No fees may be charged for distribution of the
// codes, other than a fee to cover the cost of the media and a
// reasonable handling fee.

// ***************************************************************
// ANY USE OF THIS CODE CONSTITUTES ACCEPTANCE OF THE TERMS OF THE
//                         COPYRIGHT NOTICE
// ***************************************************************

#ifndef MC__SETINV_HPP
#define MC__SETINV_HPP

#include <utility>
#include <set>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "mcop.hpp"
#include "mcfunc.hpp"

// TO DO:

namespace mc
{

template <typename T> class SetInvNode;
template <typename T> struct lt_SetInvNode;

//! @brief Pure virtual base class for set inversion
////////////////////////////////////////////////////////////////////////
//! mc::SetInv<T> is a pure virutal base C++ class implementing a
//! rigorous set inversion algorithm for nonlinear, implicit functions.
//!
//! REFERENCES:
//! - Jaulin, L., and E. Walter, <A href="http://dx.doi.org/10.1016/0005-1098(93)90106-4">Set inversion via interval analysis for nonlinear bounded-error estimation</A>, <i>Automatica</i>, <b>29</b>(4):1053--1064, 1993.
//! - Lin, Y., and M. A. Stadtherr, <A href="http://dx.doi.org/10.1021/ie0707725"> Guaranteed state and parameter estimation for nonlinear continuous-time systems with bounded-error measurements,</A> <i>Industrial & Engineering Chemistry Research</I>, <B>46</B>:7198--7207, 2007.
//! - Kieffer, M., and E. Walter, <A href="http://dx.doi.org/10.1002/acs.1194">Guaranteed estimation of the parameters of nonlinear continuous-time models: Contributions of interval analysis</A>, <i>Intermation Journal of Adaptive Control & Signal Processing</i>, <b>25</b>:191--207, 2011.
//! .
////////////////////////////////////////////////////////////////////////
template <typename T>
class SetInv
////////////////////////////////////////////////////////////////////////
{
  template <typename U> friend class SetInvNode;

public:
  typedef std::multiset< SetInvNode<T>*, lt_SetInvNode<T> > t_Nodes;
  typedef typename t_Nodes::iterator it_Nodes;
  typedef typename t_Nodes::const_iterator cit_Nodes;

  //! @brief Status of subproblems
  enum STATUS{
    INNER=0,		//!< Current node is inside the inversed set
    OUTER,		//!< Current node is outside of the inversed set
    UNDETERMINED,	//!< Current node is currently undetermined
    FAILURE,		//!< Current node remains undetermined due to failure
    ABORT		//!< Current node terminated after fatal error
  };

  //! @brief Prototype for user-supplied function
  virtual STATUS assess
    ( const unsigned int, 	// number of variables
      T*			// variable bounds - possibly contracted on return
    )=0;

  //! @brief Public constructor
  SetInv()
    : options(), _np(0), _P_root(0), _P_tmp(0)
    {}

  //! @brief Class destructor
  virtual ~SetInv()
    {
      _clean_stacks();
      delete[] _P_root;
      delete[] _P_tmp;
    }

  //! @brief mc::SetInv options
  struct Options
  {
    //! @brief Constructor
    Options():
      ABSOLUTE_TOLERANCE(1e-5), RELATIVE_TOLERANCE(1e-5),
      BRANCHING_VARIABLE_CRITERION(RGREL), DISPLAY(2), MAX_CPU_TIME(1e6),
      MAX_NODES(0)
      {}
    //! @brief Display
    void display
      ( std::ostream&out ) const;
    //! @brief Branching variable criterion
    enum CRITERION{
      RGREL=0,	//!< Relative variable range diameter
      RGABS	//!< Absolute variable range diameter
    };
    //! @brief Absolute stopping tolerance
    double ABSOLUTE_TOLERANCE;
    //! @brief Relative stopping tolerance
    double RELATIVE_TOLERANCE;
    //! @brief Variable selection criterion
    CRITERION BRANCHING_VARIABLE_CRITERION;
    //! @brief Display option
    int DISPLAY;
    //! @brief Maximum CPU time limit
    double MAX_CPU_TIME;
    //! @brief Maximum number of nodes (0 for no limit)
    unsigned int MAX_NODES;
  } options;

  //! @brief mc::SetInv exceptions
  class Exceptions
  {
  public:
    //! @brief Enumeration type for SetInv exception handling
    enum TYPE{
      BRANCH=0,		//!< Error due to an empty set of branching variables 
      INTERN=-3,	//!< SetInv internal error
      UNDEF=-33		//!< Error due to calling a function/feature not yet implemented in SetInv
    };
    //! @brief Constructor for error <a>ierr</a>
    Exceptions( TYPE ierr ) : _ierr( ierr ){}
    //! @brief Inline function returning the error flag
    int ierr()
      { return _ierr; }
    std::string what()
      {
        switch( _ierr ){
        case BRANCH:
          return "SetInv Error: empty set of branching variables";
        case INTERN:
          return "SetInv Internal Error";
        case UNDEF: default:
          return "SetInv Error: calling a feature not yet implemented";
        }
      }
  private:
    TYPE _ierr;
  };

  //! @brief Set variables
  void variables
    ( const unsigned int np, const T*P,
      const std::set<unsigned int>*exclude=0 );

  //! @brief Retreive number of variables
  unsigned int np() const
    { return _np; }

  //! @brief Apply set-inversion algorithm - on return, volumes of inner- and boundary-approximations 
  std::pair<double,double> solve
    ( std::ostream&os=std::cout );

  //! @brief Append all open nodes and inner nodes to, respectively, <a>os_open</a> and <a>os_inner</a> - Number of significant digits is set via <a>DPREC</a> (default=6)
  void output_stacks
    ( std::ostream&os_open, std::ostream&os_inner,
      const unsigned int DPREC=6 ) const;

protected:  
  //! @brief Exclusion set for branching variables (e.g. non participating)
  std::set<unsigned int> _exclude_vars;

private:  
  //! @brief Node type
  enum NODETYPE{
    ROOT=0,	//!< Root node
    LEFT=-1,	//!< Left child node
    RIGHT=1	//!< Right child node
  };

  //! @brief Number of variables in optimization problem
  unsigned int _np;
  //! @brief Variable bounds at root node
  T *_P_root;
  //! @brief Partition volume at root node
  double _volume_root;
  //! @brief Default branching variable set
  std::set<unsigned int> _branch_set;

  //! @brief Variable bounds (temporary)
  T *_P_tmp;
  //! @brief Current node index
  unsigned int _node_index;
  //! @brief Total number of nodes introduced so far
  unsigned int _node_cnt;
  //! @brief Max. number of open nodes so far
  unsigned int _node_max;
  //! @brief Set of open nodes
  t_Nodes _open_nodes;
  //! @brief Volume of open nodes
  double _open_volume;
  //! @brief Set of interior nodes
  t_Nodes _inner_nodes;
  //! @brief Volume of interior nodes
  double _inner_volume;

  //! @brief Status after subproblem assessment
  STATUS _status;

  //! @brief Starting time
  double _tstart;
  //! @brief Current time
  double _tcur;

  //! @brief maximum number of values displayed in a row
  static const unsigned int _LDISP = 4;
  //! @brief reserved space for integer variable display
  static const unsigned int _IPREC = 7;
  //! @brief reserved space for double variable display
  static const unsigned int _DPREC = 6;
  //! @brief stringstream for displaying results
  std::ostringstream _odisp;

  //! @brief Erase stored nodes in stacks
  void _clean_stacks();
  //! @brief Reinitialize variables, stacks, counters, time, etc.
  void _restart();
  //! @brief Default set of branching variables
  void _branching_variable_set();

  //! @brief Branch node and create subdomains
  STATUS _branch_node
    ( SetInvNode<T>*pNode );
  //! @brief Selects the branching variable for given node
  std::pair<unsigned int, double> _select_branching_variable
    ( SetInvNode<T>*pNode, const std::set<unsigned int>&set_branch );
  //! @brief Partition branching variable domain for given node
  std::pair<const T, const T> _partition_variable_domain
    ( SetInvNode<T>*pNode, const unsigned int ip ) const;
  //! @brief Determine whether a termination criterion is met
  bool _terminate() const;

  //! @brief Add information for display
  void _display_add
    ( const double dval );
  void _display_add
    ( const unsigned int ival );
  void _display_add
    ( const std::string &sval );

  //! @brief Initialize display
  void _display_init();
  //! @brief Final display
  void _display_final();

  //! @brief Display current buffer stream and reset it
  void _display
    ( std::ostream&os );
  //! @brief Display current open nodes in tree
  void _display_open
    ( std::ostream&os ) const;
};

//! @brief C++ base class for set-inversion nodes
////////////////////////////////////////////////////////////////////////
//! mc::SetInvNode<T> is a C++ base class for defining nodes in the
//! set-inversion algorithm.
////////////////////////////////////////////////////////////////////////
template <typename T>
class SetInvNode
////////////////////////////////////////////////////////////////////////
{
  template <typename U> friend class SetInv;
  template <typename U> friend struct lt_SetInvNode;

public:  
  //! @brief Retreive (parent) node strength
  double strength() const
    { return _strength; }
  //! @brief Retreive node index
  unsigned int index() const
    { return _index; }
  //! @brief Retreive node depth
  unsigned int depth() const
    { return _depth; }
  //! @brief Retreive node iteration
  unsigned int iter() const
    { return _iter; }

  //! @brief Compute partition volume
  double volume() const
    { double V=1.;
      for( unsigned int i=0; i<_pSetInv->_np; i++ ) V *= Op<T>::diam(_P[i]);
      return V; }
  //! @brief Retreive current variable bounds
  const T& P
    ( const unsigned int ip ) const
    { assert( ip < _pSetInv->_np ); return _P[ip]; }
  //! @brief Retreive/set bound for variable <a>ip</a> 
  T& P
    ( const unsigned int ip )
    { assert( ip < _pSetInv->_np ); return _P[ip]; }
  //! @brief Retreive pointer to variable bounds
  const T* P() const
    { return _P; }
  //! @brief Retreive/set node type (ROOT/LEFT/RIGHT)
  typename SetInv<T>::NODETYPE& type()
    { return _type; }
  //! @brief Retreive/set parent branching variable and range
  std::pair<unsigned int,T>& parent()
    { return _parent; }

  //! @brief Public constructor (root node)
  SetInvNode
    ( SetInv<T>*pSetInv, const T*P, const unsigned int index=1,
      const unsigned int iter=0 );

  //! @brief Public constructor (child node)
  SetInvNode
    ( SetInv<T>*pSetInv, const T*P, const double strength,
      const unsigned int index, const unsigned int depth,
      const unsigned int iter, const typename SetInv<T>::NODETYPE type,
      const std::pair<unsigned int,T> parent );

  //! @brief Destructor
  ~SetInvNode();

  //! @brief Assess
  typename SetInv<T>::STATUS assess();

private:
  //! @brief Private default constructor
  SetInvNode<T>()
    {}

  //! @brief Pointer to underlying set-inversion problem
  SetInv<T> *_pSetInv;
  //! @breif Stength of parent node
  double _strength;
  //! @brief Depth in SetInv tree
  unsigned int _depth;
  //! @brief Index in SetInv tree
  unsigned int _index;
  //! @brief Iteration when created in SetInv tree
  unsigned int _iter;
  //! @brief Node type (ROOT/LEFT/RIGHT)
  typename SetInv<T>::NODETYPE _type;

  //! @brief Parent node branching variable and range
  std::pair<unsigned int,T> _parent;

  //! @brief Variable bounds
  T *_P;
  //! @brief Backup variable bounds (for domain reduction)
  T *_P0;
};

//! @brief C++ structure for comparing SetInv Nodes
////////////////////////////////////////////////////////////////////////
//! mc::lt_SetInvNode is a C++ structure for comparing nodes in branch-and-
//! bound tree based on their lower bound values.
////////////////////////////////////////////////////////////////////////
template <typename T>
struct lt_SetInvNode
{
  bool operator()
    ( const SetInvNode<T>*Node1, const SetInvNode<T>*Node2 ) const
    { return( Node1->strength() > Node2->strength() ); }
};

///////////////////////////////   SetInv   ////////////////////////////////

template <typename T>
inline void
SetInv<T>::variables
( const unsigned int np, const T*P, const std::set<unsigned int>*exclude )
{
  if( !np || !P ){
    _np = 0;
    delete[] _P_root; _P_root = 0;
    delete[] _P_tmp;  _P_tmp = 0;
    return;
  }
  
  if( _np != np ){
    delete[] _P_root;
    delete[] _P_tmp;
    _np = np;
    _P_root = new T[_np];
    _P_tmp = new T[_np];
  } 
  for( unsigned int i=0; i<_np; i++ ) _P_root[i] = P[i];

  if( exclude == &_exclude_vars ) return;
  _exclude_vars.clear();
  if( exclude ) _exclude_vars = *exclude;

  return;
}

template <typename T>
inline void
SetInv<T>::output_stacks
( std::ostream&os_open, std::ostream&os_inner, const unsigned int DPREC ) const
{
  // append all open sets from _open_nodes to os_open
  if( os_open.good() ){
    os_open << std::scientific << std::setprecision(DPREC); 
    cit_Nodes cit = _open_nodes.begin();
    for( ; cit!=_open_nodes.end(); ++cit ){
      for( unsigned int ip=0; ip<_np; ip++ )
        os_open << std::setw(DPREC+9) << Op<T>::l( (*cit)->P(ip) )
                << std::setw(DPREC+9) << Op<T>::u( (*cit)->P(ip) );
      os_open << std::endl;
    }
  }
  
  // append all open sets from _inner_nodes to os_inner
  if( os_inner.good() ){
    os_inner << std::scientific << std::setprecision(DPREC); 
    cit_Nodes cit = _inner_nodes.begin();
    for( ; cit!=_inner_nodes.end(); ++cit ){
      for( unsigned int ip=0; ip<_np; ip++ )
        os_inner << std::setw(DPREC+9) << Op<T>::l( (*cit)->P(ip) )
                 << std::setw(DPREC+9) << Op<T>::u( (*cit)->P(ip) );
      os_inner << std::endl;
    }
  }
}

template <typename T>
inline std::pair< double, double >
SetInv<T>::solve
( std::ostream&os )
{
  // Create and add root node to set _open_nodes
  _restart();
  _display_init();
  _display( os );
  SetInvNode<T>* rootNode = new SetInvNode<T>( this, _P_root );
  _open_nodes.insert( rootNode );
  _open_volume = _volume_root = rootNode->volume();
  
  // Run iterative set-inversion algorithm
  for( _node_index = 1; !_open_nodes.empty() && _tcur-_tstart < options.MAX_CPU_TIME
       && ( !options.MAX_NODES || _node_index <= options.MAX_NODES );
       _node_index++, _tcur=time() ){

    // Display
    _display_open( os );
    _display_add( _node_index );
    _display_add( (*_open_nodes.begin())->strength() );
    _display_add( (unsigned int)_open_nodes.size() );
    _display_add( time()-_tstart );
    _display_add( _inner_volume );
    _display_add( _open_volume );

    // Select node and remove it for the set _open_nodes
    SetInvNode<T>* pNode = *_open_nodes.begin();
    _open_nodes.erase( _open_nodes.begin() );
    _open_volume -= pNode->volume();

    // Display
    _display_add( pNode->iter() );

    // Assess node
    _status = pNode->assess();
    switch( _status ){
    case OUTER:
      delete pNode;
      _display_add( "OUTER" );
      _display( os );
      continue;

    case INNER:
      _inner_nodes.insert( pNode );
      _inner_volume += pNode->volume();
      _display_add( "INNER" );
      _display( os );
      continue;

    case UNDETERMINED:
      // No break here - continue to branching
      
    case FAILURE:
      // Domain branching
      if( _branch_node( pNode ) == ABORT )
        _display_add( "ABORT" );
      delete pNode;
      _display( os );
      break;

    case ABORT: default:
      delete pNode;
      _display_add( "ABORT" );
      _display( os );
      break;
    }

    // Keep track of the maximum nodes in memory
    if( _open_nodes.size() > _node_max ) _node_max = _open_nodes.size();

    // Check status and termination criteria
    if( _status == ABORT || _terminate() ) break;
  }

  _display_final();
  _display( os );
  
  return std::make_pair( _inner_volume, _open_volume );
}

template <typename T>
inline bool
SetInv<T>::_terminate() const
{
  if( _open_volume < options.ABSOLUTE_TOLERANCE 
   || _open_volume < options.RELATIVE_TOLERANCE * _volume_root )
    return true;
  return false;
}

template <typename T>
inline void
SetInv<T>::_restart()
{
  _clean_stacks();
  _branching_variable_set();
  
  _node_index = _node_cnt = _node_max = 0;
  _tstart = _tcur = time();  
}

template <typename T>
inline void
SetInv<T>::_clean_stacks()
{
  // Clean open nodes stack
  it_Nodes it = _open_nodes.begin();
  for( ; it != _open_nodes.end(); it++ ) delete *it;
  _open_nodes.clear();
  _open_volume = 0.;

  // Clean interior nodes stack
  it = _inner_nodes.begin();
  for( ; it != _inner_nodes.end(); it++ ) delete *it;
  _inner_nodes.clear();
  _inner_volume = 0.;
}

template <typename T>
inline void
SetInv<T>::_branching_variable_set()
{
  _branch_set.clear();
  for( unsigned int ip=0; ip<_np; ip++ ){
    if( _exclude_vars.empty() || _exclude_vars.find(ip) == _exclude_vars.end() )
      _branch_set.insert( ip );
  }
  if( _branch_set.empty() ) throw Exceptions( Exceptions::BRANCH );
}

template <typename T>
inline std::pair<unsigned int, double>
SetInv<T>::_select_branching_variable
( SetInvNode<T>*pNode, const std::set<unsigned int>&set_branch )
{
  std::pair<unsigned int, double> branchsel( _np, -1. ); // <- Can be any negative number
  
  switch( options.BRANCHING_VARIABLE_CRITERION ){
  // Branching based on relative range diameter
  case Options::RGREL:{
    std::set<unsigned int>::iterator it = set_branch.begin();
    for( ; it != set_branch.end(); ++it ){
      double score = Op<T>::diam( pNode->P(*it) )
                   / Op<T>::diam( _P_root[*it] );
      if( score > branchsel.second )
        branchsel = std::make_pair( *it, score );
    }
    break;
   }
  // Branching based on absolute range diameter
  case Options::RGABS:{
    std::set<unsigned int>::iterator it = set_branch.begin();
    for( ; it != set_branch.end(); ++it ){
      double score = Op<T>::diam( pNode->P(*it) );
      if( score > branchsel.second )
        branchsel = std::make_pair( *it, score );
    }
    break;
   }
  }

  return branchsel;
}

template <typename T>
inline std::pair<const T, const T>
SetInv<T>::_partition_variable_domain
( SetInvNode<T>*pNode, const unsigned int ip_branch )
const
{
  // Branch at mid-point of current variable range
  std::pair<T,T> partition;
  partition.first = mc::Op<T>::l(pNode->P(ip_branch)) + Op<T>::zeroone()
    *Op<T>::diam(pNode->P(ip_branch))/2.;
  partition.second = mc::Op<T>::u(pNode->P(ip_branch)) - Op<T>::zeroone()
    *Op<T>::diam(pNode->P(ip_branch))/2.;
  return partition;
}

template <typename T>
inline typename SetInv<T>::STATUS
SetInv<T>::_branch_node
( SetInvNode<T>*pNode )
{
  // Branching variable selection
  const unsigned int ip_branch
    = _select_branching_variable( pNode, _branch_set ).first;
  if( _status == ABORT ) return ABORT;

  // Partitionning
  std::pair<const T, const T> partition
    = _partition_variable_domain( pNode, ip_branch );
  for( unsigned int ip=0; ip<_np; ip++ ) _P_tmp[ip] = pNode->P(ip);

  // Left child node
  _P_tmp[ip_branch] = partition.first;
  SetInvNode<T>*pNodeL = new SetInvNode<T>( this, _P_tmp, pNode->volume(),
    ++_node_cnt, pNode->depth()+1, _node_index, LEFT,
    std::make_pair(ip_branch, pNode->P(ip_branch)) );
  _open_nodes.insert( pNodeL );
  _open_volume += pNodeL->volume();
  
  // Right chold node
  _P_tmp[ip_branch] = partition.second;
  SetInvNode<T>*pNodeR = new SetInvNode<T>( this, _P_tmp, pNode->volume(),
    ++_node_cnt, pNode->depth()+1, _node_index, RIGHT,
    std::make_pair(ip_branch, pNode->P(ip_branch)) );
  _open_nodes.insert( pNodeR );
  _open_volume += pNodeR->volume();

  // Display
  std::ostringstream omsg;
  omsg << "BRANCH" << ip_branch;
  _display_add( omsg.str() );

  return UNDETERMINED;
}

template <typename T>
inline void
SetInv<T>::_display_open
( std::ostream&os ) const
{
  if( options.DISPLAY <= 2 ) return;
  cit_Nodes cit = _open_nodes.begin();
  for( unsigned int id=1; cit != _open_nodes.end(); ++cit, id++ ){
    std::set<unsigned int>::iterator ivar = _branch_set.begin();
    os << "Node " << id << ":" << std::scientific << std::setprecision(_DPREC)
       << std::setw(_DPREC+8) << (*cit)->strength();
    for( ; ivar != _branch_set.end(); ++ivar )
      os << (*cit)->P(*ivar) << "  ";
    os << std::endl;
  }
  return;
}

template <typename T>
inline void
SetInv<T>::_display_init()
{
  _odisp.str("");
  if( options.DISPLAY <= 1 ) return;
  _odisp << std::right
  	 << std::setw(_IPREC) << "INDEX"
  	 << std::setw(_DPREC+8) << "VOLUME  "
  	 << std::setw(_IPREC) << "STACK"
  	 << std::setw(_DPREC+8) << "CUMUL TIME "
  	 << std::setw(_DPREC+8) << "INNER   "
  	 << std::setw(_DPREC+8) << "OPEN    "
  	 << std::setw(_IPREC) << "PARENT"
  	 << std::setw(_DPREC+8) << "ACTION"
  	 << std::endl;
  
}

template <typename T>
inline void
SetInv<T>::_display_add
( const double dval )
{
  if( options.DISPLAY <= 1 ) return;
  _odisp << std::right << std::scientific << std::setprecision(_DPREC)
         << std::setw(_DPREC+8) << dval;
}

template <typename T>
inline void
SetInv<T>::_display_add
( const unsigned int ival )
{
  if( options.DISPLAY <= 1 ) return;
  _odisp << std::right << std::setw(_IPREC) << ival;
}

template <typename T>
inline void
SetInv<T>::_display_add
( const std::string &sval )
{
  if( options.DISPLAY <= 1 ) return;
  _odisp << std::right << std::setw(_DPREC+8) << sval;
}

template <typename T>
inline void
SetInv<T>::_display_final()
{
  if( options.DISPLAY <= 0 ) return;

  // Solution found within allowed time?
  if( _tcur-_tstart < options.MAX_CPU_TIME
    && ( !options.MAX_NODES || _node_index <= options.MAX_NODES ) )
    _odisp << std::endl << "#  NORMAL TERMINATION:      ";
  else
    _odisp << std::endl << "#  EXECUTION STOPPED:       ";
  _odisp << std::fixed << std::setprecision(6) << _tcur-_tstart << " CPU SEC"
         << std::endl;

  // Set-inversion results
  _odisp << "#  INNER APPROXIMATION:     "
  	 << std::scientific << std::setprecision(_DPREC)
         << "VOLUME = " << std::setw(_DPREC+6) << _inner_volume
         << ",  NODES = "  << _inner_nodes.size()
  	 << std::endl;
  _odisp << "#  BOUNDARY APPROXIMATION:  "
  	 << std::scientific << std::setprecision(_DPREC)
         << "VOLUME = " << std::setw(_DPREC+6) << _open_volume
         << ",  NODES = "  << _open_nodes.size() 
  	 << std::endl;
  _odisp << "#  TOTAL NUMBER OF NODES:   " << _node_index-1 << std::endl
         << "#  MAX OPEN NODES IN STACK: " << _node_max << std::endl;
}

template <typename T>
inline void
SetInv<T>::_display
( std::ostream &os )
{
  if( _odisp.str() == "" ) return;
  //if( options.DISPLAY > 0 ){
  os << _odisp.str() << std::endl;
  //}
  _odisp.str("");
  return;
}

template <typename T>
inline void
SetInv<T>::Options::display
( std::ostream&out ) const
{
  // Display SetInv Options
  out << std::setw(60) << "  ABSOLUTE CONVERGENCE TOLERANCE"
      << std::scientific << std::setprecision(1)
      << ABSOLUTE_TOLERANCE << std::endl;
  out << std::setw(60) << "  RELATIVE CONVERGENCE TOLERANCE"
      << std::scientific << std::setprecision(1)
      << RELATIVE_TOLERANCE << std::endl;
  out << std::setw(60) << "  MAXIMUM NUMBER OF ITERATIONS (PARTITIONS)";
  switch( MAX_NODES ){
  case 0:  out << "N/A\n"; break;
  default: out << MAX_NODES << std::endl; break;
  }
  out << std::setw(60) << "  MAXIMUM CPU TIME (SEC)"
      << std::scientific << std::setprecision(1)
      << MAX_CPU_TIME << std::endl;
  out << std::setw(60) << "  DISPLAY LEVEL"
      << DISPLAY << std::endl;
  out << std::setw(60) << "  BRANCHING STRATEGY FOR VARIABLE SELECTION";
  switch( BRANCHING_VARIABLE_CRITERION ){
  case RGREL:  out << "RGREL\n";  break;
  case RGABS:  out << "RGABS\n";  break;
  }
}

/////////////////////////////// SetInvNode ///////////////////////////////

template <typename T>
inline
SetInvNode<T>::SetInvNode
( SetInv<T>*pSetInv, const T*P, const unsigned int index,
  const unsigned int iter ):
_pSetInv( pSetInv ), _depth( 0 ), _index( index ), _iter( iter ),
_type( SetInv<T>::ROOT ), _parent( 0, T(0.) )
{
  _P0 = new T[_pSetInv->_np];
  _P  = new T[_pSetInv->_np];
  for( unsigned int i=0; i<_pSetInv->_np; i++ ) _P[i] = P[i];
  _strength = volume();
}

template <typename T>
inline
SetInvNode<T>::SetInvNode
( SetInv<T>*pSetInv, const T*P, const double strength,
  const unsigned int index, const unsigned int depth,
  const unsigned int iter, const typename SetInv<T>::NODETYPE type,
  const std::pair<unsigned int,T> parent ):
_pSetInv( pSetInv ), _strength( strength ), _depth( depth ), _index( index ),
_iter( iter ), _type( type ), _parent( parent )
{
  _P0 = new T[_pSetInv->_np];
  _P  = new T[_pSetInv->_np];
  for( unsigned int i=0; i<_pSetInv->_np; i++ ) _P[i] = P[i];
}

template <typename T>
inline
SetInvNode<T>::~SetInvNode()
{
  delete[] _P0;
  delete[] _P;
}

template <typename T>
inline typename SetInv<T>::STATUS
SetInvNode<T>::assess()
{
  for( unsigned int i=0; i<_pSetInv->_np; i++ ) _P0[i] = _P[i];
  return _pSetInv->assess( _pSetInv->_np, _P );
}

} // namespace mc

#endif
