
/*!
\page page_EA Ellipsoidal arithmetic for factorable functions
\author Mario E. Villanueva & Benoit Chachuat
\version 0.0
\date 2013
\bug No known bugs.

TBC.

*/


#include "ffunc.hpp"
#include "cpplapack.h"



namespace mc
{
//! @brief Structure defining variables in a factorable function
////////////////////////////////////////////////////////////////////////
//! mc::Ellipsoidal_Image is a C++ class that defines and computes an  
//! ellipsoidal enclosure of the image set of a factorable function.
////////////////////////////////////////////////////////////////////////
struct VFFactor
////////////////////////////////////////////////////////////////////////
{
  		 long   result_i;
  		 long loperand_i;
         long roperand_i;
  const FFOp* Oreference;		
};


//! @brief Class defining variables in a factorable function
////////////////////////////////////////////////////////////////////////
//! mc::Ellipsoidal_Image is a C++ class that defines and computes an  
//! ellipsoidal enclosure of the image set of a factorable function.
////////////////////////////////////////////////////////////////////////
class Ellipsoidal_Image
////////////////////////////////////////////////////////////////////////
{
 
	//! @Overload operator <<
	friend std::ostream& operator<< (std::ostream&, const Ellipsoidal_Image&);
		
	private:
		//! @brief Shape Matrix of the lifted Ellipsoid
		CPPL::dssmatrix _Q;
		
		//! @brief centre of the lifted Ellipsoid
		CPPL::dcovector _q;
	
		//! @brief number of variables
		long _nx;
		
		//! @brief number of factors not counting original variables
		long _nv;
		
		//! @brief number of components of the vector function
		long _nF;
		
		//! @brief number of factors
		long _nQ;
		
		//! @brief map between factorable tree and matrix index
	    std::map<long, long> _zv_index;
		
		//! @brief ... 
	    std::set<long> _Constants_Oi, _Leaves_Oi, _Auxv_Oi, _Varo_Oi;  
		
		//! @brief Projection matrix
		CPPL::dgsmatrix _P; 
		
		//! @brief Decomposition used in the matrix update
		std::list<mc::VFFactor> _vdecomp;
	
    public: 
		//! @brief Constructor
		Ellipsoidal_Image(const std::list<FFOp*>& );
		
   
    	//! @brief Computes the ellipsoidal enclosure for a given initial ellipsoidal set
		void E_enclosure(CPPL::dsymatrix& , CPPL::dcovector&);

 	//! @brief Function to retrieve the factorable decomposition used in the matrix update
		void FF2VFF_Decomposition(long ,  const std::list<mc::FFOp*>& , std::list<mc::VFFactor>&);

};

////////////////////////////////////////////////////////////////////////
///////      Function declaration for Ellipsoidal_Image       //////////
////////////////////////////////////////////////////////////////////////



///////      			Constructor					          //////////
Ellipsoidal_Image::Ellipsoidal_Image ( const std::list<FFOp*>& _Subtree_Obj )
{
	long _nAux     =  0;
//	long _nConst   =  0;
	     _nx       =  0;

	

	for (std::list<FFOp*>::const_iterator i = _Subtree_Obj.begin(); i != _Subtree_Obj.end(); ++i)
	{
		switch((*i) -> pres ->id().first)
		{   
			case    0: _Varo_Oi.insert((*i) -> pres ->id().second); break;
			case    1: _Auxv_Oi.insert((*i) -> pres ->id().second); break;
			case    2:
			case    3: _Constants_Oi.insert((*i) -> pres ->id().second); break;
			default  : break;
		}
		 
	}

	_Leaves_Oi = _Auxv_Oi; 
	_nAux      = _Auxv_Oi.size();
	_nx        = _Varo_Oi.size();
    
	FF2VFF_Decomposition(_nx,  _Subtree_Obj, _vdecomp);

	_nF   = _Leaves_Oi.size();
    _nv   =  _nAux;
	_nQ   = _nx + _nv;
    
     std::cout<<"Leaves\n";
     for (std::set<long>::iterator j = _Leaves_Oi.begin(); j != _Leaves_Oi.end(); ++j)std::cout << *j << std::endl;
    
    _Q.resize(_nQ); _Q.zero();
    _q.resize(_nQ); _q.zero();
    _P.resize(_nF,_nQ); _P.zero();
    
    
} 


///////      	V-Factorable rearrangement			          //////////
void Ellipsoidal_Image::FF2VFF_Decomposition(long nvar,  const std::list<mc::FFOp*>& _st, std::list<mc::VFFactor>& _vdc)
{

	long v_index = nvar;
	std::set<long>::iterator _it_Index;
	mc::VFFactor vfac;
	int Is_Univariate;

	for (std::list<mc::FFOp*>::const_iterator i = _st.begin(); i != _st.end(); ++i)
	{
	   switch ((*i)-> type )
	   {
	   case 3:
	   case 8:
	   case 9:
	   case 10:
	   case 11:
	   case 12:
	   case 13:
	   case 14:
	   case 15:
	   case 16:
	   case 17:
	   case 18:
	   case 19:
	   case 20:
	   case 21: Is_Univariate = 1; break;
	   default: Is_Univariate = 0; break;
	   }
	   
	   if ( (*i)-> pres->id().first == 1 )
		{ 
			_zv_index.insert(std::make_pair((*i)-> pres->id().second,v_index));
			vfac.result_i           = v_index; 
			vfac.Oreference = *i ;
			switch ( (*i)-> plop->id().first )
			{
				case   0: 	vfac.loperand_i = (*i)-> plop->id().second; break;
				case   1: 	vfac.loperand_i = _zv_index.find((*i)-> plop->id().second)->second;
							_it_Index = _Leaves_Oi.find((*i)-> plop->id().second);
							if (_it_Index != _Leaves_Oi.end() ) _Leaves_Oi.erase(_it_Index);
							break;
				default : 	break;
			}
		if ( Is_Univariate == 0 )
		{	
			switch ( (*i)-> prop->id().first )
			{
				case   0: 	vfac.roperand_i = (*i)-> prop->id().second;	break;
				case   1: 	vfac.roperand_i = _zv_index.find((*i)-> prop->id().second)->second; 
							_it_Index = _Leaves_Oi.find((*i)-> prop->id().second);
							if (_it_Index != _Leaves_Oi.end() ) _Leaves_Oi.erase(_it_Index);
							break;
				default : 	break;
		 	}		
		}
		_vdc.push_back(vfac);    	
		v_index = v_index + 1 ;
		}
	}	
}

///////      		Operator overload				          //////////
std::ostream& operator<< (std::ostream&os, const Ellipsoidal_Image&_EI)
{
os << "Ellipsoidal Image: \n\n";
os << "Number of Factors      :" << _EI._nv << "\n";
os << "Number of Variables    :" << _EI._nx << "\n";
os << "Number of Functions    :" << _EI._nF << "\n";
os << "Number of Rows/Columns :" << _EI._nQ << "\n";

os << "\n" ;

os << "Factorable representation (v-form):" << "\n" ;
for (long i = 0; i < _EI._nx; i++) os << "v_" << i << " = " << "x_" << i << std::endl;
     
for (std::list<mc::VFFactor>::const_iterator i = _EI._vdecomp.begin(); i != _EI._vdecomp.end(); ++i)
{
    os << "v_" << (*i).result_i << " = " ;
 	switch( (*i).Oreference -> type )
 	{
		case 2:   os << "v_" << (*i).loperand_i << " + " << "v_" << (*i).roperand_i; break; //PLUS
		case 3:   os << "- " << "v_"<<(*i).loperand_i ; break; 								 //NEG
		case 4:   os << "v_" << (*i).loperand_i << " - " <<"v_" << (*i).roperand_i ; break; //MINUS
		case 5:   os << "v_" << (*i).loperand_i << " * " << "v_" << (*i).roperand_i; break; //TIMES
		case 6:   

	          if((*i).Oreference -> prop -> id().first == 2)
	          {
		          os <<  (*i).Oreference -> prop -> num().n  << " * " << "v_" << (*i).loperand_i;
	          }
	          else if((*i).Oreference -> prop -> id().first == 2)
	          {
	          	  os <<  (*i).Oreference -> prop -> num().x << " * " << "v_" << (*i).loperand_i;	
	          }
	          break;																				//SCALE
		case 7:   os << "v_" << (*i).loperand_i << " / " << "v_" << (*i).roperand_i; break;//DIV 
		case 8:   os << "EXP( v_" << (*i).loperand_i << " )"; break;					    //EXP
		case 9:   os << "LOG( v_" << (*i).loperand_i << " )"; break;
		case 10:  os << "SQRT( v_" << (*i).loperand_i << " )"; break;
		case 11:  os << "SQR( v_" << (*i).loperand_i << " )"; break;
		case 12:  os << "POW( v_" << (*i).loperand_i << " , " << (*i).Oreference -> prop -> num().n << " )" ; break;
		case 13:  os << "SIN( v_" << (*i).loperand_i << " )"; break;
		case 14:  os << "COS( v_" << (*i).loperand_i << " )"; break;
		case 15:  os << "TAN( v_" << (*i).loperand_i << " )"; break;
		case 16:  os << "ASIN( v_" << (*i).loperand_i << " )"; break;
		case 17:  os << "ACOS( v_" << (*i).loperand_i << " )"; break;
		case 18:  os << "ATAN( v_" << (*i).loperand_i << " )"; break;
		case 19:  os << "FABS( v_" << (*i).loperand_i << " )"; break;
		case 20:  os << "ERF( v_" << (*i).loperand_i << " )"; break;
		case 21:  os << "FSTEP( v_" << (*i).loperand_i << " )"; break;
		case 22:  os << "MINF( v_" << (*i).loperand_i << " , " << "v_" << (*i).roperand_i << " )"; break;
		case 23:  os << "MAXF( v_" << (*i).loperand_i << " , " << "v_" << (*i).roperand_i << " )"; break;
		default:;		
  	}
	os << "\n";	
}

os << "\n";
os << "Ellipsoidal Image centre       :\n" << _EI._q << "\n";
os << "Ellipsoidal Image shape matrix :\n" << _EI._Q << "\n";
os << "Projection matrix :\n" << _EI._P << "\n";
return os;
}



} // end of namespace mc







