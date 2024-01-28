#ifndef MC__EA_HPP
#define MC__EAC_HPP

#include "ffunc.hpp"
//#include <cpplapack_dsysv.hpp>

namespace mc
{
//! @brief Class defining variables in a factorable function
////////////////////////////////////////////////////////////////////////
//! mc::Ellipsoidal_Image is a C++ class that defines and computes an
//! ellipsoidal enclosure of the image set of a factorable function.
////////////////////////////////////////////////////////////////////////
class Ellipsoidal_Image
////////////////////////////////////////////////////////////////////////
{

    private:

    public:
        //! @brief Constructor
        Ellipsoidal_Image( const FFTree*, const int, const FFVar* );


};

////////////////////////////////////////////////////////////////////////
///////      Function declaration for Ellipsoidal_Image //////////
////////////////////////////////////////////////////////////////////////


Ellipsoidal_Image::Ellipsoidal_Image
( const FFTree* _Tree_Obj, const int _nF, const FFVar* _F )
{
    long _nAux = _Tree_Obj->naux();
    std::cout << " number of auxiliary variables " << _nAux << std::endl;

}

} // end of namespace mc

#endif

