#ifndef MC__CPPLAPACK_DSYSV_HPP
#define MC__CPPLAPACK_DSYSV_HPP

#include "cpplapack.h"

namespace CPPL{ //!< namespace for CPPLapack

//=============================================================================
/*! calculate inverse.\n
  All of the arguments need not to be initialized.
  mat is not overwritten. 
*/
inline long dsysv(const dsymatrix& mat, dsymatrix& mat_inv)
{VERBOSE_REPORT;
  dsymatrix mat_cp(mat);
  mat_inv.resize(mat.n);
  mat_inv.identity(); 
  char UPLO('l');
  long NRHS(mat.n), LDA(mat.n), *IPIV(new long[mat.n]), LDB(mat.n), LWORK(-1), INFO(1);
  double *WORK( new double[1] );
  dsysv_(UPLO, mat.n, NRHS, mat_cp.array, LDA, IPIV, mat_inv.array, LDB, WORK, LWORK, INFO);

  LWORK = long(WORK[0]);
  delete [] WORK;  WORK = new double[LWORK];
  dsysv_(UPLO, mat.n, NRHS, mat_cp.array, LDA, IPIV, mat_inv.array, LDB, WORK, LWORK, INFO);
  delete [] WORK; delete [] IPIV;

  if(INFO!=0){
    WARNING_REPORT;
    std::cerr << "Serious trouble happend. INFO = " << INFO << "." << std::endl;
  } 
  return INFO;
}

}//namespace CPPL

#endif
