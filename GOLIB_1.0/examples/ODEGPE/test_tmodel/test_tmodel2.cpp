#include <fstream>
#include <iomanip>
#include "odegpe.hpp"
#include "interval.hpp"
typedef mc::Interval T;
typedef mc::TModel<T> TMT;
typedef mc::TVar<T> TVT;

const unsigned int NP = 2;
const unsigned int NF = 2;
const unsigned int NTM = 3;
const unsigned int MX = 20;
const unsigned int MY = 20;

#define SAVE_RESULTS

int main()
{
  TMT* _pTM = new TMT(NP,NTM);
  _pTM->options.BOUNDER_TYPE    = TMT::Options::LSB;
  _pTM->options.SCALE_VARIABLES = true;

#ifdef SAVE_RESULTS
  std::ofstream oorig( "TM-orig.out", std::ios_base::out );
  oorig << std::scientific << std::setprecision(5) << std::right;
  std::ofstream oupdt( "TM-updt.out", std::ios_base::out );
  oupdt << std::scientific << std::setprecision(5) << std::right;
#endif

  // Original Taylor model
  T Ip[NP] = { T(-1.,1.), T(-1.,1.) };
  TVT TVp[NP];
  for( unsigned int ip=0; ip<NP; ip++ )
    TVp[ip] = TVT( _pTM, ip, Ip[ip], 0. );

  TVT** TVf = new TVT*[1];
  TVf[0] = new TVT[NF];
  TVf[0][0] = exp(TVp[0])*TVp[1];
  TVf[0][1] = T(-1.,1.);
  mc::TMNode<T>* TVf_model = new mc::TMNode<T>( TVf, _pTM, 0 );
  std::cout << "Model: TVf[0][0]" << " =" << TVf_model->TVxk[0][0];
  std::cout << "Model: TVf[0][1]" << " =" << TVf_model->TVxk[0][1];

#ifdef SAVE_RESULTS
  for( unsigned int iX=0; iX<MX; iX++ ){ 
    for( unsigned int iY=0; iY<MY; iY++ ){ 
      double DXY[2] = { mc::Op<T>::l(Ip[0])+iX*mc::Op<T>::diam(Ip[0])/(MX-1.),
                        mc::Op<T>::l(Ip[1])+iY*mc::Op<T>::diam(Ip[1])/(MY-1.) };
      double DF = exp(DXY[0])*DXY[1];
      T Tf = TVf[0][0].P( DXY ) + TVf[0][0].R();
      oorig << std::setw(14) << DXY[0] << std::setw(14) << DXY[1] << std::setw(14) << DF
            << std::setw(14) << mc::Op<T>::l(Tf) << std::setw(14) << mc::Op<T>::u(Tf)
            << std::setw(14) << mc::Op<T>::l(TVf[0][0].B()) << std::setw(14) << mc::Op<T>::u(TVf[0][0].B())
            << std::endl;
    }
    oorig << std::endl;
  }
#endif

  // Reduced range Taylor model
  TVf = new TVT*[1];
  TVf[0] = new TVT[NF];
  Ip[0] = T(-.9,-.1); Ip[1] = T(-.8,0.);
  for( unsigned int ip=0; ip<NP; ip++ )
    TVp[ip] = TVT( _pTM, ip, Ip[ip], -0.5 );
    //TVp[ip] = TVT( _pTM, ip, Ip[ip], TVf_model->refpoint[ip]*TVf_model->scaling[ip] );

  // Update reference point and scaling in stored Taylor model
  double*coefupd = new double[_pTM->nmon()];
  for( unsigned int ix=0; ix<NF; ix++ ){
    if( !TVf_model->TVxk[0][ix].env() ){
      TVf[0][ix] = TVf_model->TVxk[0][ix];
      continue;
    }
    TVf_model->shift( 0, ix, _pTM->reference(), _pTM->scaling(), coefupd );
    TVf[0][ix].set( _pTM );
    TVf[0][ix].set( coefupd );
    TVf[0][ix].set( TVf_model->TVxk[0][ix].R() );
  }
  std::cout << "Updated: TVf[0][0]" << " =" << TVf[0][0];
  std::cout << "Updated: TVf[0][1]" << " =" << TVf[0][1];

#ifdef SAVE_RESULTS
  for( unsigned int iX=0; iX<MX; iX++ ){ 
    for( unsigned int iY=0; iY<MY; iY++ ){ 
      double DXY[2] = { mc::Op<T>::l(Ip[0])+iX*mc::Op<T>::diam(Ip[0])/(MX-1.),
                        mc::Op<T>::l(Ip[1])+iY*mc::Op<T>::diam(Ip[1])/(MY-1.) };
      double DF = exp(DXY[0])*DXY[1];
      T Tf = TVf[0][0].P( DXY ) + TVf[0][0].R();
      oupdt << std::setw(14) << DXY[0] << std::setw(14) << DXY[1] << std::setw(14) << DF
            << std::setw(14) << mc::Op<T>::l(Tf) << std::setw(14) << mc::Op<T>::u(Tf)
            << std::setw(14) << mc::Op<T>::l(TVf[0][0].B()) << std::setw(14) << mc::Op<T>::u(TVf[0][0].B())
            << std::endl;
    }
    oupdt << std::endl;
  }
#endif

  delete[] coefupd;
  delete TVf_model;
  delete[] TVf[0];
  delete[] TVf;
  delete _pTM;

}
