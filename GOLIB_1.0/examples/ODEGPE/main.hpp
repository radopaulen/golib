#ifndef MAIN_HPP
#define MAIN_HPP

#ifdef USE_PROFIL
  #include "mcprofil.hpp"
  typedef INTERVAL I;
  namespace fadbad
  {
    //! @brief C++ structure for specialization of the fadbad::Op structure to allow usage of the PROFIL interval arithmetic type INTERVAL inside the classes fadbad::T, fadbad::B and fadbad::T in FADBAD++
    template <> struct Op<INTERVAL>{
      typedef double Base;
      typedef INTERVAL I;
      static Base myInteger( const int i ) { return Base(i); }
      static Base myZero() { return myInteger(0); }
      static Base myOne() { return myInteger(1);}
      static Base myTwo() { return myInteger(2); }
      static double myPI() { return ::Constant::Pi; }
      static I myPos( const I& x ) { return  x; }
      static I myNeg( const I& x ) { return -x; }
      template <typename U> static I& myCadd( I& x, const U& y ) { return x+=y; }
      template <typename U> static I& myCsub( I& x, const U& y ) { return x-=y; }
      template <typename U> static I& myCmul( I& x, const U& y ) { return x*=y; }
      template <typename U> static I& myCdiv( I& x, const U& y ) { return x/=y; }
      static I myInv( const I& x ) { return myOne()/x; }
      static I mySqr( const I& x ) { return Sqr(x); }
      template <typename X, typename Y> static I myPow( const X& x, const Y& y ) { return ::Power( x, y ); }
      static I mySqrt( const I& x ) { return ::Sqrt( x ); }
      static I myLog( const I& x ) { return ::Log( x ); }
      static I myExp( const I& x ) { return ::Exp( x ); }
      static I mySin( const I& x ) { return ::Sin( x ); }
      static I myCos( const I& x ) { return ::Cos( x ); }
      static I myTan( const I& x ) { return ::Tan( x ); }
      static I myAsin( const I& x ) { return ::ArcSin( x ); }
      static I myAcos( const I& x ) { return ::ArcCos( x ); }
      static I myAtan( const I& x ) { return ::ArcTan( x ); }
      static bool myEq( const I& x, const I& y ) { return x==y; }
      static bool myNe( const I& x, const I& y ) { return x!=y; }
      static bool myLt( const I& x, const I& y ) { return x<y; }
      static bool myLe( const I& x, const I& y ) { return x<=y; }
      static bool myGt( const I& x, const I& y ) { return y<x; }
      static bool myGe( const I& x, const I& y ) { return y<=x; }
    };
  } // namespace fadbad
#else
  #ifdef USE_FILIB
    #include "mcfilib.hpp"
    typedef filib::interval<double,filib::native_switched,filib::i_mode_extended> I;
    //typedef filib::interval<double> I;
  namespace fadbad
  {
    //! @brief C++ structure for specialization of the fadbad::Op structure to allow usage of the PROFIL interval arithmetic type INTERVAL inside the classes fadbad::T, fadbad::B and fadbad::T in FADBAD++
    template <> struct Op<filib::interval<double,filib::native_switched,filib::i_mode_extended> >{
    //template <> struct Op<filib::interval<double> >{
      typedef double Base;
      typedef filib::interval<double,filib::native_switched,filib::i_mode_extended> I;
      //typedef filib::interval<double> I;
      static Base myInteger( const int i ) { return Base(i); }
      static Base myZero() { return myInteger(0); }
      static Base myOne() { return myInteger(1);}
      static Base myTwo() { return myInteger(2); }
      static double myPI() { return mc::PI; }
      static I myPos( const I& x ) { return  x; }
      static I myNeg( const I& x ) { return -x; }
      template <typename U> static I& myCadd( I& x, const U& y ) { return x+=y; }
      template <typename U> static I& myCsub( I& x, const U& y ) { return x-=y; }
      template <typename U> static I& myCmul( I& x, const U& y ) { return x*=y; }
      template <typename U> static I& myCdiv( I& x, const U& y ) { return x/=y; }
      static I myInv( const I& x ) { return myOne()/x; }
      static I mySqr( const I& x ) { return filib::sqr(x); }
      template <typename X> static I myPow( const X& x, const int n ) { return filib::power( x, n ); }
      template <typename X, typename Y> static I myPow( const X& x, const Y& y ) { return filib::pow( x, y ); }
      static I mySqrt( const I& x ) { return filib::sqrt( x ); }
      static I myLog( const I& x ) { return filib::log( x ); }
      static I myExp( const I& x ) { return filib::exp( x ); }
      static I mySin( const I& x ) { return filib::sin( x ); }
      static I myCos( const I& x ) { return filib::cos( x ); }
      static I myTan( const I& x ) { return filib::tan( x ); }
      static I myAsin( const I& x ) { return filib::asin( x ); }
      static I myAcos( const I& x ) { return filib::acos( x ); }
      static I myAtan( const I& x ) { return filib::atan( x ); }
      static bool myEq( const I& x, const I& y ) { return x.seq(y); }
      static bool myNe( const I& x, const I& y ) { return x.sne(y); }
      static bool myLt( const I& x, const I& y ) { return x.slt(y); }
      static bool myLe( const I& x, const I& y ) { return x.sle(y); }
      static bool myGt( const I& x, const I& y ) { return x.sgt(y); }
      static bool myGe( const I& x, const I& y ) { return x.sge(y); }
    };
  } // namespace fadbad
  #else
    #include "interval.hpp"
    typedef mc::Interval I;
  #endif
#endif
typedef mc::TModel<I> TM;

#endif
