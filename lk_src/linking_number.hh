/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  LK : Linking Number                                                     |
 |                                                                          |
 |  file          : linking_number.hh                                       |
 |  authors       : Enrico Bertolazzi                                       |
 |  affiliations  :                                                         |
 |                                                                          |
 |      Department of Industrial Engineering                                |
 |      University of Trento                                                |
 |      Via Sommarive 9, I-38123, Povo, Trento                              |
 |      email : enrico.bertolazzi@unitn.it                                  |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#ifndef LINKING_NUMBER_HH
#define LINKING_NUMBER_HH

#include <vector>
#include <iostream>
#include <algorithm>
#include <cmath>

// if C++ < C++11 define nullptr
#if __cplusplus > 199711L
  #ifndef DO_NOT_USE_CXX11
    #define LINKING_NUMBER_USE_CXX11
  #endif
#endif

#if defined(_WIN32) || defined(WIN32) || defined(_WIN64) || defined(WIN64)
  // windows architecture
  #if _MSC_VER >= 1900
    #ifndef DO_NOT_USE_CXX11
      #define LINKING_NUMBER_USE_CXX11
    #endif
  #endif
#endif

//#undef LINKING_NUMBER_USE_CXX11

#ifdef LINKING_NUMBER_USE_CXX11
  #include <thread>
  #include <mutex>
  #include <condition_variable>
  #include <atomic>
#endif

namespace LK {

  using std::hypot ;

  template <typename T>
  class Constants {
  public:
    static T             const u_epsi;
    static T             const u_epsi_x;
    static T             const u_epsi_y;
    static unsigned long const I_max;
  };

  /*\
   |               _     _____ ____
   |   _ __  _ __ | |_  |___ /|  _ \
   |  | '_ \| '_ \| __|   |_ \| | | |
   |  | |_) | | | | |_   ___) | |_| |
   |  | .__/|_| |_|\__| |____/|____/
   |  |_|
  \*/

  template <typename T>
  inline
  void
  zero3( T a[3] )
  { a[0] = a[1] = a[2] = 0 ; }

  template <typename T>
  inline
  void
  copy3( T const a[3], T b[3] )
  { b[0] = a[0] ; b[1] = a[1] ; b[2] = a[2] ;  }

  template <typename T>
  inline
  T
  norm3( T const v[3] )
  { return hypot(hypot(v[0],v[1]),v[2]) ; }

  template <typename T>
  inline
  T
  dot3( T const a[3], T const b[3] )
  { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2] ; }

  template <typename T>
  inline
  T
  dist3( T const a[3], T const b[3] )
  { return hypot(hypot(a[0]-b[0],a[1]-b[1]),a[2]-b[2]) ; }

  template <typename T>
  inline
  void
  cross( T const a[3], T const b[3], T c[3] ) {
    c[0] = a[1]*b[2] - a[2]*b[1] ;
    c[1] = a[2]*b[0] - a[0]*b[2] ;
    c[2] = a[0]*b[1] - a[1]*b[0] ;
  }

  /*\
   |     _             _
   |    /_\  _ _  __ _| |___
   |   / _ \| ' \/ _` | / -_)
   |  /_/ \_\_||_\__, |_\___|
   |             |___/
  \*/

  template <typename T> class BigAngle ;
  template <typename T> class Angle ;

  template <typename T>
  class Angle : public Constants<T> {

    using Constants<T>::u_epsi ;
    using Constants<T>::u_epsi_x ;
    using Constants<T>::u_epsi_y ;
    using Constants<T>::I_max ;

  public:
    typedef T   real_type ;
    typedef int int_type ;

  private:

    real_type x, y ;  // vector angle
    real_type e ;     // floating point error
    int_type  s ;     // sign quadrant
    int_type  sigma ; // 2*pi turn numbers

  public:

    Angle() : x(0), y(0), e(0), s(0), sigma(0) {}

    Angle( Angle const & rhs ) { *this = rhs ; }

    Angle const &
    operator = ( Angle const & rhs ) {
      x     = rhs.x ;
      y     = rhs.y ;
      e     = rhs.e ;
      sigma = rhs.sigma ;
      s     = rhs.s ;
      return *this ;
    }

    void
    build( real_type const P1[3],
           real_type const P2[3],
           real_type const Q1[3],
           real_type const Q2[3] );

    friend class BigAngle<T> ;
  };

  /*\
   |   ___ _        _             _
   |  | _ |_)__ _  /_\  _ _  __ _| |___
   |  | _ \ / _` |/ _ \| ' \/ _` | / -_)
   |  |___/_\__, /_/ \_\_||_\__, |_\___|
   |        |___/           |___/
  \*/
  template <typename T>
  class BigAngle : public Constants<T> {

    using Constants<T>::u_epsi ;
    using Constants<T>::u_epsi_x ;
    using Constants<T>::u_epsi_y ;
    using Constants<T>::I_max ;

  public:
    typedef T   real_type ;
    typedef int int_type ;

  private:
    // angle part
    real_type     X, Y ;
    int_type      SIGMA, S ;
    // error propagation
    unsigned long I ;
    real_type     R ;

  public:

    BigAngle() : X(0), Y(0), SIGMA(0), S(0), I(0), R(0) {}
    BigAngle( BigAngle<T> const & rhs ) { *this = rhs ; }

    BigAngle<T> const & operator += ( Angle<T>    const & a ) ;
    BigAngle<T> const & operator -= ( Angle<T>    const & a ) ;
    BigAngle<T> const & operator += ( BigAngle<T> const & a ) ;
    BigAngle<T> const & operator -= ( BigAngle<T> const & a ) ;

    BigAngle<T> const &
    operator = ( BigAngle<T> const & a ) {
      X     = a.X ;     Y = a.Y ;
      SIGMA = a.SIGMA ; S = a.S ;
      I     = a.I ;     R = a.R ;
      return *this ;
    }

    void
    init() {
      X     = 1  ;
      Y     = 0  ;
      S     = -1 ;
      SIGMA = 0  ;
      I     = 0  ;
      R     = 0  ;
    }

    int_type  getSigma() const { return SIGMA ; }
    real_type getError() const { return (I*u_epsi+R*u_epsi) ; }
    real_type getAngle() const ;
    real_type getFraction() const ;
    void      checkSigma() const;

    friend
    std::ostream &
    operator << ( std::ostream & stream, BigAngle<T> const & a ) {
      stream
        << "X = " << a.X << " Y = " << a.Y
        << " SIGMA = " << a.SIGMA
        << " S = " << a.S
        << " I = " << a.I
        << " R = " << a.R ;
      return stream ;
    }

  } ;

  /*\
   |   _    _      _   _           _  _            _
   |  | |  (_)_ _ | |_(_)_ _  __ _| \| |_  _ _ __ | |__  ___ _ _
   |  | |__| | ' \| / / | ' \/ _` | .` | || | '  \| '_ \/ -_) '_|
   |  |____|_|_||_|_\_\_|_||_\__, |_|\_|\_,_|_|_|_|_.__/\___|_|
   |                         |___/
  \*/

  template <typename T>
  class LinkingNumber : public Constants<T> {

    using Constants<T>::u_epsi ;
    using Constants<T>::u_epsi_x ;
    using Constants<T>::u_epsi_y ;
    using Constants<T>::I_max ;

  public:
    typedef T   real_type ;
    typedef int int_type ;
    typedef struct { real_type x, y, z ; } pnt3 ;

  private:

    #ifdef LINKING_NUMBER_USE_CXX11
    unsigned                      numThread ;
    mutable std::atomic<unsigned> _lk_computed ;
    #else
    mutable unsigned              _lk_computed ;
    #endif

    typedef struct {
      real_type P1[3];
      real_type P2[3];
      int_type  weight;
    } Segment ;

    typedef std::vector<Segment> CURVE ;
    std::vector<CURVE>           curves ;

    void
    eval_row( real_type const   P1[3],
              real_type const   P2[3],
              CURVE     const & curve,
              BigAngle<T>     & out_angle ) const ;

    void
    eval_rows( int_type      istart,
               int_type      nstep,
               CURVE const & curveA,
               CURVE const & curveB,
               BigAngle<T> & out_angle ) const ;

    void
    writhe_row( int_type          i_skip,
                real_type const   P1[3],
                real_type const   P2[3],
                CURVE     const & curve,
                BigAngle<T>     & out_angle ) const ;
    void
    writhe_rows( int_type      istart,
                 int_type      nstep,
                 CURVE const & curve,
                 BigAngle<T> & out_angle ) const ;

    // disable copy constructor and operator
    LinkingNumber( ) ;
    LinkingNumber<T> const operator = ( LinkingNumber<T> const & ) ;
    LinkingNumber( LinkingNumber<T> & ) ;

  public:

    explicit LinkingNumber( unsigned tot_curve ) ;
    ~LinkingNumber() ;

    //! initialize class for `tot_curve`
    void
    setup( unsigned tot_curve ) ;

    //! initialize curve `ncurve` and reserve memory for `reserve_pnts` points
    void
    reset( unsigned ncurve, unsigned reserve_pnts = 0 ) ;

    //! init with first point `[x,y,z]` the curve `ncurve`
    void
    init_curve( unsigned  ncurve,
                real_type x,
                real_type y,
                real_type z ) ;

    //! add point `[x,y,z]` to curve `ncurve` with weight `weight`
    void
    add_point( unsigned  ncurve,
               real_type x,
               real_type y,
               real_type z,
               int_type  weight = 1 ) ;

    //! close `ncurve` with weight `weight` for last segment
    void
    close_curve( unsigned ncurve, int_type weight = 1 );

    void
    add_segment( unsigned  ncurve,
                 real_type x1,
                 real_type y1,
                 real_type z1,
                 real_type x2,
                 real_type y2,
                 real_type z2,
                 int_type  weight );

    template <typename Tvec>
    void
    add_curve( unsigned i_curve, std::vector<Tvec> const & a ) {
      reset( i_curve, unsigned(a.size()) ) ;
      if ( !a.empty() ) {
        typename std::vector<Tvec>::const_iterator ip = a.begin() ;
        init_curve( i_curve, (*ip)[0], (*ip)[1], (*ip)[2] ) ;
        for ( ++ip ; ip != a.end() ; ++ip )
          add_point( i_curve, (*ip)[0], (*ip)[1], (*ip)[2], 1 ) ;
        close_curve( i_curve ) ;
      }
    }

    void
    add_curve( unsigned i_curve, real_type const a[][3], unsigned a_size ) ;

    void
    evaluate( unsigned   i_curve,
              unsigned   j_curve,
              int_type & ret,
              T        & fraction ) const ;

    int_type
    eval( unsigned i_curve, unsigned j_curve ) const ;

    template <typename Tvec>
    int_type
    eval( std::vector<Tvec> const & a,
          std::vector<Tvec> const & b ) {
      if ( curves.size() < 2 ) setup( 2 ) ;
      add_curve( 0, a ) ;
      add_curve( 1, b ) ;
      return eval( 0, 1 ) ;
    }

    template <typename Type>
    int_type
    eval( Type const a[], unsigned size_a,
          Type const b[], unsigned size_b ) {
      if ( curves.size() < 2 ) setup( 2 ) ;
      add_curve( 0, a, size_a ) ;
      add_curve( 1, b, size_b ) ;
      return eval( 0, 1 ) ;
    }

    void
    evals( unsigned const i_curve[], unsigned ni,
           unsigned const j_curve[], unsigned nj,
           int_type mat[] ) ;

    #ifdef LINKING_NUMBER_USE_CXX11
    int_type
    eval_mt( unsigned i_curve, unsigned j_curve ) const ;
    #else
    int_type
    eval_mt( unsigned i_curve, unsigned j_curve ) const
    { return eval( i_curve, j_curve ) ; }
    #endif

    template <typename MAT>
    void
    evals( unsigned const i_curve[], unsigned ni,
           unsigned const j_curve[], unsigned nj,
           MAT & mat ) {
      if ( ni == 0 || nj == 0 ) return ; // caso triviale
      std::vector<int> row_mat(ni*nj) ;
      evals( i_curve, ni, j_curve, nj, &row_mat.front() ) ;
      for ( unsigned i = 0 ; i < ni ; ++i )
        for ( unsigned j = 0 ; j < nj ; ++j )
          mat(i,j) = row_mat[i+j*ni] ;
    }

    unsigned lk_computed() const { return _lk_computed ; }

    real_type
    writhe( unsigned i_curve, real_type & err ) const ;

    #ifdef LINKING_NUMBER_USE_CXX11
    real_type
    writhe_mt( unsigned i_curve, real_type & err ) const ;
    #else
    real_type
    writhe_mt( unsigned i_curve, real_type & err ) const
    { return writhe( i_curve, err ) ; }
    #endif

  } ;

}

#endif

