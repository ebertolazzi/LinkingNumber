/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  LK : Linking Number                                                     |
 |                                                                          |
 |  file          : linking_number.cc                                       |
 |  authors       : Enrico Bertolazzi                                       |
 |  affiliations  :                                                         |
 |                                                                          |
 |      Department of Industrial Engineering                                |
 |      University of Trento                                                |
 |      Via Sommarive 9, I-38123, Povo, Trento                              |
 |      email : enrico.bertolazzi@unitn.it                                  |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#include "linking_number.hh"
#include <cmath>
#include <iostream>
#include <sstream>
#include <limits>

#ifdef __GCC__
#pragma GCC optimize on
#pragma GCC optimization_level 3
#endif
#ifdef __clang__
#pragma clang optimize on
#endif

#ifndef LK_ASSERT
  #include <stdexcept>
  #include <sstream>
  #define LK_ASSERT(COND,MSG)               \
    if ( !(COND) ) {                        \
      std::ostringstream ost ;              \
      ost << MSG << '\n' ;                  \
      throw std::runtime_error(ost.str()) ; \
    }
#endif

using namespace std ;

namespace LK {

  typedef long double ldouble ;

  // ------------------------------------------------------------------------

  /*\
   |              _
   |  ___ __ __ _| |___
   | (_-</ _/ _` | / -_)
   | /__/\__\__,_|_\___|
   |
  \*/

  template <typename T>
  void
  scale( T & x, T & y ) ;

  template <>
  inline
  void
  scale( float & x, float & y ) {
    int ex, ey ;
    float mx  = frexpf( x, &ex ) ;
    float my  = frexpf( y, &ey ) ;
    int   mxy = std::max( ex, ey ) ;
    x = ldexpf( mx, ex-mxy ) ;
    y = ldexpf( my, ey-mxy ) ;
  }

  template <>
  inline
  void
  scale( double & x, double & y ) {
    int ex, ey ;
    double mx  = frexp( x, &ex ) ;
    double my  = frexp( y, &ey ) ;
    int    mxy = std::max( ex, ey ) ;
    x = ldexp( mx, ex-mxy ) ;
    y = ldexp( my, ey-mxy ) ;
  }

  template <>
  inline
  void
  scale( ldouble & x, ldouble & y ) {
    int ex, ey ;
    ldouble mx  = frexpl( x, &ex ) ;
    ldouble my  = frexpl( y, &ey ) ;
    int         mxy = std::max( ex, ey ) ;
    x = ldexpl( mx, ex-mxy ) ;
    y = ldexpl( my, ey-mxy ) ;
  }

  // ------------------------------------------------------------------------

  /*\
   |                  _            _
   |   __ ___ _ _  __| |_ __ _ _ _| |_ ___
   |  / _/ _ \ ' \(_-<  _/ _` | ' \  _(_-<
   |  \__\___/_||_/__/\__\__,_|_||_\__/__/
   |
  \*/

  template <> float const
  Constants<float>::u_epsi = std::numeric_limits<float>::epsilon();

  template <> double const
  Constants<double>::u_epsi = std::numeric_limits<double>::epsilon();

  template <> ldouble const
  Constants<ldouble>::u_epsi = std::numeric_limits<ldouble>::epsilon();

  template <> float const
  Constants<float>::u_epsi_x = float(57.515)*std::numeric_limits<float>::epsilon();

  template <> double const
  Constants<double>::u_epsi_x = double(57.515)*std::numeric_limits<double>::epsilon();

  template <> ldouble const
  Constants<ldouble>::u_epsi_x = ldouble(57.515)*std::numeric_limits<ldouble>::epsilon();

  template <> float const
  Constants<float>::u_epsi_y = float(20.257)*std::numeric_limits<float>::epsilon();

  template <> double const
  Constants<double>::u_epsi_y = double(20.257)*std::numeric_limits<double>::epsilon();

  template <> ldouble const
  Constants<ldouble>::u_epsi_y = ldouble(20.257)*std::numeric_limits<ldouble>::epsilon();

  template <> unsigned long const
  Constants<float>::I_max = (unsigned long)(1.570796327/std::numeric_limits<float>::epsilon());

  template <> unsigned long const
  Constants<double>::I_max = (unsigned long)(1.570796327/std::numeric_limits<double>::epsilon());

  template <> unsigned long const
  Constants<ldouble>::I_max = (unsigned long)(1.570796327/std::numeric_limits<ldouble>::epsilon());

  // ------------------------------------------------------------------------

  /*\
   |     _             _
   |    /_\  _ _  __ _| |___
   |   / _ \| ' \/ _` | / -_)
   |  /_/ \_\_||_\__, |_\___|
   |             |___/
  \*/

  template <typename T>
  void
  Angle<T>::build( real_type const P1[3],
                   real_type const P2[3],
                   real_type const Q1[3],
                   real_type const Q2[3]) {

    real_type alpha[3] = { Q1[0]-P1[0], Q1[1]-P1[1], Q1[2]-P1[2] } ;
    real_type beta[3]  = { Q1[0]-P2[0], Q1[1]-P2[1], Q1[2]-P2[2] } ;
    real_type gamma[3] = { Q2[0]-P2[0], Q2[1]-P2[1], Q2[2]-P2[2] } ;
    real_type omega[3] = { Q2[0]-P1[0], Q2[1]-P1[1], Q2[2]-P1[2] } ;

    real_type len_alpha = norm3(alpha);
    real_type len_beta  = norm3(beta);
    real_type len_gamma = norm3(gamma);
    real_type len_omega = norm3(omega);

    if ( len_alpha > 0 ) {
      alpha[0] /= len_alpha ;
      alpha[1] /= len_alpha ;
      alpha[2] /= len_alpha ;
    }
    if ( len_beta > 0 ) {
      beta[0] /= len_beta ;
      beta[1] /= len_beta ;
      beta[2] /= len_beta ;
    }
    if ( len_gamma > 0 ) {
      gamma[0] /= len_gamma ;
      gamma[1] /= len_gamma ;
      gamma[2] /= len_gamma ;
    }
    if ( len_omega > 0 ) {
      omega[0] /= len_omega ;
      omega[1] /= len_omega ;
      omega[2] /= len_omega ;
    }

    real_type t1[3] = { alpha[0] + gamma[0], alpha[1] + gamma[1], alpha[2] + gamma[2] } ;
    real_type t2[3] ;
    cross( alpha, gamma, t2 );

    real_type tmp = 1+dot3(alpha,gamma);

    real_type x1 = tmp + dot3(beta,t1);
    real_type y1 = dot3(beta,t2);
    real_type x2 = tmp + dot3(omega,t1);
    real_type y2 = dot3(omega,t2); // WHEN omega ~ 0 y2 may have wrong sign!

    LK_ASSERT( isfinite(x1) && isfinite(y1) && isfinite(x2) && isfinite(y2),
               "found unormalized number(s):" <<
               "\nx1 = " << x1 <<
               "\ny1 = " << y1 <<
               "\nx2 = " << x2 <<
               "\ny2 = " << y2 ) ;
  
    if ( (y1 >= 0 && y2 >= 0) || (y1 <= 0 && y2 <= 0) ) {
      LK_ASSERT( x1 > u_epsi_x && x2 > u_epsi_x,
                 "Angle::build, ambiguos angle\nx1 = " << x1 <<
                 "\nx2 = " << x2 << "\ny1=y2=0 (probable intersection)" ) ;
      x     = 1 ; // x1*x2 ;
      y     = 0 ;
      s     = ( y > 0 || (y == 0 && x < 0)) ? 1 : -1 ; // s(x,y)
      sigma = 0 ;
    } else {
      // check if (x1,y1) and (x2,y2) are acceptable
      LK_ASSERT( abs(y1) > u_epsi_y || x1 > u_epsi_x,
                 "Angle::build, ambiguos angle\nx1 = " << x1 << "\ny1 = " << y1 <<
                 "\n(probable intersection)" ) ;
      LK_ASSERT( abs(y2) > u_epsi_y || x2 > u_epsi_x,
                 "Angle::build, ambiguos angle\nx2 = " << x2 << "\ny2 = " << y2 <<
                 "\n(probable intersection)" ) ;

      x = x1*x2+y1*y2 ;
      y = x1*y2-y1*x2 ;

      LK_ASSERT( isfinite(x) && isfinite(y),
                 "Angle::build, found unormalized number(s): x = " << x << "\ny = " << y );

      s     = ( y > 0 || (y == 0 && x < 0)) ? 1 : -1 ; // s(x,y)
      sigma = s*y1 > 0 ? -s : 0 ;
    }

    real_type R1 = hypot(x1,y1) ;
    real_type R2 = hypot(x2,y2) ;

    // check if point angle is far from the origin
    LK_ASSERT( 181 <= I_max*R1 && 181 <= I_max*R2,
               "Angle::build, segments too close!\nR1 = " << R1 << " R2 = " << R2 ) ;

    e = real_type(2.829 + (57.516/R1+57.516/R2)) ;

  }

  // ------------------------------------------------------------------------

  /*\
   |   ___ _        _             _
   |  | _ |_)__ _  /_\  _ _  __ _| |___
   |  | _ \ / _` |/ _ \| ' \/ _` | / -_)
   |  |___/_\__, /_/ \_\_||_\__, |_\___|
   |        |___/           |___/
  \*/

  template <typename T>
  BigAngle<T> const &
  BigAngle<T>::operator += ( BigAngle const & rhs ) {
    // add angle
    real_type nx = X*rhs.X - Y*rhs.Y ;
    real_type ny = X*rhs.Y + Y*rhs.X ;
    SIGMA += rhs.SIGMA ;

    // check sign if cross negative X axes
    int ns = ( ny > 0 || (ny == 0 && nx < 0) ) ? 1 : -1 ;
    if ( S*rhs.S > 0 && S*ns < 0 ) SIGMA -= ns ;

    // scale to avoid overflow
    scale(nx,ny) ;

    real_type E = real_type(R+rhs.R+2.829) ;

    I += rhs.I + static_cast<unsigned long>(std::floor(E)) ;
    R  = E-std::floor(E) ;
    X  = nx ;
    Y  = ny ;
    S  = ns ;
  
    LK_ASSERT( I <= I_max, "BigAngle::+=, Out of tolerance!" ) ;
    return *this ;
  }

  // ------------------------------------------------------------------------

  template <typename T>
  BigAngle<T> const &
  BigAngle<T>::operator -= ( BigAngle const & rhs ) {

    // add angle
    real_type nx = Y*rhs.Y + X*rhs.X ;
    real_type ny = Y*rhs.X - X*rhs.Y ;
    SIGMA -= rhs.SIGMA ;

    // check sign if cross negative X axes
    int ns = ( ny > 0 || (ny == 0 && nx < 0) ) ? 1 : -1 ;
    if ( S*rhs.S < 0 && S*ns < 0 ) SIGMA -= ns ;

    // scale to avoid overflow
    scale(nx,ny) ;

    real_type E = real_type(R+rhs.R+2.829) ;

    I += rhs.I + static_cast<unsigned long>(std::floor(E)) ;
    R  = E-std::floor(E) ;
    X  = nx ;
    Y  = ny ;
    S  = ns ;
  
    LK_ASSERT( I <= I_max, "BigAngle::+=, Out of tolerance!" ) ;
    return *this ;
  }

  // ------------------------------------------------------------------------

  template <typename T>
  BigAngle<T> const &
  BigAngle<T>::operator += ( Angle<T> const & rhs ) {

    // add angle X0,Y0 with x,y
    real_type nx = X*rhs.x - Y*rhs.y ;
    real_type ny = X*rhs.y + Y*rhs.x ;
    SIGMA += rhs.sigma ;

    // controllo segno e attraversamento asse X negativo
    int ns = ( ny > 0 || (ny == 0 && nx < 0)) ? 1 : -1 ;
    if ( S*rhs.s > 0 && S*ns < 0 ) SIGMA -= ns ;

    // va riscalato per evitare overflow
    scale(nx,ny) ;
      
    real_type E = real_type(2.829+rhs.e+R) ;
    I += static_cast<unsigned long>(std::floor(E)) ;
    R = E-std::floor(E) ;
    X = nx ;
    Y = ny ;
    S = ns ;

    LK_ASSERT( I <= I_max, "BigAngle::operator +=, Out of tolerance!" ) ;

    return *this ;
  }

  // ------------------------------------------------------------------------

  template <typename T>
  BigAngle<T> const &
  BigAngle<T>::operator -= ( Angle<T> const & rhs ) {

    real_type nx = Y*rhs.y + X*rhs.x;
    real_type ny = Y*rhs.x - X*rhs.y;
    SIGMA -= rhs.sigma ;

    // controllo segno e attraversamento asse X negativo
    int ns = ( ny > 0 || (ny == 0 && nx < 0)) ? 1 : -1 ;
    if ( S*rhs.s < 0 && S*ns < 0 ) SIGMA -= ns ;

    // va riscalato per evitare overflow
    scale(nx,ny) ;
      
    real_type E = real_type(2.829+rhs.e+R) ;
    I += static_cast<unsigned long>(std::floor(E)) ;
    R = E-std::floor(E) ;
    X = nx ;
    Y = ny ;
    S = ns ;

    LK_ASSERT( I <= I_max, "BigAngle::operator +=, Out of tolerance!" ) ;

    return *this ;
  }

  template <typename T>
  void
  BigAngle<T>::checkSigma() const {
    real_type m_pi  = 3.1415926535897932384626433832795028841971693993751 ;
    real_type angle = atan2( Y, X ) ;
    real_type err   = getError() ;
    LK_ASSERT( std::abs(angle) <= err*m_pi/2,
               "BigAngle::checkSigma()\nangle = " << angle <<
               "\nerror = " << err <<
               "\nX = " << X <<
               "\nY = " << Y ) ;
  }

  template <typename T>
  typename BigAngle<T>::real_type
  BigAngle<T>::getAngle() const {
    real_type m_pi  = 3.1415926535897932384626433832795028841971693993751 ;
    return atan2( Y, X )+2*SIGMA*m_pi ;
  }

  template <typename T>
  typename BigAngle<T>::real_type
  BigAngle<T>::getFraction() const {
    real_type m_pi  = 3.1415926535897932384626433832795028841971693993751 ;
    return atan2( Y, X )/(2*m_pi)+SIGMA ;
  }

  // ------------------------------------------------------------------------
  /*\
   |   _    _      _   _           _  _            _
   |  | |  (_)_ _ | |_(_)_ _  __ _| \| |_  _ _ __ | |__  ___ _ _
   |  | |__| | ' \| / / | ' \/ _` | .` | || | '  \| '_ \/ -_) '_|
   |  |____|_|_||_|_\_\_|_||_\__, |_|\_|\_,_|_|_|_|_.__/\___|_|
   |                         |___/
  \*/

  template <typename T>
  LinkingNumber<T>::LinkingNumber( unsigned tot_curve )
  : _lk_computed(0)
  {
    setup( tot_curve ) ;
  }

  template <typename T>
  LinkingNumber<T>::~LinkingNumber() { }

  // ------------------------------------------------------------------------

  template <typename T>
  void
  LinkingNumber<T>::setup( unsigned tot_curve ) {
    curves.clear() ; curves.resize( tot_curve ) ;
    #ifdef LINKING_NUMBER_USE_CXX11
    numThread = std::thread::hardware_concurrency();
    #endif
  }

  // ------------------------------------------------------------------------

  template <typename T>
  void
  LinkingNumber<T>::reset( unsigned ncurve, unsigned reserve_pnts ) {
    LK_ASSERT( ncurve < curves.size(),
               "LinkingNumber::reset( " << ncurve << ", " << reserve_pnts <<
               " ), first argument must be in [0," << curves.size()-1 << "]" ) ;
    curves[ncurve].clear() ;
    if ( reserve_pnts > 0 ) curves[ncurve].reserve(reserve_pnts) ;
  }

  // ------------------------------------------------------------------------

  template <typename T>
  void
  LinkingNumber<T>::init_curve( unsigned  ncurve,
                                real_type x,
                                real_type y,
                                real_type z ) {
    LK_ASSERT( ncurve < curves.size(),
               "LinkingNumber::init_curve( " << ncurve <<
               "...), first argument must be in [0," << curves.size()-1 << "]" ) ;
    Segment S ;
    S.P1[0] = x ;
    S.P1[1] = y ;
    S.P1[2] = z ;
    curves[ncurve].clear() ;
    curves[ncurve].push_back( S ) ;
  }

  // ------------------------------------------------------------------------

  template <typename T>
  void
  LinkingNumber<T>::add_point( unsigned  ncurve,
                               real_type x,
                               real_type y,
                               real_type z,
                               int       weight ) {

    LK_ASSERT( ncurve < curves.size(),
               "LinkingNumber::add_point( " << ncurve <<
               "...), first argument must be in [0," << curves.size()-1 << "]" ) ;

    CURVE & cur = curves[ncurve] ;

    if ( cur.empty() ) {
      Segment S ;
      S.P1[0] = x ;
      S.P1[1] = y ;
      S.P1[2] = z ;
      cur.push_back( S ) ;
    } else {
      Segment & S = cur.back() ;
      Segment S1 ;
      S1.P1[0] = S.P2[0] = x ;
      S1.P1[1] = S.P2[1] = y ;
      S1.P1[2] = S.P2[2] = z ;
      S.weight = weight ;
      cur.push_back( S1 ) ;
    }
  }

  // ------------------------------------------------------------------------

  template <typename T>
  void
  LinkingNumber<T>::close_curve( unsigned  ncurve,
                                 int       weight ) {
    LK_ASSERT( ncurve < curves.size(),
               "LinkingNumber::close_curve( " << ncurve <<
               "), argument must be in [0," << curves.size()-1 << "]" ) ;
    Segment & S0 = curves[ncurve].front() ;
    Segment & S  = curves[ncurve].back() ;

    S.P2[0]  = S0.P1[0] ;
    S.P2[1]  = S0.P1[1] ;
    S.P2[2]  = S0.P1[2] ;
    S.weight = weight ;
  }

  // ------------------------------------------------------------------------

  template <typename T>
  void
  LinkingNumber<T>::add_segment( unsigned  ncurve,
                                 real_type x1,
                                 real_type y1,
                                 real_type z1,
                                 real_type x2,
                                 real_type y2,
                                 real_type z2,
                                 int       weight ) {
    LK_ASSERT( ncurve < curves.size(),
               "LinkingNumber::add_segment( " << ncurve <<
               "...), first argument must be in [0," << curves.size()-1 << "]" ) ;
    Segment S;
    S.weight = weight ;
    S.P1[0] = x1 ; S.P1[1] = y1 ; S.P1[2] = z1 ;
    S.P2[0] = x2 ; S.P2[1] = y2 ; S.P2[2] = z2 ;
    curves[ncurve].push_back(S) ;
  }

  // ------------------------------------------------------------------------

  template <typename T>
  void
  LinkingNumber<T>::add_curve( unsigned        i_curve,
                               real_type const a[][3],
                               unsigned        a_size  ) {
    reset( i_curve, a_size ) ;
    init_curve( i_curve, a[0][0], a[0][1], a[0][2] ) ;
    for ( unsigned i = 1 ; i < a_size ; ++i )
      add_point( i_curve, a[i][0], a[i][1], a[i][2], 1 ) ;
    close_curve( i_curve ) ;
  }

  // ------------------------------------------------------------------------

  template <typename T>
  void
  LinkingNumber<T>::eval_row( real_type const   P1[3],
                              real_type const   P2[3],
                              CURVE     const & curve,
                              BigAngle<T>     & out_angle ) const {
    Angle<T> angle ;
    out_angle.init();
    for ( typename CURVE::const_iterator iQ = curve.begin() ;
          iQ != curve.end() ; ++iQ ) {
      int w = iQ->weight ;
      if ( w == 0 ) continue ;
      angle.build( P1, P2, iQ->P1, iQ->P2 ) ;
      for ( ; w > 0 ; --w ) out_angle += angle;
      for ( ; w < 0 ; ++w ) out_angle -= angle;
    }
  }

  // ------------------------------------------------------------------------

  template <typename T>
  void
  LinkingNumber<T>::eval_rows( int           istart,
                               int           nstep,
                               CURVE const & curveA,
                               CURVE const & curveB,
                               BigAngle<T> & out_angle ) const {
    out_angle.init() ;
    for ( int i = istart ; i < curveA.size() ; i += nstep ) {
      Segment const & S = curveA[i] ;
      int w = S.weight ;
      if ( w != 0 ) {
        BigAngle<T> tmp_angle ;
        eval_row( S.P1, S.P2, curveB, tmp_angle );
        for ( ; w > 0 ; --w ) out_angle += tmp_angle;
        for ( ; w < 0 ; ++w ) out_angle -= tmp_angle;
      }
    }
  }

  // ------------------------------------------------------------------------

  template <typename T>
  void
  LinkingNumber<T>::evaluate( unsigned i_curve,
                              unsigned j_curve,
                              int    & lk ) const {
    LK_ASSERT( i_curve < curves.size() && j_curve < curves.size(),
               "LinkingNumber::eval( " << i_curve << ", " << j_curve << ")\n" <<
               "arguments must be in [0," << curves.size()-1 << "]" ) ;

    CURVE const & Ci = curves[i_curve] ;
    CURVE const & Cj = curves[j_curve] ;

    if ( Ci.empty() || Cj.empty() ) { lk = 0 ; return ; } // no point

    BigAngle<T> total_angle;
    eval_rows( 0, 1, Ci, Cj, total_angle ) ;

    ++_lk_computed ;
    total_angle.checkSigma() ;
    lk = total_angle.getSigma() ;
  }

  #ifdef LINKING_NUMBER_USE_CXX11

  // ------------------------------------------------------------------------

  template <typename T>
  int
  LinkingNumber<T>::eval_mt( unsigned i_curve, unsigned j_curve ) const {

    LK_ASSERT( i_curve < curves.size() && j_curve < curves.size(),
               "LinkingNumber::eval_mt( " << i_curve << ", " << j_curve << ")\n" <<
               "arguments must be in [0," << curves.size()-1 << "]" ) ;

    CURVE const & Ci = curves[i_curve] ;
    CURVE const & Cj = curves[j_curve] ;

    if ( Ci.empty() || Cj.empty() ) return 0 ; // no point

    // launch angle computation
    std::vector<std::thread>  vec_thread(numThread) ;
    std::vector<BigAngle<T> > partial_big_angle(numThread) ;

    for ( unsigned nt = 0 ; nt < numThread ; ++nt ) {
      BigAngle<T> & angle = partial_big_angle[nt] ;
      angle.init();
      vec_thread[nt] = std::thread( &LinkingNumber<T>::eval_rows, this,
                                    nt, numThread,
                                    std::ref(Ci), std::ref(Cj),
                                    std::ref(angle) )  ;
    }

    // resume thread and add angles
    vec_thread[0].join() ;
    for ( unsigned nt = 1 ; nt < numThread ; ++nt ) {
      vec_thread[nt].join() ;
      partial_big_angle[0] += partial_big_angle[nt] ;
    }

    ++_lk_computed ;
    partial_big_angle[0].checkSigma() ;
    return partial_big_angle[0].getSigma() ;
  }

  template <typename T>
  void
  LinkingNumber<T>::evals( unsigned const i_curve[], unsigned ni,
                           unsigned const j_curve[], unsigned nj,
                           int mat[] ) {
    if ( ni == 0 || nj == 0 ) return ; // caso triviale
    // launch angle computation
    std::vector<std::thread> vec_thread(nj) ;
    for ( unsigned i = 0 ; i < ni ; ++i ) {
      unsigned const ii = i_curve[i] ;
      for ( unsigned j = 0 ; j < nj ; ++j ) {
        unsigned const jj = j_curve[j] ;
        int & mij = mat[i+j*ni];
        vec_thread[j] = std::thread( &LinkingNumber<T>::evaluate, this, ii, jj, std::ref(mij) )  ;
      }
      for ( unsigned j = 0 ; j < nj ; ++j )
        vec_thread[j].join() ;
    }
  }

  #else

  template <typename T>
  void
  LinkingNumber<T>::evals( unsigned const i_curve[], unsigned ni,
                           unsigned const j_curve[], unsigned nj,
                           int mat[] ) {
    if ( ni == 0 || nj == 0 ) return ; // caso triviale
    for ( unsigned i = 0 ; i < ni ; ++i )
      for ( unsigned j = 0 ; j < nj ; ++j )
        mat[i+j*ni] = eval(i_curve[i],j_curve[j]) ;
  }

  #endif

  // ------------------------------------------------------------------------

  template <typename T>
  void
  LinkingNumber<T>::writhe_row( int               i_skip,
                                real_type const   P1[3],
                                real_type const   P2[3],
                                CURVE     const & curve,
                                BigAngle<T>     & out_angle ) const {
    Angle<T> angle ;
    out_angle.init();
    for ( typename CURVE::const_iterator iQ = curve.begin() ;
          iQ != curve.end() ; ++iQ ) {
      if ( std::abs( (iQ-curve.begin()) - i_skip ) < 2 ) continue ;
      int w = iQ->weight ;
      if ( w == 0 ) continue ;
      angle.build( P1, P2, iQ->P1, iQ->P2 ) ;
      for ( ; w > 0 ; --w ) out_angle += angle;
      for ( ; w < 0 ; ++w ) out_angle -= angle;
    }
  }

  // ------------------------------------------------------------------------

  template <typename T>
  void
  LinkingNumber<T>::writhe_rows( int           istart,
                                 int           nstep,
                                 CURVE const & curve,
                                 BigAngle<T> & out_angle ) const {
    out_angle.init() ;
    for ( int i = istart ; i < curve.size() ; i += nstep ) {
      Segment const & S = curve[i] ;
      int w = S.weight ;
      if ( w != 0 ) {
        BigAngle<T> tmp_angle ;
        writhe_row( i, S.P1, S.P2, curve, tmp_angle );
        for ( ; w > 0 ; --w ) out_angle += tmp_angle;
        for ( ; w < 0 ; ++w ) out_angle -= tmp_angle;
      }
    }
  }

  // ------------------------------------------------------------------------

  template <typename T>
  typename LinkingNumber<T>::real_type
  LinkingNumber<T>::writhe( unsigned i_curve ) const {
    LK_ASSERT( i_curve < curves.size(),
               "LinkingNumber::writhe( " << i_curve << ")\n" <<
               "arguments must be in [0," << curves.size()-1 << "]" ) ;

    CURVE const & C = curves[i_curve] ;

    if ( C.empty() ) return 0 ; // no point

    BigAngle<T> total_angle;
    writhe_rows( 0, 1, C, total_angle ) ;

    return total_angle.getFraction();
  }

  // explicit instantiation
  template class Angle<float> ;
  template class Angle<double> ;
  template class Angle<ldouble> ;

  template class BigAngle<float> ;
  template class BigAngle<double> ;
  template class BigAngle<ldouble> ;

  template class LinkingNumber<float> ;
  template class LinkingNumber<double> ;
  template class LinkingNumber<ldouble> ;

}
