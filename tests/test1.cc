#include "linking_number.hh"
#include <cmath>
#include <limits>

using namespace std ;

//#define SCALE 1.12121e30
#define SCALE 1.12121e+3

template <typename T>
void
circle( LK::LinkingNumber<T> & lk,
        unsigned ncurve,
        int N ) {
  T dt = 8*atan(1.0)/N ;
  for ( int i = 0 ; i < N ; ++i ) {
    T theta = dt*i ;
    T x = 2.5*cos(3*theta) ;
    T y = 2.5*sin(3*theta) ;
    T z = 0 ;

    x += 1e-2*(rand() % 25) ; // perturbazione artificiale
    y += 1e-2*(rand() % 25) ;
    z += 1e-2*(rand() % 25) ;

    x *= SCALE;
    y *= SCALE;
    z *= SCALE;

    lk.add_point( ncurve, x, y, z ) ;
  }
  lk.close_curve( ncurve ) ;
}

template <typename T>
void
pnts5foil( LK::LinkingNumber<T> & lk,
           unsigned ncurve,
           int N ) {
  T dt = 8*atan(1.0)/N ;
  for ( int i = 0 ; i < N ; ++i ) {
    T theta = dt*i ;

    T x = (7.0/3.0)*sin(2*theta)-(2.0/3.0)*sin(3*theta) ;
    T y = (7.0/3.0)*cos(2*theta)+(2.0/3.0)*cos(3*theta) ;
    T z = 2*sin(5*theta) ;

    x += 1e-2*(rand() % 25) ; // perturbazione artificiale
    y += 1e-2*(rand() % 25) ;
    z += 1e-2*(rand() % 25) ;

    x *= SCALE; // perturbazione artificiale
    y *= SCALE;
    z *= SCALE;

    lk.add_point( ncurve, x, y, z ) ;
  }
  lk.close_curve( ncurve ) ;
}

//typedef LinkingNumber<long double> LK_class ;
typedef LK::LinkingNumber<double> LK_class ;
//typedef LinkingNumber<float> LK_class ;

int
main() {

  LK_class lk(2) ;

  cout << "\n\n\n\nTest N.1\n\n\n\n" ;

  for ( int nseg = 10 ; nseg < 1000 ; nseg *= 2 ) {
    cout << "\n\n\nnseg = " << nseg << '\n' ;
    lk.reset(0) ;
    lk.reset(1) ;
    circle(lk, 0, nseg) ;
    pnts5foil(lk, 1, nseg) ;
    cout << "lk = " << lk.eval(0,1) << '\n' ;
    //lk.info(cout) ;
  }
  
  cout << "\n\n\n\nTest N.2\n\n\n\n" ;
  
  double curve1[10000][3] ;
  double curve2[10000][3] ;

  int nseg = 10000;
  double m_pi = 3.1415926535897932385 ;
  for ( int i = 0 ; i < nseg ; ++i ) {
    double s = double(i)/double(nseg-1) ;

    curve1[i][0] = cos( s * 2*m_pi ) ;
    curve1[i][1] = sin( s * 2*m_pi ) ; 
    curve1[i][2] = 0 ;

    curve2[i][0] = 0 ;
    curve2[i][1] = 1+cos( -s * 2*1231*m_pi ) ;
    curve2[i][2] = sin( -s * 2*1231*m_pi ) ;
  }

  cout << "lk (1231 expected) = " << lk.eval( curve1, nseg, curve2, nseg ) << '\n' ;
  //lk.info(cout) ;

  return 0 ;
}
