#include "linking_number.hh"
#include <cmath>
#include <limits>

using namespace std ;

//typedef LK::LinkingNumber<long double> LK_class ;
//typedef LK::LinkingNumber<double> LK_class ;
typedef LK::LinkingNumber<float> LK_class ;


int
main() {

  LK_class lk(3) ;

  cout << "\n\n\n\nTest N.1\n\n\n\n" ;

  lk.reset(0) ;
  lk.reset(1) ;
  lk.reset(2) ;

  std::vector<LK_class::pnt3> c1, c2 ;

  c1.clear() ;
  c2.clear() ;

  LK_class::pnt3 p0 = {0,0,0} ;
  LK_class::pnt3 p1 = {1,0,0} ;
  LK_class::pnt3 p2 = {1,1,0} ;
  LK_class::pnt3 p3 = {0,1,0} ;

  c1.push_back(p0);
  c1.push_back(p1);
  c1.push_back(p2);
  c1.push_back(p3);

  LK_class::pnt3 q0 = {0.5,0.5,-1} ;
  LK_class::pnt3 q1 = {0.5,0.5, 1} ;
  LK_class::pnt3 q2 = {1.5,1.5, 1} ;
  LK_class::pnt3 q3 = {1.5,1.5,-1} ;

  c2.push_back(q0);
  c2.push_back(q1);
  c2.push_back(q2);
  c2.push_back(q3);
  c2.push_back(q0);
  c2.push_back(q1);
  c2.push_back(q2);
  c2.push_back(q3);

  cout << "lk = " << lk.eval( c1, c2 ) << '\n' ;
  cout << "wh(0) = " << lk.writhe( 0 ) << '\n' ;
  cout << "wh(1) = " << lk.writhe( 1 ) << '\n' ;

  int ncurve = 2 ;

  #if 0
  lk.init_curve( ncurve, -1.0, -1.0, 0.0 ) ;
  lk.add_point( ncurve, 1.0, 0.0, 0.0 ) ;
  lk.add_point( ncurve, 0.0, 0.0, 1.0 ) ;
  lk.add_point( ncurve, 0.0, 0.0,-1.0 ) ;
  lk.add_point( ncurve, 0.0, 1.0, 0.0 ) ;
  lk.close_curve( ncurve ) ;
  #else
  lk.init_curve( ncurve, 0.0, 1.0, 0.0 ) ;
  lk.add_point( ncurve, 0.0, 0.0,-1.0 ) ;
  lk.add_point( ncurve, 0.0, 0.0, 1.0 ) ;
  lk.add_point( ncurve, 1.0, 0.0, 0.0 ) ;
  lk.add_point( ncurve, -1.0, -1.0, 0.0 ) ;
  lk.close_curve( ncurve ) ;
  #endif
  cout << "wh(2) = " << lk.writhe( 2 ) << '\n' ;


  return 0 ;
}
