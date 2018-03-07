/****************************************************************************\
  Copyright (c) Enrico Bertolazzi 2018
  All Rights Reserved.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the file license.txt for more details.
\****************************************************************************/

#include "linking_number.hh"
#include "mex_utils.hh"

#define MEX_ERROR_MESSAGE \
"======================================================================\n" \
"\n" \
"lk: Compute Linking number for 2 polygonal curves, or\n" \
"    Writhe number for a polygonal curve, or\n" \
"\n" \
"USAGE:\n" \
"  [L,err] = lk(p,q) ;\n" \
"  [W,err] = lk(p) ;\n" \
"\n" \
"On input:\n" \
"\n" \
"  p = matrix n x 3 of polygon curve 1\n" \
"  q = matrix m x 3 of polygon curve 2\n" \
"\n" \
"On output:\n" \
"\n" \
"  L   = Linking number\n" \
"  W   = Writhe number\n" \
"  err = if |err| < 0.5 the linking number is certified\n" \
"\n" \
"======================================================================\n" \
"\n" \
"Autor: Enrico Bertolazzi\n" \
"       Department of Industrial Engineering\n" \
"       University of Trento\n" \
"       enrico.bertolazzi@unitn.it\n" \
"\n" \
"======================================================================\n"

using namespace std ;

namespace LK {

  extern "C"
  void
  mexFunction( int             nlhs,
               mxArray       * plhs[],
               int             nrhs,
               mxArray const * prhs[] ) {

    try {

      MEX_ASSERT( nlhs == 2, "lk(...), expected 2 output, found " << nlhs ) ;

      // check for proper number of arguments
      if ( nrhs == 1 ) {

        // compute writhe number

        mwSize number_of_dimensions = mxGetNumberOfDimensions(prhs[0]) ;
        MEX_ASSERT( number_of_dimensions == 2, "lk(p,q), p must be a matrix" ) ;
        mwSize const * dims = mxGetDimensions(prhs[0]) ;
        mwSize nr = dims[0];
        MEX_ASSERT( nr > 2 && dims[0] == 3,
                    "lk(p), p expected to be an (3 x n) matrix with n >= 3, found " <<
                    dims[0] << " x " << dims[1] ) ;
        MEX_ASSERT( mxIsDouble(prhs[0]),
                    "lk(p), p expected to be a matrix with real number" ) ;
        double * p = mxGetPr(prhs[0]) ;

        LK::LinkingNumber<double> lk(1) ;
        lk.reset(0) ;

        for ( int i = 0 ; i < nr ; ++i, p += 3 )
          lk.add_point( 0, p[0], p[1], p[2] ) ;
        lk.close_curve( 0 ) ;

        double W, err ;

        W = lk.writhe( 0, err );

        setScalarValue( plhs[0], W );
        setScalarValue( plhs[1], err );

      } else if ( nrhs == 2 ) {

        mwSize number_of_dimensions = mxGetNumberOfDimensions(prhs[0]) ;
        MEX_ASSERT( number_of_dimensions == 2, "lk(p,q), p must be a matrix" ) ;
        mwSize const * dims = mxGetDimensions(prhs[0]) ;
        mwSize nr = dims[1];

        MEX_ASSERT( nr > 2 && dims[0] == 3,
                    "lk(p,q), p expected to be an (3 x n) matrix with n >= 3, found " <<
                    dims[0] << " x " << dims[1] ) ;
        MEX_ASSERT( mxIsDouble(prhs[0]),
                    "lk(p,q), p expected to be a matrix with real number" ) ;

        double * p = mxGetPr(prhs[0]) ;

        number_of_dimensions = mxGetNumberOfDimensions(prhs[1]) ;
        MEX_ASSERT( number_of_dimensions == 2, "lk(p,q), q must be a matrix" ) ;
        dims = mxGetDimensions(prhs[1]) ;
        mwSize nc = dims[1];
        MEX_ASSERT( nc > 2 && dims[0] == 3,
                    "lk(p,q), q expected to be an (3 x n) matrix with n >= 3, found " <<
                    dims[0] << " x " << dims[1] ) ;
        MEX_ASSERT( mxIsDouble(prhs[1]),
                    "lk(p,q), q expected to be a matrix with real number" ) ;
        double * q = mxGetPr(prhs[1]) ;

        LK::LinkingNumber<double> lk(2) ;
        lk.reset(0) ;
        lk.reset(1) ;

        for ( int i = 0 ; i < nr ; ++i,  p += 3 )
          lk.add_point( 0, p[0], p[1], p[2] ) ;
        lk.close_curve( 0 ) ;

        for ( int j = 0 ; j < nc ; ++j,  q += 3 )
          lk.add_point( 1, q[0], q[1], q[2] ) ;
        lk.close_curve( 1 ) ;

        int    lknum ;
        double err;
        lk.evaluate( 0, 1, lknum, err );

        setScalarInt( plhs[0], lknum );
        setScalarValue( plhs[1], err );

      } else {
        MEX_ASSERT( false, "MEX_ERROR_MESSAGE" ) ;
      }

    } catch ( std::exception const & e ) {
      mexErrMsgTxt(e.what()) ;

    } catch (...) {
      mexErrMsgTxt("lk failed\n") ;
    }
  }
}
