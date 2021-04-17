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
"lk: Compute Linking number for n polygonal curves\n" \
"\n" \
"USAGE:\n" \
"  [L,err] = lk(p1,p2,...,pn);\n" \
"\n" \
"On input:\n" \
"\n" \
"  p1,p2,...,pn = matrices nk x 3 of polygon curve\n" \
"\n" \
"On output:\n" \
"\n" \
"  L        = Linking number matrix\n" \
"  L(i,j)   = Linking number of curve i vs curve j, L(i,i) = 0\n" \
"  err      = linking number error matrix\n" \
"  err(i,j) = L(i,j) error, if |err(i,j)| < 0.5 the L(i,j) is certified\n" \
"\n" \
"======================================================================\n" \
"\n" \
"Autor: Enrico Bertolazzi\n" \
"       Department of Industrial Engineering\n" \
"       University of Trento\n" \
"       enrico.bertolazzi@unitn.it\n" \
"\n" \
"======================================================================\n"

using namespace std;

namespace LK {

  extern "C"
  void
  mexFunction(
    int             nlhs,
    mxArray       * plhs[],
    int             nrhs,
    mxArray const * prhs[]
  ) {

    try {

      MEX_ASSERT( nlhs == 2, "lk(...), expected 2 output, found " << nlhs );

      // check for proper number of arguments
      if ( nrhs > 1 ) {

        LinkingNumber<double> lk(nrhs);

        for ( unsigned ncurve = 0; ncurve < nrhs; ++ncurve ) {

          mxArray const * & rhs = prhs[ncurve];

          mwSize number_of_dimensions = mxGetNumberOfDimensions(rhs);

          MEX_ASSERT(
            number_of_dimensions == 2,
            "lk(p1,...,pn), p" << ncurve << " must be a matrix"
          );
          mwSize const * dims = mxGetDimensions(rhs);

          unsigned nr = dims[0];
          unsigned nc = dims[1];

          MEX_ASSERT(
            nc > 2 && nr == 3,
            "lk(p1,...,pn), p" << ncurve <<
            " expected to be an (3 x n) matrix with n >= 3, found " <<
            nr << " x " << nc
          );
          MEX_ASSERT(
            mxIsDouble(rhs),
            "lk(p1,...,pn), p" << ncurve <<
            " expected to be a matrix with real number"
          );

          double * p = mxGetPr(rhs);

          lk.reset(ncurve);

          for ( int i = 0; i < nc; ++i,  p += 3 )
            lk.add_point( ncurve, p[0], p[1], p[2] );

          lk.close_curve(ncurve);
        }

        int64_t * imat = createMatrixInt64( plhs[0], nrhs, nrhs );
        double  * rmat = createMatrixValue( plhs[1], nrhs, nrhs );

        for ( unsigned i = 0; i < nrhs; ++i ) {
          imat[i*(nrhs+1)] = 0;
          rmat[i*(nrhs+1)] = 0;
          for ( unsigned j = i+1; j < nrhs; ++j ) {
            int    lknum;
            double err;
            lk.evaluate( i, j, lknum, err );
            unsigned ij = i+j*nrhs;
            unsigned ji = j+i*nrhs;
            imat[ij] = imat[ji] = lknum;
            rmat[ij] = rmat[ji] = err;
          }
        }

      } else {
        MEX_ASSERT( false, MEX_ERROR_MESSAGE );
      }

    } catch ( std::exception const & e ) {
      mexErrMsgTxt(e.what());

    } catch (...) {
      mexErrMsgTxt("lk failed\n");
    }
  }
}
