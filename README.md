Linking Number
==============

A C++ class with MATLAB interface for the computation of the linking number.
Given two parametric closed curves $\mathbf{p}(t)$
and $\mathbf{q}(t)$ where

$$
    \begin{array}{rcl}
      [0,N] & \to & \mathbb{R}^3 \cr
          t & \to & \mathbf{p}(t)
    \end{array}
    \qquad
    \begin{array}{rcl}
      [0,M] & \to & \mathbb{R}^3 \cr
          t & \to & \mathbf{q}(t)
    \end{array}
$$

and $\mathbf{p}(0)=\mathbf{p}(N)$,
$\mathbf{q}(0)=\mathbf{q}(M)$.
The linking number $L(\mathbf{p},\mathbf{q})$
of the two curve $\mathbf{q}(t)$
and $\mathbf{q}(t)$ is defined as the double integral

$$
    L(\mathbf{p},\mathbf{q})= \frac{1}{4\pi}
    \int_{0}^N
    \int_{0}^M
    \frac{(\mathbf{q}(s)-\mathbf{p}(t))\cdot(\mathbf{q}'(s)\times\mathbf{p}'(t))}
         {||\mathbf{q}(s)-\mathbf{p}(t)||^3}
    \mathrm{d}s\mathrm{d}t
$$

which can be proved its an integer (if the curve do not intersect).

The library can compute the linking number of two piecewise linear curves where the curves are defined by the support points

$$
    \{\mathbf{p}_0,\mathbf{p}_1,\ldots,\mathbf{p}_N\}
    \qquad
    \{\mathbf{q}_0,\mathbf{q}_1,\ldots,\mathbf{q}_M\}
$$

C++ interface
-------------

Load the library

```cpp
#include "linking_number.hh"
```

Instantiate the liking number manager:

```cpp
unsigned max_curve_stored = 2;
LK_class lk(max_curve_stored);
```

the class ``lk`` can store more than two curves.
Now insert the curves by points, incrementally

```cpp
unsigned ncurve=0;
lk.reset( ncurve ); // initialize curve

lk.add_point( ncurve, x0, y0, z0 );
lk.add_point( ncurve, x1, y1, z1 );

// ...

lk.add_point( ncurve, xn, yn, zn );
lk.close_curve( ncurve );
```

and compute linking number

```cpp
unsigned ncurve0=0;
unsigned ncurve1=1;

int lknumber = lk.eval(ncurve0,ncurve1);
```

A curve can be inserted by passing all the points
stored in a vector

```cpp
std::vector<Tvec> a;

// ... fill the vector

lk.add_curve( ncurve, pnts );
```

where ``Tvec`` is any structure or class such that the field ``x``, ``y`` and ``z`` are available. Another possibility is to use a  ``n x 3`` matrix as follows

```cpp
double curve[npts][3];

// fill the points

lk.add_curve( ncurve, curve, npts );
```

linking number can be computed in one shot passing
two vectors or two matrices ``n x 3`` as follows

```cpp
// mode 1

std::vector<Tvec> a;
std::vector<Tvec> b;

// ... fill the vectors

int lknumber = lk.eval( a, b );

// mode 2

double curve1[10000][3];
double curve2[10000][3];

unsigned nseg1 = xx; // number of segment of the first curve
unsigned nseg2 = xx; // number of segment of the second curve

// fill the points

int lknumber = lk.eval( curve1, nseg, curve2, nseg );
```

Having loaded ``m`` curves it is possible to compute directly a matrix of linking number where ``Mij`` is the linking number of the i-th curve vs j-th curve.
For example

```cpp
unsigned i_curve[] = {2,0};
unsigned j_curve[] = {0,2,1};

unsigned ni = 2;
unsigned nj = 3;

int mat[6];

lk.evals( i_curve, ni, j_curve, nj, mat );
```

and matrix ``mat`` contains (using Fortran addressing ``mat(i,j) = mat[i+2*j]``). the following linking numbers:

|              |              |              |
| ------------ | ------------ | ------------ |
| $L(c_2,c_0)$ | $L(c_2,c_2)$ | $L(c_2,c_1)$ |
| $L(c_0,c_0)$ | $L(c_0,c_2)$ | $L(c_0,c_1)$ |


MATLAB interface
----------------

Matlab usage is very easy, there is a unique command
``lk`` with a variable number of arguments
``lk(p1,p2,...,pn)`` where ``pk`` is a ``n x 3`` matrix
storing the points of the piecewise linear curve.

```cpp
[L,err] = lk(p1,p2,...,pn);
```

|            |                                                        |
| ---------- | ------------------------------------------------------ |
| `L`        | Linking number matrix                                  |
| `L(i,j)`   | Linking number of curve `i` vs curve `j`, `L(i,i) = 0` |
| `err`      | linking number error matrix                            |
| `err(i,j)` | if `abs(err(i,j)) < 0.5` the `L(i,j)` is certified     |

Installation
------------

Download the toolbox from <https://github.com/ebertolazzi/LinkingNumber/releases> and install as usual.
After installation compile the `mex` interface by running ``CompileLK`` on the MATAB command window.

Citation
--------

If you use this work please cite one or more work in the reference sections.
To cite specifically the software use the
[![DOI](https://zenodo.org/badge/358880071.svg)](https://zenodo.org/badge/latestdoi/358880071) of Zenodo service.

References
----------

The algorithm used in the library is detailed in

- **Enrico Bertolazzi, Riccardo Ghiloni, Ruben Specogna**,
  *Efficient computation of Linking number with certification*, 2019,
  link at arxiv <https://arxiv.org/abs/1912.13121>

This implementation of linking number was used in a various
paper here listed

- **Ana Alonso Rodríguez, Enrico Bertolazzi, Alberto Valli**,
  *The curl-div system: theory and finite element approximation*,
  chapter of the book: *Maxwell’s Equations: Analysis and Numerics*,
  de Gruyter., 20-19

- **Ana Alonso Rodríguez, Enrico Bertolazzi, Riccardo Ghiloni, Ruben Specogna**,
  *Efficient construction of 2-chains representing a basis of* $H_2(\Omega^-,\partial\Omega;\mathbb{Z})$,
  Advances in Computational Mathematics, vol.44, n.5, 2018

- **Ana Alonso Rodríguez, Enrico Bertolazzi, Riccardo Ghiloni, Ruben Specogna**,
  *Efficient construction of 2-chains with a prescribed boundary*,
  SIAM Journal on Numerical Analysis, vol.55, n.3, 2017

- **Ana Alonso Rodriguez, Enrico Bertolazzi, Riccardo Ghiloni, Alberto Valli**,
  *Finite element simulation of eddy current problems using magnetic scalar potentials*,
  Journal of Computational Physics, vol. 294, 2015

- **Ana Alonso Rodríguez, Enrico Bertolazzi, Riccardo Ghiloni, Alberto Valli**,
  *Construction of a finite element basis of the first de Rham cohomology group
  and numerical solution of 3D magnetostatic problems*,
  SIAM Journal on Numerical Analysis, vol.51, N.4, 2013
