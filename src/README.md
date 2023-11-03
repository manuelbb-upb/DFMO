# Test-Suite for Nonlinear Multi-Objective Optimization

This package provides tools to work with the test problems defined in 
[“A Derivative-Free Approach to Constrained Multiobjective Nonsmooth Optimization,” G. Liuzzi, S. Lucidi, and F. Rinaldi](http://epubs.siam.org/doi/10.1137/15M1037810)
from within Julia (1.8+).
The set of test problems *TESTMO* is given as FORTRAN90 code, and the problem interface is compatible
with the *DFMO* solver, that the authors also provide as a FORTRAN program.
Both packages can be obtained at [the Derivative-Free Library](http://www.iasi.cnr.it/~liuzzi/DFL/).
They are distributed under the [GNU GPLv3 license](https://www.gnu.org/licenses/gpl-3.0.en.html).
This enables us to also include their code in this repository (which itself is licensed under
GPLv3).


