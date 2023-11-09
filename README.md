# DFMO in Julia
## Derivative-Free Multi-Objective Optimization using Direct-Search

This project aims at providing some variants of the algorithm described in 
[“A Derivative-Free Approach to Constrained Multiobjective Nonsmooth Optimization” by G. Liuzzi, S. Lucidi, and F. Rinaldi](http://epubs.siam.org/doi/10.1137/15M1037810).
Their original *DFMO* solver is provided as a FORTRAN program.
and can be obtained at [the Derivative-Free Library](http://www.iasi.cnr.it/~liuzzi/DFL/).
It is distributed under the [GNU GPLv3 license](https://www.gnu.org/licenses/gpl-3.0.en.html).
As I took a look at it from time to time, this package itself is also licensed under GPLv3.

## Implementation Details

Just like in the Fortran version, we have a cache to store evaluation values.
It is of type `EvalCache` and the data arrays are referenced by fields `x`, `fobj`
and `viol`.
`x` is a matrix, the columns of which hold evaluation argument vectors.
The columns of `fobj` hold objective function values. 
If the problem has constraints, `fobj` holds the **penalized** objective values.
The vector `viol` stores constraint violation values.
We use the columns cyclically: the last entry has position `cache.current_position[]` 
and when that position would exceed `cache.max_store`, we start at 1 again.

In contrast to the Fortran implementation, the set of solutions and the filter are separate.
All solutions are in the filter, and we implemented this change to allow for tracking 
very few solutions or propagate a single solution.
Implementation-wise, we went with a pre-allocated and somewhat sparse design:

The filter basically is a Boolean vector of flags and an integer index vector.
If a flag with index `i` is set, 
then that slot in the filter is occupied by an element in the cache and its index in the 
cache is `filter.index[i]`.
Likewise, the solutions index into the filter (not the cache). 
That is, if `solutions.flags[i]==true`, then the `i`-th solution slot is occupied, by a filter
element with index `solutions.index[i]` and cache index `filter.index[solutions.index[i]]`.
