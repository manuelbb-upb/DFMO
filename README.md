# DFMO in Julia
## Derivative-Free Multi-Objective Optimization using Direct-Search

This project aims at providing some variants of the algorithm described in 
[“A Derivative-Free Approach to Constrained Multiobjective Nonsmooth Optimization” by G. Liuzzi, S. Lucidi, and F. Rinaldi](http://epubs.siam.org/doi/10.1137/15M1037810).
Their original *DFMO* solver is provided as a FORTRAN program
and can be obtained at [the Derivative-Free Library](http://www.iasi.cnr.it/~liuzzi/DFL/).
It is distributed under the [GNU GPLv3 license](https://www.gnu.org/licenses/gpl-3.0.en.html).
As I took a look at it from time to time, this package itself is also licensed under GPLv3.

## Examples

### Two Paraboloids
Setup problem dimensions and starting point:
```julia
n = 2
lb = fill(-4, n)
ub = fill(4, n)
x0 = lb .+ (ub .- lb) .* rand(n)
```
Define in-place objective function:
```julia
objfs! = function(y, x)
    y[1] = sum( (x .+ 1).^2 )
    y[2] = sum( (x .- 1).^2 )
end
num_objfs = 2
```

Configure direction scheme to track a single solution: 
```julia
halton_cfg = DFMO.HaltonConfig(;
    stepsize_halving = true,
    solution_selection_strategy = DFMO.FollowChildSolutions(;
        max_num_sols = 1,
        prefer_small_stepsizes = true
    )
)
```

Call optimization loop:
```julia
ret_obj = DFMO.optimize(x0, lb, ub, num_objfs, objfs!; direction_cfg=halton_cfg)
```

Read result matrices:
```julia
X_all = get_x(ret_obj.cache)
F_all = get_fobj(ret_obj.cache)

X_filter = get_x(ret_obj.cache, ret_obj.filter)
F_filter = get_fobj(ret_obj.cache, ret_obj.filter)

X_sol = get_x(ret_obj.cache, ret_obj.filter, ret_obj.dir_scheme)
F_sol = get_fobj(ret_obj.cache, ret_obj.filter, ret_obj.dir_scheme)
```

How often did we call the objective function?
```
ret.num_evals_objfs
```

### Two Paraboloids with Constraints

```julia
n = 2
lb = fill(-4, n)
ub = fill(4, n)
x0 = lb .+ (ub .- lb) .* rand(n)

objfs! = function(y, x)
    y[1] = (x[1] - 2)^2 + (x[2] - 1)^2
    y[2] = (x[1] - 2)^2 + (x[2] + 1)^2
end
num_objfs = 2

constrs! = function(y, x)
    # ‖x‖ ≥ 1 ⇔ 1 - ‖x‖ ≤ 0
    y[1] = 1 - sum(x.^2)
end
num_constrs = 1

ret_objf = DFMO.optimize(x0, lb, ub, num_objfs, objfs!, num_constrs, constrs!)
```

### Settings to Mimic Original Fortran Version
```julia
ret_objf = DFMO.optimize(
    x0, lb, ub, num_objfs, objfs!, num_constrs, constrs!;
    max_func_calls = 2_000,
    direction_cfg = DFMO.HaltonConfig(;
        exploration=true,
        stepsize_halving=true,
        solution_selection_strategy=DFMO.OriginalFiltering{Float64}()
    )    
)
```
In this case, the filter holds all solutions.

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
Then, `cache.current_position[]` stays at `cache.max_store` and `cache.pos_free[]` is the 
next column index to be filled.

In contrast to the Fortran implementation, the set of “solutions” or “iterates” 
and the filter are usually separate objects.
However, all solutions are also in the filter.
Implementation-wise, we went with a pre-allocated and somewhat sparse design:

The filter of type `Filter` basically is a Boolean vector of flags and an integer index vector.
If a flag with index `i` is set, 
then that slot in the filter is occupied by an element in the cache and its index in the 
cache is `filter.index[i]`.

The meaning of “solution” is determined by the `AbstractDirectionScheme` defined by
keyword argument `direction_cfg::AbstractDirectionSchemeConfig`.  
The default setting gives rise to an object `dir_scheme` of type `DFMO.Halton` with 
```julia
dir_scheme.solution_selection_strategy = FollowChildSolutions(; 
    max_num_sols=1, prefer_small_stepsizes=true)
```
Then `dir_scheme.solutions_flags` acts as a subset indicator for the filter slots.
That is, if `dir_scheme.solutions_flags[i] == true`, then the `i`-th element of the 
filter is a solution, has cache index `filter.index[i]` and stepsize `dir_scheme.sz_vals[i]`.
