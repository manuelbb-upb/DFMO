module DFMO

import LinearAlgebra as LA
import Primes: nextprime    # for `HaltonDirections`
using Printf: @sprintf

import UnPack: @unpack

using Logging: LogLevel, @logmsg, Info

include("CustomPrinting.jl")
using .CustomPrinting

@enum IterCode::UInt8 begin
    CONTINUE_ITERATION=0
    STOP_MAX_FUNC_CALLS=1
    STOP_MIN_STEPSIZE=2
    STOP_MAX_ITER=3
    STOP_SMALL_SPREAD_VALUES
end

include("cache.jl")
include("filter.jl")
include("utils.jl")

include("halton.jl")

struct WorkingArrays{T}
    num_vars :: Int
    num_objfs :: Int
    num_constrs :: Int
    x :: Vector{T}  # hold input values
    fx :: Vector{T} # hold objective values
    cx :: Vector{T} # hold constraint values
    zx :: Vector{T} # hold penalized objective values
    eps_iq :: Matrix{T}
end

function init_working_arrays(x0, num_vars, num_objfs, num_constrs, T)
    x = T.(x0)                  
    fx = zeros(T, num_objfs)    
    cx = zeros(T, num_constrs)  
    zx = copy(fx)
    eps_iq = Matrix{T}(undef, num_constrs, num_objfs)
    return WorkingArrays(num_vars, num_objfs, num_constrs, x, fx, cx, zx, eps_iq)
end

function optimize(
    x0, lb, ub, num_objfs, objfs!, num_constrs=0, constrs! = nothing;
    T::Type{<:AbstractFloat}=Float64,
    direction_cfg::AbstractDirectionSchemeConfig=HaltonConfig(),
    cache::Union{Nothing, EvalCache}=nothing,
    cache_cfg::Union{Nothing, EvalCacheConfig}=nothing,
    cache_max_store::Integer=DEFAULT_CACHE_MAX_STORE,
    cache_tolerance::Real=DEFAULT_CACHE_TOLERANCE(T),
    max_iter::Integer=1_000,
    max_func_calls :: Integer = 2_000,
    # stepsize_halving :: Bool = true,
    max_filter_size::Integer=cache_max_store,
    max_num_solutions::Integer=max_filter_size,
    diagonal_presampling::Bool=true,    
    stepsize_stop :: Real = 10^(-round(abs(log10(eps(T)))^0.8)), # 1e-9 for T==Float64
    coef_delta :: Real = 1,
    sort_by_spread :: Bool = false,
    log_lvl :: Integer = Info.level
)
    log_level = LogLevel(log_lvl)

    num_vars = length(x0)
    box_width = ub .- lb

    @assert num_vars > 0 "Number of variables must at least be 1."
    @assert length(box_width) == num_vars
    @assert num_objfs > 0 "Number of objective functions must at least be 1."
    @assert isnothing(constrs!) || num_constrs > 0 "Constraint function provided but `num_constrs==0`."
    @assert max_func_calls > 0 "`max_func_calls` must be positive."
    if any( lb .> x0 ) || any(ub .< x0 )
        return
    end 

    coef_delta = max(0, coef_delta)

    ## make function handles be counted
    Objfs! = CountedFunction(objfs!)
    Constrs! = CountedFunction(constrs!)

    ## initialize working arrays
    arrays = init_working_arrays(x0, num_vars, num_objfs, num_constrs, T)
    @unpack x, fx, cx, zx, eps_iq = arrays

    ## caches for values
    cache = init_cache(
        T, num_vars, num_objfs, cache, cache_cfg, cache_max_store, cache_tolerance)
    cache_index = 0

    filter = init_filter(cache, max_filter_size)
    #solutions = init_solutions(filter, max_num_solutions, sort_by_spread, T)
    solutions = fill(-1, max_num_solutions)
    previous_solutions = fill(-1, max_num_solutions)

    dir_scheme = init_direction_scheme(direction_cfg, x, zx, lb, ub, solutions, filter)

    it_code = CONTINUE_ITERATION
    ## evaluate objectives and constraints
    ## then use `cx` to compute penalizing factors 
    ## put everything in cache
    cache_index, filter_slot_index, filter_num_free = eval_cache_and_filter!(
        fx, cx, zx, eps_iq, cache, filter, x, Objfs!, Constrs!;
        initialize_eps_iq=true, check_cache=false, max_func_calls
    )
    if cache_index < 0 
        it_code = STOP_MAX_FUNC_CALLS
    end
    
    if diagonal_presampling
        box_width = ub .- lb
        i = 1
        while it_code == CONTINUE_ITERATION && i <= num_vars
            x .= lb .+ ((i-1)/(num_vars -1)) .* box_width
            
            cache_index, filter_slot_index, filter_num_free = eval_cache_and_filter!(
                fx, cx, zx, eps_iq, cache, filter, x, Objfs!, Constrs!;
                initialize_eps_iq=false, check_cache=true, max_func_calls
            )
            if cache_index < 0 
                it_code = STOP_MAX_FUNC_CALLS
            end

            i += 1
        end
    end
    spread_vals = sort_by_spread ? Vector{T}(undef, length(solutions)) : nothing

    for it_index=1:max_iter
        it_code != CONTINUE_ITERATION && break

        @logmsg log_level """\n
        #################################################
        # Iteration $(it_index) (of max $(max_iter)).
        #################################################"""

        #=
        The following condition from the Fortran code tells us when to inspect a cached element:
        ```
         condizione = (Lkappa(pn,n+q+1)>50.d0).and.   &
            (Lkappa(pn,n+q+2) > alfa_stop).and. &
            ((deltaf(pn) >= coef_delta).or..not.flag_spread).and. &
            ((Lkappa(pn,n+q+2) > 1.d-1*Lmaxalfa).or.(pn<=q))
        ```
        The condition translates to 
        "element is in filter" AND 
        "stepsize is larger than `alfa_stop`" AND
        "spread is larger than `coef_delta`" AND
        "stepsize is relatively large OR its an extremal element"

        The negation of the condition is
        "element is not in filter" OR
        "stepsize <= `alfa_stop`" OR
        "spread < `coef_delta`" OR 
        "stepsize relatively small AND not an extremal element"
        =#

        it_code = propagate_solutions!(
            dir_scheme, arrays, cache, solutions, previous_solutions, filter, 
            Objfs!, Constrs!, spread_vals,
            lb, ub, log_level, stepsize_stop, max_func_calls, it_index)
    end# for it_index=1:max_iter
    it_code = it_code == CONTINUE_ITERATION ? STOP_MAX_ITER : it_code

    return it_code, cache, filter, solutions, dir_scheme
end#function main

end#module