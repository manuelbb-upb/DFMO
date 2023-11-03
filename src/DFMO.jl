module DFMO

import LinearAlgebra as LA
import Primes: nextprime    # for `HaltonDirections`
using Printf: @sprintf

using Logging: LogLevel, @logmsg, Info

include("CustomPrinting.jl")
using .CustomPrinting

struct CountedFunction{F}
    func! :: F
    count :: Base.RefValue{Int}

    CountedFunction(func!::F) where {F} = new{F}(func!, Ref(0))
end
function (cf::CountedFunction{<:Function})(y, x)
    cf.func!(y, x)
    cf.count[] += 1
end
(cf::CountedFunction{<:Nothing})(y, x)=nothing
num_calls(cf::CountedFunction) = cf.count[]

abstract type AbstractDirectionScheme end
init_direction_scheme(::Type{<:AbstractDirectionScheme}, num_vars, lb, ub, T)=nothing
dir_type(::AbstractDirectionScheme)=Vector{<:Real}
init_stepsize_cache(::AbstractDirectionScheme, max_store)=Vector{Real}(undef, max_store)
set_stepsize_in_cache!(::AbstractDirectionScheme, sz_cache, sz, cache_index, dir_index; force_set=false, kwargs...)=nothing
get_stepsize_from_cache(::AbstractDirectionScheme, sz_cache, cache_index, dir_index)=nothing
max_stepsize_in_filter(::AbstractDirectionScheme, sz_cache, filter_flags, filter_permutation)=nothing
max_stepsize_in_cache(ds::AbstractDirectionScheme, sz_cache, cache_index)=get_stepsize_from_cache(ds, sz_cache, cache_index, 1)
prepare_for_next_iteration!(::AbstractDirectionScheme; log_level, kwargs...)=nothing
reset_current_state!(::AbstractDirectionScheme)=nothing
"Set `dir_vec` and return a direction index of type `Integer` or return `nothing` to indicate that all directions have been tried."
next_direction!(dir_vec, ::AbstractDirectionScheme)=nothing

struct HaltonDirections{T<:AbstractFloat} <: AbstractDirectionScheme 
    num_vars :: Int
    sz_init :: T
    index :: Base.RefValue{Int}
    orth_dirs :: Matrix{T}
    next_dir_index :: Base.RefValue{Int}
    last_stepsize_change :: Base.RefValue{Int}
end

# Helper for `init_direction_scheme`.
# Compute the `ih`-th Halton point in `num_vars` dimensions.
# Code very much like in Fortran source/Wikipedia, but there are better ways to do this! (HaltonSequences.jl)
function halton(num_vars, ih, T)
    @assert num_vars > 0
    dir = Vector{T}(undef, num_vars)
    if num_vars == 1
        dir[end] = 1
        return dir
    end
    base = 2
    for i=1:num_vars
        base = i==1 ? base : nextprime(base, 2)
        dir[i] = zero(T)
        f = one(eltype(dir)) / base
        j = ih
        while j > 0
            j, r = divrem(j, base)
            dir[i] += r * f
            f /= base
        end
    end
    return dir
end

function halton_base!(H, d)
    #=
    n = length(d)
    H .= 0
    H[:, 1] .= d
    
    imax = argmax( abs.(d) )
    for i=2:imax
        H[i-1, i] = 1
    end
    for i=imax+1:n
        H[i,i] = 1
    end

    _H, _ = LA.qr(H)
    H .= _H
    =#
    _H, _ = LA.qr(d)
    H .= _H
    return nothing
end

function init_direction_scheme(::Type{<:HaltonDirections}, num_vars, lb, ub, T)
    ih = 1_000 + 2 * num_vars
    d0 = halton(num_vars, ih, T)
    orth_dirs = Matrix{T}(undef, num_vars, num_vars)
    halton_base!(orth_dirs, d0)
    box_width = ub .- lb
    sz_init = min(10, maximum(box_width ./ 10))
    return HaltonDirections{T}(num_vars, sz_init, Ref(ih), orth_dirs, Ref(1), Ref(-1))
end
dir_type(::HaltonDirections{T}) where{T}=Vector{T}
function init_stepsize_cache(::HaltonDirections{T}, max_store)where{T}
    #Vector{T}(undef, max_store)
    return Vector{Union{T, Missing}}(undef, max_store)
end
function set_stepsize_in_cache!(hd::HaltonDirections, sz_cache, sz, cache_index, dir_index; force_set=false)
    if force_set || hd.last_stepsize_change[] != cache_index
        sz_cache[cache_index] = sz
        hd.last_stepsize_change[] = cache_index
    end
    return nothing
end
function get_stepsize_from_cache(hd::HaltonDirections, sz_cache, cache_index, dir_index)
    sz = sz_cache[cache_index]
    if ismissing(sz)
        sz = hd.sz_init
        sz_cache[cache_index] = sz
    end
    return sz
end
function max_stepsize_in_filter(hd::HaltonDirections, sz_cache, filter_flags, filter_permutation)
    for (i, j) in enumerate(filter_permutation)
        !filter_flags[i] && continue
        sz = sz_cache[j]
        if !ismissing(sz)
            return sz
        end
    end
    return hd.sz_init
end

function prepare_for_next_iteration!(hd::HaltonDirections{T}; log_level=Info) where T
    ih = hd.index[] += 2*hd.num_vars
    @logmsg log_level "\tComputing halton direction for halton index $(ih)."
    d0 = halton(hd.num_vars, ih, T)
    halton_base!(hd.orth_dirs, d0)
    hd.next_dir_index[] = 1
    return nothing
end

function reset_current_state!(hd::HaltonDirections)
    hd.next_dir_index[] = 1
    return nothing
end

function next_direction!(dir_vec, hd::HaltonDirections)
    i = hd.next_dir_index[]
    if i <= hd.num_vars
        dir_vec .= hd.orth_dirs[:, i]
        hd.next_dir_index[] = i + 1
        return i
    end
    return nothing
end

Base.@kwdef struct EvalCache{T<:AbstractFloat}
    ## These parameters were module parameters in Fortran    
    max_store :: Int = 100_000
    tolerance :: T = 10^(-ceil(abs(log10(eps(T)))^0.65)) # 1e-6 for T==Float64
    ncache :: Base.RefValue{Int} = Ref(0)   # cache access counter

    ## Actual caching matrices
    x :: Matrix{T}
    fobj :: Matrix{T}
    viol :: Vector{T}

    pos_free :: Base.RefValue{Int} = Ref(1)     # next free index in caches
    cur_dim :: Base.RefValue{Int} = Ref(0)      # maximum index of column holding a result
end

function setup_cache(num_vars, num_objfs, num_constrs; 
    max_store::Int=100_000, T::Type{<:AbstractFloat}=Float64
)
    x = Matrix{T}(undef, num_vars, max_store)
    fobj = Matrix{T}(undef, num_objfs, max_store)
    viol = zeros(T, max_store)
#    fconstr = Matrix{T}(undef, num_constrs, max_store)
#    return EvalCache(; max_store, x, fobj, fconstr, sz_cache)
    return EvalCache{T}(; max_store, x, fobj, viol)
end

function find_in_cache(x_query, cache)
    x_query_index = 0
    cache.ncache[] += 1
    eps = cache.tolerance
    for (j, x_c)=enumerate(eachcol(cache.x))
        if j > cache.cur_dim[]
            break
        end
        if all( abs.( x_query .- x_c ) .<= eps )
            x_query_index = j
            break
        end
    end
    return x_query_index
end

function insert_in_cache!(cache, x, fobj, constr_viol; check_cache=true)
    new_cache_index = if check_cache
        find_in_cache(x, cache)
    else
        -1
    end
    if new_cache_index != 0
        new_cache_index = j = cache.pos_free[]
        cache.x[:, j] .= x
        cache.fobj[:, j] .= fobj
        cache.viol[j] = constr_viol
        j = j < cache.max_store ? j+1 : 1
        cache.pos_free[] = j
        if cache.cur_dim[] < cache.max_store
            cache.cur_dim[] += 1
        end
    else
        @warn "`insert_in_cache!`: Point is in cache already, at pos $(new_cache_index)."
    end
    return new_cache_index
end

function main(
    x0, lb, ub, num_objfs, objfs!, num_constrs=0, constrs! =nothing;
    T::Type{<:AbstractFloat}=Float64,
    direction_cfg::Type{<:AbstractDirectionScheme}=HaltonDirections,
    max_store::Int=100_000,
    max_iter::Int=1_000,
    max_func_calls :: Integer = 2_000,
    # stepsize_halving :: Bool = true,
    max_filter_size::Int=max_store,
    diagonal_presampling::Bool=true,    
    stepsize_stop :: Real = 10^(-round(abs(log10(eps(T)))^0.8)), # 1e-9 for T==Float64
    eta :: Real = 10^(-ceil(abs(log10(eps(T)))^0.65)), # 1e-6 for T==Float64
    gamma :: Real = eta,
    coef_delta :: Real = 1,
    sort_by_spread :: Bool = false,
    log_lvl :: Int = Info.level
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

    ## make function handles be counted
    Objfs! = CountedFunction(objfs!)
    Constrs! = CountedFunction(constrs!)

    ## initialize working arrays
    x = T.(x0)                  # hold input values
    step = zero(x)
    x_trial = copy(x)
    x_beta = copy(x)
    step_linesearch = copy(step)
    fx = zeros(T, num_objfs)    # hold objective values, `fob` in Fortran
    cx = zeros(T, num_constrs)  # hold constraint values, `ciq` in Fortran
    zx = copy(fx)               # hold penalized objective values, `fpen` or `finiz` in Fortran
    zx_trial = copy(zx)
    zx_beta = copy(zx)
    d = similar(x)              # direction vector
    cache_index = 0

    ## initialize an object defining how directions are computed
    dir_scheme = init_direction_scheme(direction_cfg, num_vars, lb, ub, T)
    ## and caches for values …
    cache = setup_cache(num_vars, num_objfs, num_constrs; max_store, T)
    ## … and stepsizes 
    sz_cache = init_stepsize_cache(dir_scheme, max_store)

    ## initialize filter as a set of two indices
    ## if `filter_flags[i] == true`, then `j=filter_permutation[i]`
    ## is the index of a point in the filter belonging to `cache` and `sz_cache`.
    ## most of the time we loop through `filter_permutation` and then check `filter_flags`
    filter_permutation = fill(-1, max_filter_size)
    filter_flags = zeros(Bool, max_filter_size)
    filter_spread_vals = Vector{T}(undef, max_filter_size)

    ## evaluate objectives and constraints
    ## then use `cx` to compute penalizing factors 
    ## put everything in cache
    eps_iq = Matrix{T}(undef, num_constrs, num_objfs)
    cache_index = eval_cache_and_filter!(
        fx, cx, zx, eps_iq, cache, filter_permutation, filter_flags, x, Objfs!, Constrs!;
        initialize_eps_iq=true, check_cache=false, max_func_calls
    )
    
    if diagonal_presampling
        box_width = ub .- lb
        for i=1:num_vars
            x .= lb .+ ((i-1)/(num_vars -1)) .* box_width
            
            cache_index = eval_cache_and_filter!(
                fx, cx, zx, eps_iq, cache, filter_permutation, filter_flags, x, Objfs!, Constrs!;
                initialize_eps_iq=false, check_cache=true, max_func_calls
            )
        end
    end

    tol_sgn_switch = 10^(-round(abs(log10(eps(Float64)))^0.76)) # 1e-8 for T == Float64
    ## tol_success = 10^(-ceil(abs(log10(eps(Float64)))^0.9)) 1e-12
    next_filter_flags = deepcopy(filter_flags)
    next_filter_permutation = deepcopy(filter_permutation)
    for it_index=1:max_iter
        @logmsg log_level """\n
        #################################################
        # Iteration $(it_index) (of max $(max_iter)).
        #################################################"""
        if sort_by_spread
            compute_spread_and_sort!(filter_spread_vals, filter_permutation, filter_flags, cache, coef_delta)
        end
        sz_filter_max = max_stepsize_in_filter(dir_scheme, sz_cache, filter_flags, filter_permutation)

        next_filter_flags .= filter_flags
        next_filter_permutation .= filter_permutation

        num_points_found = 0
        early_exit = false

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

        all_stepsizes_too_small = true
        for (i, ci) in enumerate(filter_permutation)
            !filter_flags[i] && continue
            @logmsg log_level "Cached element $(ci) is in filter." 

            ## inverting `condizione` (DFMO.f90, ll. 184):
            ## a) consider only filter elements with stepsize big enough
            sz_max = max_stepsize_in_cache(dir_scheme, sz_cache, ci)
            if sz_max <= stepsize_stop 
                @logmsg log_level "\tStepsize is $(sz_max) <= $(stepsize_stop), skipping."
                continue
            end
            all_stepsizes_too_small = false

            ## b) favor larger spread values
            spread = filter_spread_vals[i]
            if spread < coef_delta 
                @logmsg log_level "Spread value is $(spread) < $(coef_delta), skipping."
                continue
            end
            
            ## c) if an element is not filter-extremal, it has to have large *relative* stepsize
            sz_max <= sz_filter_max/10 && i > num_objfs && continue 

            num_points_found += 1

            x .= cache.x[:, ci]
            zx .= cache.fobj[:, ci]
            @logmsg log_level "Spread value is $(spread), vectors are\n$(CustomPrinting.pretty_str(CustomPrinting.VectorTable("x"=>x, "zx" => zx)))"
            #src viol = cache.viol[ci]

            dir_index = next_direction!(d, dir_scheme)
            while !isnothing(dir_index) && !early_exit
                alpha = get_stepsize_from_cache(dir_scheme, sz_cache, ci, dir_index)
                for sgn_index=1:2
                    @logmsg log_level "Testing direction $(dir_index) and α=$(alpha), try $(sgn_index); d=$(d)"

                    @. x_trial = x + alpha * d
                    project_into_box!(x_trial, lb, ub)
                    @. step = x_trial - x

                    reduce_alpha = false
                    if LA.norm(step) < tol_sgn_switch 
                        if sgn_index == 1
                            d .*= -1
                            continue
                        else
                            reduce_alpha = true
                        end
                    end

                    cache_index_trial = eval_and_cache!(
                        fx, cx, zx_trial, eps_iq, cache, x_trial, Objfs!, Constrs!;
                        initialize_eps_iq=false, check_cache=true, max_func_calls
                    )

                    if cache_index_trial < 1
                        early_exit = true
                        break
                    end
                    
                    offset_alpha = gamma * alpha^2

                    trial_point_is_dominated = filter_strictly_dominates_point_at_index(
                        cache, next_filter_permutation, next_filter_flags, 
                        cache_index_trial, offset_alpha
                    )
                    if trial_point_is_dominated 
                        @logmsg log_level @sprintf("\t Trial point %d is dominated by the filter with offset %.2e.", cache_index_trial, offset_alpha)
                    end

                    if reduce_alpha || trial_point_is_dominated
                        ## failure step
                        @logmsg log_level "\t!Failure step!"
                        for (_i, _ci) in enumerate(next_filter_permutation)
                            !next_filter_flags[_i] && continue
                            if ci == _ci
                                set_stepsize_in_cache!(dir_scheme, sz_cache, alpha/2, ci, dir_index)
                                break
                            end
                        end
                    else
                        ## do line-search
                        while true
                            beta = alpha * 2

                            @. x_beta = x + beta * d
                            project_into_box!(x_beta, lb, ub)
                            @. step_linesearch = x_trial - x_beta

                            if LA.norm(step_linesearch) < tol_sgn_switch 
                                if sgn_index == 1
                                    @logmsg log_level "\tSwitching sign..."
                                    d .*= -1
                                    alpha /= 2
                                else
                                    @logmsg log_level "Failure in other direction?"
                                end
                                break
                            end#if
                            cache_index_beta = eval_and_cache!(
                                fx, cx, zx_beta, eps_iq, cache, x_beta, Objfs!, Constrs!;
                                initialize_eps_iq=false, check_cache=true, max_func_calls
                            )

                            if cache_index_beta < 1
                                early_exit = true
                                break
                            end

                            offset_beta = gamma * beta^2
                            offset_beta_alpha = offset_beta - gamma * alpha^2
                            if any(zx_beta[ℓ] >= zx_trial[ℓ] - offset_beta_alpha for ℓ=eachindex(zx_beta))
                                # add_point x_trial to filter
                                @logmsg log_level "\t* Adding trial point $(cache_index_trial) with α=$(alpha)."
                                update_filter!(next_filter_permutation, next_filter_flags, cache, cache_index_trial; force_add=i)
                                set_stepsize_in_cache!(dir_scheme, sz_cache, alpha, cache_index_trial, dir_index)
                            end
                            if filter_strictly_dominates_point_at_index(
                                cache, next_filter_permutation, next_filter_flags, 
                                cache_index_beta, offset_beta
                            )
                                break # from while-loop doing line search
                            end#if

                            @. x_trial = x_beta
                            @. zx_trial = zx_beta
                            cache_index_trial = cache_index_beta
                            alpha = beta
                        end#while (line search / projected expansion)

                    end#if (failure step or linesearch)
                end#for sgn_index=1:2
                
                early_exit && break
                
                if stepsize_halving
                    old_alpha = get_stepsize_from_cache(dir_scheme, sz_cache, ci, dir_index)
                    if old_alpha > stepsize_stop
                        new_alpha = old_alpha/2
                        @logmsg log_level @sprintf("\tReducing α[%d] from %.2e to %.2e for %d.", dir_index, old_alpha, new_alpha, ci)
                        set_stepsize_in_cache!(dir_scheme, sz_cache, new_alpha, ci, dir_index)
                    end
                end
                    
                dir_index = next_direction!(d, dir_scheme)
            end#while !isnothing(dir_index)
            early_exit && break
            reset_current_state!(dir_scheme)
        end#for element in filter
        early_exit && break
        all_stepsizes_too_small && break
        
        filter_permutation .= next_filter_permutation
        filter_flags .= next_filter_flags
 
        prepare_for_next_iteration!(dir_scheme; log_level)
        if num_points_found == 0
            coef_delta *= 0.95
            @logmsg log_level "Did not use any points in this round, reducing `coef_delta` to $(coef_delta)."
            if coef_delta < 10^round((log10(eps(Float64))))
                break
            end
        else
            @logmsg log_level "Investigated $(num_points_found) points."
        end
    end# for it_index=1:max_iter

    return cache, filter_permutation, filter_flags, dir_scheme, filter_spread_vals
end#function main

function project_into_box!(x, lb, ub)
    for i = eachindex(x)
        x[i] = max(lb[i], min(ub[i], x[i]))
    end
end

@views function filter_strictly_dominates_point_at_index(
    cache, filter_permutation, filter_flags, cache_index, filter_offset=0
)
    point_is_strictly_dominated = false
    fx = cache.fobj[:, cache_index]
    for (i, i_infilter)=enumerate(filter_flags)
        !i_infilter && continue
        ci = filter_permutation[i]
        fx_i = cache.fobj[:, ci]
        if all( fx_i .- filter_offset .< fx )
            point_is_strictly_dominated = true
            break
        end
    end
    return point_is_strictly_dominated
end

function eval_and_cache!(
    fx, cx, zx, eps_iq, cache, x, 
    @nospecialize(Objfs!), @nospecialize(Constrs!);
    initialize_eps_iq::Bool=false,
    max_func_calls::Integer=typemax(Int),
    check_cache::Bool=true
)
    cache_index = if check_cache
        find_in_cache(x, cache)
    else
        -1
    end
    
    cache_index > 0 && return cache_index

    if num_calls(Objfs!) >= max_func_calls || num_calls(Constrs!) >= max_func_calls
        @warn "Maximum number of function calls exhausted."
        return -1
    end

    Objfs!(fx, x)
    Constrs!(cx, x)
    if initialize_eps_iq
        for i=axes(eps_iq, 1)
            eps_iq[i, :] = max(0, cx[i]) < 1 ? 1e-3 : 1e-1
        end
    end
    viol = penalize!(zx, fx, cx, eps_iq)
    ## put result into cache
    cache_index = insert_in_cache!(cache, x, zx, viol; check_cache=false)

    return cache_index  
end

function eval_cache_and_filter!(
    fx, cx, zx, eps_iq, cache, filter_permutation, filter_flags, x, 
    @nospecialize(Objfs!), @nospecialize(Constrs!);
    initialize_eps_iq::Bool=false,
    max_func_calls::Integer=typemax(Int),
    check_cache::Bool=true
)
    cache_index = eval_and_cache!(
        fx, cx, zx, eps_iq, cache, x, Objfs!, Constrs!;
        initialize_eps_iq, max_func_calls, check_cache
    )
    
    ## update the filter
    update_filter!(filter_permutation, filter_flags, cache, cache_index)

    return cache_index  
end

function update_filter!(
    filter_permutation, filter_flags, cache, cache_index;
    force_add::Integer=0
)
    @assert length(filter_permutation) == length(filter_flags)

    fx = @view(cache.fobj[:, cache_index])

    sub_ind = 0
    for (i, is_in_filter)=enumerate(filter_flags)
        if !is_in_filter
            if sub_ind < 1
                sub_ind = i
            end
            continue
        end

        j = filter_permutation[i]

        φx = @view(cache.fobj[:, j])
        if all(φx .<= fx)
            sub_ind = -1
            break
        end

        if all(fx .<= φx) && any( fx .< φx )
            filter_flags[i] = false
            if sub_ind < 1
                sub_ind = i
            end
        end
    end#for
    if sub_ind == 0
        # Point `fx` is not dominated by filter, and no filter point is dominated by `fx`
        # and there is no free spot in the filter
        filter_size = length(filter_flags)
        if iszero(force_add) || filter_size == 0
            @warn "`update_filter!`: Filter is full, cannot add point."
        else
            sub_ind = min(force_add, filter_size)
        end
    end
    if sub_ind > 0
        filter_flags[sub_ind] = true
        filter_permutation[sub_ind] = cache_index
    end
    return sub_ind
end

@views function compute_spread_and_sort!(
    filter_spread_vals, filter_permutation, filter_flags, cache, coef_delta)

    current_size = sum(filter_flags)
    if current_size <= 1
        filter_spread_vals[filter_flags] .= coef_delta
        return nothing
    end

    num_objfs = size(cache.fobj, 1)
    filter_index = filter_permutation[filter_flags]
    sort_ind = collect(eachindex(filter_index))
    
    ## reset spread value array
    T = eltype(filter_spread_vals)
    Tinf = T(Inf)
    filter_spread_vals .= -Tinf
    delta = filter_spread_vals[filter_flags]
    delta .= 0
    for i=1:num_objfs
        Fi = cache.fobj[i, filter_index]     # TODO not optimal in column major...
        sortperm!(sort_ind, Fi; rev=true)
        max_ind = first(sort_ind)
        min_ind = last(sort_ind)
        delta[ max_ind ] += Tinf
        delta[ min_ind ] += Tinf
        divisor_i = Fi[max_ind] - Fi[min_ind]
        for j=2:current_size-1
            delta[sort_ind[j]] += abs( Fi[sort_ind[j-1]] - Fi[sort_ind[j + 1]] ) / divisor_i
        end
        #= 
        for (src_ind, trgt_ind)=enumerate(sort_ind[2:end-1])
            delta[trgt_ind] += abs( Φi[src_ind] - Φi[src_ind+2] ) / divisor_i
        end
        =#
    end
    delta ./= current_size

    sortperm!(sort_ind, delta; rev=true)
    delta .= delta[sort_ind]
    filter_index .= filter_index[sort_ind]

    return nothing
end

function penalize!(zx, fx, cx, eps_iq) # `functs_pen` in Fortran
    T = eltype(zx)
    viol = zero(T)

    zx .= fx
    for i = eachindex(cx)
        viol_i = max(0, cx[i])
        viol += viol_i
        for j = eachindex(fx)
            zx[j] += viol_i / eps_iq[i, j]      # TODO this is a bit slower than if we swapped dimensions, but allows for computation of viol
        end
    end

    return viol
end

end#module