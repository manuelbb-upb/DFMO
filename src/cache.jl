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

float_type(cache::EvalCache{T}) where{T}=T

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
    ## modified:
    fx, cx, zx, eps_iq, 
    cache, filter, solutions,
    ## not modified:
    x, 
    @nospecialize(Objfs!), @nospecialize(Constrs!);
    ## kwargs:
    initialize_eps_iq::Bool=false,
    max_func_calls::Integer=typemax(Int),
    check_cache::Bool=true
)
    cache_index = eval_and_cache!(
        fx, cx, zx, eps_iq, cache, x, Objfs!, Constrs!;
        initialize_eps_iq, max_func_calls, check_cache
    )
    
    ## update the filter
    filter_slot_index, filter_num_free = update_filter!(filter, cache, cache_index)
    sync!(solutions, filter)

    return cache_index, filter_slot_index, filter_num_free
end