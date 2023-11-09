const DEFAULT_CACHE_MAX_STORE = 100_000
function DEFAULT_CACHE_TOLERANCE(T::Type{<:AbstractFloat})
    return 10^(-round(abs(log10(eps(T)))^0.65)) # 1e-6 if `T==Float64`
end

Base.@kwdef struct EvalCacheConfig
    max_store :: Int = DEFAULT_CACHE_MAX_STORE
    tolerance :: Union{Real, Nothing} = nothing
end

Base.@kwdef struct EvalCache{T<:AbstractFloat}
    ## These parameters were module parameters in Fortran    
    max_store :: Int = DEFAULT_CACHE_MAX_STORE
    tolerance :: T = DEFAULT_CACHE_TOLERANCE(T)
    ncache :: Base.RefValue{Int} = Ref(0)   # cache access counter

    ## Actual caching matrices
    x :: Matrix{T}
    fobj :: Matrix{T}
    viol :: Vector{T}

    pos_free :: Base.RefValue{Int} = Ref(1)     # next free index in caches
    current_position :: Base.RefValue{Int} = Ref(0)      # maximum index of column holding a result
end

float_type(cache::EvalCache{T}) where{T}=T

Base.convert(::Type{EvalCache{T}}, cache::EvalCache{T}) where T<:AbstractFloat=cache
function Base.convert(::Type{EvalCache{F}}, cache::EvalCache{T}) where {F<:AbstractFloat, T<:AbstractFloat}
    return EvalCache(
        cache.max_store,
        F(cache.tolerance),
        copy(cache.ncache),
        Matrix{F}(cache.x),
        Matrix{F}(cache.fobj),
        Vector{F}(cache.viol),
        copy(cache.pos_free),
        copy(cache.current_position)
    )
end

function init_cache(
    T, num_vars, num_objfs, old_cache, cache_cfg, 
    cache_max_store, cache_tolerance
)
    if isa(old_cache, EvalCache)
        if float_type(old_cache) != T
            @warn "Provided cache has precision different from $(T), converting (and copying!)."
            return convert(EvalCache{T}, old_cache)
        else
            return old_cache
        end
    end
    if isa(cache_cfg, EvalCacheConfig)
        max_store = max(0, cache_cfg.max_store)
        tolerance = isnothing(cache_cfg.tolerance) ? DEFAULT_CACHE_TOLERANCE(T) : T(cache_cfg.tolerance)
    else
        max_store = max(0, cache_max_store)
        tolerance = isnothing(cache_tolerance) ? DEFAULT_CACHE_TOLERANCE(T) : T(cache_tolerance)
    end
    x = Matrix{T}(undef, num_vars, max_store)
    fobj = Matrix{T}(undef, num_objfs, max_store)
    viol = zeros(T, max_store)
    return EvalCache{T}(; tolerance, max_store, x, fobj, viol)
end

function find_in_cache(x_query, cache)
    x_query_index = 0
    cache.ncache[] += 1
    eps = cache.tolerance
    for (j, x_c)=enumerate(eachcol(cache.x))
        if j > cache.current_position[]
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
        if cache.current_position[] < cache.max_store
            cache.current_position[] += 1
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
    cache, filter, 
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

    return cache_index, filter_slot_index, filter_num_free
end