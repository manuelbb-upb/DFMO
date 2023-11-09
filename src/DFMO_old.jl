module DFMO

using Primes: nextprime
import LinearAlgebra as LA

struct FilterSlot{T<:AbstractFloat}
    x :: Vector{T}
    f :: Vector{T}
    nondom_flag :: Vector{Bool}
    alph_d :: Vector{T}
    alphs_orth :: Vector{T}
end

struct Filter{T<:AbstractFloat}
    size :: Int                     # `ndim` in Fortran
    X :: Matrix{T}                  # columns `1:n` in Fortran matrices
    F :: Matrix{T}                  # columns `n+1:n+q` in Fortran matrices
    nondom_flags :: Vector{Bool}    # encoded by column `n+q+1`
    alph_d :: Vector{T}            # stepsize for dense direction, column `n+q+2`
    alphs_orth :: Matrix{T}        # stepsizes for orthogonal directions, columns `n+q+2+1:n+q+2+n`
end

function sort_with_index!(filter, index)
    @assert length(index) == filter.size
    filter.X .= filter.X[:, index]
    filter.F .= filter.F[:, index]
    filter.nondom_flags .= filter.nondom_flags[index]
    filter.alph_d .= filter.alph_d[index]
    filter.alphs_orth .= filter.alphs_orth[:, index]
    return nothing
end

function sort_by_spread!(delta, filter::Filter{T}) where T
    #delta = zeros(T, filter.size)

    #=
    if !any(nondom_flags)
        return
    end
    =#
    ncomp = sum(filter.nondom_flags)
    if ncomp <= 1
        return false
    end
    
    num_objfs = size(filter.F, 1)
    for i = 1:num_objfs
        Fi = @view(filter.F[i, filter.nondom_flags])
        ipt = sortperm(Fi; rev=true)

        Φi = @view(Fi[ipt])
        _delta = @view(delta[filter.nondom_flags])
        _delta[first(ipt)] += T(Inf)
        _delta[last(ipt)] += T(Inf)

        divisor_i = first(Φi) - last(Φi)
        for j=2:length(ipt)-1
            _delta[ipt[j]] += abs( Φi[j-1] - Φi[j+1] ) / divisor_i
        end
    end

    delta ./= ncomp
    for j = 1:filter.size
        if !filter.nondom_flags[j]
            delta[j] = -100
        end
    end

    ipt = sortperm(delta; rev=true)
    delta .= delta[ipt]
    sort_with_index!(filter, ipt)

    return true
end

function Filter(
    num_vars::Int, num_objfs::Int; size::Int=100_000, T::Type{<:AbstractFloat}=Float64
)
    X = Matrix{T}(undef, num_vars, size)
    F = Matrix{T}(undef, num_objfs, size)
    nondom_flags = zeros(Bool, size)
    alph_d = Vector{T}(undef, size)
    alphs_orth = Matrix{T}(undef, num_vars, size)
    return Filter(size, X, F, nondom_flags, alph_d, alphs_orth)
end

function slot(filter, j=1)
    x = copy( filter.X[:, j] )
    f = copy( filter.F[:, j] )
    nondom_flag = [filter.nondom_flags[j], ]
    alph_d = [filter.alph_d[j],]
    alphs_orth = filter.alph_d[:, j]
    return FilterSlot(x, f, nondom_flag, alph_d, alphs_orth)
end

function overwrite_slot_with_other!(filter, trgt, src)
    filter.X[:, trgt] .= filter.X[:, src]   
    filter.F[:, trgt] .= filter.F[:, src]   
    filter.alphs_orth[:, trgt] .= filter.alphs_orth[:, src]   
    filter.alph_d[trgt] = filter.alph_d[src]   
    filter.nondom_flags[trgt] = filter.nondom_flags[src]   
    return nothing
end

function extract_slot!(dummy, filter, j)
    dummy.x .= filter.X[:, j]
    dummy.f .= filter.F[:, j]
    dummy.nondom_flag[end] = filter.nondom_flag[j]
    dummy.alph_d[end] = filter.alph_d[j]
    dummy.alphs_orth .= filter.alphs_orth[:, j]
    return nothing
end

function push_slot!(filter, dummy, j=1)
    filter.X[:, j] .= dummy.x
    filter.F[:, j] .= dummy.f
    filter.nondom_flags[j] = dummy.nondom_flag[end]
    filter.nondom_alph_d[j] = dummy.alph_d[end]
    filter.nondom_alphs_orth[:, j] .= dummy.alphs_orth
    return nothing
end

function update_filter!(filter, x, fx, alph_d, alphs_orth)
    fx_pos = 0
    for j=1:filter.size
        if !filter.nondom_flags[j]
            if fx_pos < 1
                fx_pos = j
            end
            continue
        end
        φx = filter.F[:, j]
        if all(φx .<= fx)
            fx_pos = -1
            break
        end

        if all(fx .<= φx) && any( fx .< φx )
           filter.nondom_flags[j] = 0
           if fx_pos < 1 fx_pos = j end
        end
    end#for

    if fx_pos == 0
        # Point `fx` is not dominated by filter, and no filter point is dominated by `fx`
        # and there is no free spot in the filter
        # (or filter has zero size)
        @warn "`update_filter!`: Filter is full, cannot add point."
    elseif fx_pos > 0
        filter.nondom_flags[fx_pos] = true
        filter.X[:, fx_pos] .= x
        filter.F[:, fx_pos] .= fx
        filter.alph_d[fx_pos] = alph_d
        filter.alphs_orth[:, fx_pos] .= alphs_orth
        return true
    end
    return false
end

function is_outside_of_filter(fx, filter; delta=1e-2) # `intorno`
    is_outside = true
    for j=1:filter.size
        !filter.nondom_flags[j] && continue
        φx = filter.F[:, j]
        if maximum( abs.(fx .- φx) ) < min(1, delta)
            is_outside = false
            break
        end
    end
    return is_outside
end

Base.@kwdef struct EvalCache{T<:AbstractFloat}
    ## These parameters were module parameters in Fortran
    maxstore :: Int = 100_000
    tolerance :: T = 10^(-ceil(abs(log10(eps(T)))^0.65)) # 1e-6 for T==Float64
    ncache :: Base.RefValue{Int} = Ref(0)   # cache counter

    ## Actual caching matrices
    x :: Matrix{T}
    fobj :: Matrix{T}
    fconstr :: Matrix{T}

    pos_free :: Base.RefValue{Int} = Ref(1)     # next free index in caches
    current_position :: Base.RefValue{Int} = Ref(0)      # maximum index of column holding a result
end

function setup_cache(num_vars, num_objfs, num_constrs; max_store::Int=100_000, T::Type{<:AbstractFloat}=Float64)
    x = Matrix{T}(undef, num_vars, max_store)
    fobj = Matrix{T}(undef, num_objfs, max_store)
    fconstr = Matrix{T}(undef, num_constrs, max_store)
    return EvalCache(; max_store, x, fobj, fconstr)
end

"""
    find_in_cache!(fx, cx, cache, x_query)

Look for objective and constraint value vectors corresponding to 
`x_query`. If `cache` contains suitable values (up to a Inf-norm tolerance
of `cache.tolerance` in decision space), then `fx` and `cx` are modified in-place
to hold those values.
"""
function find_in_cache!(fobj, fconstr, cache, x_query)
    x_query_is_in_cache = false
    cache.ncache[] += 1
    eps = cache.tolerance
    for (j, x_c)=enumerate(eachcol(cache.x))
        if all( abs.( x_query .- x_c ) .<= eps )
            x_query_is_in_cache = true
            fobj .= cache.fobj[:, j]
            fconstr .= cache.fconstr[:, j]
        end
    end
    return x_query_is_in_cache
end

function insert_in_cache!(cache, x, fobj, fconstr; check_cache=true)
    if !check_cache || !find_in_cache!(fobj, fconstr, cache, x)
        j = cache.pos_free[]
        cache.x[:, j] .= x
        cache.fobj[:, j] .= fobj
        cache.fconstr[:, j] .= fconstr
        j = j < cache.max_store ? j+1 : 1
        cache.pos_free[] = j
        if cache.current_position[] < cache.max_store
            cache.current_position[] += 1
        end
    else
        @warn "`insert_in_cache!`: Point is in cache already."
    end
    return nothing
end

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

"Fill columns of `points` with samples along the diagonal of box 
with bottom left corner `lb` and diagonal vector `box_width`."
function diagonal!(points, lb, box_width) # `diagonale` in Fortran
    num_vars = length(lb)
    for j=eachindex(lb)
        points[:, j] = lb + ((i-1)/(num_vars-1)) * box_width
    end
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
            zx[j] += viol_i / eps_iq[i, j]
        end
    end

    return viol
end

function main(x0, lb, ub, num_objfs, objfs!, num_constrs=0, constrs!=nothing; T::Type{<:AbstractFloat}=Float64)
    num_vars = length(x0)
    box_width = ub .- lb

    @assert num_vars > 0 "Number of variables must at least be 1."
    @assert length(box_width) == num_vars
    @assert num_objfs > 0 "Number of objective functions must at least be 1."
    @assert isnothing(constrs!) || num_constrs > 0 "Constraint function provided but `num_constrs==0`."

    if any( lb .> x0 ) || any(ub .< x0 )
        return
    end 

    x = T.(x0)
    fx = zeros(T, num_objfs)    # `fob` in Fortran 
    cx = zeros(T, num_constrs)  # `ciq` in Fortran
    zx = deepcopy(fx)           # `fpen` or `finiz` in Fortran

    cache = setup_cache(num_vars, num_objfs, num_constrs; T)

    L_kappa = Filter(; num_vars, num_objfs, T)
    L_new = Filter(; num_vars, num_objfs, T)
    L_tilde = Filter(; num_vars, num_objfs, T)
 
    alpha_d_init = min(10, maximum(box_width/10)) |> T
    alphas_orth_init = [ T(min(10), box_width[i]/10) for i=1:num_vars ]

    Objfs! = CountedFunction(objfs!)
    Constrs! = CountedFunction(constrs!)
   
    if !find_in_cache!(fx, cx, cache, x)    
        Objfs!(fx, x)
        Constrs!(cx, x)
        insert_in_cache!(cache, x, fx, cx; check_cache=false)
    end
   
    eps_iq = Matrix{T}(undef, num_constrs, num_objfs)
    for i=1:num_constrs
        eps_iq[i, :] = max(0, cx[i]) < 1 ? 1e-3 : 1e-1
    end

    #=
    eps_iq = Matrix{T}(undef, num_objfs, num_constrs)
    for k=1:num_objfs
        for i=1:num_constrs
            eps_iq[k, i] = max(0, cx[i]) < 1 ? 1e-3 : 1e-1
        end
    end
    =#

    # NOTE in the Fortran version, here is another call to
    # `startp` and the functions are evaluated again (lines 140 ff) 
    points = Matrix{T}(undef, num_vars, num_vars)
    diagonal!(points, lb, box_width)

    viol = penalize!(zx, fx, cx, eps_iq)
    filter_changed = update_filter!(L_kappa, x, zx, alpha_d_init, alphas_orth_init)

    for j=1:num_vars
        x .= points[:, j]
        if !find_in_cache!(fx, cx, cache, x)    
            Objfs!(fx, x)
            Constrs!(cx, x)
            insert_in_cache!(cache, x, fx, cx; check_cache=false)
        end
        viol = penalize!(zx, fx, cx, eps_iq)
        if is_outside_of_filter(zx, L_kappa) || is_outside_of_filter(zx, L_new) || is_outside_of_filter(zx, L_tilde) 
            update_filter!(L_kappa, x, zx, alpha_d_init, alphas_orth_init)
        end
    end

    eps_iq_init = copy(eps_iq)
    zx_init = copy(fx) # TODO check why it is not `zx`
    viol_init = viol

    alpha_stop = 10 ^ floor(abs(log10(eps(T)))^0.8) # 1e-9 for T==Float64
    max_func_calls = 2_000

    eta = 10^(-ceil(abs(log10(eps(T)))^0.65)) # 1e-6 for T==Float64
    gamma = 10^(-ceil(abs(log10(eps(T)))^0.65))

    d_dense = similar(x)
    index_halton = 1_000 + 2 * num_vars

    if num_vars <= 1
        d_dense[end] = 1
    else
        halton!(d_dense, index_halton)
    end

    _H, _ = LA.qr(d_dense)
    H = Matrix(_H[:,:])

    ## Sort entries in L_kappa according to objective values. If `q == num_objfs`,
    ## then afterwards the first entry will have lowest value w.r.t. q-th objective.
    ## The next entry has lowest value with respect to objective `q-1` etc. ...
    dummy = slot(L_kappa)
    for i=1:num_objfs
        fbest = T(Inf)
        ibest = 0
        for pn=1:L_kappa.size
            !L_kappa.nondom_flags[pn] && continue
            φ = L_kappa.F[i, pn]    # TODO currently slow in column-major arrays
            if fbest > φ
                fbest = φ
                ibest = pn
            end
        end
        ## Shift entries 2:ibest one position to the right and put ibest first
        extract_slot!(dummy, L_kappa, ibest)
        for pn=ibest:-1:2
            overwrite_slot_with_other!(L_kappa, pn, pn-1)
        end
        push_slot!(L_kappa, dummy, 1)
    end

    # sort by spread metric/crowding distance
    # TODO why sort twice??
    deltaf = zeros(T, L_kappa.size)
    spread_flag = sort_by_spread!(deltaf, L_kappa)  # `spread_ord` in Fortran

    L_maxalpha = L_maxalpha2 = zero(T)
    for pn=1:L_kappa.size
        !L_kappa.nondom_flags[pn] && continue
        if L_maxalpha < L_kappa.alph_d[pn]
            L_maxalpha = L_kappa.alph_d[pn]
        end
        if L_maxalpha2 < L_kappa.alph_d[pn]
            L_maxalpha2 = L_kappa.alph_d[pn]
        end
        for j = 1:num_vars
            if L_maxalpha2 < L_kappa.alphs_orth[j, pn]
                L_maxalpha2 = L_kappa.alphs_orth[j, pn]
            end
        end
    end

    # Refine L_kappa (the filter)

    for pn=1:L_kappa.size
        if L_kappa.nondom_flags[pn]
            # compute constraint violation value
            viol = zero(T)
            if num_objfs >= 1 
                x = L_kappa.X[:, pn]
                if !find_in_cache!(fx, cx, cache, x)
                    Constrs!(cx, x)                 # Todo: why not store these values in filter too?
                end
                viol = max(0, maximum(cx))
            end
            alpha_max = zero(T)
            for i=1:num_vars
                if alpha_max < L_kappa.alphs_orth[i, pn]
                    alpha_max = L_kappa.alphs_orth[i, pn]
                end
            end#for
            if alpha_max < L_kappa.alph_d[pn]
                alpha_max = L_kappa.alph_d[pn]
            end
        end

        condizione = (
            L_kappa.nondom_flags[pn] &&
            L_kappa.alph_d[pn] > alpha_stop &&
            (!spread_flag || deltaf[pn] >= coef_delta) &&
            (L_kappa.alph_d[pn] > L_maxalpha / 10 || pn <= q)       # TODO why this condition??
        )

        if condizione
            x = L_kappa.X[:, pn]        # TODO `.=` to change arrays instead?
            f = L_kappa.F[:, pn]

            alphas_orth = L_kappa.alphs_orth[:, pn]
            alhpa_d = L_kappa.alph_d[pn]

            z = copy(x)   # TODO why? # TODO preallocate z
            dz = zero(z)

            if num_vars==1
                @warn "`dir_coord` not implemented..."
            elseif num_vars >= 1
                for ii=1:num_vars
                    d_dense = H[:, ii]
                    for _=1:2
                        z .= x
                        dz .= alhpa_d .* d_dense
                        z .+= dz
                        z .= max.(lb, min.(ub, z))

                        if LA.norm(dz) <= 10^(-ceil(abs(log10(eps(T)))^0.75)) # 1e-8 if T==Float64
                            d_dense .*= -1
                            continue
                        end
                        if !find_in_cache!(fx, cx, cache, z)
                            Objfs!(fx, z)
                            Constrs!(cx, z)
                            insert_in_cache!(cache, z, fx, cx; check_cache=false)
                        end
                        viol = penalize!(zx, fx, cx, eps_iq)
                    end
                end#for ii
            end#dir_coord vs dir_dense
        end#condizione
    end

    # for each slot in filter
    #  if it is nondominated and if there are constraints
    #    try to find constraint values in cache or compute them
    #    compute maximum constraint violation
    #    compute largest alpha
    #  if spreading distance has been calculated and Δ(slot) == 0, print it
    # `condizione`: slot is nondominated and stepsize_d is > alfa_stop and
    #               spread is >= coef_delta or not computed and
    #               stepsize > 0.1*`L_maxalpha` || pn <= q ??????
    #  read stepsizese, set `trovate=true`
end

function halton!(dir, ih)
    num_vars = length(dir)
    if num_vars <= 1
        dir[end] = 1
        return nothing
    end
    base = 2
    for i=1:num_vars
        base = i==1 ? base : nextprime(base, 2)
        dir[i] = 0
        f = one(eltype(dir)) / base
        j = ih
        while j > 0
            j, r = divrem(j, base)
            dir[i] += r * f
            f /= base
        end
    end
    return 
end

end