struct Filter
    size :: Int
    flags :: Vector{Bool}
    index :: Vector{Int}
    generations :: Vector{Int}
    parents :: Vector{Int}
    directions :: Vector{Int}
end

function copy_filter!(dest, src)
    dest.flags .= src.flags
    dest.index .= src.index
    dest.generations .= src.generations
    dest.parents .= src.parents
    dest.directions .= src.directions
end

function init_filter(cache, max_filter_size)
    ## initialize filter as a set of two indices
    ## if `filter_flags[i] == true`, then `j=filter_index[i]`
    ## is the index of a point in the filter belonging to `cache`.
    filter_size = max(0, min(max_filter_size, cache.max_store))

    filter_flags = zeros(Bool, filter_size)
    filter_index = fill(-1, filter_size)

    filter_generations = copy(filter_index)
    filter_parents = copy(filter_index)
    filter_directions = copy(filter_index)

    #return Filter(filter_size, filter_flags, tmp_flags, filter_index)
    return Filter(filter_size, filter_flags, filter_index, filter_generations, filter_parents, filter_directions)
end

"""
    update_filter!(filter, cache, ci; force_add=0)

Try to add cached object to filter based on values `cache.fobj[:, ci]`.
`force_add` can be an index to a filter slot. If no other slots are available,
then this slot is overwritten to add the cached objects index.

Return the filter slot index and the number of available free slots.
"""
function update_filter!(
    ## modified
    filter,
    ## not modified
    cache, cache_index;
    force_add::Integer=0
)
    filter_index = filter.index
    filter_flags = filter.flags
    #src tmp_flags = filter.tmp_flags
    #src tmp_flags .= false

    fx = @view(cache.fobj[:, cache_index])

    num_free = 0
    sub_ind = 0
    for (i, is_occupied_slot)=enumerate(filter_flags)
        if !is_occupied_slot
            num_free += 1
            if sub_ind < 1
                sub_ind = i
            end
            continue
        end

        j = filter_index[i]

        φx = @view(cache.fobj[:, j])
        if all(φx .<= fx)
            ## fx is dominated (or in the filter already)
            if force_add > 0
                ## if we want to add it anyways, we have to remove the dominating element
                num_free += 1
                filter_flags[i] = false
                #src tmp_flags[i] = true
                if sub_ind < 1
                    sub_ind = i
                end
                continue
            else
                ## the filter should not include dominated elements, break here
                sub_ind = -1
                break
            end
        end

        if all(fx .<= φx) && any( fx .< φx )
            ## fx dominates an element that is in the filter, so we have to remove it
            filter_flags[i] = false
            #src tmp_flags[i] = true
            num_free += 1
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
            num_free = 1
        end
    end
    if sub_ind > 0
        filter_flags[sub_ind] = true
        filter_index[sub_ind] = cache_index
        #src tmp_flags[sub_ind] = true           # TODO probably redundant
        num_free -= 1
    end
    return sub_ind, num_free
end

function compute_spread!(svals::AbstractVector, F_vals::AbstractMatrix)
    num_objfs, num_vecs = size(F_vals)
    sort_ind = collect(1:num_vecs)
    svals .= 0
    for Fi = eachrow(F_vals)
        sortperm!(sort_ind, Fi; rev=true)
        max_ind = first(sort_ind)
        min_ind = last(sort_ind)
        svals[ max_ind ] += Inf
        svals[ min_ind ] += Inf
        divisor_i = Fi[max_ind] - Fi[min_ind]
        for j=2:num_vecs-1
            svals[sort_ind[j]] += abs( Fi[sort_ind[j-1]] - Fi[sort_ind[j+1]] ) / divisor_i
        end
    end
    svals ./= num_vecs
    return sort_ind
end

function compute_spread!(svals::Nothing, solutions, filter, cache; sort_solutions=true)
    return nothing
end

@views function compute_spread!(svals, solutions, filter, cache; sort_solutions=true)
    sol_flags = solutions .> 0
    sols = solutions[sol_flags]
    F_index = filter.index[sols]
    if length(F_index) <= 1
        svals[sol_flags] .= Inf
        return nothing
    end

    svals .= nextfloat(-Inf)
    
    F_vals = cache.fobj[:, F_index]
    _svals = svals[sol_flags]
    s_ind = compute_spread!(_svals, F_vals)
    if sort_solutions
        sortperm!(s_ind, _svals; rev=true)
        _svals .= _svals[s_ind]
        sols .= sols[s_ind]
    end
    return nothing
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