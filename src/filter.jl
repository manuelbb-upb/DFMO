struct Filter
    size :: Int
    flags :: Vector{Bool}
    tmp_flags :: Vector{Bool}
    index :: Vector{Int}
end

function init_filter(cache, max_filter_size)
    ## initialize filter as a set of two indices
    ## if `filter_flags[i] == true`, then `j=filter_index[i]`
    ## is the index of a point in the filter belonging to `cache`.
    filter_size = max(0, min(max_filter_size, cache.max_store))

    filter_flags = zeros(Bool, filter_size)
    tmp_flags = deepcopy(filter_flags)
    filter_index = fill(-1, filter_size)

    return Filter(filter_size, filter_flags, tmp_flags, filter_index)
end

struct Solutions{S}
    flags :: Vector{Bool}
    index :: Vector{Int}
    spread_vals :: S
end

function init_solutions(filter, max_num_sols, has_spread, T)
    num_sols = max(0, min(filter.size, max_num_sols))
    sols_flags = zeros(Bool, num_sols)
    sols_index = fill(-1, num_sols)
    spread_vals = has_spread ? Vector{T}(undef, num_sols) : nothing
    return Solutions(sols_flags, sols_index, spread_vals)
end

function set_flag!(solutions, filter_slot_index; force_add=-1, soft_force=true)
    filter_slot_index < 1 && return filter_slot_index
    forced_slot_index = min(force_add, length(solutions.flags))
    solutions_slot_index = soft_force ? -1 : forced_slot_index
    if solutions_slot_index < 0
        for (i, is_occupied_sol_slot) = enumerate(solutions.flags)
            if is_occupied_sol_slot 
                if solutions.index[i] == filter_slot_index
                    return i
                else
                    continue
                end
            end
            solutions_slot_index = i
            break
        end
    end
    if solutions_slot_index < 0
        solutions_slot_index = forced_slot_index
    end
    solutions_slot_index < 1 && return solutions_slot_index

    solutions.flags[solutions_slot_index] = true
    solutions.index[solutions_slot_index] = filter_slot_index
    svals = solutions.spread_vals
    if !isnothing(svals) 
        svals[solutions_slot_index] = NaN
    end
    return solutions_slot_index
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
    tmp_flags = filter.tmp_flags
    tmp_flags .= false

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
                tmp_flags[i] = true
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
            tmp_flags[i] = true
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
        tmp_flags[sub_ind] = true           # TODO probably redundant
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

@views function compute_spread!(solutions, filter, cache; sort_solutions=true)
    isnothing(solutions.spread_vals) && return nothing

    F_index = filter.index[solutions.index[solutions.flags]]

    if length(F_index) <= 1
        solutions.spread_vals[solutions.flags] .= Inf
        return nothing
    end

    solutions.spread_vals .= -Inf
    F_vals = cache.fobj[:, F_index]
    svals = solutions.spread_vals[solutions.flags]
    s_ind = compute_spread!(svals, F_vals)
    if sort_solutions
        sort_by_spread!(solutions, s_ind)
    end
    return nothing
end

@views function sort_by_spread!(solutions, sort_ind = nothing)
    isnothing(solutions.spread_vals) && return nothing
    if isnothing(sort_ind)
        sort_ind = sortperm(solutions.spread_vals[solutions.flags]; rev=true)
    else
        sortperm!(sort_ind, solutions.spread_vals[solutions.flags]; rev=true)
    end
    solutions.index[solutions.flags] .= solutions.index[sort_ind]
    solutions.spread_vals[solutions.flags] .= solutions.spread_vals[sort_ind]
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

function sync!(solutions, filter; set_new_solutions=false)
    for (si, sol_slot_is_occupied) = enumerate(solutions.flags)
        !sol_slot_is_occupied && continue
        fi = solutions.index[si]
        !filter.tmp_flags[fi] && continue
        solutions.flags[si] = false
        filter.tmp_flags[fi] = false
    end
    if set_new_solutions
        for (fi, filter_slot_has_changed) = enumerate(filter.tmp_flags)
            !filter_slot_has_changed && continue
            set_flag!(solutions, fi)
        end
    end
end