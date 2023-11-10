abstract type AbstractDirectionScheme end
abstract type AbstractDirectionSchemeConfig end

"""
    init_direction_scheme(
        cfg::AbstractDirectionSchemeConfig, x, zx, lb, ub, solutions, filter)

Return an object of type `AbstractDirectionScheme` corresponding to configuration `cfg`.
"""
function init_direction_scheme(
    ::AbstractDirectionScheme, x, zx, lb, ub, solutions, filter)::AbstractDirectionScheme
    return nothing
end

"Return the maximum stepsize currently associated with any solution and direction."
max_stepsize(::AbstractDirectionScheme, solutions)=missing

function propagate_solutions!(
    hd::DS, arrays, cache, solutions_int_indices, previous_solutions_int_indices, filter, 
    Objfs!, Constrs!,
    spread_vals, lb, ub, log_level, stepsize_stop, max_func_calls, it_index
) where{DS <: AbstractDirectionScheme}
    return error("`propagate_solutions!` not implemented for $(DS).")
end

"""
    HaltonConfig(;
        tol_sgn_switch::Union{Real, Nothing}=nothing,
        coef_delta::Union{Real, Nothing}=nothing,
        do_subiterations::Bool=true
    )

Return a configuration object of type `HaltonConfig <: AbstractDirectionSchemeConfig`
that is used to initialize an `AbstractDirectionScheme` of type `HaltonDirections`.
"""
Base.@kwdef struct HaltonConfig <: AbstractDirectionSchemeConfig
    tol_sgn_switch :: Union{Real, Nothing} = nothing
    coef_delta :: Union{Real, Nothing} = nothing
    gamma :: Union{Real, Nothing} = nothing
    do_subiterations :: Bool = false
    stepsize_halving :: Bool = true
end

struct HaltonDirections{T<:AbstractFloat} <: AbstractDirectionScheme 
    num_vars :: Int
    sz_init :: T
    sz_vals :: Vector{T}
    index :: Base.RefValue{Int}
    orth_dirs :: Matrix{T}

    x_trial :: Vector{T}
    x_beta :: Vector{T}
    zx_trial :: Vector{T}
    zx_beta :: Vector{T}
    step_trial :: Vector{T}
    step_beta :: Vector{T}

    tol_sgn_switch :: T
    coef_delta :: Base.RefValue{T}
    gamma :: T

    do_subiterations :: Bool
    stepsize_halving :: Bool 

    next_filter :: Filter
end

function init_direction_scheme(
    cfg::HaltonConfig, x, zx, lb, ub, solutions, filter
)
    num_vars = length(x)
    T = eltype(x)

    ih = 1_000 + 2 * num_vars
    d0 = halton(num_vars, ih, T)
    orth_dirs = Matrix{T}(undef, num_vars, num_vars)
    
    halton_base!(orth_dirs, d0)
    box_width = ub .- lb
    sz_init = min(10, maximum(box_width ./ 10))

    #src sz_vals = fill(sz_init, num_sols)
    sz_vals = fill(sz_init, filter.size)

    x_trial = copy(x)
    x_beta = copy(x)
    zx_trial = copy(zx)
    zx_beta = copy(zx)
    step_trial = zero(x)
    step_beta = copy(step_trial)

    tol_sgn_switch = if isnothing(cfg.tol_sgn_switch)
        10^(-round(abs(log10(eps(T)))^0.76)) # 1e-8 for T == Float64
    else
        T(cfg.tol_sgn_switch)
    end

    gamma = if isnothing(cfg.gamma) || cfg.gamma <= 0
        10^(-round(abs(log10(eps(T)))^0.65)) # 1e-6 for T==Float64
    else
        T(cfg.gamma)
    end

    coef_delta_ref = if isnothing(cfg.coef_delta)
        Ref(one(T))
    else
        Ref(T(cfg.coef_delta))
    end

    return HaltonDirections{T}(
        num_vars, sz_init, sz_vals, Ref(ih), orth_dirs,
        x_trial, x_beta, zx_trial, zx_beta, step_trial, step_beta,
        tol_sgn_switch, coef_delta_ref, gamma, cfg.do_subiterations, cfg.stepsize_halving,
        deepcopy(filter),
        #src deepcopy(solutions)
    )
end

@views function max_stepsize(hd::HaltonDirections, solutions)
    return maximum(hd.sz_vals[si] for si in solutions if si > 0; init=eltype(hd.sz_vals)(-Inf))
end

function prepare_for_next_iteration!(hd::HaltonDirections{T}; log_level=Info) where T
    ih = hd.index[] += 2*hd.num_vars
    @logmsg log_level "\tComputing halton direction for halton index $(ih)."
    d0 = halton(hd.num_vars, ih, T)
    halton_base!(hd.orth_dirs, d0)
    return nothing
end

function preselect_solutions!(
    hd::HaltonDirections,
    next_solution_int_indices, previous_solution_int_indices, filter, it_index,
    log_level;
    prefer_large_stepsize::Bool=false
)
    @logmsg log_level "##### SOLUTION SELECTION #####"
    max_num_sols = length(next_solution_int_indices)
    next_solution_int_indices .= -1
    if it_index == 1
        num_sols = 0
        for (fi, fi_slot_occupied) = enumerate(filter.flags)
            !fi_slot_occupied && continue
            num_sols == max_num_sols && break
            num_sols += 1
            next_solution_int_indices[num_sols] = fi
        end
        @logmsg log_level "\tSOLUTION SELECTION: Chose $(num_sols) first elements in filter: $([next_solution_int_indices[1:num_sols]])."
    else
        ## 0) a bit of cheeting... if the filter has few enough elements, we select all
        
        num_filter = sum(filter.flags)
        if num_filter <= max_num_sols
            next_solution_int_indices[1:num_filter] .= findall(filter.flags)
            return nothing
        end
        
        ## 1) for each previous solution, try to find a child solution first
        ##    if there are multiple elements in the filter, that where generated from a 
        #     previous solution, then choose that with the largest stepsize
        sol_pos = 0
        last_gen = it_index - 1
        for pfi in previous_solution_int_indices
            pfi <  1 && continue
            pfi_child = -1
            pfi_child_sz = prefer_large_stepsize ? -Inf : Inf
            pfi_is_in_filter = false
            for (fi, fi_slot_occupied) = enumerate(filter.flags)
                !fi_slot_occupied && continue
                if pfi == fi
                    pfi_is_in_filter = true
                end
                filter.parents[fi] != pfi && continue
                filter.generations[fi] != last_gen && continue
                fi_sz = hd.sz_vals[fi]
                if (prefer_large_stepsize && fi_sz > pfi_child_sz) ||
                    (!prefer_large_stepsize && fi_sz < pfi_child_sz)
                    pfi_child_sz = fi_sz
                    pfi_child = fi
                end#if
            end#for (fi, fi_slot_occupied) = enumerate(filter.flags)
            if pfi_child > 0
                sol_pos += 1
                next_solution_int_indices[sol_pos] = pfi_child
                @logmsg log_level "\tSOLUTION SELECTION [$(sol_pos)]: Previous solution $(pfi) has filter descendant $(pfi_child)."
            else
                @logmsg log_level "\tSOLUTION SELECTION: Previous solution $(pfi) has no descendant in filter."
                if pfi_is_in_filter && !(pfi in next_solution_int_indices[1:sol_pos+1])
                    sol_pos += 1
                    next_solution_int_indices[sol_pos] = pfi
                    @logmsg log_level "\tSOLUTION SELECTION [$(sol_pos)]: Previous solution $(pfi) is still in filter."
                end
            end
        end#fi in previous_solution_int_indices

        @logmsg log_level "\tSOLUTION SELECTION: Found $(sol_pos) (descendant) solution(s)."
        ## 2)  if there are still unused solution slots, try to find child solutions 
        ##     stemming from different directions
        num_sols = sol_pos
        found_more_elements = num_sols < max_num_sols
        while found_more_elements
            found_more_elements = false
            for si in next_solution_int_indices[1:sol_pos]
                si_parent = filter.parents[si]
                new_si = -1
                new_si_sz = prefer_large_stepsize ? -Inf : Inf
                found_new_direction_sol = false
                for (fi, fi_slot_occupied) = enumerate(filter.flags)
                    !fi_slot_occupied && continue
                    filter.parents[fi] != si_parent && continue
                    filter.generations[fi] != last_gen && continue
                    fi_direction = filter.directions[fi]
                    fi_already_a_solution = false
                    fi_has_new_direction = true
                    for oi in next_solution_int_indices[1:num_sols]
                        if fi == oi
                            fi_already_a_solution = true
                            break
                        end
                        oi_direction = filter.directions[oi]
                        if fi_direction == oi_direction
                            fi_has_new_direction = false
                        end
                    end#for oi in next_solution_int_indices
                    fi_already_a_solution && continue
                    found_more_elements = true
                    fi_sz = hd.sz_vals[fi]
                    found_new_direction_sol && !fi_has_new_direction && continue
                    if fi_has_new_direction && !found_new_direction_sol
                        new_si = fi
                        new_si_sz = fi_sz
                        found_new_direction_sol = true
                        continue
                    end
                    if (prefer_large_stepsize && fi_sz > new_si_sz) ||
                        (!prefer_large_stepsize && fi_sz < new_si_sz)
                        new_si_sz = fi_sz
                        new_si = fi
                    end
                end#for (fi, fi_slot_occupied) = enumerate(filter.flags)
                if new_si > 1
                    num_sols += 1
                    next_solution_int_indices[num_sols] = new_si
                    @logmsg log_level "\tSOLUTION SELECTION [$(num_sols)]: Previous solution $(si_parent) has filter descendant $(new_si)."
                end
                num_sols == max_num_sols && break
            end#for si in next_solution_int_indices[1: sol_pos]
            num_sols == max_num_sols && break
        end#while

        for si in previous_solution_int_indices
            num_sols >= max_num_sols && break
            for (fi, fi_slot_occupied) = enumerate(filter.flags)
                !fi_slot_occupied && continue
                if fi == si && !(si in next_solution_int_indices[1:num_sols])
                    num_sols += 1
                    next_solution_int_indices[num_sols] = si
                    @logmsg log_level "\tSOLUTION SELECTION [$(num_sols)]: Previous solution $(si) is still in filter."
                    break
                end
            end
        end

        @logmsg log_level "\tSOLUTION SELECTION: Chose $(num_sols) elements from filter: $([next_solution_int_indices[1:num_sols]])."
    end#if else

    return nothing
end

function propagate_solutions!(
    hd::HaltonDirections, arrays, cache, solutions_int_indices, previous_solutions_int_indices, filter, 
    Objfs!, Constrs!,
    spread_vals, lb, ub, log_level, stepsize_stop, max_func_calls, it_index
)

    @unpack x, cx, fx, zx, eps_iq, num_vars, num_objfs = arrays
    @unpack x_trial, x_beta, zx_trial, zx_beta, step_trial, step_beta, tol_sgn_switch, 
        next_filter, gamma = hd;

    preselect_solutions!(hd, solutions_int_indices, previous_solutions_int_indices, filter, it_index, log_level)
    #src compute_spread!(spread_vals, solutions_int_indices, filter, cache; sort_solutions=true)

    @logmsg log_level "#### PROPAGATION" 

    sz_max = max_stepsize(hd, solutions_int_indices)
    @logmsg log_level "  Maximum stepsize is $(sz_max)." 

    if sz_max < stepsize_stop
        return STOP_MIN_STEPSIZE
    end

    copy_filter!(next_filter, filter)
    it_code = CONTINUE_ITERATION
    coef_delta = hd.coef_delta[]
   
    num_points_found = 0


    for (i, si) in enumerate(solutions_int_indices)
        si < 1 && continue
        ci = filter.index[si]
        @logmsg log_level "  Considering (num, filter, cache) $((i, si, ci))." 
        
        alpha_start = alpha = hd.sz_vals[si]
        ## inverting `condizione` (DFMO.f90, ll. 184):
        ## a) consider only filter elements with stepsize big enough
        if alpha <= stepsize_stop
            @logmsg log_level "\t\tStepsize is $(alpha) <= $(stepsize_stop), SKIPPING."
            continue
        end
        #=
        ## b) favor larger spread values
        if !isnothing(spread_vals) 
            spread = spread_vals[i]
            if spread < coef_delta
                @logmsg log_level "\t\tSpread value is $(spread) < $(coef_delta), SKIPPING."
                continue
            end
        else
            spread = nothing
        end
        =#
        spread = nothing
        ## c) if an element is not filter-extremal, it has to have large *relative* stepsize
        alpha <= sz_max/10 && i > num_objfs && continue

        num_points_found += 1

        ci_x = ci
        si_x = si
        si_filter_slot_overwritten = false
        x .= cache.x[:, ci_x]
        zx .= cache.fobj[:, ci_x]
        #@logmsg log_level "\tSpread is $(spread) >= $(coef_delta), vectors are\n$(CustomPrinting.pretty_str(CustomPrinting.VectorTable("x"=>x, "zx" => zx); line_prefix="\t\t"))"
        @logmsg log_level "\tVectors are\n$(CustomPrinting.pretty_str(CustomPrinting.VectorTable("x"=>x, "zx" => zx); line_prefix="\t\t"))"

        any_dir_success = false
        for (dir_index, d) = enumerate(eachcol(hd.orth_dirs))
            x .= cache.x[:, ci_x]
            ci_trial = ci_x
            si_trial = si_x
            #src alpha_start = alpha = hd.sz_vals[i_x]
            alpha = alpha_start
            num_dir_sols = 0
            for sgn in (1, -1)
                num_dir_sols > 0 && continue
                @logmsg log_level "\tLineSearch (filter, cache, direction) $((si_x, ci_x, dir_index)), α=$(sgn*alpha)."

                @. x_trial = x + sgn * alpha * d
                project_into_box!(x_trial, lb, ub)
                @. step_trial = x_trial - x
                if LA.norm(step_trial) < tol_sgn_switch
                    continue
                end

                ci_trial = eval_and_cache!(
                    fx, cx, zx_trial, eps_iq, cache, x_trial, Objfs!, Constrs!;
                    initialize_eps_iq=false, check_cache=true, max_func_calls
                )
                if ci_trial < 1 
                    it_code = STOP_MAX_FUNC_CALLS
                    break#for sgn in (1, -1)
                end

                offset_alpha = gamma * alpha^2
                trial_point_is_dominated = filter_strictly_dominates_point_at_index(
                    cache, next_filter.index, next_filter.flags, 
                    ci_trial, offset_alpha
                )
                if trial_point_is_dominated
                    @logmsg log_level @sprintf("\t\tTrial point %d is dominated by the filter with offset %.2e.", ci_trial, offset_alpha)
                    continue
                else
                    ## do line-search
                    while true
                        beta = alpha * 2

                        @. x_beta = x + sgn * beta * d
                        project_into_box!(x_beta, lb, ub)
                        @. step_beta = x_trial - x_beta

                        if LA.norm(step_beta) < tol_sgn_switch 
                            @logmsg log_level "\t\tHitting problem boundary."
                            break#while true
                        end#if

                        ci_beta = eval_and_cache!(
                            fx, cx, zx_beta, eps_iq, cache, x_beta, Objfs!, Constrs!;
                            initialize_eps_iq=false, check_cache=true, max_func_calls
                        )
                        if ci_beta < 1
                            it_code = STOP_MAX_FUNC_CALLS
                            break#while true
                        end

                        offset_beta = gamma * beta^2
                        offset_beta_alpha = offset_beta - gamma * alpha^2
                        
                        if any(zx_beta[ℓ] >= zx_trial[ℓ] - offset_beta_alpha for ℓ=eachindex(zx_beta))
                            ## add_point x_trial to filter
                            si_trial, num_filter_free = update_filter!(next_filter, cache, ci_trial; force_add=si_x)
                            if si_trial == si_x
                                si_filter_slot_overwritten = true
                            end
                            ## set stepsize for next (sub)iteration
                            hd.sz_vals[si_trial] = alpha
                            ## set metadata
                            next_filter.directions[si_trial]=dir_index
                            next_filter.generations[si_trial]=it_index
                            next_filter.parents[si_trial]=si

                            @logmsg log_level "\t\t* Added (filter, cache) $((si_trial, ci_trial)), α=$(sgn*alpha)."
                            num_dir_sols += 1
                        end
                        if filter_strictly_dominates_point_at_index(
                            cache, next_filter.index, next_filter.flags, 
                            ci_beta, offset_beta
                        )
                            break # from while-loop doing line search
                        end#if

                        @. x_trial = x_beta
                        @. zx_trial = zx_beta
                        ci_trial = ci_beta
                        alpha = beta
                    end#while (line search / projected expansion)
                end#if trial_point_is_dominated
                it_code != CONTINUE_ITERATION && break#for sgn in (1, -1)
                if num_dir_sols > 0
                    if hd.do_subiterations
                        #=
                        if si_x != si_trial && alpha_start > stepsize_stop
                            alpha_new = alpha_start/2
                            @logmsg log_level @sprintf("*Reducing α[%d] from %.2e to %.2e.", ci_x, alpha_start, alpha_new)
                            hd.sz_vals[si_x] = alpha_new
                        end
                        =#
                        ci_x = ci_trial
                        si_x = si_trial
                    end
                    any_dir_success = true
                else
                    ## failure step with direction `dir_index`
                    if alpha_start > stepsize_stop && hd.sz_vals[si] != alpha_start && !si_filter_slot_overwritten
                        alpha_new = alpha_start/2
                        @logmsg log_level @sprintf("\tReducing α[%d] from %.2e to %.2e.", si, alpha_start, alpha_new)
                        hd.sz_vals[si] = alpha_new
                    end
                end#if dir_success
            end#for sgn in (1, -1)
            it_code != CONTINUE_ITERATION && break#for (dir_index, d) = enumerate(eachcol(hd.orth_dirs))
        end#for (dir_index, d) = enumerate(eachcol(hd.orth_dirs))
        it_code != CONTINUE_ITERATION && break#for (dir_index, d) = enumerate(eachcol(hd.orth_dirs))
        if hd.stepsize_halving && alpha_start > stepsize_stop && !si_filter_slot_overwritten
            _alpha = hd.sz_vals[si]
            alpha_new = _alpha/2
            @logmsg log_level @sprintf("\tReducing α[%d] from %.2e to %.2e.", si, _alpha, alpha_new)
            hd.sz_vals[si] = alpha_new
        end
    end#for si in solutions_int_indices

    @logmsg log_level "\tInvestigated $(num_points_found) solutions."

    copy_filter!(filter, next_filter)
    copyto!(previous_solutions_int_indices, solutions_int_indices)
    
    it_code != CONTINUE_ITERATION && return it_code
    
    prepare_for_next_iteration!(hd; log_level)

    #=
    if num_points_found <= num_objfs && !isnothing(spread_vals)
        coef_delta *= 0.95
        @logmsg log_level "\tReducing coef_delta to $(coef_delta)."
        if coef_delta < 10^round((log10(eps(Float64))))
            return STOP_SMALL_SPREAD_VALUES
        end
        hd.coef_delta[] = coef_delta
    end
    =#
     
    return it_code
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
