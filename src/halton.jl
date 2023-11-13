abstract type AbstractDirectionScheme end
abstract type AbstractDirectionSchemeConfig end

"""
    init_direction_scheme(
        cfg::AbstractDirectionSchemeConfig, x, zx, lb, ub, solutions, filter)

Return an object of type `AbstractDirectionScheme` corresponding to configuration `cfg`.
"""
function init_direction_scheme(
    ::AbstractDirectionScheme, x, zx, lb, ub, filter)::AbstractDirectionScheme
    return nothing
end

function propagate_solutions!(
    hd::DS, arrays, cache, filter, 
    Objfs!, Constrs!,
    lb, ub, log_level, stepsize_stop, max_func_calls, it_index
) where{DS <: AbstractDirectionScheme}
    return error("`propagate_solutions!` not implemented for $(DS).")
end

abstract type AbstractHaltonStrategy end
Base.@kwdef struct FollowChildSolutions <: AbstractHaltonStrategy
    max_num_sols :: Int = 1
    prefer_small_stepsizes :: Bool = true
end

Base.@kwdef struct OriginalFiltering{T<:AbstractFloat} <: AbstractHaltonStrategy
    coef_delta_ref :: Base.RefValue{T} = Ref(1.0)
    spread_vals :: Vector{T} = T[]
    spread_hash :: Base.RefValue{UInt64} = Ref(zero(UInt64))
    min_coef_delta :: T = 10^(round(log10(eps(Float64))))
end

Base.@kwdef struct HaltonConfig{S<:AbstractHaltonStrategy} <: AbstractDirectionSchemeConfig
    max_dirs_per_it :: Union{Int, Nothing} = nothing
    tol_sgn_switch :: Union{Real, Nothing} = nothing
    gamma :: Union{Real, Nothing} = nothing
    do_subiterations :: Bool = false
    stepsize_halving :: Bool = true
    exploration :: Bool = false
    solution_selection_strategy :: S = FollowChildSolutions()
end

struct Halton{T<:AbstractFloat, S} <: AbstractDirectionScheme
    num_vars :: Int
    
    sz_init :: T
    sz_vals :: Vector{T}
    
    seed :: Base.RefValue{Int}
    orth_dirs :: Matrix{T}

    x_trial :: Vector{T}
    x_beta :: Vector{T}
    zx_trial :: Vector{T}
    zx_beta :: Vector{T}
    step_trial :: Vector{T}
    step_beta :: Vector{T}

    tol_sgn_switch :: T
    gamma :: T

    do_subiterations :: Bool
    stepsize_halving :: Bool
    exploration :: Bool 
    solution_selection_strategy :: S

    next_filter :: Filter

    solutions_flags :: Vector{Bool}
    tmp_solutions_flags :: Vector{Bool}
    solutions_choice_index :: Base.RefValue{Int}
end

function get_x(cache::EvalCache, filter::Filter, hd::Halton)
    return get_x(cache, filter, hd.solutions_flags)
end
function get_fobj(cache::EvalCache, filter::Filter, hd::Halton)
    return get_fobj(cache, filter, hd.solutions_flags)
end
function get_viol(cache::EvalCache, filter::Filter, hd::Halton)
    return get_viol(cache, filter, hd.solutions_flags)
end

function init_direction_scheme(
    cfg::HaltonConfig, x, zx, lb, ub, filter
)
    T = eltype(x)
    num_vars = length(x)
    max_dirs_per_it = isnothing(cfg.max_dirs_per_it) ? num_vars : min(max(1, cfg.max_dirs_per_it), num_vars)

    seed = 1_000 + 2 * num_vars
    d0 = halton(num_vars, seed, T)
    orth_dirs = Matrix{T}(undef, num_vars, max_dirs_per_it)
    halton_base!(orth_dirs, d0)
    
    box_width = ub .- lb
    sz_init = min(10, maximum(box_width ./ 10))

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

    next_filter = deepcopy(filter)
    solutions_flags = zeros(Bool, filter.size)
    tmp_solutions_flags = deepcopy(solutions_flags)
    solutions_choice_index = Ref(-1)
    return Halton{T, typeof(cfg.solution_selection_strategy)}(
        num_vars,
        sz_init, sz_vals,
        Ref(seed), orth_dirs, 
        x_trial, x_beta, zx_trial, zx_beta, step_trial, step_beta,
        tol_sgn_switch, gamma,
        cfg.do_subiterations, cfg.stepsize_halving, cfg.exploration, 
        cfg.solution_selection_strategy,
        next_filter, solutions_flags, tmp_solutions_flags, solutions_choice_index
    )
end#init_direction_scheme

function choose_child_solutions!(
    strategy::OriginalFiltering,
    next_solutions_flags, last_solutions_flags, sz_vals, filter, last_filter, 
    cache, stepsize_stop, it_index, log_level)

    @logmsg log_level "##### SOLUTION SELECTION #####"
    next_solutions_flags .= false

    alpha_max = maximum(sz_vals[filter.flags]; init=0)
    if alpha_max <= stepsize_stop 
        return STOP_MIN_STEPSIZE
    end
    alpha_max /= 10

    @unpack spread_vals, coef_delta_ref = strategy
    if isempty(spread_vals)
        append!(spread_vals, fill(NaN, filter.size))
    end

    if strategy.spread_hash[] != filter.mod_hash[]
        @logmsg log_level " recomputing SPREAD VALUES"
        sz_index = compute_spread!(spread_vals, filter, cache)
        #src sz_vals[filter.flags] .= sz_vals[filter.flags][sz_index]
        strategy.spread_hash[] = filter.mod_hash[]
    end
    
    coef_delta = coef_delta_ref[]
    _i = 0
    num_considered = 0
    num_log_msgs = 0
    for (i, i_is_filter_position) = enumerate(filter.flags)
        !i_is_filter_position && continue
        alpha = sz_vals[i]
        if alpha <= stepsize_stop 
            if num_log_msgs <= 10
                @logmsg log_level "  Filter element $(i) has too small a stepsize ($(alpha))."
                num_log_msgs += 1
            end
            continue
        end
        
        if spread_vals[i] <= coef_delta 
            if num_log_msgs <= 10
                @logmsg log_level "  Filter element $(i) has too small a spread ($(spread_vals[i]) vs $(coef_delta))."
                num_log_msgs += 1
            end
            continue
        end

        _i += 1
        if _i > size(cache.fobj, 1) && alpha <= alpha_max 
            if num_log_msgs <= 10
                @logmsg log_level "  Filter element $(i) has too small a relative stepsize ($(alpha) vs $(alpha_max))."
                num_log_msgs += 1
            end
            continue
        end

        num_considered += 1
        next_solutions_flags[i] = true
        @logmsg log_level "  Choosing filter element $(i) with spread $(spread_vals[i]) for next iteration!"
    end

    if num_considered == 0
        coef_delta_ref[] = coef_delta * 0.95
        if coef_delta_ref[] < strategy.min_coef_delta
            return STOP_SMALL_SPREAD_VALUES
        end
    end
    return CONTINUE_ITERATION
end

function choose_child_solutions!(
    strategy::FollowChildSolutions,
    next_solutions_flags, last_solutions_flags, sz_vals, filter, last_filter,
    cache, stepsize_stop, it_index, log_level)
    
    @logmsg log_level "##### SOLUTION SELECTION #####"
    next_solutions_flags .= false
    num_sols = 0
    @unpack max_num_sols = strategy
    
    if it_index == 1
        for (fi, fi_slot_occupied) = enumerate(filter.flags)
            !fi_slot_occupied && continue
            sz_vals[fi] <= stepsize_stop && continue
            num_sols == max_num_sols && break
            num_sols += 1
            next_solutions_flags[fi]=true
        end
    else
        found_more_candidates = true
        last_gen = it_index - 1
        while num_sols < max_num_sols && found_more_candidates
            _num_sols = num_sols
            for (pfi, pfi_was_solution) = enumerate(last_solutions_flags)
                num_sols >= max_num_sols && break
                !pfi_was_solution && continue

                if filter.flags[pfi] && filter.index[pfi] == last_filter.index[pfi]
                    if !next_solutions_flags[pfi] && sz_vals[pfi] > stepsize_stop
                        num_sols += 1
                        next_solutions_flags[pfi] = true
                        @logmsg log_level "\tSOLUTION SELECTION [$(num_sols)]: Previous solution is still in filter."
                    
                        continue
                    end
                end

                si = -1
                si_sz = strategy.prefer_small_stepsizes ? Inf : -Inf
                for (fi, fi_slot_occupied) = enumerate(filter.flags)
                    !fi_slot_occupied && continue
                    filter.parents[fi] != pfi && continue
                    filter.generations[fi] != last_gen && continue
                    next_solutions_flags[fi] && continue
                    fi_sz = sz_vals[fi]
                    fi_sz <= stepsize_stop && continue
                    if (
                        (strategy.prefer_small_stepsizes && fi_sz < si_sz) ||
                        (!strategy.prefer_small_stepsizes && fi_sz > si_sz)
                    )
                        si_sz = fi_sz
                        si = fi
                    end
                end
                if si > 0
                    num_sols += 1
                    next_solutions_flags[si] = true                
                    @logmsg log_level "\tSOLUTION SELECTION [$(num_sols)]: Previous solution $(pfi) has filter descendant $(si)."
                end
            end
            found_more_candidates = (num_sols != _num_sols)
        end
    end
    @logmsg log_level "\tSOLUTION SELECTION: Chose $(num_sols) elements in filter."
    if num_sols == 0
        return STOP_MIN_STEPSIZE
    else
        return CONTINUE_ITERATION
    end
end

function propagate_solutions!(
    hd::Halton, 
    arrays, cache, filter, 
    Objfs!, Constrs!,
    lb, ub, log_level, stepsize_stop, max_func_calls, it_index
)
    @unpack x, cx, fx, zx, eps_iq, num_vars, num_objfs = arrays
    @unpack x_trial, x_beta, zx_trial, zx_beta, step_trial, step_beta, sz_vals,
        tol_sgn_switch, gamma, next_filter, solutions_flags, solutions_choice_index, 
        tmp_solutions_flags, stepsize_halving = hd

    @logmsg log_level "#### PROPAGATION"

    copy_filter!(next_filter, filter)
    
    if solutions_choice_index[] != it_index - 1
        ## select those element from the filter that have the smallest stepsize associated with them
        #src choose_minimal_solutions!(solutions_flags, sz_vals, filter, max_num_sols)
        choose_child_solutions!(
            hd.solution_selection_strategy, solutions_flags, tmp_solutions_flags, sz_vals, 
            next_filter, filter, cache, stepsize_stop, it_index, log_level)
        solutions_choice_index[] = it_index - 1
    end
    tmp_solutions_flags .= solutions_flags
    
    it_code = CONTINUE_ITERATION
   
    for (si, si_is_solution_index) = enumerate(solutions_flags)
        !si_is_solution_index && continue
        @assert filter.flags[si]
        alpha_start = sz_vals[si]
        
        ci = filter.index[si]
        @logmsg log_level "  Considering $(si) ($(ci)) with α=$(alpha_start)."

        alpha_fail = alpha_start/2
        si_filter_slot_overwritten = false
        si_trial = si_x = si 
        success_step = false
        dir_success = false

        # Print vectors
        ci_x = filter.index[si_x]
        x .= cache.x[:, ci_x]
        zx .= cache.fobj[:, ci_x]
        @logmsg log_level "\tConstraint Violation is $(cache.viol[ci_x]) and vectors are\n$(CustomPrinting.pretty_str(CustomPrinting.VectorTable("x"=>x, "zx" => zx); line_prefix="\t\t"))"
 
        for (dir_index, d) = enumerate(eachcol(hd.orth_dirs))
            (success_step && !hd.do_subiterations && !hd.exploration) && break
            
            if dir_success && hd.do_subiterations
                si_x = si_trial
                @logmsg log_level "\t (do_subiterations) Changing base to $(si_x)."
                ci_x = next_filter.index[si_x]
                x .= cache.x[:, ci_x]
                zx .= cache.fobj[:, ci_x]
            end
            dir_success = false
        
            alpha = alpha_start
            for sgn in (1, -1)
                dir_success && break
                @logmsg log_level "\t  Testing direction $(dir_index) with α=$(sgn*alpha)."
                
                @. x_trial = x + sgn * alpha * d
                project_into_box!(x_trial, lb, ub)
                @. step_trial = x_trial - x
                
                if LA.norm(step_trial) < tol_sgn_switch
                    @logmsg log_level "\t\tInitial step to small."
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
                    @logmsg log_level @sprintf("\t\tci=%d is dominated by the filter with offset %.2e.", ci_trial, offset_alpha)
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
                        offset_beta = gamma * beta^2
                        offset_beta_alpha = offset_beta - gamma * alpha^2
                        
                        if ci_beta < 1 || any(zx_beta[ℓ] >= zx_trial[ℓ] - offset_beta_alpha for ℓ=eachindex(zx_beta))
                            ## add_point x_trial to filter
                            si_trial, num_filter_free = update_filter!(next_filter, cache, ci_trial; force_add=si_x)
                            if si_trial == si && next_filter.index[si_trial] != ci
                                si_filter_slot_overwritten = true
                            end
                            ## set stepsize for next (sub)iteration
                            sz_vals[si_trial] = alpha
                            
                            ## set metadata (not necessary for Halton)
                            next_filter.directions[si_trial]=dir_index
                            next_filter.generations[si_trial]=it_index
                            next_filter.parents[si_trial]=si

                            @logmsg log_level "\t\t* Added (filter, cache) $((si_trial, ci_trial)), α=$(sgn*alpha)."
                            dir_success = success_step = true
                        end
                        if ci_beta < 1
                            it_code = STOP_MAX_FUNC_CALLS
                            break#while true
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
                it_code != CONTINUE_ITERATION && break
            end#for sgn in (1, -1)
            it_code != CONTINUE_ITERATION && break
            if !dir_success && !si_filter_slot_overwritten
                if sz_vals[si] != alpha_fail
                    @logmsg log_level @sprintf("\t(failure step) Reducing α[%d] from %.2e to %.2e.", si, alpha_start, alpha_fail)
                    sz_vals[si] = alpha_fail
                end
            end
        end#for (dir_index, dir) = enumerate ...
        it_code != CONTINUE_ITERATION && break

        if stepsize_halving && !si_filter_slot_overwritten
            alpha_start = min(sz_vals[si], alpha_start)
            alpha_fail = alpha_start / 2
            @logmsg log_level @sprintf("\t(stepsize_halving) Reducing α[%d] from %.2e to %.2e.", si, alpha_start, alpha_fail)
            sz_vals[si] = alpha_fail
        end
    end#for si....
    
    #src choose_minimal_solutions!(solutions_flags, sz_vals, filter, max_num_sols)
    _it_code = choose_child_solutions!(
        hd.solution_selection_strategy, tmp_solutions_flags, solutions_flags, sz_vals, 
        next_filter, filter, cache, stepsize_stop, it_index+1, log_level)
    solutions_flags .= tmp_solutions_flags
    solutions_choice_index[] = it_index
    
    copy_filter!(filter, next_filter)   # after (!!!) `choose_child_solutions`

    prepare_next_halton_directions!(hd; log_level)
    it_code = it_code != CONTINUE_ITERATION ? it_code : _it_code
    return it_code
end

function prepare_next_halton_directions!(hd::Halton{T}; log_level) where T
    seed = hd.seed[] += 2*hd.num_vars
    @logmsg log_level "\tComputing halton direction for halton index $(seed)."
    d0 = halton(hd.num_vars, seed, T)
    halton_base!(hd.orth_dirs, d0)
    return nothing
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
    dir .-= .5
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
    
    ncols = size(H, 2)
    if ncols > 1
        _H, _ = LA.qr(d)
        H[:, 1:ncols] .= _H[:, 1:ncols]
    else
        H[:] .= d ./ LA.norm(d)
    end
    
    return nothing
end
