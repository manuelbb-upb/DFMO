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

function query_solutions_filter_index(::AbstractDirectionScheme, filter)
    return Int[]
end

function query_solutions_cache_index(ds::AbstractDirectionScheme, filter)
    return filter.index[query_solutions_filter_index(ds, filter)]
end

Base.@kwdef struct HaltonDriveToZeroConfig <: AbstractDirectionSchemeConfig
    max_num_sols :: Int = 1
    max_dirs_per_it :: Union{Int, Nothing} = 1
    tol_sgn_switch :: Union{Real, Nothing} = nothing
    gamma :: Union{Real, Nothing} = nothing
    do_subiterations :: Bool = false
    stepsize_halving :: Bool = false
end

struct HaltonDriveToZero{T<:AbstractFloat} <: AbstractDirectionScheme
    num_vars :: Int
    max_num_sols :: Int
    
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

    next_filter :: Filter

    solutions_flags :: Vector{Bool}
    solutions_choice_index :: Base.RefValue{Int}
end

function init_direction_scheme(
    cfg::HaltonDriveToZeroConfig, x, zx, lb, ub, filter
)
    T = eltype(x)
    num_vars = length(x)
    max_num_sols = cfg.max_num_sols
    max_dirs_per_it = isnothing(cfg.max_num_sols) ? num_vars : min(max(1, cfg.max_dirs_per_it), num_vars)

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
    solutions_choice_index = Ref(-1)
    return HaltonDriveToZero{T}(
        num_vars, max_num_sols,
        sz_init, sz_vals,
        Ref(seed), orth_dirs, 
        x_trial, x_beta, zx_trial, zx_beta, step_trial, step_beta,
        tol_sgn_switch, gamma,
        cfg.do_subiterations, cfg.stepsize_halving,
        next_filter, solutions_flags, solutions_choice_index
    )
end#init_direction_scheme

function choose_minimal_solutions!(solutions_flags, sz_vals, filter, max_num_sols)
    solutions_flags .= false
    filter_sort_index = sortperm(sz_vals[filter.flags])
    num_sols = min(max_num_sols, length(filter_sort_index))
    @views solutions_flags[filter.flags][filter_sort_index[1:num_sols]] .= true
    return nothing
end

function propagate_solutions!(
    hd::HaltonDriveToZero, 
    arrays, cache, filter, 
    Objfs!, Constrs!,
    lb, ub, log_level, stepsize_stop, max_func_calls, it_index
)
    @unpack x, cx, fx, zx, eps_iq, num_vars, num_objfs = arrays
    @unpack x_trial, x_beta, zx_trial, zx_beta, step_trial, step_beta, sz_vals,
        tol_sgn_switch, gamma, next_filter, solutions_flags, solutions_choice_index, 
        stepsize_halving, max_num_sols = hd

    @logmsg log_level "#### PROPAGATION"

    if solutions_choice_index[] != it_index - 1
        ## select those element from the filter that have the smallest stepsize associated with them
        choose_minimal_solutions!(solutions_flags, sz_vals, filter, max_num_sols)
        solutions_choice_index[] = it_index - 1
    end
    
    copy_filter!(next_filter, filter)
    it_code = CONTINUE_ITERATION
   
    for (si, si_is_solution_index) = enumerate(solutions_flags)
        !si_is_solution_index && continue
        @assert filter.flags[si]
        alpha_start = sz_vals[si]
        @logmsg log_level "  Considering $(si) with α=$(alpha_start)."
        if alpha_start < stepsize_stop
            @logmsg log_level "\tStepsize is $(alpha_start) <= $(stepsize_stop), SKIPPING."
            continue
        end

        alpha_fail = alpha_start/2
        si_filter_slot_overwritten = false
        si_trial = si_x = si 
        success_step = false
        dir_success = false

        # Print vectors
        ci_x = filter.index[si_x]
        x .= cache.x[:, ci_x]
        zx .= cache.fobj[:, ci_x]
        @logmsg log_level "\tVectors are\n$(CustomPrinting.pretty_str(CustomPrinting.VectorTable("x"=>x, "zx" => zx); line_prefix="\t\t"))"
 
        for (dir_index, d) = enumerate(eachcol(hd.orth_dirs))
            (success_step && !hd.do_subiterations) && break

            if dir_success
                si_x = si_trial
                ci_x = filter.index[si_x]
                x .= cache.x[:, ci_x]
                zx .= cache.fobj[:, ci_x]
            end
            dir_success = false
        
            alpha = alpha_start
            for sgn in (1, -1)
                dir_success && break

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
                if trial_point_is_dominated && sgn < 0
                    @logmsg log_level @sprintf("\t\tci=%d is dominated by the filter with offset %.2e.", ci_trial, offset_alpha)
                    if !si_filter_slot_overwritten && sz_vals[si] != alpha_fail
                        @logmsg log_level @sprintf("\t\tReducing α[%d] from %.2e to %.2e.", si, alpha_start, alpha_fail)
                        sz_vals[si] = alpha_fail
                    end
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
                            if si_trial == si
                                si_filter_slot_overwritten = true
                            end
                            ## set stepsize for next (sub)iteration
                            sz_vals[si_trial] = alpha
                            ## set metadata (not necessary for HaltonDriveToZero)
                            next_filter.directions[si_trial]=dir_index
                            next_filter.generations[si_trial]=it_index
                            next_filter.parents[si_trial]=si

                            @logmsg log_level "\t\t* Added (filter, cache) $((si_trial, ci_trial)), α=$(sgn*alpha)."
                            dir_success = success_step = true
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
        end#for (dir_index, dir) = enumerate ...
        it_code != CONTINUE_ITERATION && break
        if stepsize_halving &&!si_filter_slot_overwritten
            alpha_start = sz_vals[si]
            alpha_fail = alpha_start / 2
            @logmsg log_level @sprintf("\t(stepsize_halving) Reducing α[%d] from %.2e to %.2e.", si, alpha_start, alpha_fail)
            sz_vals[si] = alpha_fail
        end
    end#for si....
    copy_filter!(filter, next_filter)
    choose_minimal_solutions!(solutions_flags, sz_vals, filter, max_num_sols)
    solutions_choice_index[] = it_index
    prepare_next_halton_directions!(hd; log_level)
    return it_code
end

function prepare_next_halton_directions!(hd::HaltonDriveToZero{T}; log_level) where T
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
