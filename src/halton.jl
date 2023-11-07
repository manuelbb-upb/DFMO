
abstract type AbstractDirectionScheme end
init_direction_scheme(::Type{<:AbstractDirectionScheme}, num_vars, lb, ub, solutions, T)=nothing
max_stepsize(::AbstractDirectionScheme, solutions)=missing
#=
set_stepsize_in_cache!(::AbstractDirectionScheme, sz_cache, sz, cache_index, dir_index)=nothing
get_stepsize_from_cache(::AbstractDirectionScheme, sz_cache, cache_index, dir_index)=nothing
max_stepsize_in_filter(::AbstractDirectionScheme, sz_cache, filter_flags, filter_permutation)=nothing
max_stepsize_in_cache(ds::AbstractDirectionScheme, sz_cache, cache_index)=get_stepsize_from_cache(ds, sz_cache, cache_index, 1)
prepare_for_next_iteration!(::AbstractDirectionScheme; log_level, kwargs...)=nothing
reset_current_state!(::AbstractDirectionScheme)=nothing
"Set `dir_vec` and return a direction index of type `Integer` or return `nothing` to indicate that all directions have been tried."
next_direction!(dir_vec, ::AbstractDirectionScheme)=nothing
=#

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

    do_subiterations :: Bool

    next_filter :: Filter
    next_solutions :: Solutions
end

function init_direction_scheme(::Type{<:HaltonDirections}, x, zx, lb, ub, solutions, filter)
    num_vars = length(x)
    num_sols = length(solutions.flags)
    T = eltype(x)

    ih = 1_000 + 2 * num_vars
    d0 = halton(num_vars, ih, T)
    orth_dirs = Matrix{T}(undef, num_vars, num_vars)
    
    halton_base!(orth_dirs, d0)
    box_width = ub .- lb
    sz_init = min(10, maximum(box_width ./ 10))

    sz_vals = fill(sz_init, num_sols)

    x_trial = copy(x)
    x_beta = copy(x)
    zx_trial = copy(zx)
    zx_beta = copy(zx)
    step_trial = zero(x)
    step_beta = copy(step_trial)

    tol_sgn_switch = 10^(-round(abs(log10(eps(T)))^0.76)) # 1e-8 for T == Float64

    return HaltonDirections{T}(
        num_vars, sz_init, sz_vals, Ref(ih), orth_dirs,
        x_trial, x_beta, zx_trial, zx_beta, step_trial, step_beta,
        tol_sgn_switch, Ref(one(T)), true,
        deepcopy(filter), deepcopy(solutions)
    )
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

@views function max_stepsize(hd::HaltonDirections, solutions)
    return maximum(hd.sz_vals[solutions.flags]; init=eltype(hd.sz_vals)(-Inf))
end

function prepare_for_next_iteration!(hd::HaltonDirections{T}; log_level=Info) where T
    ih = hd.index[] += 2*hd.num_vars
    @logmsg log_level "\tComputing halton direction for halton index $(ih)."
    d0 = halton(hd.num_vars, ih, T)
    halton_base!(hd.orth_dirs, d0)
    return nothing
end

function propagate_solutions!(
    hd::HaltonDirections, arrays, cache, solutions, filter, 
    Objfs!, Constrs!,
    lb, ub, log_level, stepsize_stop, max_func_calls, gamma
)

    @unpack x, cx, fx, zx, eps_iq, num_vars, num_objfs = arrays
    @unpack x_trial, x_beta, zx_trial, zx_beta, step_trial, step_beta, tol_sgn_switch, 
        next_solutions, next_filter = hd;

    sz_max = max_stepsize(hd, solutions)
    if sz_max < stepsize_stop
        return STOP_MIN_STEPSIZE
    end

    num_points_found = 0

    num_sols_total = length(solutions.flags)
    num_free_solution_slots = num_sols_total - sum(solutions.flags)
    max_num_sols_per_dir = max(1, div(num_free_solution_slots, num_sols_total * num_vars))

    next_filter.index .= filter.index
    next_filter.flags .= filter.flags
    next_filter.tmp_flags .= filter.tmp_flags

    next_solutions.flags .= solutions.flags
    next_solutions.index .= solutions.index
    #src next_solutions.spread_vals .= solutions.spread_vals

    it_code = CONTINUE_ITERATION
    coef_delta = hd.coef_delta[]

    for (i, solution_slot_is_occupied) = enumerate(solutions.flags)
        it_code != CONTINUE_ITERATION && break
        !solution_slot_is_occupied && continue
        filter_position = solutions.index[i]
        ci = filter.index[filter_position]

        @logmsg log_level "Considering solution $(i) with cache index $(ci)." 

        alpha_start = alpha = hd.sz_vals[i]
        ## inverting `condizione` (DFMO.f90, ll. 184):
        ## a) consider only filter elements with stepsize big enough
        if alpha <= stepsize_stop
            @logmsg log_level "\tStepsize is $(alpha) <= $(stepsize_stop), skipping."
            continue
        end

        ## b) favor larger spread values
        if !isnothing(solutions.spread_vals) 
            spread = solutions.spread_vals[i]
            if spread < coef_delta
                @logmsg log_level "Spread value is $(spread) < $(coef_delta), skipping."
                continue
            end
        else
            spread = nothing
        end

        ## c) if an element is not filter-extremal, it has to have large *relative* stepsize
        alpha <= sz_max/10 && i > num_objfs && continue

        num_points_found += 1

        ci_x = ci
        i_x = i
        x .= cache.x[:, ci_x]
        zx .= cache.fobj[:, ci_x]
        @logmsg log_level "Spread is $(spread) >= $(coef_delta), vectors are\n$(CustomPrinting.pretty_str(CustomPrinting.VectorTable("x"=>x, "zx" => zx); line_prefix="\t\t"))"

        any_dir_success = false
        for (dir_index, d) = enumerate(eachcol(hd.orth_dirs))
            x .= cache.x[:, ci_x]
            ci_trial = ci_x
            i_trial = i_x
            #src alpha_start = alpha = hd.sz_vals[i_x]
            alpha = alpha_start
            num_dir_sols = 0
            for sgn in (1, -1)
                (num_dir_sols > 0 && (
                    #src num_free_solution_slots == 0 || 
                    hd.do_subiterations || 
                    (dir_index < num_vars && (num_dir_sols >= max_num_sols_per_dir)))
                ) && continue
                @logmsg log_level "LineSearch (cache, solution, direction) $((ci_x, i_x, dir_index)), α=$(sgn*alpha); d=$(d)"

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
                    @logmsg log_level @sprintf("\t Trial point %d is dominated by the filter with offset %.2e.", ci_trial, offset_alpha)
                    continue
                else
                    ## do line-search
                    while true
                        beta = alpha * 2

                        @. x_beta = x + sgn * beta * d
                        project_into_box!(x_beta, lb, ub)
                        @. step_beta = x_trial - x_beta

                        if LA.norm(step_beta) < tol_sgn_switch 
                            @logmsg log_level "\tHitting problem boundary."
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
                        #src if reduce(max, zx_beta[ℓ] >= zx_trial[ℓ] - offset_beta_alpha for ℓ=eachindex(zx_beta))
                            ## add_point x_trial to filter
                            filter_slot_index, num_filter_free = update_filter!(next_filter, cache, ci_trial; force_add=filter_position)
                            sync!(next_solutions, next_filter)
                            soft_force = num_dir_sols < max_num_sols_per_dir || dir_index == num_vars
                            i_trial = set_flag!(next_solutions, filter_slot_index; force_add=i_x, soft_force)
                            hd.sz_vals[i_trial] = alpha

                            @logmsg log_level "\t* Added (cache, filter, solution) $((ci_trial, filter_slot_index, i_trial)), α=$(sgn*alpha)."
                            
                            num_free_solution_slots = num_sols_total - sum(solutions.flags)
                            max_num_sols_per_dir = div(num_free_solution_slots, num_sols_total * num_vars)
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
                        if i_x != i_trial && alpha_start > stepsize_stop
                            alpha_new = alpha_start/2
                            @logmsg log_level @sprintf("*Reducing α[%d] from %.2e to %.2e.", ci_x, alpha_start, alpha_new)
                            hd.sz_vals[i_x] = alpha_new
                        end
                        ci_x = ci_trial
                        i_x = i_trial
                    end
                    any_dir_success = true
                end#if dir_success
            end#for sgn in (1, -1)
            it_code != CONTINUE_ITERATION && break#for (dir_index, d) = enumerate(eachcol(hd.orth_dirs))
        end#for (dir_index, d) = enumerate(eachcol(hd.orth_dirs))
        
        if (!hd.do_subiterations || !any_dir_success) && alpha_start > stepsize_stop
            alpha_new = alpha_start/2
            @logmsg log_level @sprintf("*Reducing α[%d] from %.2e to %.2e.", ci_x, alpha_start, alpha_new)
            hd.sz_vals[i_x] = alpha_new
        end
    end

    @logmsg log_level "Investigated $(num_points_found) points."
    
    filter.flags .= next_filter.flags
    filter.index .= next_filter.index
    #src filter.tmp_flags .= next_filter.tmp_flags
    solutions.flags .= next_solutions.flags
    solutions.index .= next_solutions.index

    it_code != CONTINUE_ITERATION && return it_code
    
    prepare_for_next_iteration!(hd; log_level)

    if num_points_found == 0 && !isnothing(solutions.spread_vals)
        coef_delta *= 0.95
        @logmsg log_level "Reducing coef_delta to $(coef_delta)."
        if coef_delta < 10^round((log10(eps(Float64))))
            return STOP_SMALL_SPREAD_VALUES
        end
        hd.coef_delta[] = coef_delta
    end
     
    return it_code
end