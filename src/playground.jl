using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using DFMO
using CairoMakie
#%%

n = 2
lb = fill(-4, n)
ub = fill(4, n)

objfs! = function(y, x)
    y[1] = sum( (x .+ 1).^2 )
    y[2] = sum( (x .- 1).^2 )
end
num_objfs = 2

x0 = lb .+ (ub .- lb) .* rand(n)
#%%

halton_cfg = DFMO.HaltonConfig(;
    max_dirs_per_it = 1,
    exploration = false,
    stepsize_halving = true,
    do_subiterations = false,
    solution_selection_strategy = DFMO.FollowChildSolutions(;
        max_num_sols = 1,
        prefer_small_stepsizes = true
    )
)

ret = DFMO.optimize(
    x0, lb, ub, num_objfs, objfs!;
    max_func_calls = 100,
    max_iter = 20,
    direction_cfg = halton_cfg
)
#%%
halton_cfg2 = DFMO.HaltonConfig(;
    max_dirs_per_it = nothing,
    exploration = true,
    do_subiterations = false,
    solution_selection_strategy = DFMO.OriginalFiltering{Float64}(;)
)

ret2 = DFMO.optimize(
    x0, lb, ub, num_objfs, objfs!;
    max_func_calls = 2_000,
    max_iter = typemax(Int),
    direction_cfg = halton_cfg2
)

#%%

function plot_ret(ret; getter_func = get_fobjf, fig=nothing, ax=nothing)
    if isnothing(fig)
        fig = Figure()
    end
    if isnothing(ax)
        ax = Axis(fig[1,1])
    end
    
    F_eval = getter_func(ret.cache)
    F_filter = getter_func(ret.cache, ret.filter)
    F_fin = getter_func(ret.cache, ret.filter, ret.dir_scheme)
    
    sz_filter = ret.dir_scheme.sz_vals[ret.filter.flags]
    sz_min, sz_max = extrema(sz_filter)
    sz_div = sz_max != sz_min ? abs(sz_max - sz_min) : 1
    scale_sz(vals) = 10 .+ 10 .* (vals .- sz_min) ./ sz_div

    @show ms_filter = scale_sz(sz_filter)
    ms_fin = scale_sz(ret.dir_scheme.sz_vals[ret.dir_scheme.solutions_flags])

    scatter!(ax, F_eval; markersize=10, color=:red)
    scatter!(ax, F_filter; 
        markersize=ms_filter,
        color=:green,
        strokecolor=:white,
        strokewidth=1,
    )
    scatter!(ax, F_fin;
        markersize=ms_fin,
        color=:blue
    )

    return fig, ax
end
#%%
fig, ax1 = plot_ret(ret; getter_func=get_x)
lines!(ax1, [(-1, -1), (1, 1)]; color=:orange)
ax2 = Axis(fig[1,2])
plot_ret(ret; getter_func=get_fobj, fig, ax=ax2)
fig
#%%
fig, ax1 = plot_ret(ret2; getter_func=get_x)
lines!(ax1, [(-1, -1), (1, 1)]; color=:orange)
ax2 = Axis(fig[1,2])
plot_ret(ret2; getter_func=get_fobj, fig, ax=ax2)
fig

#%%

objfs!2 = function(y, x)
    y[1] = (x[1] - 2)^2 + (x[2] - 1)^2
    y[2] = (x[1] - 2)^2 + (x[2] + 1)^2
end
num_objfs = 2

constrs! = function(y, x)
    # ‖x‖ ≥ 1 ⇔ 1 - ‖x‖ ≤ 0
    y[1] = 1 - sum(x.^2)
end
num_constrs = 1
#%%
ret3 = DFMO.optimize(
    x0, lb, ub, num_objfs, objfs!2, num_constrs, constrs!;
    max_func_calls = 20,
    max_iter = 6,
    direction_cfg = halton_cfg,
)

fig, ax1 = plot_ret(ret3; getter_func=get_x)
lines!(ax1, [(2, 1), (2,-1)]; color=:orange)
lines!(ax1, cos.(0:0.1:2*π), sin.(0:0.1:2*π); color=:orange)
ax2 = Axis(fig[1,2])
ax2.limits = ((0,3), (0,3))
plot_ret(ret3; getter_func=get_fobj, fig, ax=ax2)
fig