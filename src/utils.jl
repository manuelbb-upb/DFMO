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

function penalize!(zx, fx, cx, eps_iq) # `functs_pen` in Fortran
    T = eltype(zx)
    viol = zero(T)

    zx .= fx
    for i = eachindex(cx)
        viol_i = max(0, cx[i])
        viol += viol_i
        for j = eachindex(fx)
            zx[j] += viol_i / eps_iq[i, j]      # TODO this is a bit slower than if we swapped dimensions, but allows for computation of viol
        end
    end

    return viol
end

function project_into_box!(x, lb, ub)
    for i = eachindex(x)
        x[i] = max(lb[i], min(ub[i], x[i]))
    end
end

get_x(cache)=cache.x[:, 1:cache.current_position[]]
get_fobj(cache)=cache.fobj[:, 1:cache.current_position[]]
get_viol(cache)=cache.viol[1:cache.current_position[]]

function get_x(cache::EvalCache, filter::Filter)
    return get_x(cache)[:, filter.index[filter.flags]]
end
function get_fobj(cache::EvalCache, filter::Filter)
    return get_fobj(cache)[:, filter.index[filter.flags]]
end
function get_viol(cache::EvalCache, filter::Filter)
    return get_viol(cache)[filter.index[filter.flags]]
end

function get_x(cache::EvalCache, filter::Filter, solutions)
    return get_x(cache)[:, filter.index[solutions]]
end
function get_fobj(cache::EvalCache, filter::Filter, solutions)
    return get_fobj(cache)[:, filter.index[solutions]]
end
function get_viol(cache::EvalCache, filter::Filter, solutions)
    return get_viol(cache)[filter.index[solutions]]
end

export get_x, get_fobj, get_viol