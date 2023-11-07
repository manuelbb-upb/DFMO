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