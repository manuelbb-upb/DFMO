module CustomPrinting

import Formatting: sprintf1

function mantissa_and_exponent(num)
    if !(num isa AbstractFloat)
        num = Float64(num)
    end
    if iszero(num)
        return num, 0
    else
        exponent = floor(log10(abs(num)))
        mantissa = num / 10^(exponent)
        return mantissa, exponent
    end
end

const EXP_STRS = Base.ImmutableDict(
    1 => "¹",
    2 => "²",
    3 => "³",
    4 => "⁴",
    5 => "⁵", 
    6 => "⁶",
    7 => "⁷",
    8 => "⁸",
    9 => "⁹", 
    0 => "⁰"
)

function pretty_vector_str(
    vec::AbstractVector{T}; ndigits=3, oneline::Bool=false, indent::Int=0) where T
    if T <: AbstractFloat
        F = T
        mants = zero(vec)
        exps = zero(vec)
    else
        F = Float64
        mants = zeros(Float64, length(vec))
        exps = zeros(Float64, length(vec))
    end
    indent_str = lpad("", indent)
    s = indent_str * "{$(F)}["
    if !oneline s *= "\n$(indent_str)" end
    e_min = typemax(Int)
    e_max = typemin(Int)
    for (i, x) = enumerate(vec)
        xm, xe = mantissa_and_exponent(x)
        mants[i] = xm
        exps[i] = xe
        if xe <= e_min
            e_min = xe
        end
        if e_max <= xe
            e_max = xe
        end
    end

    str_len = Int(3 + (e_max - e_min) + ndigits)
    for i = eachindex(mants)
        xe = exps[i]
        if xe > e_min
            de = xe - e_min
            mants[i] *= 10^de
        end
        s *= sprintf1("% $(str_len).$(ndigits)f, ", mants[i])
        if !oneline s *= "\n$(indent_str)" end
    end

    s *= "] ⋅ 10"
    
    if e_min < 0 
        s*= "⁻" 
        e_min *= -1
    end

    for edig in digits(Int(abs(e_min)))
        s *= EXP_STRS[edig]
    end

    return s
end

struct VectorTableEntry{T}
    mantissas :: Vector{T}
    exponents :: Vector{T}

    min_exponent :: T
    max_exponent :: T
    n_digits :: Int
end

function VectorTableEntry(vec::AbstractVector{F}) where F<:Real
    T = F <: AbstractFloat ? F : Float64
    len = length(vec)
    mantissas = zeros(T, len)
    exponents = zeros(T, len)
    min_exponent = typemax(T)
    max_exponent = typemin(T)
    for (i, x) = enumerate(vec)
        xm, xe = mantissa_and_exponent(x)
        mantissas[i] = xm
        exponents[i] = xe
        if xe <= min_exponent
            min_exponent = xe
        end
        if max_exponent <= xe
            max_exponent = xe
        end
    end 
    return VectorTableEntry{T}(mantissas, exponents, min_exponent, max_exponent, 1)
end

function denormalize_mantissas(vec_table_entry)
    e_max = vec_table_entry.max_exponent
    e_min = vec_table_entry.min_exponent

    e_max == e_min && return vec_table_entry
    mants = copy(vec_table_entry.mantissas)
    exps = copy(vec_table_entry.exponents)
    n_digits = 1 + round(Int, e_max - e_min)
    for i = eachindex(mants)
        exps_i = exps[i]
        if exps_i > e_min
            diff = exps_i - e_min
            mants[i] *= 10^diff
            exps[i] = e_min
        end
    end
    return VectorTableEntry(mants, exps, e_min, e_min, n_digits)
end

function pretty_str(
    vec_table_entry::VectorTableEntry; 
    one_line::Bool=false, line_prefix::String="",
    digits_after_comma::Integer=3,
)
    v = denormalize_mantissas(vec_table_entry)
    s = line_prefix * "[" 

    if !one_line s *= "\n" end
    len = 1 + v.n_digits + 1 + digits_after_comma
    num_entries = length(v.mantissas)
    for (i, m) in enumerate(v.mantissas)
        s *= line_prefix * sprintf1("% $(len).$(digits_after_comma)f", m)        
        if i < num_entries s *= ", " end
        if !one_line s *= "\n" end
    end

    s *= line_prefix * "] ⋅ 10"
    if v.min_exponent < 0
        s *= "⁻"
        e = Int(-v.min_exponent)
    else
        e = Int(v.min_exponent)
    end
    s *= EXP_STRS[e]
    return s
end
 
struct VectorTable{VTE}
    names :: Vector{String}
    entries :: Vector{VTE}
    vec_lens :: Vector{Int}
    max_vec_len :: Int
end

function VectorTable(name_vec_pairs::Vararg{Pair{String, <:AbstractVector{<:Real}}})
    names = String[]
    max_vec_len = 0
    vec_lens = Int[]

    if isempty(name_vec_pairs)
        entries = VectorTableEntry[]
    else
    
        T = mapreduce(p -> eltype(last(p)), Base.promote_type, name_vec_pairs)
        entries = VectorTableEntry{T}[]
        
        for (n, v) in name_vec_pairs
            push!(names, n)
            push!(entries, VectorTableEntry(T.(v)))

            vec_len_entry = length(v)
            push!(vec_lens, vec_len_entry)
            if vec_len_entry >= max_vec_len
                max_vec_len = vec_len_entry
            end
        end
    end

    return VectorTable(names, entries, vec_lens, max_vec_len)
end

function pretty_str(
    vec_table::VectorTable; 
    line_prefix::String="",
    digits_after_comma::Integer=3,
    index_letter::String="i",
    index_column::Bool=true,
)

    first_col_width = index_column ? max(1, length(index_letter), ndigits(vec_table.max_vec_len)) : 0
    entries = [denormalize_mantissas(e) for e in vec_table.entries]
    col_widths = zeros(Int, length(entries))
    
    if index_column
        s = lpad(index_letter, first_col_width) * " |"
    else
        s = ""
    end
    for (j, name) in enumerate(vec_table.names)
        entry = entries[j]
        cw = max(length(name), ndigits(round(Int, entry.max_exponent)) + 2, 2 + entry.n_digits + digits_after_comma)
        col_widths[j] = cw
        s *= " " * lpad("$(name)", cw) * " |" 
    end

    s *= "\n"

    if index_column
        s *= lpad("", first_col_width) * " |"
    end

    for (j, entry) = enumerate(entries)
        cw = col_widths[j]
        int_exp = round(Int, entry.max_exponent)
        if int_exp < 0
            s *= " " * lpad("E$(int_exp)", cw) * " |"
        else
            s *= " " * lpad("E+$(int_exp)", cw) * " |"
        end
    end

    s *= "\n" * repeat("=", first_col_width + 2 + sum(3 + cw for cw in col_widths)) * "\n"

    for i=1:vec_table.max_vec_len
        if index_column
            s *= lpad("$i", first_col_width) * " |"
        end
        for (j, entry) = enumerate(entries)
            cw = col_widths[j]
            if i <= length(entry.mantissas)
                s *= " " * sprintf1("% $(cw).$(digits_after_comma)f", entry.mantissas[i]) * " |"
            else
                s *= " " * lpad("", cw) * " |"
            end
        end
        if i != vec_table.max_vec_len
            s *= "\n"
        end
    end
    return s
end

end#module