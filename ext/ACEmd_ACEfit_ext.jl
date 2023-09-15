module ACEmd_ACEfit_ext

using ACEfit
using ACEmd
using Folds
using StaticArrays

function ACEfit.feature_matrix(data, basis; energy=true, force=true, virial=true, kwargs...)
    # Basis functions are on different collumns.
    # Energy is on fist row.
    # Force is flattened on several rows, so that each basis function is on same collumns.
    # Virial is in practice triangular matrix, so first remove double values and then flatten.
    # This is equal to only flattening lower triangular matrix.

    blocks = []
    if energy && haskey(data, :energy)
        e = ace_energy(basis, data; kwargs...)
        push!(blocks, e')
    end
    if force && haskey(data, :force)
        f = ace_forces(basis, data; kwargs...)
        tf = reinterpret.(Float64, f)
        f_bock = reduce(hcat, tf)
        push!(blocks, f_bock)
    end
    if virial && haskey(data, :virial)
        v = ace_virial(basis, data; kwargs...)
        tv = map( v ) do m
            m[SVector(1,5,9,6,3,2)]
        end
        v_block = reduce(hcat, tv)
        push!(blocks, v_block)
    end
    return reduce(vcat, blocks)
end


function ACEfit.target_vector(
    data; 
    energy=true, 
    force=true, 
    virial=true,
    energy_ref=nothing,
    kwargs...
)
    extract_virial(m::AbstractMatrix) = m[SVector(1,5,9,6,3,2)]
    extract_virial(v::AbstractVector) = extract_virial( reduce(hcat, v) ) # this could be wrong, might need transpose
    
    blocks = []
    if energy && haskey(data, :energy)
        e = data[:energy]
        if !isnothing(energy_ref)
            e -= ace_energy(energy_ref, data) |> ustrip
        end
        push!(blocks, [e])
    end
    if force && haskey(data, :force)
        tf = reduce(vcat, data[:force]) # flatten force
        push!(blocks, tf)
    end
    if virial && haskey(data, :virial)
        tv = extract_virial( data[:virial] )
        push!(blocks, tv)
    end
   return reduce(vcat, blocks)
end

function ACEfit.weight_vector(
    data;
    energy=true,
    force=true,
    virial=true, 
    energy_default_weight=1, 
    force_default_weight=1, 
    virial_default_weight=1,
    kwargs...
)
    blocks = []
    if energy && haskey(data, :energy)
        e = haskey(data, :energy_weight) ? data[:energy] : energy_default_weight
        we = e / (sqrt ∘ length)(data)
        push!(blocks, [we])
    end
    if force && haskey(data, :force)
        f = haskey(data, :force_weight) ? data[:force_weight] : force_default_weight
        wf = f * ones(3*length(data))
        push!(blocks, wf)
    end
    if virial && haskey(data, :virial)
        v = haskey(data, :virial_weight) ? data[:virial_weight] : virial_default_weight
        wv = ( v / (sqrt ∘ length)(data) ) * ones(6)
        push!(blocks, wv)
    end
   W = reduce(vcat, blocks)
   if haskey(data, :weight)
        W *= data[:weight]
   end
   return W
end


function ACEfit.assemble_weights(data::AbstractArray; kwargs...)
    w = map( data ) do d
        ACEfit.weight_vector(d; kwargs...)
    end
    return reduce(vcat, w)
end


function ACEfit.assemble(data::AbstractArray, basis; kwargs...)
    W = Threads.@spawn ACEfit.assemble_weights(data; kwargs...)
    raw_data = Folds.map( data ) do d
        A = ACEfit.feature_matrix(d, basis; kwargs...)
        Y = ACEfit.target_vector(d; kwargs...)
        (A, Y)
    end
    A = [ a[1] for a in raw_data ]
    Y = [ a[2] for a in raw_data ]

    A_final = reduce(vcat, A)
    Y_final = reduce(vcat, Y)
    return A_final, Y_final, fetch(W)
end

end