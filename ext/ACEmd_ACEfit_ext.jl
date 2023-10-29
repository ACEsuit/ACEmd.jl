module ACEmd_ACEfit_ext

using ACEfit
using ACEmd
using Distributed
using StaticArrays

function ACEfit.feature_matrix(
    data,
    basis; 
    energy=true, 
    force=true, 
    virial=true,
    energy_key=:energy,
    force_key=:force,
    virial_key=:virial,
    kwargs...)
    # Basis functions are on different collumns.
    # Energy is on fist row.
    # Force is flattened on several rows, so that each basis function is on same collumns.
    # Virial is in practice triangular matrix, so first remove double values and then flatten.
    # This is equal to only flattening lower triangular matrix.

    blocks = []
    if energy && haskey(data, energy_key)
        e = ace_energy(basis, data; kwargs...)
        push!(blocks, e')
    end

    # check if both force and virial are to be calculated and if so use 
    # optimized call for them
    cal_f = force && ( haskey(data, force_key) || hasatomkey(data, force_key) )
    cal_v = virial && haskey(data, virial_key)
    if cal_f && cal_v
        F_and_V = ace_forces_virial(basis, data; kwargs...)

        f = F_and_V[:forces]
        tf = reinterpret.(Float64, f)
        f_bock = reduce(hcat, tf)
        push!(blocks, f_bock)

        v = F_and_V[:virial]
        tv = map( v ) do m
            m[SVector(1,5,9,6,3,2)]
        end
        v_block = reduce(hcat, tv)
        push!(blocks, v_block)
        cal_f = false
        cal_v = false
    end

    if cal_f
        f = ace_forces(basis, data; kwargs...)
        tf = reinterpret.(Float64, f)
        f_bock = reduce(hcat, tf)
        push!(blocks, f_bock)
    end
    if cal_v
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
    energy_key=:energy,
    force_key=:force,
    virial_key=:virial,
    energy_ref=nothing,
    kwargs...
)
    extract_virial(m::AbstractMatrix) = m[SVector(1,5,9,6,3,2)]
    extract_virial(v::AbstractVector) = extract_virial( reduce(hcat, v)' ) # this could be wrong, might need transpose
    
    blocks = []
    if energy && haskey(data, energy_key)
        e = data[energy_key] |> ustrip
        if !isnothing(energy_ref)
            e -= ace_energy(energy_ref, data) |> ustrip
        end
        push!(blocks, [e])
    end
    if force
        if haskey(data, force_key)
            tf = reduce(vcat, data[force_key]) # flatten force
            push!(blocks, tf)
        elseif  hasatomkey(data, force_key)
            tf = reduce(vcat, [ x[force_key] for x in data ] )
            push!(blocks, tf)
        end
    end
    if virial && haskey(data, virial_key)
        tv = extract_virial( data[virial_key] )
        push!(blocks, tv)
    end
   return reduce(vcat, blocks)
end

function ACEfit.weight_vector(
    data;
    energy=true,
    force=true,
    virial=true,
    energy_key=:energy,
    force_key=:force,
    virial_key=:virial, 
    energy_default_weight=1, 
    force_default_weight=1, 
    virial_default_weight=1,
    kwargs...
)
    blocks = []
    if energy && haskey(data, energy_key)
        e = haskey(data, :energy_weight) ? data[:energy_weight] : energy_default_weight
        we = e / (sqrt ∘ length)(data)
        push!(blocks, [we])
    end
    if force && ( haskey(data, force_key) || hasatomkey(data, force_key) )
        f = haskey(data, :force_weight) ? data[:force_weight] : force_default_weight
        wf = f * ones(3*length(data))
        push!(blocks, wf)
    end
    if virial && haskey(data, virial_key)
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
    #raw_data = Folds.map( data ) do d # this will bug out
    raw_data = pmap( data ) do d
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