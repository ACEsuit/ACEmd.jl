
function ace_energy(calc, at; domain=1:length(at), executor=ThreadedEx(), energy_unit=default_energy, kwargs...)
    nlist = neighborlist(at, get_cutoff(calc); storelist=false)
    Etot = Folds.sum( domain, executor ) do i
        _, R, Z = neigsz(nlist, at, i)
        ace_evaluate(calc, R, Z, _atomic_number(at,i))
    end
    return Etot * energy_unit
end

function ace_energy(V::OneBody, at::Atoms; domain=1:length(at), energy_unit=default_energy, kwargs...)
    E = sum( domain ) do i
        ACE1.evaluate(V, chemical_symbol(at.Z[i]) )
    end
    return E * energy_unit
end

function ace_energy(V::OneBody, as::AbstractSystem; domain=1:length(as), energy_unit=default_energy, kwargs...)
    E = sum( domain ) do i
        ACE1.evaluate(V, atomic_symbol(as, i) )
    end
    return E * energy_unit
end

# Generate interface for array potentials
for ace_method in [ :ace_energy, :ace_forces, :ace_virial ]
    @eval begin
        function $ace_method(calc::AbstractArray, at;
                domain=1:length(at),
                executor=ThreadedEx(),
                energy_unit=default_energy,
                length_unit=default_length,
                cutoff_unit=default_length,
                kwargs...
            )
            tmp = asyncmap( calc ) do V
                $ace_method(V, at;
                    domain=domain,
                    executor=executor,
                    energy_unit=energy_unit,
                    length_unit=length_unit,
                    cutoff_unit=cutoff_unit,
                    kwargs...
                )
            end
            return sum(tmp)
        end
    end
end

# Generate interface for ACEpotential type
for ace_method in [ :ace_energy, :ace_forces, :ace_virial ]
    @eval begin
        function $ace_method(calc::ACEpotential, at;
                domain=1:length(at),
                executor=ThreadedEx(),
                energy_unit=calc.energy_unit,
                length_unit=calc.length_unit,
                cutoff_unit=calc.cutoff_unit,
                kwargs...
            )
            tmp = asyncmap( calc ) do V
                $ace_method(V, at;
                    domain=domain,
                    executor=executor,
                    energy_unit=energy_unit,
                    length_unit=length_unit,
                    cutoff_unit=cutoff_unit,
                    kwargs...
                )
            end
            return sum(tmp)
        end
    end
end




## forces

function ace_forces(V, at;
        domain=1:length(at),
        executor=ThreadedEx(),
        energy_unit=default_energy,
        length_unit=default_length,
        cutoff_unit=default_length,
        kwargs...
    )
    nlist = neighborlist(at, get_cutoff(V; cutoff_unit=cutoff_unit) )
    F = Folds.sum( domain, executor ) do i
        j, R, Z = neigsz(nlist, at, i)
        _, tmp = ace_evaluate_d(V, R, Z, _atomic_number(at,i))

        #TODO make this faster
        f = zeros(eltype(tmp.dV), length(at))
        for k in eachindex(j)
            f[j[k]] -= tmp.dV[k]
            f[i]    += tmp.dV[k]
        end
        f
    end
    return F * (energy_unit / length_unit)
end


function ace_forces(::OneBody, at::Atoms; energy_unit=default_energy, length_unit=default_length, kwargs...)
    T = (eltype ∘ eltype)(at.X)
    F = [ SVector{3}( zeros(T, 3) ) * (energy_unit / length_unit) for _ in 1:length(at) ]
    return F
end

function ace_forces(::OneBody, as::AbstractSystem; energy_unit=default_energy, length_unit=default_length, kwargs...)
    T = eltype( ustrip.( position(as, 1) )  )
    F = [ SVector{3}( zeros(T, 3) ) * (energy_unit / length_unit) for _ in 1:length(as) ]
    return F
end


## virial

function ace_virial(V, at;
        domain=1:length(at),
        executor=ThreadedEx(),
        energy_unit=default_energy,
        length_unit=default_length,
        cutoff_unit=default_length,
        kwargs...
    )
    nlist = neighborlist(at, get_cutoff(V; cutoff_unit=cutoff_unit) )
    vir = Folds.sum( domain, executor ) do i
        j, R, Z = neigsz(nlist, at, i)
        _, tmp = ace_evaluate_d(V, R, Z, _atomic_number(at,i))
        site_virial = -sum( zip(R, tmp.dV) ) do (Rⱼ, dVⱼ)
            dVⱼ * Rⱼ'
        end
        site_virial
    end
    return vir * (energy_unit * length_unit)
end

function ace_virial(::OneBody, at::Atoms; energy_unit=default_energy, length_unit=default_length, kwargs...)
    T = (eltype ∘ eltype)(at.X)
    return SMatrix{3,3}(zeros(T, 3,3)) * (energy_unit * length_unit)
end

function ace_virial(::OneBody, as::AbstractSystem; energy_unit=default_energy, length_unit=default_length, kwargs...)
    T = eltype( ustrip.( position( as[begin] ) )  )
    return SMatrix{3,3}(zeros(T, 3,3)) * (energy_unit * length_unit)
end


## Combinations
# these will be optimized later

function ace_energy_forces(pot, data; kwargs...)
    E = ace_energy(pot, data; kwargs...)
    F = ace_forces(pot, data; kwargs...)
    return Dict("energy"=>E, "forces"=>F)
end


function ace_energy_forces_virial(pot, data; kwargs...)
    E = ace_energy(pot, data; kwargs...)
    F = ace_forces(pot, data; kwargs...)
    V = ace_virial(pot, data; kwargs...)
    return Dict("energy"=>E, "forces"=>F, "virial"=>V)
end

function ace_forces_virial(pot, data; kwargs...)
    F = ace_forces(pot, data; kwargs...)
    V = ace_virial(pot, data; kwargs...)
    return Dict("forces"=>F, "virial"=>V)
end