"""
    ace_energy(potential, atoms, Kwargs)

Calculates ACE potential energy for atomic system.
The `atoms` object needs to be in `AtomsBase` compatable format.
The returned energy has a unit as defined by `Unitful`.

Parallel execution is done with Transducers.jl and there is an option to
use different executors. Look for `ThreadedEx` or other executors for more details on how to control it. 

# Kwargs
- `domain=1:length(atoms)`  :  choose subset of atoms to which energy is calculated
- `executor=ThreadedEx()`   :  used to control multithreading using Transducers.jl
- `energy_unit`  :   used to override energy unit for the calculation
- `length_unit`  :   used to override lenght unit for the calculation
- `cutoff_unit`  :   used to override unit that cutoff radius is defined
"""
function ace_energy(calc, at; domain=1:length(at), executor=ThreadedEx(), energy_unit=default_energy, kwargs...)
    nlist = neighborlist(at, get_cutoff(calc); storelist=false)
    Etot = Folds.sum( domain, executor ) do i
        _, R, Z = neigsz(nlist, at, i)
        ace_evaluate(calc, R, Z, _atomic_number(at,i))
    end
    return Etot * energy_unit
end

function ace_energy(V::ACE1.OneBody, at::ACE1.Atoms; domain=1:length(at), energy_unit=default_energy, kwargs...)
    E = sum( domain ) do i
        ACE1.evaluate(V, ACE1.chemical_symbol(at.Z[i]) )
    end
    return E * energy_unit
end

function ace_energy(V::ACE1.OneBody, as::AbstractSystem; domain=1:length(as), energy_unit=default_energy, kwargs...)
    E = sum( domain ) do i
        ACE1.evaluate(V, atomic_symbol(as, i) )
    end
    return E * energy_unit
end

# Generate interface for array potentials
for ace_method in [ :ace_energy, :ace_forces, :ace_virial, :ace_atom_energies ]
    @eval begin
        function $ace_method(calc::AbstractArray, at;
                domain=1:length(at),
                executor=ThreadedEx(),
                ntasks=Threads.nthreads(),
                energy_unit=default_energy,
                length_unit=default_length,
                cutoff_unit=default_length,
                kwargs...
            )
            tmp = map( calc ) do V
                Threads.@spawn $ace_method(V, at;
                    domain=domain,
                    executor=executor,
                    ntasks=ntasks,
                    energy_unit=energy_unit,
                    length_unit=length_unit,
                    cutoff_unit=cutoff_unit,
                    kwargs...
                )
            end
            return sum(fetch, tmp)
        end
    end
end

# Generate interface for ACEpotential type
for ace_method in [ :ace_energy, :ace_forces, :ace_virial, :ace_atom_energies ]
    @eval begin
        function $ace_method(calc::ACEpotential, at;
                domain=1:length(at),
                executor=ThreadedEx(),
                ntasks=Threads.nthreads(),
                energy_unit=calc.energy_unit,
                length_unit=calc.length_unit,
                cutoff_unit=calc.cutoff_unit,
                kwargs...
        )
            tmp = asyncmap( calc ) do V
                $ace_method(V, at;
                    domain=domain,
                    executor=executor,
                    ntasks=ntasks,
                    energy_unit=energy_unit,
                    length_unit=length_unit,
                    cutoff_unit=cutoff_unit,
                    kwargs...
                )
            end
            return sum( tmp )
        end
    end
end




## forces

"""
    ace_forces(potential, atoms, Kwargs)

Calculates forces for ACE potential for given atomic system.
The `atoms` object needs to be in `AtomsBase` compatable format.
The returned energy has a unit as defined by `Unitful`.

Parallel execution is done with Transducers.jl and there is an option to
use different executors. Look for `ThreadedEx` or other executors for more details on how to control it.
`ntasks` parameter is used to define number of task that the calculation is divided into. 

# Kwargs
- `domain=1:length(atoms)`          :  choose subset of atoms to which energy is calculated
- `executor=ThreadedEx()`           :  used to control multithreading using Transducers.jl
- `ntasks=Threads.nthreads()`       :  how many tasks are used in the calculation
- `energy_unit`  :   used to override energy unit for the calculation
- `length_unit`  :   used to override lenght unit for the calculation
- `cutoff_unit`  :   used to override unit that cutoff radius is defined
"""
function ace_forces(
    V,
    at;
    domain=1:length(at),
    executor=ThreadedEx(),
    ntasks=Threads.nthreads(),
    energy_unit=default_energy,
    length_unit=default_length,
    cutoff_unit=default_length,
    kwargs...
)   
    nlist = neighborlist(at, get_cutoff(V; cutoff_unit=cutoff_unit) )
    F = Folds.sum( collect(chunks(domain, ntasks)), executor ) do (d, _)
        ace_forces(V, at, nlist; domain=d)
    end
    return F * (energy_unit / length_unit)
end


function ace_forces(
    V, at, nlist;
    domain=1:length(at),
    kwargs...
)   
    f = zeros(SVector{3, Float64}, length(at))
    for i in domain
        j, R, Z = neigsz(nlist, at, i)
        _, tmp = ace_evaluate_d(V, R, Z, _atomic_number(at,i))

        for k in eachindex(j)
            f[j[k]] -= tmp.dV[k]
        end
        f[i] += sum(tmp.dV)
    end
    return f
end


function ace_forces(::ACE1.OneBody, at::ACE1.Atoms; energy_unit=default_energy, length_unit=default_length, kwargs...)
    T = (eltype ∘ eltype)(at.X)
    F = [ SVector{3}( zeros(T, 3) ) * (energy_unit / length_unit) for _ in 1:length(at) ]
    return F
end

function ace_forces(::ACE1.OneBody, as::AbstractSystem; energy_unit=default_energy, length_unit=default_length, kwargs...)
    T = eltype( ustrip.( position(as, 1) )  )
    F = [ SVector{3}( zeros(T, 3) ) * (energy_unit / length_unit) for _ in 1:length(as) ]
    return F
end


## virial

"""
    ace_virial(potential, atoms, Kwargs)

Calculates virial for ACE potential for given atomic system.
The `atoms` object needs to be in `AtomsBase` compatable format.
The returned energy has a unit as defined by `Unitful`.

Parallel execution is done with Transducers.jl and there is an option to
use different executors. Look for `ThreadedEx` for more details on how to control it. 

# Kwargs
- `domain=1:length(atoms)`  :  choose subset of atoms to which energy is calculated
- `executor=ThreadedEx()`   :  used to control multithreading using Transducers.jl
- `energy_unit`  :   used to override energy unit for the calculation
- `length_unit`  :   used to override lenght unit for the calculation
- `cutoff_unit`  :   used to override unit that cutoff radius is defined
"""
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

function ace_virial(::ACE1.OneBody, at::ACE1.Atoms; energy_unit=default_energy, length_unit=default_length, kwargs...)
    T = (eltype ∘ eltype)(at.X)
    return SMatrix{3,3}(zeros(T, 3,3)) * (energy_unit * length_unit)
end

function ace_virial(::ACE1.OneBody, as::AbstractSystem; energy_unit=default_energy, length_unit=default_length, kwargs...)
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


## Individual atom energies

"""
    ace_atom_energies(potential, atoms; kwargs)

Calculates ACE potential energy for each atom.
The `atoms` object needs to be in `AtomsBase` compatable format.
The returned energy has a unit as defined by `Unitful`.

Parallel execution is done with Transducers.jl and there is an option to
use different executors. Look for `ThreadedEx` or other executors for more details on how to control it. 

# Kwargs
- `domain=1:length(atoms)`  :  choose subset of atoms to which energy is calculated
- `executor=ThreadedEx()`   :  used to control multithreading using Transducers.jl
- `energy_unit`  :   used to override energy unit for the calculation
- `length_unit`  :   used to override lenght unit for the calculation
- `cutoff_unit`  :   used to override unit that cutoff radius is defined
"""
function ace_atom_energies(calc, at; domain=1:length(at), executor=ThreadedEx(), energy_unit=default_energy, kwargs...)
    nlist = neighborlist(at, get_cutoff(calc); storelist=false)
    Etot = Folds.map( domain, executor ) do i
        _, R, Z = neigsz(nlist, at, i)
        ace_evaluate(calc, R, Z, _atomic_number(at,i))
    end
    return Etot * energy_unit
end

function ace_atom_energies(V::ACE1.OneBody, at::ACE1.Atoms; domain=1:length(at), energy_unit=default_energy, kwargs...)
    E = map( domain ) do i
        ACE1.evaluate(V, ACE1.chemical_symbol(at.Z[i]) )
    end
    return E * energy_unit
end

function ace_atom_energies(V::ACE1.OneBody, as::AbstractSystem; domain=1:length(as), energy_unit=default_energy, kwargs...)
    E = map( domain ) do i
        ACE1.evaluate(V, atomic_symbol(as, i) )
    end
    return E * energy_unit
end


## Basis calls

function ace_energy(basis::ACE1.IPBasis, at; domain=1:length(at), cutoff_unit=default_length, executor=ThreadedEx(), kwargs...)
    nlist = neighborlist(at, get_cutoff(basis; cutoff_unit=cutoff_unit) )
    E = Folds.sum( domain, executor ) do i
        j, R, Z = neigsz(nlist, at, i)
        ace_evaluate(basis, R, Z, _atomic_number(at,i))
    end
    return E
end

function ace_forces(
    basis::ACE1.IPBasis,
    at;
    domain=1:length(at),
    executor=ThreadedEx(),
    ntasks=Threads.nthreads(),
    cutoff_unit=default_length,
    kwargs...
)   
    nlist = neighborlist(at, get_cutoff(basis; cutoff_unit=cutoff_unit) )
    F = Folds.sum( collect(chunks(domain, ntasks)), executor ) do (d, _)
        ace_forces(basis, at, nlist; domain=d)
    end
    return [ Vector(f) for f in eachrow(F)]
end


function ace_forces(
    shipB::ACE1.IPBasis,
    at,
    nlist;
    domain=1:length(at)
)
    f = zeros(SVector{3, Float64}, length(shipB), length(at))
    for i in domain
        j, R, Z = neigsz(nlist, at, i)
        dB = ace_evaluate_d(shipB, R, Z, _atomic_number(at,i))

        for k in eachindex(j)
            f[:, j[k]] .-= @view dB[:, k]
        end
        for k in eachindex(j)
            f[:,i] .+= @view dB[:, k]
        end
    end
    return f
end


function ace_forces(
    pair_basis::ACE1.PolyPairBasis,
    at;
    domain=1:length(at),
    executor=ThreadedEx(),
    ntasks=Threads.nthreads(),
    cutoff_unit=default_length,
    kwargs...
)   
    nlist = neighborlist(at, get_cutoff(pair_basis; cutoff_unit=cutoff_unit) )
    pair_data = [ (i, j, R)  for (i,j,R) in pairs(nlist) if i in domain ]
    F = Folds.sum( collect(chunks(pair_data, ntasks)), executor ) do (d, _)
        ace_forces(pair_basis, at, pair_data, d)
    end
    return [ Vector(f) for f in eachrow(F)]
end


function ace_forces(
    pair_basis::ACE1.PolyPairBasis,
    at,
    pair_data,
    domain
)
    f = zeros(SVector{3, Float64}, length(pair_basis), length(at))
    tmp = ACE1.alloc_temp_d(pair_basis)
    for I in domain
        # Pairpotential evaluator from JuLIP.
        # This is temporary implementation and will be changed in the future
        i, j, R = pair_data[I]
        r = norm(R)
        Zi = _atomic_number(at, i)
        Zj = _atomic_number(at, j)
        Ii = ACE1.z2i(pair_basis, Zi)
        Ij = ACE1.z2i(pair_basis, Zj)
        dJ = tmp.dJ[Ii, Ij]
        ACE1.evaluate_d!(tmp.J[Ii, Ij], dJ, tmp.tmpd_J[Ii, Ij], pair_basis.J[Ii, Ij], r, Zi, Zj)
        idx0 = pair_basis.bidx0[Ii, Ij] 
        r_tmp = R/(2r)
        for n = 1:length(pair_basis.J[Ii, Ij])
            f[idx0 + n, i] +=  dJ[n] * r_tmp
        end
        for n = 1:length(pair_basis.J[Ii, Ij])
            f[idx0 + n, j] -= dJ[n] * r_tmp
        end
    end
    return f
end


function ace_virial(
    basis::ACE1.IPBasis,
    at;
    domain=1:length(at),
    executor=ThreadedEx(),
    cutoff_unit=default_length,
    kwargs...
)
    nlist = neighborlist(at, get_cutoff(basis; cutoff_unit=cutoff_unit) )
    vir_tot = Folds.sum( domain, executor ) do i
        vir = zeros(SMatrix{3,3,Float64}, length(basis))
        j, R, Z = neigsz(nlist, at, i)
        dB = ace_evaluate_d(basis, R, Z, _atomic_number(at,i))
        for ib in 1:length(basis)
            site_virial = -sum( zip( R, view(dB, ib,:) ) ) do (Rⱼ, dVⱼ)
                dVⱼ * Rⱼ'
            end
            vir[ib] += site_virial
        end
        vir
    end
    return vir_tot
 end

function ace_virial(
    pair_basis::ACE1.PolyPairBasis,
    at;
    domain=1:length(at),
    executor=ThreadedEx(),
    ntasks=Threads.nthreads(),
    cutoff_unit=default_length,
    kwargs...
)   
    nlist = neighborlist(at, get_cutoff(pair_basis; cutoff_unit=cutoff_unit) )
    pair_data = [ (i, j, R)  for (i,j,R) in pairs(nlist) if i in domain ]
    vir = Folds.sum( collect(chunks(pair_data, ntasks)), executor ) do (d, _)
        ace_virial(pair_basis, at, pair_data, d)
    end
    return vir
end


function ace_virial(
    pair_basis::ACE1.PolyPairBasis,
    at,
    pair_data,
    domain
)
    vir = zeros(SMatrix{3,3,Float64}, length(pair_basis))
    tmp = ACE1.alloc_temp_d(pair_basis)
    for I in domain
        # Pairpotential evaluator from JuLIP.
        # This is temporary implementation and will be changed in the future
        i, j, R = pair_data[I]
        r = norm(R)
        Zi = _atomic_number(at, i)
        Zj = _atomic_number(at, j)
        Ii = ACE1.z2i(pair_basis, Zi)
        Ij = ACE1.z2i(pair_basis, Zj)
        dJ = tmp.dJ[Ii, Ij]
        ACE1.evaluate_d!(tmp.J[Ii, Ij], dJ, tmp.tmpd_J[Ii, Ij], pair_basis.J[Ii, Ij], r, Zi, Zj)
        idx0 = pair_basis.bidx0[Ii, Ij] 
        r_tmp =  R * R'
        for n = 1:length(pair_basis.J[Ii, Ij])
            vir[idx0 + n] -= (dJ[n]/(2r)) * r_tmp
        end
    end
    return vir
end