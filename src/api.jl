
function ace_energy(calc, at; domain=1:length(at), executor=ThreadedEx())
    nlist = neighborlist(at, cutoff(calc))
    Etot = Folds.sum( domain, executor ) do i
        _, R, Z = neigsz(nlist, at, i)
        ace_evaluate(calc, R, Z, _atomic_number(at,i))
    end
    return Etot
end

function ace_energy(V::OneBody, at::Atoms; domain=1:length(at), executor=nothing)
    E = sum( domain ) do i
        ACE1.evaluate(V, chemical_symbol(at.Z[i]) )
    end
    return E
end


for ace_method in [ :ace_energy, :ace_forces, :ace_virial ]
    @eval begin
        function $ace_method(calc::AbstractArray, at::Atoms; domain=1:length(at), executor=ThreadedEx())
            tmp = asyncmap( calc ) do V
                $ace_method(V, at; domain=domain, executor=executor)
            end
            return sum(tmp)
        end
    end
end




## forces

function ace_forces(V, at::Atoms; domain=1:length(at), executor=ThreadedEx())
    nlist = neighbourlist(at, cutoff(V))
    F = Folds.sum( domain, executor ) do i
        j, R, Z = neigsz(nlist, at, i)
        _, tmp = ace_evaluate_d(V, R, Z, at.Z[i])

        #TODO make this faster
        f = zeros(eltype(tmp.dV), length(at))
        for k in eachindex(j)
            f[j[k]] -= tmp.dV[k]
            f[i]    += tmp.dV[k]
        end
        f
    end
    return F
end


function ace_forces(::OneBody, at::Atoms; kwargs...)
    T = (eltype ∘ eltype)(at.X)
    F = [ SVector{3}( zeros(T, 3) ) for i in 1:length(at) ]
    return F
end


## virial

function ace_virial(V, at::Atoms; domain=1:length(at), executor=ThreadedEx())
    nlist = neighbourlist(at, cutoff(V))
    vir = Folds.sum( domain, executor ) do i
        j, R, Z = neigsz(nlist, at, i)
        _, tmp = ace_evaluate_d(V, R, Z, at.Z[i])
        site_virial = -sum( zip(R, tmp.dV) ) do (Rⱼ, dVⱼ)
            dVⱼ * Rⱼ'
        end
        site_virial
    end
    return vir
end

function ace_virial(::OneBody, at::Atoms; kwargs...)
    T = (eltype ∘ eltype)(at.X)
    return SMatrix{3,3}(zeros(T, 3,3))
end