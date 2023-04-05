
function ace_energy(calc, at::Atoms; domain=1:length(at), executor=ThreadedEx())
    nlist = neighbourlist(at, cutoff(calc))
    Etot = Folds.sum( domain, executor ) do i
        _, R, Z = neigsz(nlist, at, i)
        ace_evaluate(calc, R, Z, at.Z[i])
    end
    return Etot
end

function ace_energy(calc::AbstractArray, at::Atoms; domain=1:length(at), executor=ThreadedEx())
    E = asyncmap( calc ) do V
        ace_energy(V, at; domain=domain, executor=executor)
    end
    return sum(E)
end

function ace_energy(V::OneBody, at::Atoms; domain=1:length(at), executor=nothing)
    E = sum( domain ) do i
        ACE1.evaluate(V, chemical_symbol(at.Z[i]) )
    end
    return E
 end





## forces

function ace_forces(V, at::Atoms; domain=1:length(at), executor=ThreadedEx())
    nlist = neighbourlist(at, cutoff(V))
    F = Folds.sum( domain, executor ) do i
        j, R, Z = neigsz(nlist, at, i)
        _, tmp = ace_evaluate_d(V, R, Z, at.Z[i])
        f = sum(tmp.dV) # = F[i]
        s = SparseVector( length(at), collect(j), tmp.dV )
        ss = SparseVector( length(at), [i], [f] )
        ss - s
    end
    return Vector( F )
end


function ace_virial(V, at::Atoms; domain=1:length(at), executor=ThreadedEx())
    nlist = neighbourlist(at, cutoff(V))
    vir = Folds.sum( domain, executor ) do i
        j, R, Z = neigsz(nlist, at, i)
        _, tmp = ace_evaluate_d(V, R, Z, at.Z[i])
    end
end