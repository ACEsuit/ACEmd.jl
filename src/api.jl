import Base.Threads.@spawn


@inline function ace_evaluate!(tmp, calc, R::AbstractVector, species::AbstractVector, species0 )
    return ACE1.evaluate!(tmp, calc, R, species, species0)
end

function ace_evaluate(calc, R::AbstractVector, species::AbstractVector, species0)
    tmp = ACE1.alloc_temp(calc, length(R))
    return ace_evaluate!(tmp, calc, R, species, species0)
end


function energy_nonthreaded!(tmp, calc, at::Atoms; domain=1:length(at))
    # tmp = ACE1.alloc_temp(calc, at)
    nlist = neighbourlist(at, cutoff(calc))
    Etot = sum( domain ) do i
        _, R, Z = neigsz!(tmp, nlist, at, i)
        ace_evaluate!(tmp, calc, R, Z, at.Z[i]) 
    end
    return Etot
end 


function energy_nonthreaded(calc, at::Atoms; domain=1:length(at))
    nlist = neighbourlist(at, cutoff(calc))
    Etot = sum( domain ) do i
        _, R, Z = neigsz(nlist, at, i)
        ace_evaluate(calc, R, Z, at.Z[i])
    end
    return Etot
end

function energy_nonthreaded(calc, at::Atoms, nlist; domain=1:length(at))
    #nlist = neighbourlist(at, cutoff(calc))
    Etot = sum( domain ) do i
        _, R, Z = neigsz(nlist, at, i)
        ace_evaluate(calc, R, Z, at.Z[i])
    end
    return Etot
end


function energy(calc, at::Atoms; domain=1:length(at), executor=ThreadedEx())
    nlist = neighbourlist(at, cutoff(calc))
    Etot = Folds.sum( domain, executor ) do i
        _, R, Z = neigsz(nlist, at, i)
        ace_evaluate(calc, R, Z, at.Z[i])
    end
    return Etot
end

function energy_floops(calc, at::Atoms; domain=1:length(at), executor=ThreadedEx())
    nlist = neighbourlist(at, cutoff(calc))
    @floop executor for i in domain
        _, R, Z = neigsz(nlist, at, i)
        @reduce (Etot += ace_evaluate(calc, R, Z, at.Z[i]) )
    end
    return Etot
end


function energy_tasks(calc, at::Atoms; ntasks=1)
    nlist = neighbourlist(at, cutoff(calc))
    Δ = (Int ∘ floor)( length(at) / ntasks )
    tasks = map( 1:ntasks ) do i
        s = 1+(i-1)*Δ : i*Δ
        @spawn energy_nonthreaded(calc, at::Atoms, nlist; domain=s)
    end
    Etot = sum(tasks) do t
        fetch(t)
    end
    return Etot
end


