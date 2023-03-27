

function energy_nonthreaded!(tmp, calc, at::Atoms; domain=1:length(at))
    nlist = neighbourlist(at, cutoff(calc))
    Etot = sum( domain ) do i
        _, R, Z = neigsz!(tmp, nlist, at, i)
        evaluate!(tmp, calc, R, Z, at.Z[i])
    end
    return Etot
end