
function neigsz!(tmp, nlist::PairList, at::Atoms, i::Integer)
    # from JuLIP
    j, R = neigs!(tmp.R, nlist, i)
    Z = tmp.Z
    for n in eachindex(j)
       Z[n] = at.Z[j[n]]
    end
    return j, R, (@view Z[1:length(j)])
 end
 
 function neigsz(nlist::PairList, at::Atoms, i::Integer)
    # from JuLIP
    j, R = NeighbourLists.neigs(nlist, i)
    return j, R, at.Z[j]
 end


 function load_ace_model(fname; old_format=false)
    pot_tmp = load_dict(fname)["IP"]
    pot = read_dict(pot_tmp)
    if old_format
        return pot
    else
        return pot.components
    end
 end


 ## CellListMap

neighborlist(at::Atoms, cutoff) = ACE1.neighbourlist(at, cutoff)

function neighborlist(ab::AbstractSystem, cutoff)
    function push_pair!(i, j, x, y, d2, pairs, cutoff) 
        d = sqrt(d2)
        if d < cutoff
            push!(pairs, (i, j, ustrip(y-x)))
        end
        return pairs
    end
    function reduce_pairs(pairs, pairs_threaded)
        for i in eachindex(pairs_threaded)
            append!(pairs, pairs_threaded[i])
        end
        return pairs
    end
    ucutoff = cutoff*u"Å"  # ACE cutoff is in Å
    tmp = bounding_box(ab)
    cell = [ tmp[i][i] for i in eachindex(tmp) ]

    box = Box(cell, ucutoff)
    
    #TODO allow dynamic types here
    pairs = Tuple{Int,Int, SVector{3, Float64}}[]

    cl = CellList(position(ab), box, parallel=true)

    map_pairwise!(
        (x, y, i, j, d2, pairs) -> push_pair!(i, j, x, y, d2, pairs, ucutoff),
        pairs,box,cl,
        reduce=reduce_pairs,
        parallel=true
    )
    return pairs
end


function neigsz(list, ab::AbstractSystem, i)
    tmp = filter( x-> x[1] == i || x[2] == i, list )

    R = map( tmp ) do x
        x[1] == i ? x[3] : -x[3]
    end

    j = map( tmp ) do x
        x[1] == i ? x[2] : x[1]
    end

    Z = map( j ) do k
        _atomic_number(ab,k)
    end
    return j, R, Z
end

_atomic_number(at::Atoms, i)          = at.Z[i]
_atomic_number(ab::AbstractSystem, i) = ACE1.AtomicNumber(AtomsBase.atomic_number(ab,i))