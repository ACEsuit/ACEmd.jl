
function neigsz!(tmp, nlist::PairList, at::Atoms, i::Integer)
    # from JuLIP
    j, R = neigs!(tmp.R, nlist, i)
    Z = tmp.Z
    for n in eachindex(j)
       Z[n] = at.Z[j[n]]
    end
    return j, R, (@view Z[1:length(j)])
end
 
function neigsz(nlist::PairList, at, i::Integer)
    j, R = NeighbourLists.neigs(nlist, i)
    Z = map(j) do k
        _atomic_number(at, k)
    end
    return j, R, Z 
end


function load_ace_model(fname; old_format=false)
    pot_tmp = load_dict(fname)["IP"]
    pot = read_dict(pot_tmp)
    if old_format
        return pot
    else
        return ACEpotential(pot.components)
    end
end


function neighborlist(ab::AbstractSystem, cutoff; kwargs...)
    cell = ustrip.( u"Å", hcat( bounding_box(ab)... ) )
    pbc = map( boundary_conditions(ab) ) do x
        x == Periodic()
    end
    r = map( 1:length(ab)) do i 
        ustrip.(u"Å", position(ab,i))
    end
    nlist = PairList(r, cutoff, cell, pbc; int_type=Int)
    return nlist
end



 ## CellListMap

neighborlist(at::Atoms, cutoff; kwargs...) = ACE1.neighbourlist(at, cutoff; kwargs...)


function neighborlist_clm(ab::AbstractSystem, cutoff; kwargs...)
    function push_pair!(i, j, x, y, d2, pairs, cutoff) 
        d = sqrt(d2)
        if d < cutoff
            push!(pairs, (i, j, ustrip.(y-x)))
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


function neigsz(list, ab, i)
    T = typeof( list[begin][3] )
    R = T[]
    j = Int[]
    Z = JuLIP.AtomicNumber[]

    for x in list
        if x[1] == i
            push!(R, x[3])
            push!(j, x[2])
            push!(Z, _atomic_number(ab,x[2]))
        elseif x[2] == i
            push!(R, -x[3])
            push!(j, x[1])
            push!(Z, _atomic_number(ab,x[1]))
        end
    end
    return j, R, Z
end



_atomic_number(at::Atoms, i)          = at.Z[i]
_atomic_number(ab::AbstractSystem, i) = ACE1.AtomicNumber(AtomsBase.atomic_number(ab,i))