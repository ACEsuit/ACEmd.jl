
struct ACEpotential
    potentials::Vector{AbstractCalculator}
end


function Base.iterate(acep::ACEpotential, state::Int=1)
    if 0 < state <= length( acep )
        return acep.potentials[state], state + 1
    else
        return nothing
    end
end

Base.length(acep::ACEpotential) = length(acep.potentials)

Base.getindex(acep::ACEpotential, i) = acep.potentials[i]

function Base.show(io::IO, ::MIME"text/plain", acep::ACEpotential)
    print(io, "ACE potential with ", length(acep), " subpotentials")
end