
struct ACEpotential{TE,TL}
    potentials::Vector{AbstractCalculator}
    energy_unit::TE
    length_unit::TL
    function ACEpotential(potentials; energy_unit=u"hartree", lenght_unit=u"Å")
        new{typeof(energy_unit), typeof(lenght_unit)}(potentials, energy_unit, lenght_unit)
    end
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