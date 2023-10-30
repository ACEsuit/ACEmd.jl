

AtomsCalculators.@generate_interface function AtomsCalculators.potential_energy(
    system,
    calculator::ACEmd.ACEpotential;
    kwargs...
)
    return ace_energy(calculator, system; kwargs...)
end



AtomsCalculators.@generate_interface function AtomsCalculators.forces(
    system,
    calculator::ACEmd.ACEpotential;
    kwargs...
)
    return ace_forces(calculator, system; kwargs...)
end



AtomsCalculators.@generate_interface function AtomsCalculators.virial(
    system,
    calculator::ACEmd.ACEpotential;
    kwargs...
)
    return ace_virial(calculator, system; kwargs...)
end