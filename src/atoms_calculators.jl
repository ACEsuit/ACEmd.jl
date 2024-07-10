
AtomsCalculators.energy_unit(ace::ACEpotential) = ace.energy_unit
AtomsCalculators.length_unit(ace::ACEpotential) = ace.length_unit


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


function AtomsCalculators.energy_forces(system, calculator::ACEpotential; kwargs...)
    return ace_energy_forces(calculator, system;  kwargs...)
end


function AtomsCalculators.energy_forces_virial(system, calculator::ACEpotential; kwargs...)
    return ace_energy_forces_virial(calculator, system;  kwargs...)
end