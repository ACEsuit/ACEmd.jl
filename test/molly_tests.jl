# This is not run at the moment due to Molly having only AtomsCalculators v0.1 support.
# The plan is also to remove the extension in the future, but saving test here for now.

using ACEmd
using ExtXYZ
using Molly
using Test

fname_ace = joinpath(pkgdir(ACEmd), "data", "TiAl.json")
fname_xyz = joinpath(pkgdir(ACEmd), "data", "TiAl-big.xyz")
fname_train = joinpath(pkgdir(ACEmd), "data", "TiAl-train.xyz")

@testset "Molly support" begin
    pot = load_ace_model(fname_ace)
    data = ExtXYZ.Atoms(read_frame(fname_xyz))

    sys = Molly.System(data, pot)
    
    @test ace_energy(pot, data)  ≈ Molly.potential_energy(sys)
    @test all( ace_forces(pot, data) .≈ Molly.forces(sys) )

    simulator = VelocityVerlet(
        dt=1.0u"fs",
        coupling=AndersenThermostat(300u"K", 1.0u"ps"),
    )
    simulate!(sys, simulator, 10)
end 