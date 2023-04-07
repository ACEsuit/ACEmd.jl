using ACEapi
using ACE1
#using AtomsBase
#using ExtXYZ
using Test

fname_ace = joinpath(pkgdir(ACEapi), "data", "TiAl.json")
fname_xyz = joinpath(pkgdir(ACEapi), "data", "TiAl.xyz")

@testset "ACEapi.jl" begin
    # Write your tests here.
    pot = load_ace_model(fname_ace)
    pot_old = load_ace_model(fname_ace; old_format=true)
    data = read_extxyz(fname_xyz)[end]


    @test ace_energy(pot, data) ≈ ACE1.energy(pot_old, data)
    @test ace_forces(pot, data) ≈ ACE1.forces(pot_old, data)
    @test ace_virial(pot, data) ≈ ACE1.virial(pot_old, data)
end
