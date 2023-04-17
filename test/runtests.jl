using ACEapi
using ACE1
using AtomsBase
using ExtXYZ
using Test

fname_ace = joinpath(pkgdir(ACEapi), "data", "TiAl.json")
fname_xyz = joinpath(pkgdir(ACEapi), "data", "TiAl-big.xyz")


@testset "JuLIP comparison" begin
    pot = load_ace_model(fname_ace)
    pot_old = load_ace_model(fname_ace; old_format=true)
    data = read_extxyz(fname_xyz)[end]

    @test ace_energy(pot, data) ≈ ACE1.energy(pot_old, data)
    @test ace_forces(pot, data) ≈ ACE1.forces(pot_old, data)
    @test ace_virial(pot, data) ≈ ACE1.virial(pot_old, data)
end


@testset "AtomsBase" begin
    pot = load_ace_model(fname_ace)
    julip_data = read_extxyz(fname_xyz)[end]
    ab_data = FastSystem(ExtXYZ.Atoms(read_frame(fname_xyz)))

    @test ace_energy(pot, julip_data) ≈ ace_energy(pot, ab_data)
    @test ace_forces(pot, julip_data) ≈ ace_forces(pot, ab_data)
    @test ace_virial(pot, julip_data) ≈ ace_virial(pot, ab_data)
end

@testset "Combination interface" begin
    pot = load_ace_model(fname_ace)
    data = FastSystem(ExtXYZ.Atoms(read_frame(fname_xyz)))

    E = ace_energy(pot, data)
    F = ace_forces(pot, data)
    V = ace_virial(pot, data)

    ace_ef = ace_energy_forces(pot, data)
    ace_efv = ace_energy_forces_virial(pot, data)
    ace_fv = ace_forces_virial(pot, data)

    @test ace_ef["energy"] ≈ E
    @test ace_ef["forces"] ≈ F
    
    @test ace_efv["energy"] ≈ E
    @test ace_efv["forces"] ≈ F
    @test ace_efv["virial"] ≈ V

    @test ace_fv["forces"] ≈ F
    @test ace_fv["virial"] ≈ V
end