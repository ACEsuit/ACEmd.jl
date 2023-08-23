using ACEmd
using ACE1
using AtomsBase
using ExtXYZ
using Molly
using Unitful
using UnitfulAtomic
using Test

fname_ace = joinpath(pkgdir(ACEmd), "data", "TiAl.json")
fname_xyz = joinpath(pkgdir(ACEmd), "data", "TiAl-big.xyz")

const u_energy = ACEmd.default_energy
const u_length = ACEmd.default_length


@testset "JuLIP comparison" begin
    pot_old = load_ace_model(fname_ace; old_format=true)
    pot = load_ace_model(pot_old)
    data = read_extxyz(fname_xyz)[end]

    @test ace_energy(pot, data) ≈ ACE1.energy(pot_old, data) * u_energy
    @test all( ace_forces(pot, data) .≈ ACE1.forces(pot_old, data) * (u_energy / u_length) )
    @test ace_virial(pot, data) ≈ ACE1.virial(pot_old, data) * (u_energy * u_length)
end


@testset "AtomsBase" begin
    pot = load_ace_model(fname_ace)
    julip_data = read_extxyz(fname_xyz)[end]
    ab_data = FastSystem(ExtXYZ.Atoms(read_frame(fname_xyz)))

    @test ace_energy(pot, julip_data) ≈ ace_energy(pot, ab_data)
    @test all( ace_forces(pot, julip_data) .≈ ace_forces(pot, ab_data) )
    @test ace_virial(pot, julip_data) ≈ ace_virial(pot, ab_data)
end

@testset "Units" begin
    pot = load_ace_model(fname_ace; energy_unit=u"eV", length_unit=u"pm")
    data = FastSystem(ExtXYZ.Atoms(read_frame(fname_xyz)))

    @test unit( ace_energy(pot, data) ) == u"eV"
    @test unit( ace_forces(pot, data)[1][1] ) == u"eV" / u"pm"
    @test unit( ace_virial(pot, data)[1,1] ) == u"eV" * u"pm"

    @test unit( ace_energy(pot, data; energy_unit=u_energy) ) == u_energy
    @test unit( ace_forces(pot, data;  energy_unit=u_energy, length_unit=u_length)[1][1]) == u_energy / u_length
    @test unit( ace_virial(pot, data; energy_unit=u_energy, length_unit=u_length)[1,1] ) == u_energy * u_length

    @test unit( ACEmd.get_cutoff(pot; cutoff_unit=u"m") ) == u_length
    @test unit( ACEmd.get_cutoff(pot[2]; cutoff_unit=u"pm") ) == u"pm"
    
    pot1 = load_ace_model(fname_ace; energy_unit=u"eV", length_unit=u"pm", cutoff_unit=u"pm")
    @test unit( ACEmd.get_cutoff(pot1) ) == u"pm"

    @test_throws AssertionError load_ace_model(fname_ace; energy_unit=u"s")
    @test_throws AssertionError load_ace_model(fname_ace; length_unit=u"s")
    @test_throws AssertionError load_ace_model(fname_ace; cutoff_unit=u"s")
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
    @test all( ace_ef["forces"] .≈ F )
    
    @test ace_efv["energy"] ≈ E
    @test all( ace_efv["forces"] .≈ F )
    @test ace_efv["virial"] ≈ V

    @test all( ace_fv["forces"] .≈ F )
    @test ace_fv["virial"] ≈ V
end

@testset "Molly support" begin
    pot = load_ace_model(fname_ace)
    data = ExtXYZ.Atoms(read_frame(fname_xyz))

    sys = Molly.System(data, pot)
    
    @test ace_energy(pot, data)  ≈ Molly.potential_energy(sys)
    @test all( ace_forces(pot, data) .≈ Molly.forces(sys) )
end 
