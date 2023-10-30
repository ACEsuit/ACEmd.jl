using ACEmd
using ACE1
using ACE1x
using ACEfit
using AtomsBase
using AtomsCalculators.AtomsCalculatorsTesting
using ExtXYZ
using Molly
using Unitful
using UnitfulAtomic
using Test

fname_ace = joinpath(pkgdir(ACEmd), "data", "TiAl.json")
fname_xyz = joinpath(pkgdir(ACEmd), "data", "TiAl-big.xyz")
fname_train = joinpath(pkgdir(ACEmd), "data", "TiAl-train.xyz")


const u_energy = ACEmd.default_energy
const u_length = ACEmd.default_length


@testset "JuLIP comparison" begin
    pot_old = load_ace_model(fname_ace; old_format=true)
    pot = load_ace_model(pot_old)
    data = read_extxyz(fname_xyz)[end]

    @test ace_energy(pot, data) ≈ ACE1.energy(pot_old, data) * u_energy
    @test all( ace_forces(pot, data) .≈ ACE1.forces(pot_old, data) * (u_energy / u_length) )
    @test ace_virial(pot, data) ≈ ACE1.virial(pot_old, data) * u_energy

    @testset "Basis evaluations" begin
        model = acemodel(
            elements = [:Ti, :Al],
			order = 3,
			totaldegree = 6,
			rcut = 5.5,
			Eref = [:Ti => -1586.0195, :Al => -105.5954]
        )
        basis = model.basis
        # ACE basis
        @test all( ace_energy(basis.BB[2], data) .≈ ACE1.energy(basis.BB[2], data)  )
        @test (all∘map)( ace_forces(basis.BB[2], data), ACE1.forces(basis.BB[2], data)) do a,b
            all( a .≈ b  )
        end
        @test (all∘map)( ace_virial(basis.BB[2], data), ACE1.virial(basis.BB[2], data)) do a,b
            all( a .≈ b  )
        end

        # pair potential basis
        @test all( ace_energy(basis.BB[1], data) .≈ ACE1.energy(basis.BB[1], data)  )
        @test (all∘map)( ace_forces(basis.BB[1], data), ACE1.forces(basis.BB[1], data)) do a,b
            all( a .≈ b  )
        end
        @test (all∘map)( ace_virial(basis.BB[1], data), ACE1.virial(basis.BB[1], data)) do a,b
            all( a .≈ b  )
        end

    end
end


@testset "AtomsBase" begin
    pot = load_ace_model(fname_ace)
    julip_data = read_extxyz(fname_xyz)[end]
    ab_data = FastSystem(ExtXYZ.Atoms(read_frame(fname_xyz)))

    @test ace_energy(pot, julip_data) ≈ ace_energy(pot, ab_data)
    @test all( ace_forces(pot, julip_data) .≈ ace_forces(pot, ab_data) )
    @test ace_virial(pot, julip_data) ≈ ace_virial(pot, ab_data)
    F = ace_forces(pot, ab_data)
    @test typeof(F) <: Array
    @test ace_energy(pot, ab_data) ≈ sum( ace_atom_energies(pot, ab_data) )

    @test all( ace_energy(pot[3].pibasis, ab_data) .≈ ace_energy(pot[3].pibasis, julip_data)  )

    @testset "Basis evaluations" begin
        model = acemodel(
            elements = [:Ti, :Al],
			order = 3,
			totaldegree = 6,
			rcut = 5.5,
			Eref = [:Ti => -1586.0195, :Al => -105.5954]
        )
        basis = model.basis
        # ACE basis
        @test all( ace_energy(basis.BB[2], ab_data) .≈ ace_energy(basis.BB[2], julip_data)  )
        @test (all∘map)( ace_forces(basis.BB[2], ab_data), ace_forces(basis.BB[2], julip_data) ) do a,b
            all( a .≈ b  )
        end
        @test (all∘map)( ace_virial(basis.BB[2], ab_data), ace_virial(basis.BB[2], julip_data) ) do a,b
            all( a .≈ b  )
        end

        # pair potential basis
        @test all( ace_energy(basis.BB[1], ab_data) .≈ ace_energy(basis.BB[1], julip_data)  )
        @test (all∘map)( ace_forces(basis.BB[1], ab_data), ace_forces(basis.BB[1], julip_data) ) do a,b
            all( a .≈ b  )
        end
        @test (all∘map)( ace_virial(basis.BB[1], ab_data), ace_virial(basis.BB[1], julip_data) ) do a,b
            all( a .≈ b  )
        end

    end
end

@testset "Units" begin
    pot = load_ace_model(fname_ace; energy_unit=u"eV", length_unit=u"pm")
    data = FastSystem(ExtXYZ.Atoms(read_frame(fname_xyz)))

    @test unit( ace_energy(pot, data) ) == u"eV"
    @test unit( ace_forces(pot, data)[1][1] ) == u"eV" / u"pm"
    @test unit( ace_virial(pot, data)[1,1] ) == u"eV"

    @test unit( ace_energy(pot, data; energy_unit=u_energy) ) == u_energy
    @test unit( ace_forces(pot, data;  energy_unit=u_energy, length_unit=u_length)[1][1]) == u_energy / u_length
    @test unit( ace_virial(pot, data; energy_unit=u_energy, length_unit=u_length)[1,1] ) == u_energy

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

    @test ace_ef[:energy] ≈ E
    @test all( ace_ef[:forces] .≈ F )
    
    @test ace_efv[:energy] ≈ E
    @test all( ace_efv[:forces] .≈ F )
    @test ace_efv[:virial] ≈ V

    @test all( ace_fv[:forces] .≈ F )
    @test ace_fv[:virial] ≈ V
end

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


@testset "ACEfit extension" begin
    data = ExtXYZ.load(fname_train)
    #data_julip = read_extxyz( fname_train )

    model = acemodel(
            elements = [:Ti, :Al],
			order = 3,
			totaldegree = 6,
			rcut = 5.5,
			Eref = [:Ti => -1586.0195, :Al => -105.5954]
    )
    basis = model.basis

    weights = Dict( "default" => Dict("E" => 5.0, "F" => 2.0 , "V" => 3.0 ) )
    datakeys = (energy_key = "energy", force_key = "force", virial_key = "virial")
    # AtomsData needs ACEpotentials and ACEpotentials uses ACEmd this can lead to issues in testing
    # so skipping test here.
    # Use test-fit.jl tests this use that to check that tests are functioning.
    # include( joinpath(pkgdir(ACEmd), "test", "test-fit.jl") )
    #train_julip = [ACEpotentials.AtomsData(t; weights=weights, v_ref=model.Vref, datakeys...) for t in data_julip]
    
    #A, Y, W = ACEfit.assemble(train_julip, basis)
    a, y, w = ACEfit.assemble(data, basis; energy_default_weight=5, force_default_weight=2, virial_default_weight=3, energy_ref=model.Vref)

    #tol = 1e-15 # tolerance to to test assembly errors
    #@test maximum(abs2,  A - a  ) < tol
    #@test maximum(abs2,  Y - y  ) < tol
    #@test maximum(abs2,  W - w  ) < tol
    @test size(a,1) == length(y) == length(w)
    @test size(a,2) == length(basis)

    # change virial key, which leads to virials being skipped
    a1, y1, w1 = ACEfit.assemble(
        data, 
        basis; 
        energy_default_weight=5, 
        force_default_weight=2, 
        virial_default_weight=3,
        virial_key=:false_virial, 
        energy_ref=model.Vref
    )
    @test size(a1,1) == length(y1) == length(w1)
    @test size(a,1) != size(a1,1)
    @test length(y) != length(y1)
    @test length(w) != length(w1)
end


@testset "AtomsCalculators interface" begin
    pot = load_ace_model(fname_ace)
    data = ExtXYZ.load(fname_xyz)

    test_potential_energy(data, pot)
    test_forces(data, pot)
    test_virial(data, pot)
end