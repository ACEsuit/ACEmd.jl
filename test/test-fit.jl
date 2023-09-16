using ACE1x
using ACEfit
using ACEmd
using ACEpotentials
using ExtXYZ
using Test

fname_train = joinpath(pkgdir(ACEmd), "data", "TiAl-train.xyz")

@testset "ACEfit extension" begin
    #data = ExtXYZ.Atoms.( ExtXYZ.read_frames(fname_train) )
    data_julip = read_extxyz( fname_train )
    data = FlexibleSystem.(data_julip)

    model = acemodel(
            elements = [:Ti, :Al],
			order = 3,
			totaldegree = 6,
			rcut = 5.5,
			Eref = [:Ti => -1586.0195, :Al => -105.5954]
    )
    basis = model.basis.BB[2]

    weights = Dict( "default" => Dict("E" => 5.0, "F" => 2.0 , "V" => 3.0 ) )
    datakeys = (energy_key = "energy", force_key = "force", virial_key = "virial")
    train_julip = [ACEpotentials.AtomsData(t; weights=weights, v_ref=model.Vref, datakeys...) for t in data_julip]
    
    A, Y, W = ACEfit.assemble(train_julip, basis)
    a, y, w = ACEfit.assemble(data, basis; energy_default_weight=5, force_default_weight=2, virial_default_weight=3, energy_ref=model.Vref)

    tol = 1e-15 # tolerance to to test assembly errors
    @test maximum(abs2,  A - a  ) < tol
    @test maximum(abs2,  Y - y  ) < tol
    @test maximum(abs2,  W - w  ) < tol
end