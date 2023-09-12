using BenchmarkTools
using ACEmd
using ACE1
using ACE1x
using ExtXYZ
using AtomsBase

SUITE = BenchmarkGroup()

fname_ace = joinpath(pkgdir(ACEmd), "data", "TiAl.json")
fname_xyz = joinpath(pkgdir(ACEmd), "data", "TiAl-big.xyz")
fname_xyz_huge = joinpath(pkgdir(ACEmd), "data", "TiAl-huge.xyz")

pot = load_ace_model(fname_ace)
pot_julip = load_ace_model(fname_ace; old_format=true)
data_julip = read_extxyz(fname_xyz)[end]
data_julip_huge = read_extxyz(fname_xyz_huge)[end]
data = FastSystem(ExtXYZ.Atoms(read_frame(fname_xyz)))
data_huge = FastSystem(ExtXYZ.Atoms(read_frame(fname_xyz_huge)))

model = acemodel(
    elements = [:Ti, :Al],
	order = 3,
	totaldegree = 6,
	rcut = 5.5,
	Eref = [:Ti => -1586.0195, :Al => -105.5954]
)

basis = model.basis


SUITE["JuLIP"] = BenchmarkGroup([],
        "energy" => BenchmarkGroup(),
        "forces" => BenchmarkGroup(),
        "virial" => BenchmarkGroup()
)
SUITE["AtomsBase"] = BenchmarkGroup([],
        "energy" => BenchmarkGroup(),
        "forces" => BenchmarkGroup(),
        "virial" => BenchmarkGroup()
)

SUITE["JuLIP"]["energy"]["big"] = @benchmarkable ACE1.energy($pot_julip, $data_julip)
SUITE["JuLIP"]["energy"]["huge"] = @benchmarkable ACE1.energy($pot_julip, $data_julip_huge)
SUITE["JuLIP"]["energy"]["basis"] = @benchmarkable ACE1.energy($(basis.BB[2]), $data_julip)
SUITE["AtomsBase"]["energy"]["big"] = @benchmarkable ace_energy($pot, $data)
SUITE["AtomsBase"]["energy"]["huge"] = @benchmarkable ace_energy($pot, $data_huge)
SUITE["AtomsBase"]["energy"]["basis"] = @benchmarkable ace_energy($(basis.BB[2]), $data)

SUITE["JuLIP"]["forces"]["big"] = @benchmarkable ACE1.forces($pot_julip, $data_julip)
SUITE["JuLIP"]["forces"]["huge"] = @benchmarkable ACE1.forces($pot_julip, $data_julip_huge)
SUITE["JuLIP"]["forces"]["basis"] = @benchmarkable ACE1.forces($(basis.BB[2]), $data_julip)
SUITE["AtomsBase"]["forces"]["big"] = @benchmarkable ace_forces($pot, $data)
SUITE["AtomsBase"]["forces"]["huge"] = @benchmarkable ace_forces($pot, $data_huge)
SUITE["AtomsBase"]["forces"]["basis"] = @benchmarkable ace_forces($(basis.BB[2]), $data)

SUITE["JuLIP"]["virial"]["big"] = @benchmarkable ACE1.virial($pot_julip, $data_julip)
SUITE["JuLIP"]["virial"]["huge"] = @benchmarkable ACE1.virial($pot_julip, $data_julip_huge)
SUITE["JuLIP"]["virial"]["basis"] = @benchmarkable ACE1.virial($(basis.BB[2]), $data_julip)
SUITE["AtomsBase"]["virial"]["big"] = @benchmarkable ace_virial($pot, $data)
SUITE["AtomsBase"]["virial"]["huge"] = @benchmarkable ace_virial($pot, $data_huge)
SUITE["AtomsBase"]["virial"]["basis"] = @benchmarkable ace_virial($(basis.BB[2]), $data)



## These test serial calculations compared to default async
## They should be always slower so skipping them here
# function serial_energy(pot, data)
#     return sum(  pot  ) do p
#         ace_energy(p, data)
#     end
# end

# function serial_forces(pot, data)
#     return sum(  pot  ) do p
#         ace_forces(p, data)
#     end
# end

#SUITE["AtomsBase", "energy serial", "big"] = @benchmarkable serial_energy($pot, $data)
#SUITE["AtomsBase", "forces serial", "big"] = @benchmarkable serial_forces($pot, $data)
