using BenchmarkTools
using ACEapi
using ACE1
using ExtXYZ
using AtomsBase

const SUITE = BenchmarkGroup()

fname_ace = joinpath(pkgdir(ACEapi), "data", "TiAl.json")
fname_xyz = joinpath(pkgdir(ACEapi), "data", "TiAl-big.xyz")
fname_xyz_huge = joinpath(pkgdir(ACEapi), "data", "TiAl-huge.xyz")

pot = load_ace_model(fname_ace)
pot_julip = load_ace_model(fname_ace; old_format=true)
data_julip = read_extxyz(fname_xyz)[end]
data_julip_huge = read_extxyz(fname_xyz_huge)[end]
data = FastSystem(ExtXYZ.Atoms(read_frame(fname_xyz)))
data_huge = FastSystem(ExtXYZ.Atoms(read_frame(fname_xyz_huge)))



SUITE["JuLIP"] = BenchmarkGroup(["JuLIP","Modes", "Size"])
SUITE["AtomsBase"] = BenchmarkGroup(["AtomsBase", "Modes", "Size"])

SUITE["JuLIP", "energy", "big"] = @benchmarkable ACE1.energy($pot_julip, $data_julip)
SUITE["JuLIP", "energy", "huge"] = @benchmarkable ACE1.energy($pot_julip, $data_julip_huge)
SUITE["AtomsBase", "energy", "big"] = @benchmarkable ace_energy($pot, $data)
SUITE["AtomsBase", "energy", "huge"] = @benchmarkable ace_energy($pot, $data_huge)

SUITE["JuLIP", "forces", "big"] = @benchmarkable ACE1.forces($pot_julip, $data_julip)
SUITE["JuLIP", "forces", "huge"] = @benchmarkable ACE1.forces($pot_julip, $data_julip_huge)
SUITE["AtomsBase", "forces", "big"] = @benchmarkable ace_forces($pot, $data)
SUITE["AtomsBase", "forces", "huge"] = @benchmarkable ace_forces($pot, $data_huge)

SUITE["JuLIP", "virial", "big"] = @benchmarkable ACE1.virial($pot_julip, $data_julip)
SUITE["JuLIP", "virial", "huge"] = @benchmarkable ACE1.virial($pot_julip, $data_julip_huge)
SUITE["AtomsBase", "virial", "big"] = @benchmarkable ace_virial($pot, $data)
SUITE["AtomsBase", "virial", "huge"] = @benchmarkable ace_virial($pot, $data_huge)



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