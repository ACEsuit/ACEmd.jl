module ACE_ASEext

using ACEapi
using AtomsBase
using ASEconvert
using Unitful
using UnitfulAtomic
using PythonCall

calculator = pyimport("ase.calculators.calculator")
Calculator = calculator.Calculator

ACEcalculator = pytype("ACEcalculator", (Calculator,),[
    "implemented_properties" => ["energy", "forces"],

    pyfunc(
        name  = "__init__",
        function (self, potential, atoms=nothing)
            calculator.Calculator.__init__(self, atoms=atoms)

            self.potential = potential
            return
        end
    ),

    pyfunc(
        name = "calculate",
        function (self, atoms=nothing, properties=["energy"], system_changes=nothing)
            if "energy" in properties
                E = ace_energy(pyconvert(ACEpotential,self.potential), pyconvert(AbstractSystem, atoms))
                self.results["energy"] = ustrip(u"eV", E*u"hartree")
            end

            if "forces" in properties
                F = ace_forces(pyconvert(ACEpotential,self.potential), pyconvert(AbstractSystem, atoms))
                conv = ustrip(u"eV/Å", 1u"hartree/Å")
                self.results["forces"] = F * conv
            end
        end

    )
])

ACEapi.ase_ace_calculator(pot::ACEpotential) = ACEcalculator(pot)

end