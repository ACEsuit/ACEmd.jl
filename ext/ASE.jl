module ASE

using ACEapi
using AtomsBase
using ASEconvert
using Unitful
using UnitfulAtomic
using PythonCall


export ACEcalculator

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
                F = ace_energy(pyconvert(ACEpotential,self.potential), pyconvert(AbstractSystem, atoms))
                self.results["forces"] = F
            end
        end

    )
])

end