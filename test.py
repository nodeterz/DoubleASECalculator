import numpy as np
from ase import units
from ase.build import bulk
from ase.md import MDLogger
from ase.calculators.emt import EMT
from ase.md.verlet import VelocityVerlet
from ase.calculators.lj import LennardJones
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from force_energy_check_calculator import ForceEnergyCheckCalculator


atoms = bulk("Cu", "fcc", a=3.61, cubic=True)
atoms = atoms * (3, 3, 3)
atoms.write('posinp.extxyz')
calc1 = LennardJones(epsilon=3.837914780509462, sigma=1.840593221635568)
calc2 = LennardJones(epsilon=3.837914780509462, sigma=1.940593221635568)
atoms.calc=calc1
reference_energy_main=atoms.get_potential_energy()
reference_energies_main=atoms.get_potential_energies()
atoms.calc=calc2
reference_energy_secondary=atoms.get_potential_energy()
reference_energies_secondary=atoms.get_potential_energies()

atoms.calc = ForceEnergyCheckCalculator(
    main=calc1,
    secondary=calc2,
    reference_energy_main=reference_energy_main,
    reference_energy_secondary=reference_energy_secondary,
    reference_energies_main=reference_energies_main,
    reference_energies_secondary=reference_energies_secondary,
    force_angle_threshold=0.95,
    energy_difference_per_atom=0.001,
)

#forces = atoms.get_forces()
#energy = atoms.get_potential_energy()
#stress = atoms.get_stress()

temp_k = 600 + 273.15
MaxwellBoltzmannDistribution(atoms, temperature_K=temp_k)
dyn = VelocityVerlet(atoms, timestep=2.0 * units.fs)
dyn.attach(MDLogger(dyn, atoms, 'nve.log', mode='w', peratom=True), interval=10)
print(f"Starting NVE MD at {temp_k} K...")
dyn.run(10000)
print("Simulation complete.")
