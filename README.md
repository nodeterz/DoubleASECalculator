# ForceEnergyCheckCalculator

`ForceEnergyCheckCalculator` is an ASE wrapper calculator designed to estimate the reliability (or certainty) of a primary atomistic calculator by continuously comparing its predictions against a secondary reference calculator.

The wrapper computes forces and energies from both calculators and automatically stores suspicious structures whenever the disagreement exceeds user-defined thresholds.

This is useful for:

- active learning workflows
- uncertainty detection
- ML potential validation
- out-of-domain structure detection
- hybrid fast/accurate simulation pipelines
- online dataset generation

---

# Features

## Force agreement monitoring

The calculator computes:

- `forces_1` from the main calculator
- `forces_2` from the secondary calculator

and evaluates their cosine similarity:

```math
\cos(\theta)=\frac{F_1 \cdot F_2}{|F_1||F_2|}
```

If the similarity exceeds a threshold, the structure is automatically saved.

---

## Zero-force handling

If either force norm is zero, cosine similarity becomes unstable.

In this case the calculator instead computes:

```math
|F_1 - F_2|
```

and saves the structure if the disagreement exceeds:

```text
5 meV / Å
```

(default: `0.005 eV/Å`)

---

## Energy disagreement detection

The wrapper also compares total energies:

```math
\frac{E_1 - E_2}{N_\mathrm{atoms}}
```

If the energy disagreement per atom exceeds a threshold, the structure is saved.

---

## Transparent ASE integration

The wrapper behaves like a normal ASE calculator and returns all requested properties directly from the main calculator.

Supported properties are automatically mirrored from the main calculator:

- energy
- forces
- stress
- per-atom energies
- free energy
- any additional implemented ASE property

---

# Installation

Simply copy the calculator into your project:

```bash
force_energy_check_calculator.py
```

Dependencies:

```bash
pip install ase numpy
```

---

# Example

```python
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
```

---

# Output

Suspicious structures are automatically appended to:

```text
force_triggered.xyz
energy_triggered.xyz
```

Each saved structure contains metadata in `atoms.info`, for example:

```python
atoms.info["trigger_reason"]
atoms.info["force_cosine_similarity"]
atoms.info["energy_difference_per_atom"]
```

---

# Intended Usage

A common workflow is:

- `main` = fast ML potential
- `secondary` = expensive DFT calculator

The wrapper identifies configurations where the ML model may become unreliable and stores them for later retraining or analysis.

---

# Parameters

| Parameter | Description |
|---|---|
| `main` | Primary calculator whose outputs are returned |
| `secondary` | Reference calculator used for comparison |
| `force_angle_threshold` | Cosine similarity threshold |
| `energy_difference_per_atom` | Energy disagreement threshold |
| `reference_energy_main` | Reference energy from the main calculator |
| `reference_energy_secondary` | Reference energy from the second calculator |
| `reference_energies_main` | Reference atomic energies from the main calculator (Optional)|
| `reference_energies_secondary` | Reference atomic energies from the second calculator (Optional) |
| `force_difference_threshold` | Force disagreement threshold for zero-norm cases |
| `force_dump_file` | File for force-triggered structures |
| `energy_dump_file` | File for energy-triggered structures |
---

# Notes

- The wrapper always returns results from the **main** calculator.
- The secondary calculator is only used for consistency checks.
- Structures are appended incrementally to XYZ files for later inspection.
- The calculator is fully compatible with ASE optimizers and MD workflows.
- The reference energies must be calculated from a similar structure, and this is required to cancel out the energy bias.
---

# License

MIT
