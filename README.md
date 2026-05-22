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
from ase.build import bulk
from ase.calculators.emt import EMT

atoms = bulk("Cu")

calc = ForceEnergyCheckCalculator(
    main=EMT(),
    secondary=EMT(),
    force_angle_threshold=0.95,
    energy_difference_per_atom=0.01,
)

atoms.calc = calc

energy = atoms.get_potential_energy()
forces = atoms.get_forces()
stress = atoms.get_stress()
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
| `force_difference_threshold` | Force disagreement threshold for zero-norm cases |
| `force_dump_file` | File for force-triggered structures |
| `energy_dump_file` | File for energy-triggered structures |

---

# Notes

- The wrapper always returns results from the **main** calculator.
- The secondary calculator is only used for consistency checks.
- Structures are appended incrementally to XYZ files for later inspection.
- The calculator is fully compatible with ASE optimizers and MD workflows.

---

# License

MIT
