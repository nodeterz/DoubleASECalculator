import numpy as np
import sys
from ase.calculators.calculator import Calculator, all_changes
from ase.io import write


class ForceEnergyCheckCalculator(Calculator):
    """
    ASE wrapper calculator that compares two calculators and
    automatically saves suspicious structures.

    Structures are saved when:
      - force cosine similarity exceeds threshold
      - force difference norm exceeds 5 meV/Å in zero-norm cases
      - energy difference per atom exceeds threshold

    All requested properties are returned from the main calculator.
    """

    implemented_properties = []

    def __init__(
        self,
        main,
        secondary,
        force_angle_threshold,
        energy_difference_per_atom,
        referance_energy_main=None,
        referance_energy_secondary=None,
        referance_energies_main=None,
        referance_energies_secondary=None,
        force_difference_threshold=0.005,  # eV/Å = 5 meV/Å
        force_dump_file="force_triggered.xyz",
        energy_dump_file="energy_triggered.xyz",
        **kwargs,
    ):
        super().__init__(**kwargs)

        self.main = main
        self.secondary = secondary

        self.referance_energy_main = referance_energy_main
        self.referance_energy_secondary = referance_energy_secondary
        
        self.referance_energies_main = referance_energies_main
        self.referance_energies_secondary = referance_energies_secondary

        self.force_angle_threshold = force_angle_threshold
        self.energy_difference_per_atom = (energy_difference_per_atom)

        self.force_difference_threshold = (force_difference_threshold)

        self.force_dump_file = force_dump_file
        self.energy_dump_file = energy_dump_file

        # Mirror supported properties from main calculator
        self.implemented_properties = list(getattr(main, "implemented_properties", []))

        self.last_positions = None

    def calculate( self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        
        current_pos = atoms.get_positions()
        if self.last_positions is not None and np.allclose(current_pos, self.last_positions):
            # ============================================================
            # Return requested properties from MAIN calculator
            # ============================================================
            self.results = {}
            for prop in properties:
                if prop == "forces":
                    self.results["forces"] = self.main.get_forces(atoms)
                elif prop == "energy":
                    self.results["energy"] = self.main.get_potential_energy(atoms)
                elif prop == "free_energy":
                    try:
                        self.results["free_energy"] = self.main.get_potential_energy(atoms,force_consitent=True)
                    except:
                        self.results["free_energy"] = self.main.get_potential_energy(atoms)
                elif prop == "stress":
                    self.results["stress"] = self.main.get_stress(atoms,voigt=False)
                elif prop == "energies":
                    self.results["energies"] = self.main.get_potential_energies(atoms)
                else:
                    # Generic fallback to main calculator
                    if hasattr(self.main, "get_property"):
                        self.results[prop] = (self.main.get_property(prop, atoms))
            self.last_positions = atoms.get_positions().copy()
        else:
            # ============================================================
            # Forces
            # ============================================================
            forces_1 = self.main.get_forces(atoms)
            forces_2 = self.secondary.get_forces(atoms)

            flat_f1 = forces_1.ravel()
            flat_f2 = forces_2.ravel()

            norm1 = np.linalg.norm(flat_f1)
            norm2 = np.linalg.norm(flat_f2)

            # ------------------------------------------------------------
            # Normal cosine similarity check
            # ------------------------------------------------------------
            if norm1 > 0.0 and norm2 > 0.0:
                cosine = (np.dot(flat_f1, flat_f2) / (norm1 * norm2))
                if np.abs(cosine) < self.force_angle_threshold:
                    atoms_to_save = atoms.copy()
                    atoms_to_save.info[ "trigger_reason" ] = "force_alignment"
                    atoms_to_save.info[ "force_cosine_similarity" ] = float(cosine)
                    atoms_to_save.new_array('force_diff_norm',np.linalg.norm(forces_1-forces_2,axis=1))
                    atoms_to_save.new_array('force_diff',(forces_1-forces_2))
                    atoms_to_save.new_array('force_1',(forces_1))
                    atoms_to_save.new_array('force_2',(forces_2))
                    write( self.force_dump_file, atoms_to_save, append=True)
            # ------------------------------------------------------------
            # Handle zero-norm case
            # ------------------------------------------------------------
            else:
                # Difference norm in eV/Å
                diff_norm = np.linalg.norm(flat_f1 - flat_f2)
                if ( diff_norm > self.force_difference_threshold):
                    atoms_to_save = atoms.copy()
                    atoms_to_save.info[ "trigger_reason" ] = "force_difference"
                    atoms_to_save.info[ "force_difference_norm" ] = float(diff_norm)
                    write( self.force_dump_file, atoms_to_save, append=True)
            # ============================================================
            # Energies
            # ============================================================
            ener1 = self.main.get_potential_energy(atoms)
            ener2 = self.secondary.get_potential_energy(atoms)
            n_atoms = atoms.get_global_number_of_atoms()
            energy_diff_per_atom = ((ener1-self.referance_energy_main) - (ener2-self.referance_energy_secondary))/n_atoms
            if (self.referance_energies_main is not None) and (self.referance_energies_secondary is not None):
                eners1=self.main.get_potential_energies(atoms)
                eners2=self.secondary.get_potential_energies(atoms)
                energies_diff = (eners1-self.referance_energies_main)-(eners2-self.referance_energies_secondary)
            if ( np.abs(energy_diff_per_atom) > self.energy_difference_per_atom):
                atoms_to_save = atoms.copy()
                atoms_to_save.info[ "trigger_reason" ] = "energy_difference"
                atoms_to_save.info["energy_main"] = float( ener1)
                atoms_to_save.info[ "energy_secondary" ] = float(ener2)
                atoms_to_save.info[ "energy_difference_per_atom" ] = float(energy_diff_per_atom)
                atoms_to_save.info[ "energy_difference" ] = float(energy_diff_per_atom*n_atoms)
                atoms_to_save.new_array('energies_difference',np.array(energies_diff))
                write( self.energy_dump_file, atoms_to_save, append=True)
            # ============================================================
            # Return requested properties from MAIN calculator
            # ============================================================
            self.results = {}
            for prop in properties:
                if prop == "forces":
                    self.results["forces"] = forces_1
                elif prop == "energy":
                    self.results["energy"] = ener1
                elif prop == "free_energy":
                    self.results["free_energy"] = ener1
                elif prop == "stress":
                    self.results["stress"] = (self.main.get_stress(atoms))
                elif prop == "energies":
                    self.results["energies"] = (self.main.get_potential_energies(atoms))
                else:
                    # Generic fallback to main calculator
                    if hasattr(self.main, "get_property"):
                        self.results[prop] = (self.main.get_property(prop, atoms))
            self.last_positions = atoms.get_positions().copy()
