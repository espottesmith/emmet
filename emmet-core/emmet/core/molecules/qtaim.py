import re

from typing import List, Optional, Dict, Any
from hashlib import blake2b

from pydantic import Field

from monty.json import MSONable

from emmet.core.mpid import MPculeID
from emmet.core.material import PropertyOrigin
from emmet.core.qchem.task import TaskDocument
from emmet.core.molecules.molecule_property import PropertyDoc


__author__ = "Evan Spotte-Smith <ewcspottesmith@cmu.edu>"


class CriticalPoint(MSONable):
    def __init__(
        self,
        name: str,
        position: List[float],
        total_electron_density: float,
        alpha_electron_density: float,
        beta_electron_density: float,
        spin_density: float,
        laplacian_electron_density: float,
        localized_orbital_locator: float,
        energy_density: float,
        lagrangian_kinetic_energy: float,
        hamiltonian_kinetic_energy: float,
        electron_localization_function: float,
        average_local_ionization_energy: float,
        delta_g_promolecular: float,
        delta_g_hirsh: float,
        esp_nuclear: float,
        esp_electron: float,
        esp_total: float,
        gradient_norm: float,
        laplacian_norm: float,
        determinant_hessian: float,
        ellipticity_electron_density: float,
        eta_index: float,
        atom_index: Optional[int] = None,
        atom_indices: Optional[List[int]] = None,
        bond_names: Optional[List[str]] = None,
        ring_names: Optional[List[str]] = None,
        connected_bond_paths: Optional[List[int]] = None
    ):
        """
        A critical point (CP) parsed from Multiwfn quantum theory of atoms in molecules (QTAIM) analysis.

        All values are given in atomic units

        :param name (str): Name of this CP. This is mainly important for bonds and rings, as bond names and ring
            names are associated with ring and cage CPs (see below).
        :param position (List[float]): Cartesian (x, y, z) location of this CP
        :param total_electron_density (float): Value of the total electron density at this CP
        :param alpha_electron_density (float): Value of the electron density of alpha electrons at this CP 
        :param beta_electron_density (float): Value of the electron density of beta electrons at this CP
        :param spin_density (float): Value of the spin density at this CP
        :param laplacian_electron_density (float): Laplacian of the electron density at this CP
        :param localized_orbital_locator (float): Localized orbital locator function value at this CP. This value
            is somewhat related to the electron localization function (see below), in that it seeks to quantify
            high-localization points.
            See:  J. Mol. Struct. (THEOCHEM), 527, 51 (2000); DOI: 10.1016/S0166-1280(00)00477-2
        :param energy_density (float): Energy density at this CP
        :param lagrangian_kinetic_energy (float): Lagrangian (or positive definite) kinetic energy density at this CP
        :param hamiltonian_kinetic_energy (float): Hamiltonian kinetic energy density at this CP
        :param electron_localization_function (float): Electron localization function (ELF) at this CP, which seeks to
            quantify the relative extent of electron confinement based on the Pauli kinetic energy density and the
            Thomas-Fermi kinetic energy density of the uniform electron gas.
        :param average_local_ionization_energy (float): Average local ionization energy at this CP, indicating how
            weakly bound electrons are locally
        :param delta_g_promolecular (float): The dg function from the independent gradient model (IGM), obtained using
            the promolecular approximation that atomic densities are the same in a molecule as in isolation
            See Phys. Chem. Chem. Phys., 19, 17928 (2017); DOI: 10.1039/C7CP02110K
        :param delta_g_hirsh (float): The dg function obtained using the Hirshfeld partition of the electron density
        :param esp_nuclear (float): Electronstatic potential for nuclear charges at this CP
        :param esp_electron (float): Electrostatic potential for electronic charges at this CP
        :param esp_total (float): Total electrostatic potential at this CP
        :param gradient_norm (float): Gradient norm of the electron density at this CP
        :param determinant_hessian (float): Determinant of the Hessian matrix for this CP
        :param ellipticity_electron_density (float): The ellipticity of this CP, defined based on the Hessian
            eigenvalues of this CP
        :param eta_index (float): Eta index of this CP, which may indicate the covalent character of an interaction.
            See Angew. Chem. Int. Ed., 53, 2766 (2014); DOI: 10.1002/anie.201308609
        :param atom_index (Optional[int]): For an atom CP, the index of the atom associated with this CP
        :param atom_indices (Optional[List[int]]): For bond, ring, and cage CPs, the indices of the atoms surrounding
            this CP
        :param bond_names (Optional[List[str]]): For ring and cage CPs, the names of the bond CPs surrounding this CP
        :param ring_names (Optional[List[str]]): For cage CPs, the names of the ring CPs surrounding this CP
        """
        
        self.name = name
        self.position = position
        self.total_electron_density = total_electron_density
        self.alpha_electron_density = alpha_electron_density
        self.beta_electron_density = beta_electron_density
        self.spin_density = spin_density
        self.laplacian_electron_density = laplacian_electron_density
        self.localized_orbital_locator = localized_orbital_locator
        self.energy_density = energy_density
        self.lagrangian_kinetic_energy = lagrangian_kinetic_energy
        self.hamiltonian_kinetic_energy = hamiltonian_kinetic_energy
        self.electron_localization_function = electron_localization_function
        self.average_local_ionization_energy = average_local_ionization_energy
        self.delta_g_promolecular = delta_g_promolecular
        self.delta_g_hirsh = delta_g_hirsh
        self.esp_nuclear = esp_electron
        self.esp_electron = esp_electron
        self.esp_total = esp_total
        self.gradient_norm = gradient_norm
        self.determinant_hessian = determinant_hessian
        self.ellipticity_electron_density = ellipticity_electron_density
        self.eta_index = eta_index
        self.atom_index = atom_index
        self.atom_indices = atom_indices
        self.bond_names = bond_names
        self.ring_names = ring_names


class QTAIMDoc(PropertyDoc):
    property_name: str = "quantum theory of atoms in molecules"

    atoms: List[CriticalPoint] = Field(
        ..., description="Atom or nuclear critical points"
    )

    bonds: List[CriticalPoint] = Field(
        [], description="Bond critical points"
    )

    rings: List[CriticalPoint] = Field(
        [], description="Ring critical points"
    )

    cages: List[CriticalPoint] = Field(
        [], description="Cage critical points"
    )

    @classmethod
    def from_task(
        cls,
        task: TaskDocument,
        molecule_id: MPculeID,
        deprecated: bool = False,
        **kwargs
    ):  # type: ignore[override]
        """
        Construct a QTAIM document with electron critical points from a task

        :param task: document from which QTAIm properties can be extracted
        :param molecule_id: MPculeID
        :param deprecated: bool. Is this document deprecated?
        :param kwargs: to pass to PropertyDoc
        :return:
        """
        
        if task.output.qtaim is None:
            raise ValueError("No QTAIM info in task {}!".format(task.task_id))
        
        qtaim = task.output.qtaim

        if task.output.optimized_molecule is not None:
            mol = task.output.optimized_molecule
        else:
            mol = task.output.initial_molecule

        atoms = list()
        bonds = list()
        rings = list()
        cages = list()

        mapping = {
            "pos_ang": "position",
            "density_total": "total_electron_density",
            "density_alpha": "alpha_electron_density",
            "density_beta": "beta_electron_density",
            "lap_e_density": "laplacian_electron_density",
            "lol": "localized_orbital_locator",
            "Lagrangian_K": "lagrangian_kinetic_energy",
            "Hamiltonian_K": "hamiltonian_kinetic_energy",
            "e_loc_func": "electron_localization_function",
            "ave_loc_ion_E": "average_local_ionization_energy",
            "esp_nuc": "esp_nuclear",
            "esp_e": "esp_electron",
            "grad_norm": "gradient_norm",
            "det_hessian": "determinant_hessian",
            "ellip_e_dens": "ellipticity_electron_density",
            "eta": "eta_index",
            "atom_inds": "atom_indices",
        }

        for group, grouplist in [
            ("atom", atoms),
            ("bond", bonds),
            ("ring", rings),
            ("cage", cages)
        ]:
            for name, desc in qtaim[group].items():
                this_data = dict()
                for k, v in desc.items():
                    if k in mapping:
                        this_data[mapping[k]] = v
                    elif k in ["lap_norm"]:
                        continue
                    else:
                        this_data[k] = v
                    
                    if group != "atom":
                        this_data["name"] = name
                
                grouplist.append(this_data)

        id_string = (
            f"qtaim-{molecule_id}-{task.task_id}-{task.lot_solvent}"
        )
        h = blake2b()
        h.update(id_string.encode("utf-8"))
        property_id = h.hexdigest()

        return super().from_molecule(
            meta_molecule=mol,
            property_id=property_id,
            molecule_id=molecule_id,
            level_of_theory=task.level_of_theory,
            solvent=task.solvent,
            lot_solvent=task.lot_solvent,
            atoms=[CriticalPoint(**acp) for acp in atoms],
            bonds=[CriticalPoint(**bcp) for bcp in bonds],
            rings=[CriticalPoint(**rcp) for rcp in rings],
            cages=[CriticalPoint(**ccp) for ccp in cages],
            origins=[
                PropertyOrigin(
                    name="qtaim", task_id=task.task_id
                )
            ],
            warnings=list(),
            deprecated=deprecated,
            **kwargs,
        )
