import re

from typing import List, Optional, Dict, Any
from hashlib import blake2b

from pydantic import Field

from monty.json import MSONable

from emmet.core.mpid import MPculeID
from emmet.core.material import PropertyOrigin
from emmet.core.qchem.task import TaskDocument
from emmet.core.molecules.molecule_property import PropertyDoc


# TODO: change my e-mail on everything
__author__ = "Evan Spotte-Smith <ewcspottesmith@cmu.edu>"


class CriticalPoint(MSONable):
    def __init__(
        self,
        critical_point_index: int,
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
        connected_bond_paths: Optional[List[int]] = None
    ):
        """

        """
        
        # TODO
        pass


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

        """
        
        #TODO
        pass