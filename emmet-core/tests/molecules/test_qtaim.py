import pytest

from monty.serialization import loadfn

from emmet.core.qchem.task import TaskDocument
from emmet.core.molecules.qtaim import QTAIMDoc


@pytest.fixture(scope="session")
def simple(test_dir):
    task = TaskDocument(**loadfn((test_dir / "simple_qtaim.json.gz")))
    return task


def test_qtaim(simple):
    # Test QTAIM parsing
    doc = QTAIMDoc.from_task(
        simple, "b9ba54febc77d2a9177accf4605767db-C1Li2O3-1-2", deprecated=False
    )

    assert doc.property_name == "quantum theory of atoms in molecules"

    assert len(doc.atoms) == 3
    assert len(doc.bonds) == 2
    assert len(doc.rings) == 0
    assert len(doc.cages) == 0

    assert doc.atoms[0].name == "1_H"
    assert doc.bonds[0].atom_indices == [1, 2]
    assert doc.bonds[0].position == pytest.approx([-0.603118228751, -0.0, 0.327664239395])
