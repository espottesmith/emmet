import pytest
from maggma.stores import JSONStore, MemoryStore

from emmet.builders.qchem.molecules import MoleculesAssociationBuilder, MoleculesBuilder
from emmet.builders.molecules.qtaim import QTAIMBuilder


__author__ = "Evan Spotte-Smith <ewcspottesmith@cmu.edu>"


@pytest.fixture(scope="session")
def tasks_store(test_dir):
    return JSONStore(test_dir / "C3H6.json.gz")


@pytest.fixture(scope="session")
def assoc_store(tasks_store):
    assoc_store = MemoryStore(key="molecule_id")
    stage_one = MoleculesAssociationBuilder(tasks=tasks_store, assoc=assoc_store)
    stage_one.run()

    return assoc_store


@pytest.fixture(scope="session")
def qtaim_store():
    return MemoryStore()


def test_qtaim_builder(tasks_store, assoc_store, qtaim_store):
    builder = QTAIMBuilder(tasks_store, assoc_store, qtaim_store)
    builder.run()

    assert qtaim_store.count() == 8
    assert qtaim_store.count({"deprecated": True}) == 0
