"""The shipped reference data must not drift. These values are pinned."""

from __future__ import annotations

import numpy as np
import pytest

import tricyclepy as tp


def test_neuro_ref_shape_and_first_rows():
    ref = tp.load_neuro_ref()
    assert ref.rotation.shape == (500, 2)
    assert len(ref) == 500
    assert ref.ensembl[0] == "ENSMUSG00000040204"
    assert ref.symbol[0] == "Pclaf"
    assert ref.SYMBOL[0] == "PCLAF"
    assert ref.rotation[0, 0] == pytest.approx(-0.1989885, abs=1e-7)
    assert ref.rotation[0, 1] == pytest.approx(0.11245580, abs=1e-7)


def test_neuro_ref_identifiers_are_unique():
    ref = tp.load_neuro_ref()
    for names in (ref.ensembl, ref.symbol, ref.SYMBOL):
        assert len(set(names)) == 500


def test_upper_case_symbol_column_is_the_upper_cased_symbol():
    ref = tp.load_neuro_ref()
    assert [s.upper() for s in ref.symbol] == list(ref.SYMBOL)


@pytest.mark.parametrize(
    "species,gname_type,expected_first",
    [
        ("mouse", "ENSEMBL", "ENSMUSG00000040204"),
        ("mouse", "SYMBOL", "Pclaf"),
        ("human", "SYMBOL", "PCLAF"),
        # Human ignores gname_type: the reference is always keyed on SYMBOL.
        ("human", "ENSEMBL", "PCLAF"),
    ],
)
def test_get_rotation_keys(species, gname_type, expected_first):
    names, rotation = tp.get_rotation(gname_type=gname_type, species=species)
    assert names[0] == expected_first
    assert rotation.shape == (500, 2)


def test_get_rotation_rejects_bad_arguments():
    with pytest.raises(ValueError):
        tp.get_rotation(species="rat")
    with pytest.raises(ValueError):
        tp.get_rotation(gname_type="entrez")


def test_revelio_stage_order_and_sizes():
    stages = tp.load_revelio_gene_list()
    # Order matters: the Schwabe rule treats these five as a cycle.
    assert list(stages) == ["G1.S", "S", "G2", "G2.M", "M.G1"]
    assert [len(v) for v in stages.values()] == [147, 162, 182, 216, 148]
    assert stages["G1.S"][0] == "ABCA7"


def test_ensembl_to_symbol_map():
    human = tp.load_ensembl_to_symbol("human")
    assert human["ENSG00000141510"] == "TP53"
    assert len(human) > 30_000
    mouse = tp.load_ensembl_to_symbol("mouse")
    assert mouse["ENSMUSG00000059552"] == "Trp53"


def test_go_cell_cycle_gene_sets_are_populated():
    for species in ("human", "mouse"):
        for gname_type in ("SYMBOL", "ENSEMBL"):
            genes = tp.load_go_cell_cycle_genes(species=species,
                                                gname_type=gname_type)
            assert genes.size > 2000
    assert "TOP2A" in set(tp.load_go_cell_cycle_genes("human", "SYMBOL"))


def test_provenance_records_the_r_versions_used():
    info = tp.provenance()
    assert info["tricycle"] == "1.12.0"
    assert info["org.Hs.eg.db"].startswith("3.")
    assert tp.TRICYCLE_R_VERSION == info["tricycle"]


def test_reference_data_loaders_are_cached_not_copied():
    # The loaders are lru_cached, so callers must not mutate what they return.
    assert tp.load_neuro_ref() is tp.load_neuro_ref()
    assert np.shares_memory(tp.get_rotation()[1], tp.load_neuro_ref().rotation)
