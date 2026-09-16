"""Reference data shipped with the package.

Every object here is carried over verbatim from the R package ``tricycle``
(version 1.12.0, Bioconductor 3.19) or from the Bioconductor annotation
packages that its ``ENSEMBL`` code paths query. Nothing is regenerated. A
regenerated reference would not reproduce the R numbers.

Provenance is recorded in ``data/PROVENANCE.tsv``.
"""

from __future__ import annotations

import csv
import gzip
from dataclasses import dataclass
from functools import lru_cache
from importlib import resources
from typing import Dict, List, Sequence

import numpy as np

__all__ = [
    "NeuroRef",
    "load_neuro_ref",
    "load_revelio_gene_list",
    "load_ensembl_to_symbol",
    "load_go_cell_cycle_genes",
    "get_rotation",
    "provenance",
]

SPECIES = ("mouse", "human")
GNAME_TYPES = ("ENSEMBL", "SYMBOL")


def _data_path(name: str):
    return resources.files(__package__).joinpath("data", name)


@dataclass(frozen=True)
class NeuroRef:
    """The pre-learned Neurosphere reference, as the R ``neuroRef`` data frame.

    Attributes
    ----------
    rotation
        ``(500, 2)`` float array. Columns are ``pc1.rot`` and ``pc2.rot``.
    ensembl
        Mouse Ensembl gene identifiers, unversioned.
    symbol
        Mouse gene symbols, as capitalised by MGI.
    SYMBOL
        The upper-cased mouse symbols. R uses these as the human gene names.
    """

    rotation: np.ndarray
    ensembl: np.ndarray
    symbol: np.ndarray
    SYMBOL: np.ndarray

    def __len__(self) -> int:
        return self.rotation.shape[0]


@lru_cache(maxsize=1)
def load_neuro_ref() -> NeuroRef:
    """Return the Neurosphere reference projection matrix and its gene names."""
    pc1: List[float] = []
    pc2: List[float] = []
    ensembl: List[str] = []
    symbol: List[str] = []
    SYMBOL: List[str] = []
    with _data_path("neuroRef.tsv").open("r", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            pc1.append(float(row["pc1.rot"]))
            pc2.append(float(row["pc2.rot"]))
            ensembl.append(row["ensembl"])
            symbol.append(row["symbol"])
            SYMBOL.append(row["SYMBOL"])
    return NeuroRef(
        rotation=np.column_stack([np.asarray(pc1, dtype=np.float64),
                                  np.asarray(pc2, dtype=np.float64)]),
        ensembl=np.asarray(ensembl, dtype=object),
        symbol=np.asarray(symbol, dtype=object),
        SYMBOL=np.asarray(SYMBOL, dtype=object),
    )


def get_rotation(gname_type: str = "ENSEMBL", species: str = "mouse"):
    """Return ``(gene_names, rotation)`` for the internal Neurosphere reference.

    This is the Python form of the R private function ``.getRotation``. For
    ``species="human"`` the reference gene names are the upper-cased mouse
    symbols, whatever ``gname_type`` says -- that is what R does, and it is why
    a human Ensembl input has to be mapped to symbols before projection.

    Parameters
    ----------
    gname_type
        ``"ENSEMBL"`` or ``"SYMBOL"``. Ignored when ``species="human"``.
    species
        ``"mouse"`` or ``"human"``.

    Returns
    -------
    names : numpy.ndarray
        Object array of gene names, length 500.
    rotation : numpy.ndarray
        ``(500, 2)`` float array with columns PC1 and PC2.
    """
    if species not in SPECIES:
        raise ValueError(f"species must be one of {SPECIES}, got {species!r}")
    if gname_type not in GNAME_TYPES:
        raise ValueError(f"gname_type must be one of {GNAME_TYPES}, got {gname_type!r}")
    ref = load_neuro_ref()
    if species == "human":
        names = ref.SYMBOL
    else:
        names = ref.ensembl if gname_type == "ENSEMBL" else ref.symbol
    return names, ref.rotation


@lru_cache(maxsize=1)
def load_revelio_gene_list() -> Dict[str, List[str]]:
    """Return the five-stage Revelio marker list, keyed by stage.

    Stage order is ``G1.S, S, G2, G2.M, M.G1``. Order matters: the Schwabe
    assignment treats the stages as a cycle and rejects a cell whose top two
    stages are not neighbours in this order.
    """
    stages: Dict[str, List[str]] = {}
    with _data_path("RevelioGeneList.tsv").open("r", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            stages.setdefault(row["stage"], []).append(row["gene"])
    return stages


@lru_cache(maxsize=2)
def load_ensembl_to_symbol(species: str = "human") -> Dict[str, str]:
    """Return the Ensembl-to-symbol map used by the R ``ENSEMBL`` code paths.

    The map reproduces ``AnnotationDbi::mapIds(db, keys, "SYMBOL", "ENSEMBL",
    multiVals = "first")`` for ``org.Hs.eg.db`` or ``org.Mm.eg.db``. An
    identifier absent from the map is unmapped, exactly as R returns ``NA``.
    """
    if species not in SPECIES:
        raise ValueError(f"species must be one of {SPECIES}, got {species!r}")
    out: Dict[str, str] = {}
    path = _data_path(f"ensembl2symbol_{species}.tsv.gz")
    with path.open("rb") as raw, gzip.open(raw, "rt", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        next(reader)  # header
        for ensembl, symbol in reader:
            out[ensembl] = symbol
    return out


@lru_cache(maxsize=4)
def load_go_cell_cycle_genes(species: str = "human",
                             gname_type: str = "SYMBOL") -> np.ndarray:
    """Return the GO:0007049 (cell cycle) gene set, with duplicates kept.

    R's ``AnnotationDbi::select(db, keytype = "GOALL", keys = "GO:0007049",
    columns = gname_type)`` returns one row per gene/evidence pair, so the
    vector holds repeats and ``NA``. The repeats do not change the result of
    ``run_pca_cc_genes`` -- it only asks which input genes are members -- but
    they do change the reported denominator, so they are kept.
    """
    if species not in SPECIES:
        raise ValueError(f"species must be one of {SPECIES}, got {species!r}")
    if gname_type not in GNAME_TYPES:
        raise ValueError(f"gname_type must be one of {GNAME_TYPES}, got {gname_type!r}")
    path = _data_path(f"go_cc_{species}_{gname_type.lower()}.txt")
    text = path.read_text(encoding="utf-8").splitlines()
    return np.asarray(text, dtype=object)


@lru_cache(maxsize=1)
def provenance() -> Dict[str, str]:
    """Return the software versions the shipped reference data came from."""
    out: Dict[str, str] = {}
    for line in _data_path("PROVENANCE.tsv").read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        key, _, value = line.partition("\t")
        out[key] = value
    return out


def map_ensembl_to_symbol(gname: Sequence[str], species: str = "human",
                          mapping: Dict[str, str] | None = None) -> np.ndarray:
    """Map Ensembl identifiers to symbols, leaving unmapped entries as ``None``.

    Version suffixes are not stripped. R does not strip them either, so an
    input of ``ENSG00000141510.16`` fails to map on both sides.
    """
    table = load_ensembl_to_symbol(species) if mapping is None else mapping
    return np.asarray([table.get(g) for g in gname], dtype=object)
