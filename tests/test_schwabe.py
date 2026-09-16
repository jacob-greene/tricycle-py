"""Five-stage assignment: the selection rules, the cycle rule, the tolerance."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest

import tricyclepy as tp
from tricyclepy.schwabe import _assign, _scale_columns, _stage_score


def _called(adata):
    """Stage calls as an object array, with unassigned cells as None.

    pandas turns an unassigned categorical into NaN, which is not None, so a
    naive `!= None` test silently counts it as a call.
    """
    values = adata.obs["CCStage"].astype(object).to_numpy()
    return np.asarray([None if pd.isna(v) else v for v in values], dtype=object)

RNG = np.random.default_rng(4242)

STAGES = ["G1.S", "S", "G2", "G2.M", "M.G1"]


def _synthetic(n_per_stage=40, n_markers=10, noise=0.05):
    """Cells on a latent cycle, with one marker block peaking per stage.

    The structure has to be cyclic, not five isolated clusters. The Schwabe
    rule rejects a cell whose two best stages are not neighbours on the cycle,
    so isolated clusters leave the runner-up stage random and almost every cell
    unassigned. A latent angle puts the runner-up next door, which is the
    situation the rule was designed for.
    """
    genes = [f"{s}_{k}" for s in STAGES for k in range(n_markers)]
    n_cells = n_per_stage * len(STAGES)
    n_stages = len(STAGES)
    centres = 2 * np.pi * np.arange(n_stages) / n_stages
    angle = np.repeat(centres, n_per_stage)
    angle = angle + RNG.uniform(-0.25, 0.25, n_cells) * (2 * np.pi / n_stages)

    x = RNG.normal(3.0, noise, (n_cells, len(genes)))
    for si in range(n_stages):
        cols = slice(si * n_markers, (si + 1) * n_markers)
        x[:, cols] += 5.0 * np.cos(angle - centres[si])[:, None]
    truth = np.repeat(np.asarray(STAGES, dtype=object), n_per_stage)
    cycle_genes = {s: [f"{s}_{k}" for k in range(n_markers)] for s in STAGES}
    adata = ad.AnnData(X=x, var=pd.DataFrame(index=genes),
                       obs=pd.DataFrame(index=[f"c{i}" for i in range(n_cells)]))
    return adata, cycle_genes, truth


def test_clean_stage_structure_is_recovered():
    adata, cycle_genes, truth = _synthetic()
    tp.estimate_schwabe_stage(adata, cycle_genes=cycle_genes)
    called = _called(adata)
    assigned = called != None                                     # noqa: E711
    assert assigned.mean() > 0.9
    assert (called[assigned] == truth[assigned]).all()


def test_result_is_categorical_with_the_stage_order_preserved():
    adata, cycle_genes, _ = _synthetic()
    tp.estimate_schwabe_stage(adata, cycle_genes=cycle_genes)
    assert list(adata.obs["CCStage"].cat.categories) == STAGES


def test_batches_are_scored_separately():
    adata, cycle_genes, truth = _synthetic()
    batch = np.where(np.arange(adata.n_obs) % 2 == 0, "a", "b")
    adata.obs["batch"] = batch
    # Shift one batch wholesale. Per-batch z-scoring must absorb it.
    adata.X[batch == "b"] += 4.0
    tp.estimate_schwabe_stage(adata, cycle_genes=cycle_genes, batch_key="batch")
    called = _called(adata)
    assigned = called != None                                     # noqa: E711
    assert assigned.mean() > 0.9
    assert (called[assigned] == truth[assigned]).all()


def test_batch_and_batch_key_are_mutually_exclusive():
    adata, cycle_genes, _ = _synthetic()
    adata.obs["batch"] = "a"
    with pytest.raises(ValueError, match="not both"):
        tp.estimate_schwabe_stage(adata, cycle_genes=cycle_genes,
                                  batch=np.zeros(adata.n_obs), batch_key="batch")


def test_too_few_cell_cycle_genes_is_refused():
    adata, cycle_genes, _ = _synthetic(n_markers=4)
    with pytest.raises(ValueError, match="Less than 30"):
        tp.estimate_schwabe_stage(adata, cycle_genes=cycle_genes)


def test_missing_batch_labels_are_refused():
    adata, cycle_genes, _ = _synthetic()
    batch = np.asarray(["a"] * adata.n_obs, dtype=object)
    batch[3] = None
    with pytest.raises(ValueError, match="missing values"):
        tp.estimate_schwabe_stage(adata, cycle_genes=cycle_genes, batch=batch)


# --- the two rejection rules, tested on the decision function directly ------


def _scores(order, gap=1.0):
    """Five z-scores whose ranking is `order`, separated by `gap`."""
    s = np.zeros(5)
    for rank, stage in enumerate(order):
        s[stage] = -rank * gap
    return s[None, :]


@pytest.mark.parametrize("top_two,assigned", [
    ((0, 1), True),     # neighbours
    ((1, 0), True),
    ((3, 4), True),
    ((4, 0), True),     # neighbours across the wrap: |4 - 0| == 4 == K - 1
    ((0, 4), True),
    ((0, 2), False),    # two apart
    ((1, 4), False),    # three apart
    ((0, 3), False),
])
def test_only_neighbouring_stages_are_assignable(top_two, assigned):
    rest = [i for i in range(5) if i not in top_two]
    scores = _scores(list(top_two) + rest, gap=1.0)
    out = _assign(scores, STAGES, tolerance=0.3)
    assert (out[0] is not None) is assigned
    if assigned:
        assert out[0] == STAGES[top_two[0]]


def test_a_gap_below_the_tolerance_leaves_the_cell_unassigned():
    close = np.array([[1.0, 0.9, 0.0, -1.0, -2.0]])   # neighbours, gap 0.1
    assert _assign(close, STAGES, tolerance=0.3)[0] is None
    assert _assign(close, STAGES, tolerance=0.05)[0] == "G1.S"


def test_ties_are_broken_by_stage_order_as_r_does():
    """R's order() is stable, so among equal scores the earlier stage wins.

    One tie pattern does not pin this. NumPy's quicksort agrees with a stable
    sort on many small inputs by accident, so a single case passes whether or
    not stability is asked for. Every tie position is swept instead, at the
    five shipped stages and at a longer custom cycle where an unstable sort
    reorders almost every draw.
    """
    for k in range(4):
        scores = np.full((1, 5), -5.0)
        scores[0, k] = scores[0, k + 1] = 1.0        # adjacent pair, tied
        assert _assign(scores, STAGES, tolerance=-1.0)[0] == STAGES[k]

    long_stages = [f"S{i}" for i in range(30)]
    rng = np.random.default_rng(3)
    for _ in range(200):
        scores = rng.integers(0, 3, (1, 30)).astype(float)
        top = float(scores.max())
        first = int(np.flatnonzero(scores[0] == top)[0])
        out = _assign(scores, long_stages, tolerance=-1.0)[0]
        if out is not None:
            assert out == long_stages[first]


def test_the_correlation_floor_is_strict_not_inclusive():
    """R keeps a marker when `cor > corThres`, so equality is excluded."""
    rng = np.random.default_rng(7)
    block = rng.normal(0.0, 1.0, (6, 40)) + 5.0       # genes by cells
    rows = np.arange(6)

    sub = block[rows]
    mean_v = sub.mean(axis=0)
    centred = sub - sub.mean(axis=1, keepdims=True)
    mc = mean_v - mean_v.mean()
    cor = (centred @ mc) / (np.sqrt((centred ** 2).sum(axis=1))
                            * np.sqrt((mc ** 2).sum()))
    cut = float(np.sort(cor)[2])                      # one gene sits exactly here

    at_threshold = _stage_score(block, rows, cut, "S", 1)
    just_below = _stage_score(block, rows, np.nextafter(cut, -np.inf), "S", 1)
    # The gene whose correlation equals the cut must be excluded, so the two
    # scores must differ. An inclusive comparison makes them identical.
    assert not np.allclose(at_threshold, just_below)
    assert np.allclose(at_threshold, block[rows][cor > cut].mean(axis=0))


def test_nan_scores_raise_rather_than_assign_silently():
    """One NaN stage poisons every stage for that cell, and R errors there.

    R's double `scale()` takes a mean across the five stages per cell, so a
    single NA spreads to all five. The decision then has no answer, and R
    stops with "missing value where TRUE/FALSE needed".
    """
    bad = np.full((1, 5), np.nan)
    with pytest.raises(ValueError, match="NaN"):
        _assign(bad, STAGES, tolerance=0.3)


def test_one_nan_stage_spreads_to_the_whole_cell():
    raw = np.array([[1.0, 2.0, np.nan, 4.0, 5.0],
                    [2.0, 1.0, np.nan, 3.0, 6.0],
                    [3.0, 5.0, np.nan, 2.0, 1.0]])
    scaled = _scale_columns(_scale_columns(raw).T).T
    assert np.isnan(scaled).all()


def test_scale_columns_matches_r_scale_with_n_minus_one_denominator():
    m = RNG.normal(size=(20, 3)) * 5 + 2
    got = _scale_columns(m)
    want = (m - m.mean(axis=0)) / m.std(axis=0, ddof=1)
    assert np.allclose(got, want, atol=1e-14)


def test_marker_correlation_floor_drops_uncorrelated_genes():
    """A marker that does not track its stage mean must not raise the score."""
    adata, cycle_genes, truth = _synthetic()
    # Replace one G2 marker with pure noise.
    adata.X[:, 20] = RNG.normal(0.0, 1.0, adata.n_obs)
    tp.estimate_schwabe_stage(adata, cycle_genes=cycle_genes, cor_thres=0.2)
    called = _called(adata)
    assigned = called != None                                     # noqa: E711
    assert assigned.mean() > 0.9
    assert (called[assigned] == truth[assigned]).all()


def test_a_correlation_floor_above_one_warns_and_then_has_no_answer():
    """Nothing can correlate above 1, so every stage score becomes NaN."""
    adata, cycle_genes, _ = _synthetic()
    with pytest.warns(RuntimeWarning, match="too few genes"), \
            pytest.raises(ValueError, match="NaN"):
        tp.estimate_schwabe_stage(adata, cycle_genes=cycle_genes, cor_thres=1.5)
