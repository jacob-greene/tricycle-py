# Identifier code paths, scored against R

**Verdict: PASS.**
The main equivalence run exercises `species="human"` with `gname_type="SYMBOL"`. These arms exercise the others, including both paths where R queries `AnnotationDbi` at run time and this package reads a frozen extract instead.

| Arm | Cells | Projection genes, Python / R | Median circular difference | Maximum | Maximum arc error |
|---|---:|---|---|---|---|
| mouse_ensembl | 400 | 500 / 500 | 0.000e+00 | 3.997e-15 | 1.446e-14 |
| mouse_symbol | 400 | 500 / 500 | 0.000e+00 | 3.997e-15 | 1.446e-14 |
| human_ensembl | 8627 | 437 / 437 | 0.000e+00 | 6.173e-14 | 1.263e-14 |

Bounds are the ones `compare.py` fixes: median 1e-09 rad, maximum 1e-06 rad, arc 1e-09.

