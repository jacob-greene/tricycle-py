# Attribution and licensing

## What this package is

`tricycle-py` is a Python port of the R/Bioconductor package
[`tricycle`](https://github.com/hansenlab/tricycle). It is a derivative work.
The method, the reference embedding and the marker lists are not original to
this package.

The port was written against `tricycle` 1.12.0 (Bioconductor 3.19), function by
function, from the R sources. Data objects were extracted from the installed R
package rather than regenerated, because a regenerated reference would not
reproduce the R numbers.

## Original authors

| Work | Authors | Reference |
|---|---|---|
| `tricycle` R package, as credited by the package itself | Shijie C. Zheng, sole author and maintainer | Upstream `DESCRIPTION`, `Authors@R` |
| The method the package implements | Shijie C. Zheng, Genevieve Stein-O'Brien, Jonathan J. Augustin, Jared Slosberg, Giovanni A. Carosso, Briana Winer, Gloria Shin, Hans T. Bjornsson, Loyal A. Goff, Kasper D. Hansen | Zheng SC et al. *Universal prediction of cell cycle position using transfer learning.* Genome Biology 23:41 (2022). doi:10.1186/s13059-021-02581-y |
| Five-stage assignment method | Daniel Schwabe et al. | Schwabe D et al. *The transcriptome dynamics of single cells during the cell cycle.* Molecular Systems Biology 16:e9946 (2020). doi:10.15252/msb.20209946 |
| `RevelioGeneList` marker genes | Michael L. Whitfield et al., by way of the `Revelio` package | Whitfield ML et al. *Identification of genes periodically expressed in the human cell cycle and their expression in tumors.* Molecular Biology of the Cell 13:1977-2000 (2002). doi:10.1091/mbc.02-02-0030 |
| Local regression | William S. Cleveland et al. | Cleveland WS, Grosse E, Shyu WM. *Local regression models.* In: Statistical Models in S (1992) |

## Copyright

    tricycle-py, a Python port of the R package tricycle.
    Copyright (C) 2026 Jacob Greene.

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the Free
    Software Foundation, version 3.

    This program is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
    more details.

    You should have received a copy of the GNU General Public License along
    with this program. If not, see <https://www.gnu.org/licenses/>.

Copyright in the original `tricycle` R package remains with its authors. This
notice covers only the port.

## Licence

`tricycle` is licensed **GPL-3**. Its `DESCRIPTION` file states `License: GPL-3`.

A port is a derivative work under copyright, so the GPL-3 terms carry over to
this package. Two consequences follow, and neither is optional:

1. **This package is licensed GPL-3.** A permissive licence such as MIT or
   BSD-3 is not available for it. The full licence text is in `LICENSE`.
2. **Redistribution must carry the same terms.** Anyone who distributes this
   package, modified or not, must do so under GPL-3, must keep these notices,
   and must offer the corresponding source.

Publishing this repository on GitHub is distribution, so the licence applies
from the first push.

### What GPL-3 does not restrict

Running the package on your own data, and publishing results computed with it,
carry no licence obligation. The GPL governs distribution of the software, not
the output of running it.

### Third-party data shipped in `src/tricyclepy/data/`

| File | Origin | Licence of origin |
|---|---|---|
| `neuroRef.tsv` | `tricycle` R package data object `neuroRef` | GPL-3 |
| `RevelioGeneList.tsv` | `tricycle` R package data object `RevelioGeneList` | GPL-3 |
| `neurosphere_example.npz` | `tricycle` R package data object `neurosphere_example`, with the `tricyclePosition` R computes on it | GPL-3 |
| `ensembl2symbol_human.tsv.gz` | `org.Hs.eg.db` 3.19.1, `mapIds(..., "SYMBOL", "ENSEMBL", multiVals = "first")` | Artistic-2.0 |
| `ensembl2symbol_mouse.tsv.gz` | `org.Mm.eg.db` 3.19.1, same query | Artistic-2.0 |
| `go_cc_{human,mouse}_{symbol,ensembl}.txt` | GO:0007049 membership via `org.*.eg.db` `GOALL` | Gene Ontology content is CC BY 4.0 |

Artistic-2.0 and CC BY 4.0 are both compatible with redistribution inside a
GPL-3 work. Attribution for all three is this file.

`PROVENANCE.tsv` in the same directory records the exact package versions the
extraction ran against. `equivalence/extract_reference_data.R` produced every file except
`neurosphere_example.npz`, which `equivalence/extract_example_data.R` and
`equivalence/build_example_data.py` produced.

## Relationship to the upstream project

This is an independent port. It is not endorsed by, affiliated with, or
supported by the `tricycle` authors or the Hansen lab. Report bugs in this port
here. Report method questions upstream.
