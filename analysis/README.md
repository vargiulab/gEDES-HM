# `analysis/`

The computed data behind the numbers reported in the revision, and the scripts
that produce it. No figure-generation code and no rendered figures: the figures
are in the paper and its Supporting Information and can be drawn from these
tables in any plotting software.

No absolute paths. Each script finds its inputs through `common/adk_paths.py`,
which derives the archive root from its own position; set `ADK_ARCHIVE_ROOT`
only if you separate `analysis/` from the rest of the archive. Python 3.9 or
newer with `numpy`, `scipy` and `pandas`; `cpptraj` for `tableS2_rmsd_check/`.

## What is here

| Directory | Key files | Manuscript item |
| --- | --- | --- |
| `common/` | `adk_geom.py` (the coordinate definitions), `adk_paths.py` | |
| `structures/` | input coordinates in E. coli numbering | |
| `collective_coordinates/` | `adk_coordinates.tsv`, `adk_angles.tsv`, `nearest_representative.tsv` | Tables S3, S4, Figure 2 |
| `ensemble_transfer/` | `model_coordinates.tsv`, `transfer_stats.txt` | Figure S2 |
| `tableS2_rmsd_check/` | `tableS2_published.tsv` vs `tableS2_recomputed.tsv` | Table S2 |
| `tableS6_check/` | `check_tableS6.sh` | Table S6 |
| `conformer_analysis/` | `ensemble_size_analysis.py` | Table S7 |
| `statistics/` | `gEDES_docking_data.csv`, `gEDES_SI_statistics.py` | Main text, Results |

## The coordinates

Both established descriptions of the ADK transition, in E. coli (1AKE/4AKE)
numbering, defined once in `common/adk_geom.py`:

- **Distances** `dLID` (d1) and `dNMP` (d2), mass-weighted centres of mass over
  all atoms, LID 118-160, NMP 30-67, CORE 1-29, 68-117, 161-214.
- **Angles** `theta1` (LID-CORE) and `theta2` (NMP-CORE), backbone atoms only,
  LID 123-155, NMP 50-59, hinge 161-165, CORE 1-8, 79-85, 104-110, 190-198,
  following Li, Liu and Ji, *Biophys. J.* **2015**, *109*, 647.

Every row also carries `cov_*`, the fraction of each domain present. This is
what separates the six short-form targets (systems 15 to 20), whose LID is
largely absent from the sequence, from long-form structures with a few
unresolved residues.

## The data tables

`adk_coordinates.tsv` (Table S3), 132 rows: the two E. coli references, seven
experimental intermediates, 22 homologue holo structures, 101 gEDES
representatives. Columns `label`, `role`, `dLID`, `dNMP`, `theta1`, `theta2`,
`cov_*`.

`model_coordinates.tsv` (Figure S2), 2,222 rows: the same coordinates for all
101 homology models of each of the 22 targets. Columns `target`, `model`, then
as above. `transfer_stats.txt` is the per-target summary and the source of the
median 0.16 Å in the Figure S2 caption.

`nearest_representative.tsv` (Table S4), 16 rows: for each long-form homologue,
the member of the 101 closest to it over the LID. Columns `system`, `tau` (the
closure coordinate, 0 at 4AKE and 1 at 1AKE), `nearest_representative`,
`LID_CA_RMSD`, `delta_tau`. The comparison was made after placing all
structures in a common reference frame built from the experimental structures
alone, not by superposing each pair individually. Table S4 prints the middle
two columns; `tau` and `delta_tau` are kept here as provenance.

`gEDES_docking_data.csv` (the binding-site and sequence-identity correlations
quoted in the Results), 26 rows: `rmsd_BM`, `rmsd_lig`,
`Fnat`, `rmsd_BS`, `rmsd_prot` and the two sequence identities per complex.

## Regenerating

```bash
cd collective_coordinates && python3 compute_adk_coordinates.py
cd ../ensemble_transfer   && python3 compute_models.py
cd ../statistics          && python3 gEDES_SI_statistics.py
cd ../tableS6_check       && bash check_tableS6.sh
```

The first two reproduce their shipped tables byte for byte. `check_tableS6.sh`
gives 24/26 (92.3%), 2/26 (7.7%), 24/26 (92.3%) and 25/26 (96.2%) for the four
receptor conditions, with complexes 12 and 14 unrecovered for the first and
third and 12 alone for the redocking control.

`tableS2_rmsd_check/` is a comparison rather than a byte-for-byte
reproduction: it holds Table S2 as printed in the SI next to a recomputation
from the archived ensembles (`run_all.sh`, requires cpptraj). Every
minimum-RMSD value agrees exactly. The percentage columns differ by up to six
percentage points, and the whole-protein columns of the second ligand of the
dual-ligand systems 2, 9, 15 and 21 are left empty in the published table,
which counts whole-protein quantities once per target, but are filled in the
recomputation, which counts them once per complex. Neither difference changes
a reported result: both give 84.6% of complexes with a binding-site model
within 2.5 Angstrom for the 101 DOPE-selected models and 92.3% for all 1010.
