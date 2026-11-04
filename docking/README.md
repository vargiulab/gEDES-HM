# `docking/`

Final HADDOCK3 docking results for every complex, under all four receptor
conditions compared in the manuscript. All files sit directly in this
directory; the receptor condition is carried by the filename prefix.

## Filename prefixes

| Prefix | Receptor condition | Receptors per target |
| --- | --- | --- |
| *(none)* | **HM_gEDES**, homology models built from the gEDES representatives of the reference protein, DOPE-selected | 101 |
| `APO_` | **HM_4AKE(apo)**, homology models built on the single experimental apo template, PDB 4AKE | 10 |
| `HOLO-1ake_` | **HM_1AKE(holo)**, homology models built on the single experimental holo template, PDB 1AKE | 10 |
| `redocking_HOLO-expt_` | **HOLO_expt**, the target's own experimental holo receptor with the ligand removed, redocked | 1 |

`HM_1AKE(holo)` is not a condition a prospective application would have, since
it presupposes an experimental bound structure of the reference protein. It is
included as the demanding upper comparison at the modelling stage.
`HOLO_expt` is the redocking control: it defines the practical ceiling of the
docking protocol on this ligand set and the noise floor of the RMSD metrics.

## Filename pattern

```
[<prefix>]<system>_cluster_ranks[_<ligand>].tsv
[<prefix>]<system>_best_docking_complex[_<ligand>].pdb
```

`<system>` is 1 to 22 and matches the folders of
`proteins_structural_models/`. Four targets carry two ligands each and were
docked as separate complexes, so the ligand suffix distinguishes them and the
22 systems give 26 complexes:

| System | Ligand suffixes |
| --- | --- |
| 2 | `_ADP`, `_AMP` |
| 9 | `_ADP`, `_AMP` |
| 15 | `_ADP182`, `_ADP183` |
| 21 | `_ADP214`, `_ADP215` |

Each condition therefore has 26 `*_cluster_ranks*.tsv` files. Best-complex PDBs
are provided for the `HM_gEDES` condition only.

## `*_cluster_ranks*.tsv`

One row per pose cluster, ranked by HADDOCK score. The clusters are **cpptraj
(Amber) clusters** of the HADDOCK `emref` poses, not HADDOCK clusters: the
docking run ends at `emref` and no HADDOCK clustering module is executed. See
`scripts/clusters.cpptraj` (`hieragglo`, `epsilon 1.5`, on binding-site plus
ligand heavy atoms). Within each cluster the poses are ranked by HADDOCK score,
and the cluster is ranked by the mean of its three best scores.

| Column | Meaning |
| --- | --- |
| `rank` | cluster rank by HADDOCK score |
| `cluster_id` | cpptraj cluster identifier |
| `size` | number of members in the cluster |
| `avg_top3_score` | mean HADDOCK score of the three best members |
| `pose_score` | HADDOCK score of the representative pose |
| `pose_combined_BS-lig` | binding-mode RMSD of the representative pose, RMSD_BM, over binding-site plus ligand heavy atoms, symmetry-corrected |
| `pose_lig_rmsd` | ligand heavy-atom RMSD, symmetry-corrected |
| `pose_bs_rmsd` | binding-site heavy-atom RMSD |
| `pose_wholepro_rmsd` | whole-protein Cα RMSD |
| `pose_id` | pose identifier, e.g. `emref_1304` |
| `pose_fnat` | fraction of native contacts |

All RMSD values are against the experimental holo reference of that complex,
the same `ref_clean*.pdb.gz` used in `proteins_structural_models/`.

The number of clusters scales with the size of the receptor ensemble, from a
few thousand for the 101-conformer `HM_gEDES` runs down to a few tens for the
single-receptor redocking control.

## `*_best_docking_complex*.pdb`

The best docked receptor-ligand complex for that complex under the `HM_gEDES`
condition. The `REMARK FILENAME` line gives the HADDOCK pose id (`emref_<N>`),
which is the `pose_id` of that pose in the matching `*_cluster_ranks*.tsv`. The
ligand is written with residue name `UNK`, the HADDOCK convention.

## Reproducing Table S6

`analysis/tableS6_check/check_tableS6.sh` reads these tables and reports, for
each condition, the success rate at RMSD_BM <= 2.5 Å, the median RMSD_BM and
which complexes are unrecovered.

Note that Table S6 of the SI prints the representative pose of the cluster
containing the best native-like solution, chosen on all metrics together rather
than by minimising any single one. A plain minimum over `pose_combined_BS-lig`
therefore returns values at or below the tabulated ones; the success counts and
the identity of the unrecovered complexes are unaffected.
