# Scripts: homology modeling and HADDOCK3 docking workflow

This folder holds the scripts used for homology modeling, docking, and post-processing
in this work. Each script carries detailed usage notes and rationale in its own header;
the summary below is the high-level map of how they fit together.

Binding-site (BS) residue definitions are the same for a given target across all stages:
homology-modeling RMSD, docking, and docking post-processing. The authoritative
per-system definition is the `BS_RESIDUES` value inside each system's
`run_calc_rmsd-cpptraj.sh` (under `proteins_structural_models/`).

## Script inventory

| Script | Stage | Purpose |
| --- | --- | --- |
| `DOPE_score_select.sh` | Homology modeling | Pick the lowest-DOPE model per `results_*` directory. |
| `copy_best_structures.sh` | Homology modeling | Copy the selected best models into `selected_best_structures/`. |
| `cpptraj_rmsd_HM.sh` | Homology modeling | WholePro (Cα) and BS (heavy-atom) RMSD of models vs reference. |
| `haddock_run_redocking_bs.sh` | Docking | HADDOCK3 redocking / validation against experimental holo structure. |
| `haddock_run_ensemble_crossdocking_bs.sh` | Docking | HADDOCK3 ensemble cross-docking of all homology models in one run (ends at `emref`; no HADDOCK clustering module). |
| `clusters.cpptraj` | Post-processing | Post-run cpptraj (Amber) clustering of docked poses (per-system BS definition). |
| `ap5_symswap_haddock.py` | Post-processing | Build a C2-symmetry-swapped dual reference for the AP5 ligand. |
| `g5p_symswap_haddock.py` | Post-processing | Build a C2-symmetry-swapped dual reference for the G5P ligand. |
| `cpptraj_rmsd.sh` | Post-processing | Symmetry-corrected RMSD (BS, ligand, whole-protein, combined) of docked poses. |

## Workflow: homology modeling

1. **Generate models.** Each system folder under
   `proteins_structural_models/<setup>/<n>/` provides a MODELLER alignment
   (`<n>_adk_<ID>.pir`) and driver (`modeller_run.py`). Running `modeller_run.py`
   generates the homology models for that system: for each template it builds
   `MODELS_PER_TMPL = 10` models into a `results_<i>/` folder. The two setups differ only
   in the template list at the top of the driver:

   - **gEDES**: `TEMPLATE_DIR = '../reference_gEDES_templates/'`, `N_TEMPLATES = 101`
     (`0.pdb` ... `100.pdb`) → 101 x 10 = 1010 models per target, from which the
     DOPE-selected ensemble of 101 models is taken (steps 2-3).
   - **apo**: `TEMPLATES = ['../reference_apo_template_4ake.pdb']`, the single apo
     template (PDB 4AKE) → 10 models per target.

   Both paths are relative to the system folder the driver sits in, so the drivers run
   straight from the archive. The templates ship gzipped to keep the deposit small;
   decompress them in place first:

   ```bash
   gunzip proteins_structural_models/gEDES/reference_gEDES_templates/*.gz
   gunzip proteins_structural_models/apo/reference_apo_template_4ake.pdb.gz
   ```

2. **Select by DOPE.** `DOPE_score_select.sh` scans the `results_*` model directories and
   writes `lowest_dope_scores.txt`, the lowest-DOPE model per directory (best 1 in 10).

3. **Collect best models.** `copy_best_structures.sh` copies the lowest-DOPE structures
   listed in `lowest_dope_scores.txt` into a `selected_best_structures/` folder, renamed
   by their source index.

4. **RMSD analysis.** `cpptraj_rmsd_HM.sh` builds the concat files and computes WholePro
   (Cα) and BS (heavy-atom) RMSD versus the reference. Run it per system through the
   provided wrapper, which supplies that system's BS definition and protein range, e.g.:

   ```bash
   # proteins_structural_models/gEDES/1/run_calc_rmsd-cpptraj.sh
   BS_RESIDUES="9-15,31,32,36,57,59,85,88,92,119,123,132,133,156,167,200,202" \
     PROT_FIRST=1 PROT_LAST=214 bash cpptraj_rmsd_HM.sh
   ```

## Workflow: docking (HADDOCK3)

1. **Redocking / validation.** `haddock_run_redocking_bs.sh` builds the ligand topologies
   and the HADDOCK3 inputs, and validates docking against the experimental holo structure.
   It is the generic per-model docking script.

2. **Production docking.** The DOPE-selected best structures are used as input for the
   gEDES ensembles; the same protocol is applied to the apo setups.
   `haddock_run_redocking_bs.sh` is the generic driver;
   `haddock_run_ensemble_crossdocking_bs.sh` docks all homology models of a target together
   as a single ensemble. For flipped-ligand (symmetric) targets, the post-processing dual
   references are built with `ap5_symswap_haddock.py` (AP5) and `g5p_symswap_haddock.py`
   (Ap4A/B4P).

   The ensemble run executes the modules `topoaa` → `rigidbody` → `seletop` → `flexref`
   → `emref` → `caprieval`, with rigid-body sampling scaled so that each ensemble
   conformer is sampled ~250 times. Energy-minimisation refinement (`emref`) is the final
   refinement stage.

3. **Clustering.** Pose clustering is done afterwards with Amber/cpptraj:
   `clusters.cpptraj` clusters the `emref` poses on the binding-site + ligand
   heavy atoms (`hieragglo`, `epsilon 1.5`, best-fit on the binding site) and writes the
   cluster members, representatives, and summary used to build the
   `docking/*_cluster_ranks.tsv` tables. The file ships with the BS definition of
   system 22 as an example; update the BS residue selection (both the `rms BSfit` and
   the `cluster ... rms` masks) for each system before running.

4. **RMSD analysis.** `cpptraj_rmsd.sh` computes, per docked frame, the BS heavy-atom,
   ligand heavy-atom, whole-protein Cα, and binding-mode (combined BS+ligand) RMSD versus 
   the reference, with symmetry correction.

## Requirements

Exact versions and options are documented in each script header. In broad terms the
workflow uses:

- **MODELLER** for homology modeling (`modeller_run.py`).
- **HADDOCK3** for docking.
- **AmberTools `cpptraj`** for clustering and RMSD.
- A POSIX shell (`bash`) and standard Unix utilities.
