"""Path resolution for the analysis scripts of this archive.

No absolute path is hard-coded anywhere in `analysis/`. Every script that needs
a file outside its own directory asks this module for it, and this module
resolves the location in two steps:

  1. an environment variable, if the user sets one;
  2. otherwise a location derived from the position of this file, which is
     correct as long as the archive is unpacked with its layout intact.

Roots
-----
ARCHIVE_ROOT   the top of the archive, i.e. the directory containing
               `proteins_structural_models/`, `docking/`, `scripts/` and
               `analysis/`.                       env: ADK_ARCHIVE_ROOT
DOCKING_ROOT   the directory holding the per-complex HADDOCK3 output used by
               the pose-level analyses. Defaults to <ARCHIVE_ROOT>/docking; set
               it to the raw run directory if you re-run those analyses on the
               full HADDOCK3 output rather than on the archived subset.
                                                  env: ADK_DOCKING_ROOT

Everything is built with os.path.join, so the module works unchanged on Linux,
macOS and Windows.
"""
import os

_HERE = os.path.dirname(os.path.abspath(__file__))          # analysis/common
ANALYSIS_ROOT = os.path.abspath(os.path.join(_HERE, os.pardir))
ARCHIVE_ROOT = os.path.abspath(
    os.environ.get("ADK_ARCHIVE_ROOT")
    or os.path.join(ANALYSIS_ROOT, os.pardir))

MODELS_ROOT = os.path.join(ARCHIVE_ROOT, "proteins_structural_models")
GEDES_MODELS = os.path.join(MODELS_ROOT, "gEDES")
APO_MODELS = os.path.join(MODELS_ROOT, "apo")
GEDES_TEMPLATES = os.path.join(GEDES_MODELS, "reference_gEDES_templates")

DOCKING_ROOT = os.path.abspath(
    os.environ.get("ADK_DOCKING_ROOT") or os.path.join(ARCHIVE_ROOT, "docking"))

STRUCTURES = os.path.join(ANALYSIS_ROOT, "structures")
ECOLI = os.path.join(STRUCTURES, "ecoli")
HOMOLOGUES = os.path.join(STRUCTURES, "homologues_ecnum")
INTERMEDIATES = os.path.join(STRUCTURES, "intermediates_ecnum")


def add_common_to_path():
    """Convenience for scripts: put analysis/common on sys.path."""
    import sys
    if _HERE not in sys.path:
        sys.path.insert(0, _HERE)


def require(path, what):
    """Fail with a message that says which root to set, not a bare traceback."""
    if not os.path.exists(path):
        raise SystemExit(
            f"missing {what}: {path}\n"
            "The archive layout appears to be incomplete. Either unpack the "
            "whole archive, or point ADK_ARCHIVE_ROOT (and, for docking "
            "analyses, ADK_DOCKING_ROOT) at the right directories."
        )
    return path
