"""Copy every figure from reports/ into figures/, grouped by analysis part.

The analysis scripts write their figures next to their tables in reports/. This script gathers them in one place
(figures/<part>/<same file name>) so they can be browsed on GitHub; see figures/README.md for what each shows.
Run it after re-running an analysis to refresh the copies.

Usage: python scripts/collect_figures.py [project_root]
"""
import shutil
import sys
from pathlib import Path

ROOT = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(".")

PARTS = {  # figures/<part>  <-  reports/<folder>
    "01_pca_embeddings": "ops0077_cell_cycle_clock",
    "02_score_markers_dna": "ops0077_cell_cycle_score",
    "02_score_markers_only": "ops0077_cell_cycle_score_only_markers",
    "03_phase_contrast_prediction": "ops0077_phase_predicts_cc_phase",
    "04_knockout_evaluation": "ops0077_cc_transfer/ko3_evaluation",
    "04_knockout_guide_analysis": "ops0077_cc_transfer/ko3_guide_analysis",
    "04_knockout_density_profiles": "ops0077_cc_transfer/ko3_density_profiles",
}


def main():
    n = 0
    for part, folder in PARTS.items():
        src, dst = ROOT / "reports" / folder, ROOT / "figures" / part
        dst.mkdir(parents=True, exist_ok=True)
        for png in sorted(src.glob("*.png")):
            shutil.copy2(png, dst / png.name)
            n += 1
        print(f"{part:32s} <- reports/{folder}")
    print(f"{n} figures copied to {ROOT / 'figures'}")


if __name__ == "__main__":
    main()
