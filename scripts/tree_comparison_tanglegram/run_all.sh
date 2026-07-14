#!/usr/bin/env bash
set -euo pipefail

# v4 pipeline. Colors for all tanglegrams/heatmaps are defined manually in
# tanglegram_colors.py — edit that file to recolor, then re-run steps 3–5.

python3 step1_reduce_16s.py
python3 step2_compare_all_rootings.py
python3 step3_tanglegrams_all_rootings.py
python3 step4_tanglegram_metadata_strips_all_rootings.py
python3 step5_metric_summary_figures.py
