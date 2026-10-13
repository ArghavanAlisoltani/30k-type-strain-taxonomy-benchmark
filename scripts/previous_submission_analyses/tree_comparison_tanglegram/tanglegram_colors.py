#!/usr/bin/env python3
"""
=============================================================================
tanglegram_colors.py  —  MANUAL color definitions for all tanglegram figures
=============================================================================

Single source of truth for every color used by:
    step3_tanglegrams_all_rootings.py
    step4_tanglegram_metadata_strips_all_rootings.py
    step5_metric_summary_figures.py

UPDATE (this revision)
----------------------
1. Colors below are the user-supplied palettes (domain / kingdom / phylum),
   applied verbatim.
2. NEW fourth metadata strip: "W" = Carl Woese's historical group
   (column `carl_woese_historical_group` in the new metadata CSV). Colored by
   WOESE_COLORS below.
3. Tip labels are now BLACK (TIP_LABEL_COLOR), not colored by domain.
4. Metadata now comes from metadata_with_official_phyla_and_woese_groups.csv,
   which has no missing kingdoms — so there is no longer an "Unknown" kingdom.
   ("Unknown" entries are kept in the maps only as a safety fallback.)

HOW TO RECOLOR
--------------
Edit any hex below and re-run step3 / step4 / step5. Any taxon not listed,
or with blank/NA metadata, is drawn in MISSING_COLOR.
=============================================================================
"""

# Fallback for unlisted / blank / NA values.
MISSING_COLOR = "#BDBDBD"

# Tip-label text color. Black per request (no domain-based tip coloring).
TIP_LABEL_COLOR = "#000000"

# Which column colors the connecting lines between the two trees.
# Options: "domain" | "kingdom" | "phylum" | "woese". Lines stay informative;
# only the tip-label TEXT is forced black (above).
LINE_COLOR_BY = "domain"


# ---------------------------------------------------------------------------
# DOMAIN   (lowercase aliases included so harmonized/raw spellings both work)
# ---------------------------------------------------------------------------
DOMAIN_COLORS = {
    "Archaea":  "#FF7F00",
    "Archea":   "#FF7F00",
    "Bacteria": "#2171B5",
    "bacteria": "#2171B5",
    "Unknown":  "#BDBDBD",
}


# ---------------------------------------------------------------------------
# KINGDOM
# ---------------------------------------------------------------------------
KINGDOM_COLORS = {
    "Bacillati":         "#08519C",
    "Fusobacteriati":    "#C6DBEF",
    "Methanobacteriati": "#FFA600",
    "Nanobdellati":      "#FB5607",
    "Promethearchaeati": "#FFCC99",
    "Pseudomonadati":    "#6BAED6",
    "Thermoproteati":    "#fdbf6f",
    "Thermotogati":      "#3A86FF",
    "Unknown":           "#BDBDBD",
}


# ---------------------------------------------------------------------------
# PHYLUM   (phylum_order zipped with phylum_colors, in declaration order)
# ---------------------------------------------------------------------------
_PHYLUM_ORDER = [
    "Abditibacteriota", "Acidobacteriota", "Actinomycetota", "Aquificota",
    "Armatimonadota", "Atribacterota", "Bacillota", "Bacteroidota",
    "Balneolota", "Bdellovibrionota", "Caldisericota", "Calditrichota",
    "Campylobacterota", "Chlamydiota", "Chlorobiota", "Chloroflexota",
    "Chrysiogenota", "Coprothermobacterota", "Cyanobacteriota", "Deferribacterota",
    "Deinococcota", "Desulfobacterota", "Dictyoglomerota", "Elusimicrobiota",
    "Fibrobacterota", "Fidelibacterota", "Fusobacteriota", "Gemmatimonadota",
    "Halobacteriota",
    "Ignavibacteriota", "Kiritimatiellota", "Lentisphaerota", "Methanobacteriota",
    "Microcaldota", "Minisyncoccota", "Mycoplasmatota", "Myxococcota",
    "Nanobdellota", "Nitrososphaerota", "Nitrospinota", "Nitrospirota",
    "Planctomycetota", "Promethearchaeota", "Pseudomonadota", "Rhodothermota",
    "Spirochaetota", "Synergistota", "Thermodesulfobacteriota", "Thermodesulfobiota",
    "Thermomicrobiota", "Thermoplasmatota", "Thermoproteota", "Thermosulfidibacterota",
    "Thermotogota", "Verrucomicrobiota", "Vulcanimicrobiota",
]
_PHYLUM_COLORS_LIST = [
    "#08306B", "#08519C", "#2171B5", "#4292C6", "#6BAED6", "#9ECAE1", "#C6DBEF", "#807DBA",
    "#2F4B7C", "#264653", "#005F73", "#74A9CF", "#A6BDDB", "#D0D1E6", "#00441B", "#006D2C",
    "#238B45", "#41AB5D", "#74C476", "#A1D99B", "#C7E9C0", "#7F7F7F", "#000000", "#06D6A0",
    "#66C2A4", "#99D8C9", "#CCECE6", "#3F007D", "#FFF3B0", "#54278F", "#6A51A3", "#EF476F",
    "#FFE66D", "#FFD166", "#9E9AC8", "#BCBDDC", "#DADAEB", "#FFCA3A", "#FFA600", "#FF9DA7",
    "#E377C2", "#9D4EDD", "#F4A261", "#D45087", "#C77DFF", "#E0AAFF", "#253494", "#2C7FB8",
    "#B07AA1", "#A05195", "#FF7F0E", "#FB5607", "#2BB3A3", "#1D91C0", "#756BB1", "#17BECF",
]
PHYLUM_COLORS = dict(zip(_PHYLUM_ORDER, _PHYLUM_COLORS_LIST))
PHYLUM_COLORS["Unknown"] = "#BDBDBD"


# ---------------------------------------------------------------------------
# CARL WOESE historical group  (the new 4th strip, "W")
# Values come from `carl_woese_historical_group`. "12. Archaea" is orange to
# match the Archaea domain color; the catch-all group is light grey so the
# 12 historically-named groups stand out. Edit freely.
# ---------------------------------------------------------------------------
WOESE_COLORS = {
    "1. Purple Bacteria":                   "#9D4EDD",
    "2. Gram-positive High G+C":            "#1F78B4",
    "3. Gram-positive Low G+C":             "#A6CEE3",
    "4. Cyanobacteria":                     "#33A02C",
    "5. Spirochetes":                       "#B2DF8A",
    "6. Green Sulfur Bacteria":             "#006D2C",
    "7. Bacteroides-Flavobacteria":         "#FF7F0E",
    "8. Planctomyces":                      "#E31A1C",
    "9. Chlamydiae":                        "#FB9A99",
    "10. Radio-resistant Micrococci":       "#6A3D9A",
    "11. Green Non-sulfur Bacteria":        "#B15928",
    "12. Archaea":                          "#FF7F00",
    "Outside Carl Woese 12-group mapping":  "#E8E8E8",
    "Unknown":                              "#BDBDBD",
}


# ---------------------------------------------------------------------------
# STEP 5 kingdom-monophyly heatmap status colors.
# ---------------------------------------------------------------------------
STATUS_COLORS = {
    "neither": "#F0F0F0",
    "16S":     "#80B1D3",
    "FastAAI": "#FDB462",
    "both":    "#66C2A5",
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
_BLANKS = {"", "nan", "NaN", "None", "NA", "na", "N/A", "Unknown"}


def resolve(value, mapping):
    """Manual color for `value` from `mapping`, else MISSING_COLOR."""
    if value is None:
        return MISSING_COLOR
    key = str(value).strip()
    return mapping.get(key, MISSING_COLOR)


def color_map_for_column(column):
    """Return the manual color map for a metadata column name."""
    table = {
        "domain":  DOMAIN_COLORS,
        "kingdom": KINGDOM_COLORS,
        "phylum":  PHYLUM_COLORS,
        "woese":   WOESE_COLORS,
    }
    if column not in table:
        raise ValueError(f"no manual color map for column {column!r}")
    return table[column], MISSING_COLOR


def present_color_map(values, mapping):
    """Subset `mapping` to categories present in `values`, preserving the
    declaration order of `mapping`. Present-but-unmapped values are appended."""
    seen = []
    for v in values:
        if v is None:
            continue
        k = str(v).strip()
        if k in _BLANKS or k in seen:
            continue
        seen.append(k)
    out = {}
    for k in mapping:
        if k in seen:
            out[k] = mapping[k]
    for k in seen:
        if k not in out:
            out[k] = MISSING_COLOR
    return out
