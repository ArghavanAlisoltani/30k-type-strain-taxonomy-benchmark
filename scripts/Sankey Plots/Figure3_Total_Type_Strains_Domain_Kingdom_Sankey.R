# ============================================================
# REVISED FIGURE 3: Total Type Strains -> Domains -> Kingdoms
# Reviewer fixes:
# 1. Thicker thin lines
# 2. Higher contrast / more reader-friendly palette
# 3. Better visibility for small groups
# 4. Publication-quality output
# ============================================================

rm(list = ls())

library(tidyverse)
library(stringr)
library(ggplot2)
library(scales)
library(svglite)
library(grid)

# ============================================================
# INPUT FILE
# ============================================================

taxonomy_file <- "/Users/akhimomtaz/Desktop/figure_fix/final_clean_taxonomy_30495.csv"

# ============================================================
# OUTPUT DIRECTORY
# ============================================================

outdir <- "FIGURE3_REVISED"

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

# ============================================================
# COLORS
# Cleaner, higher contrast, print-friendly
# ============================================================

col_total_box    <- "#4B5563"
col_archaea_box  <- "#B23A48"
col_bacteria_box <- "#3B6FB6"
col_kingdom_box  <- "#6B7280"

flow_colors <- c(
  "Archaea"               = "#D77A84",
  "Bacteria"              = "#7FA6D9",
  "Methanobacteriati"     = "#C97A2B",
  "Nanobdellati"          = "#D9A404",
  "Promethearchaeati"     = "#B65AA3",
  "Thermoproteati"        = "#D94C4C",
  "Bacillati"             = "#2FA57A",
  "Pseudomonadati"        = "#6A5ACD",
  "Fusobacteriati"        = "#1FA9C9",
  "Thermotogati"          = "#D64FA3",
  "Unclassified Bacteria" = "#7DAA3C",
  "Unclassified Archaea"  = "#9C6644"
)

# ============================================================
# HELPER FUNCTIONS
# ============================================================

clean_names_lower <- function(df) {
  names(df) <- names(df) |>
    str_trim() |>
    str_to_lower() |>
    str_replace_all("\\s+", "_")
  df
}

curve_fun <- function(x0, y0, x1, y1, id, val, col) {
  t <- seq(0, 1, length.out = 220)

  tibble(
    x = (1 - t)^3 * x0 +
      3 * (1 - t)^2 * t * (x0 + 0.68 * (x1 - x0)) +
      3 * (1 - t) * t^2 * (x1 - 0.68 * (x1 - x0)) +
      t^3 * x1,

    y = (1 - t)^3 * y0 +
      3 * (1 - t)^2 * t * y0 +
      3 * (1 - t) * t^2 * y1 +
      t^3 * y1,

    id = id,
    val = val,
    col = col
  )
}

# ============================================================
# READ TAXONOMY
# ============================================================

taxonomy_df <- read_csv(taxonomy_file, show_col_types = FALSE)
taxonomy_df <- clean_names_lower(taxonomy_df)

required_cols <- c("accession", "domain", "kingdom")
missing_cols <- setdiff(required_cols, names(taxonomy_df))

if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

taxonomy_df <- taxonomy_df %>%
  mutate(
    accession = str_trim(as.character(accession)),
    domain    = str_trim(as.character(domain)),
    kingdom   = str_trim(as.character(kingdom)),

    kingdom = case_when(
      domain == "Bacteria" &
        kingdom %in% c(
          "Bacteria incertae sedis",
          "Unclassified_Deltaproteobacteria"
        ) ~ "Unclassified Bacteria",

      domain == "Archaea" &
        !kingdom %in% c(
          "Methanobacteriati",
          "Nanobdellati",
          "Promethearchaeati",
          "Thermoproteati"
        ) ~ "Unclassified Archaea",

      TRUE ~ kingdom
    )
  ) %>%
  filter(
    accession != "",
    domain != "",
    kingdom != ""
  ) %>%
  distinct(accession, .keep_all = TRUE)

cat("\nTotal rows used:", nrow(taxonomy_df), "\n")

# ============================================================
# COUNTS
# ============================================================

total_count <- n_distinct(taxonomy_df$accession)

domain_df <- taxonomy_df %>%
  count(domain, name = "total_genomes") %>%
  mutate(
    y = case_when(
      domain == "Archaea"  ~ 5.25,
      domain == "Bacteria" ~ 0.35,
      TRUE ~ 0
    ),
    box_fill = case_when(
      domain == "Archaea"  ~ col_archaea_box,
      domain == "Bacteria" ~ col_bacteria_box,
      TRUE ~ "#999999"
    ),
    flow_col = flow_colors[domain]
  )

kingdom_df <- taxonomy_df %>%
  count(domain, kingdom, name = "total_genomes") %>%
  mutate(
    y = case_when(
      kingdom == "Methanobacteriati"     ~ 8.05,
      kingdom == "Nanobdellati"          ~ 6.80,
      kingdom == "Promethearchaeati"     ~ 5.20,
      kingdom == "Thermoproteati"        ~ 3.70,

      kingdom == "Bacillati"             ~ 2.15,
      kingdom == "Pseudomonadati"        ~ 0.98,
      kingdom == "Fusobacteriati"        ~ -0.38,
      kingdom == "Thermotogati"          ~ -1.78,
      kingdom == "Unclassified Bacteria" ~ -3.28,
      kingdom == "Unclassified Archaea"  ~ 2.20,

      TRUE ~ 0
    ),
    flow_col = flow_colors[kingdom]
  ) %>%
  arrange(desc(y))

# ============================================================
# MANUAL LAYOUT
# ============================================================

x_total <- -1.10
y_total <- 2.45

bridge_xmin <- 0.18
bridge_xmax <- 0.50
bridge_ymin <- 1.68
bridge_ymax <- 3.20

x_from_total <- bridge_xmax

x_domain  <- 5.08
x_kingdom <- 9.18

title_x <- -2.45
title_y <- 9.75

# ============================================================
# BUILD CURVES
# ============================================================

curves <- list()

# Total -> Domain
for (i in seq_len(nrow(domain_df))) {
  curves[[length(curves) + 1]] <- curve_fun(
    x0  = x_from_total,
    y0  = y_total,
    x1  = x_domain,
    y1  = domain_df$y[i],
    id  = paste0("total_domain_", i),
    val = domain_df$total_genomes[i],
    col = domain_df$flow_col[i]
  )
}

# Domain -> Kingdom
for (i in seq_len(nrow(kingdom_df))) {
  y_dom <- domain_df$y[domain_df$domain == kingdom_df$domain[i]]

  curves[[length(curves) + 1]] <- curve_fun(
    x0  = x_domain,
    y0  = y_dom,
    x1  = x_kingdom,
    y1  = kingdom_df$y[i],
    id  = paste0("domain_kingdom_", i),
    val = kingdom_df$total_genomes[i],
    col = kingdom_df$flow_col[i]
  )
}

curve_df <- bind_rows(curves)

# ============================================================
# PLOT
# ============================================================

p <- ggplot() +

  # bridge block
  geom_rect(
    aes(
      xmin = bridge_xmin,
      xmax = bridge_xmax,
      ymin = bridge_ymin,
      ymax = bridge_ymax
    ),
    fill = col_total_box,
    color = NA,
    alpha = 0.98
  ) +

  # flow curves
  geom_path(
    data = curve_df,
    aes(
      x = x,
      y = y,
      group = id,
      color = col,
      linewidth = val
    ),
    alpha = 0.80,              # increased from faint look
    lineend = "round"
  ) +

  # total box
  geom_label(
    aes(
      x = x_total,
      y = y_total,
      label = paste0("Total type strains\n", comma(total_count))
    ),
    fill = col_total_box,
    color = "white",
    size = 6.2,
    fontface = "bold",
    label.size = 0,
    label.padding = unit(0.34, "lines"),
    lineheight = 0.92
  ) +

  # domain boxes
  geom_label(
    data = domain_df,
    aes(
      x = x_domain,
      y = y,
      label = paste0(domain, "\n", comma(total_genomes)),
      fill = box_fill
    ),
    color = "white",
    size = 5.7,
    fontface = "bold",
    label.size = 0,
    label.padding = unit(0.32, "lines"),
    lineheight = 0.92
  ) +

  # kingdom boxes
  geom_label(
    data = kingdom_df,
    aes(
      x = x_kingdom,
      y = y,
      label = paste0(kingdom, "\n", comma(total_genomes))
    ),
    fill = col_kingdom_box,
    color = "white",
    size = 4.3,
    fontface = "bold",
    label.size = 0,
    label.padding = unit(0.28, "lines"),
    lineheight = 0.92
  ) +

  # title
  annotate(
    "text",
    x = title_x,
    y = title_y,
    label = "Type Strain Genomes → Domains → Kingdoms",
    hjust = 0,
    size = 8.2,
    fontface = "bold",
    color = "black"
  ) +

  scale_color_identity() +
  scale_fill_identity() +

  # IMPORTANT: thicker minimum linewidth
  scale_linewidth(
    range = c(0.60, 28),   # increased minimum from very thin
    guide = "none"
  ) +

  coord_cartesian(
    xlim = c(-2.8, 11.2),
    ylim = c(-4.6, 10.5),
    clip = "off"
  ) +

  theme_void(base_size = 18) +
  theme(
    legend.position = "none",
    plot.margin = margin(55, 90, 55, 90),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    text = element_text(family = "Helvetica", face = "bold")
  )

print(p)

# ============================================================
# SAVE OUTPUTS
# ============================================================

ggsave(
  file.path(outdir, "Figure3_REVISED.png"),
  p,
  width = 16,
  height = 10,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(outdir, "Figure3_REVISED.pdf"),
  p,
  width = 16,
  height = 10,
  bg = "white"
)

ggsave(
  file.path(outdir, "Figure3_REVISED.svg"),
  p,
  width = 16,
  height = 10,
  bg = "white"
)

cat("\nSaved files to:", normalizePath(outdir), "\n")
