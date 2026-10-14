# ============================================================
# FIGURE 5 - ORIGINAL SANKEY STYLE, EXACT KINGDOM SPLIT
#
# Figure 5A: Archaea (original domain -> kingdom -> phylum view)
# Figure 5B: Bacillati, Fusobacteriati, Thermotogati,
#            Unclassified Bacteria (kingdom -> phylum)
# Figure 5C: all Pseudomonadati phyla (kingdom -> phylum)
#
# Visual design based directly on the supplied
# Figure5_centered_domains_no_all_genomes(3).R.
# Keep the original color-coded kingdom nodes, Bezier curves,
# status markers, four horizontal metric bars and right-side legend.
#
# Reviewer improvements:
#  - single-line phylum name + genome count
#  - readable text, strong vector Sankey connections
#  - "N/A" only when no same-species pairs are possible
#  - "0%" when pair comparisons exist and failure rate is zero
#  - PDF/SVG vector graphics and optional high-DPI PNG
# ============================================================

# ============================================================
# 1. PACKAGES
# ============================================================
library(tidyverse)
library(stringr)
library(ggplot2)
library(scales)
library(svglite)
library(grid)

# ============================================================
# 2. INPUT FILES
# ============================================================
# PORTABLE FILE DETECTION (Mac, Windows, Linux).
# Recommended: keep this R script and the four CSV files in one folder.
# Optional override: Sys.setenv(FIGURE5_DATA_DIR = "/full/path/to/CSVs")
#
# Works with both canonical CSV names and downloaded names such as:
#   final_clean_taxonomy_30495(3).csv
#   aai_mash_ani_16s(8).csv

get_script_dir <- function() {
  arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(arg)) {
    return(dirname(normalizePath(sub("^--file=", "", arg[1]), mustWork = FALSE)))
  }
  frames <- sys.frames()
  ofiles <- vapply(frames, function(f) {
    if (!is.null(f$ofile)) as.character(f$ofile)[1] else ""
  }, character(1))
  ofiles <- ofiles[nzchar(ofiles)]
  if (length(ofiles)) {
    return(dirname(normalizePath(tail(ofiles, 1), mustWork = FALSE)))
  }
  getwd()
}

script_dir <- get_script_dir()
env_dir <- Sys.getenv("FIGURE5_DATA_DIR", unset = "")
data_folders <- unique(c(
  if (nzchar(env_dir)) env_dir,
  script_dir,
  getwd(),
  file.path(script_dir, "plot_ready_16s_outputs"),
  file.path(getwd(), "plot_ready_16s_outputs")
))
data_folders <- data_folders[dir.exists(data_folders)]

output_dir <- file.path(script_dir, "Figure5_revised")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

find_csv <- function(pattern, friendly_name) {
  found <- unlist(lapply(data_folders, function(folder) {
    list.files(folder, pattern = pattern, full.names = TRUE, ignore.case = TRUE)
  }), use.names = FALSE)
  found <- unique(found)
  if (!length(found)) {
    stop(paste0(
      "\nMISSING INPUT FILE: ", friendly_name,
      "\nNo CSV filename matched: ", pattern,
      "\nSearched folders:\n  ", paste(data_folders, collapse = "\n  "),
      "\n\nFIX: Copy the four CSV files into the SAME folder as this R script,",
      "\nor download and unzip Figure5_R_ready_to_run.zip and run its script.",
      "\nAlternatively, set FIGURE5_DATA_DIR to your CSV folder."
    ), call. = FALSE)
  }
  if (length(found) > 1) {
    message("Multiple files match ", friendly_name, "; choosing: ", basename(found[1]))
  }
  normalizePath(found[1])
}

taxonomy_file <- find_csv("^final_clean_taxonomy.*\\.csv$", "taxonomy")
archaea_file <- find_csv("^Archaea_16s_phylum_summary_WITH_LOWER_RANKS_CORRECTED.*\\.csv$", "Archaea 16S summary")
bacteria_file <- find_csv("^Bacteria_16s_phylum_summary_WITH_LOWER_RANKS_CORRECTED.*\\.csv$", "Bacteria 16S summary")
metric_file <- find_csv("^aai_mash_ani_16s.*\\.csv$", "Mash/FastAAI/FastANI metrics")

message("\nInput files detected successfully:")
message("  Taxonomy: ", basename(taxonomy_file))
message("  Archaea:  ", basename(archaea_file))
message("  Bacteria: ", basename(bacteria_file))
message("  Metrics:  ", basename(metric_file), "\n")

png_dpi <- 400

# ============================================================
# 3. OFFICIALLY RECOGNIZED PHYLA
# ============================================================
official_phyla <- c(
  "Pseudomonadota", "Campylobacterota", "Myxococcota", "Bdellovibrionota",
  "Actinomycetota", "Bacillota", "Mycoplasmatota", "Cyanobacteriota",
  "Spirochaetota", "Chlorobiota", "Bacteroidota", "Planctomycetota",
  "Chlamydiota", "Deinococcota", "Chloroflexota", "Methanobacteriota",
  "Thermoproteota", "Abditibacteriota", "Acidobacteriota", "Aquificota",
  "Armatimonadota", "Atribacterota", "Balneolota", "Caldisericota",
  "Calditrichota", "Chrysiogenota", "Coprothermobacterota", "Deferribacterota",
  "Elusimicrobiota", "Fibrobacterota", "Fidelibacterota", "Fusobacteriota",
  "Gemmatimonadota", "Kiritimatiellota", "Lentisphaerota", "Minisyncoccota",
  "Nitrospinota", "Nitrospirota", "Rhodothermota", "Synergistota",
  "Thermodesulfobacteriota", "Thermodesulfobiota", "Thermomicrobiota",
  "Thermosulfidibacterota", "Thermotogota", "Verrucomicrobiota",
  "Vulcanimicrobiota", "Nanobdellota", "Nitrososphaerota", "Thermoplasmatota",
  "Bipolaricaulota", "Desulfobacterota", "Dictyoglomerota", "Dictyoglomota",
  "Halobacteriota", "Ignavibacteriota", "Poseidoniota", "Microcaldota",
  "Promethearchaeota"
)

# ============================================================
# 4. CARL WOESE PHYLA
# ============================================================
woese_phyla <- c(
  "Pseudomonadota", "Campylobacterota", "Myxococcota", "Bdellovibrionota",
  "Desulfobacterota", "Thermodesulfobacteriota", "Actinomycetota", "Bacillota",
  "Mycoplasmatota", "Cyanobacteriota", "Spirochaetota", "Chlorobiota",
  "Bacteroidota", "Planctomycetota", "Chlamydiota", "Deinococcota",
  "Chloroflexota", "Methanobacteriota", "Thermoproteota"
)
extra_official_phyla <- setdiff(official_phyla, woese_phyla)

# ============================================================
# 5. COLORS
# ============================================================
domain_colors <- c("Archaea" = "#C44E52", "Bacteria" = "#4C72B0")
kingdom_colors <- c(
  "Methanobacteriati" = "#D55E00",
  "Nanobdellati" = "#A66B00",
  "Promethearchaeati" = "#8B4D81",
  "Thermoproteati" = "#BC3046",
  "Unclassified Archaea" = "#9C6644",
  "Bacillati" = "#1B9E77",
  "Pseudomonadati" = "#7570B3",
  "Fusobacteriati" = "#008D91",
  "Thermotogati" = "#B82674",
  "Unclassified Bacteria" = "#567F19"
)
phylum_status_colors <- c(
  "Carl Woese phylum" = "#1976D2",
  "Additional official" = "#B87500",
  "Unofficial / other" = "#575B60"
)

# ============================================================
# 6. HELPER FUNCTIONS
# ============================================================
safe_first_numeric <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_real_)
  x[1]
}

make_percent_label <- function(x) {
  ifelse(is.na(x), "N/A", paste0(round(x, 1), "%"))
}

curve_fun <- function(x0, y0, x1, y1, id, val, col) {
  t <- seq(0, 1, length.out = 150)
  tibble(
    x = (1 - t)^3 * x0 +
      3 * (1 - t)^2 * t * (x0 + 0.65 * (x1 - x0)) +
      3 * (1 - t) * t^2 * (x1 - 0.65 * (x1 - x0)) +
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
# 7. READ TAXONOMY DATA
# ============================================================
taxonomy_df <- read_csv(taxonomy_file, show_col_types = FALSE)
names(taxonomy_df) <- names(taxonomy_df) |> str_trim() |> str_to_lower()

required_taxonomy_columns <- c("accession", "domain", "kingdom", "phylum", "species")
missing_taxonomy_columns <- setdiff(required_taxonomy_columns, names(taxonomy_df))
if (length(missing_taxonomy_columns) > 0) {
  stop(paste0("Missing columns in taxonomy file: ", paste(missing_taxonomy_columns, collapse = ", ")))
}

taxonomy_df <- taxonomy_df %>%
  mutate(
    accession = str_trim(as.character(accession)),
    domain = str_trim(as.character(domain)),
    kingdom = str_trim(as.character(kingdom)),
    phylum = str_trim(as.character(phylum)),
    species = str_trim(as.character(species)),
    kingdom = case_when(
      domain == "Bacteria" & (is.na(kingdom) | kingdom == "" | kingdom %in% c(
        "Bacteria incertae sedis", "Unclassified_Deltaproteobacteria", "Unclassified Deltaproteobacteria"
      )) ~ "Unclassified Bacteria",
      domain == "Archaea" & (is.na(kingdom) | kingdom == "" | !kingdom %in% c(
        "Methanobacteriati", "Nanobdellati", "Promethearchaeati", "Thermoproteati"
      )) ~ "Unclassified Archaea",
      TRUE ~ kingdom
    ),
    phylum_status = case_when(
      phylum %in% woese_phyla ~ "Carl Woese phylum",
      phylum %in% extra_official_phyla ~ "Additional official",
      TRUE ~ "Unofficial / other"
    ),
    phylum_status = factor(phylum_status, levels = c(
      "Carl Woese phylum", "Additional official", "Unofficial / other"
    ))
  ) %>%
  filter(
    !is.na(accession), accession != "",
    domain %in% c("Archaea", "Bacteria"),
    !is.na(phylum), phylum != ""
  ) %>%
  distinct(accession, .keep_all = TRUE)

overall_genome_count <- n_distinct(taxonomy_df$accession)
cat("\nTotal unique genomes:", comma(overall_genome_count), "\n")

# ============================================================
# 8. READ 16S SUMMARY DATA
# ============================================================
archaea_16s <- read_csv(archaea_file, show_col_types = FALSE)
bacteria_16s <- read_csv(bacteria_file, show_col_types = FALSE)
names(archaea_16s) <- names(archaea_16s) |> str_trim() |> str_to_lower()
names(bacteria_16s) <- names(bacteria_16s) |> str_trim() |> str_to_lower()
archaea_16s <- archaea_16s %>% mutate(domain = "Archaea")
bacteria_16s <- bacteria_16s %>% mutate(domain = "Bacteria")

sixteen_s_df <- bind_rows(archaea_16s, bacteria_16s) %>%
  mutate(
    domain = str_trim(as.character(domain)),
    phylum = str_trim(as.character(phylum)),
    total_genomes_16s = suppressWarnings(as.numeric(total_genomes)),
    missing_16s = suppressWarnings(as.numeric(missing_16s))
  ) %>%
  group_by(domain, phylum) %>%
  summarise(
    total_genomes_16s = sum(total_genomes_16s, na.rm = TRUE),
    missing_16s = sum(missing_16s, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(percent_missing = if_else(total_genomes_16s > 0, missing_16s / total_genomes_16s * 100, NA_real_))

# ============================================================
# 9. MASH / FASTAAI / FASTANI + PAIR ELIGIBILITY
# ============================================================
pair_df <- taxonomy_df %>%
  filter(!is.na(species), species != "") %>%
  count(domain, phylum, species, name = "n_species_genomes") %>%
  group_by(domain, phylum) %>%
  summarise(
    same_species_pairs = sum(n_species_genomes * (n_species_genomes - 1) / 2),
    .groups = "drop"
  )
cat("Total eligible same-species pairs: ", comma(sum(pair_df$same_species_pairs)), "\n", sep = "")

raw_metric <- read_csv(metric_file, show_col_types = FALSE)
names(raw_metric) <- str_trim(names(raw_metric))
required_metric_columns <- c("Phylum", "Total #", "Mash Failed, %", "Jaccard Failed,%", "ANI Failed, %")
missing_metric_columns <- setdiff(required_metric_columns, names(raw_metric))
if (length(missing_metric_columns) > 0) {
  stop("Missing metric columns: ", paste(missing_metric_columns, collapse = ", "))
}

metric_df <- raw_metric %>%
  transmute(
    phylum = str_trim(as.character(Phylum)),
    source_pairs = suppressWarnings(as.numeric(`Total #`)),
    mash_failed_percent = suppressWarnings(as.numeric(`Mash Failed, %`)),
    fastaai_failed_percent = suppressWarnings(as.numeric(`Jaccard Failed,%`)),
    fastani_failed_percent = suppressWarnings(as.numeric(`ANI Failed, %`))
  ) %>%
  filter(!is.na(source_pairs), !is.na(phylum), phylum != "") %>%
  mutate(phylum = if_else(phylum == "Unclassified Deltaproteobacteria", "Unclassified_Deltaproteobacteria", phylum)) %>%
  group_by(phylum) %>%
  summarise(
    source_pairs = safe_first_numeric(source_pairs),
    mash_failed_percent = safe_first_numeric(mash_failed_percent),
    fastaai_failed_percent = safe_first_numeric(fastaai_failed_percent),
    fastani_failed_percent = safe_first_numeric(fastani_failed_percent),
    .groups = "drop"
  )

pair_check <- pair_df %>%
  left_join(metric_df, by = "phylum") %>%
  filter(!is.na(source_pairs), source_pairs != same_species_pairs)
if (nrow(pair_check) > 0) {
  warning("Source denominators disagree with taxonomy counts for: ", paste(pair_check$phylum, collapse = ", "))
}

metric_check <- pair_df %>%
  filter(same_species_pairs > 0) %>%
  left_join(metric_df, by = "phylum") %>%
  filter(is.na(mash_failed_percent) | is.na(fastaai_failed_percent) | is.na(fastani_failed_percent))
if (nrow(metric_check) > 0) {
  stop(
    "Some phyla have eligible comparisons but missing metric values: ",
    paste(metric_check$phylum, collapse = ", "),
    "\nPlease resolve these in the source CSV."
  )
}

# ============================================================
# 10. KINGDOM ORDER
# ============================================================
kingdom_order_table <- tribble(
  ~domain, ~kingdom, ~kingdom_order,
  "Archaea", "Methanobacteriati", 1,
  "Archaea", "Nanobdellati", 2,
  "Archaea", "Promethearchaeati", 3,
  "Archaea", "Thermoproteati", 4,
  "Archaea", "Unclassified Archaea", 5,
  "Bacteria", "Bacillati", 1,
  "Bacteria", "Pseudomonadati", 2,
  "Bacteria", "Fusobacteriati", 3,
  "Bacteria", "Thermotogati", 4,
  "Bacteria", "Unclassified Bacteria", 5
)

# ============================================================
# 11. COMPLETE COMBINED PLOT FUNCTION
# ============================================================
make_domain_taxonomy_plot <- function(
    domain_to_plot,
    output_prefix,
    kingdoms_to_plot = NULL,
    show_domain = TRUE,
    page_width = 16.54,
    page_height = 11.69
) {
  stopifnot(domain_to_plot %in% c("Archaea", "Bacteria"))

  phyla_df <- taxonomy_df %>%
    filter(domain == domain_to_plot, is.null(kingdoms_to_plot) | kingdom %in% kingdoms_to_plot) %>%
    count(domain, kingdom, phylum, phylum_status, name = "total_genomes") %>%
    left_join(kingdom_order_table, by = c("domain", "kingdom")) %>%
    mutate(
      domain_order = case_when(domain == "Archaea" ~ 1, domain == "Bacteria" ~ 2, TRUE ~ 99),
      kingdom_order = replace_na(kingdom_order, 99)
    ) %>%
    left_join(sixteen_s_df, by = c("domain", "phylum")) %>%
    left_join(pair_df, by = c("domain", "phylum")) %>%
    mutate(same_species_pairs = replace_na(same_species_pairs, 0)) %>%
    left_join(metric_df, by = "phylum") %>%
    mutate(
      mash_failed_percent = if_else(same_species_pairs == 0, NA_real_, mash_failed_percent),
      fastaai_failed_percent = if_else(same_species_pairs == 0, NA_real_, fastaai_failed_percent),
      fastani_failed_percent = if_else(same_species_pairs == 0, NA_real_, fastani_failed_percent)
    ) %>%
    mutate(
      domain_color = unname(domain_colors[domain]),
      kingdom_color = unname(kingdom_colors[kingdom]),
      kingdom_color = if_else(is.na(kingdom_color), "#6B7280", kingdom_color),
      status_color = unname(phylum_status_colors[as.character(phylum_status)])
    ) %>%
    arrange(domain_order, kingdom_order, desc(total_genomes), phylum)

  if (nrow(phyla_df) == 0) stop("No taxonomy rows are available for plotting.")

  phyla_df <- phyla_df %>%
    mutate(
      row_number_global = row_number(),
      domain_change = case_when(row_number_global == 1 ~ FALSE, TRUE ~ domain != lag(domain)),
      kingdom_change = case_when(
        row_number_global == 1 ~ FALSE,
        TRUE ~ domain != lag(domain) | kingdom != lag(kingdom)
      ),
      domain_gap_increment = if_else(domain_change, 2.8, 0),
      kingdom_gap_increment = if_else(kingdom_change & !domain_change, 0.75, 0),
      accumulated_gap = cumsum(domain_gap_increment + kingdom_gap_increment),
      y_ascending = row_number_global * 1.38 + accumulated_gap,
      y = max(y_ascending) - y_ascending + 1.38,
      phylum_label = paste0(phylum, " (", comma(total_genomes), ")"),
      missing_label = make_percent_label(percent_missing),
      mash_label = make_percent_label(mash_failed_percent),
      fastaai_label = make_percent_label(fastaai_failed_percent),
      fastani_label = make_percent_label(fastani_failed_percent)
    )

  kingdom_df <- phyla_df %>%
    group_by(domain, domain_order, kingdom, kingdom_order, kingdom_color) %>%
    summarise(
      kingdom_genome_total = sum(total_genomes, na.rm = TRUE),
      min_phylum_y = min(y),
      max_phylum_y = max(y),
      .groups = "drop"
    ) %>%
    mutate(
      y = (min_phylum_y + max_phylum_y) / 2,
      total_genomes = kingdom_genome_total,
      kingdom_label = case_when(
        kingdom == "Unclassified Bacteria" ~ paste0("Unclassified\nBacteria\n", comma(total_genomes)),
        TRUE ~ paste0(kingdom, "\n", comma(total_genomes))
      )
    ) %>%
    select(domain, domain_order, kingdom, kingdom_order, kingdom_color, total_genomes, kingdom_label, y) %>%
    arrange(domain_order, kingdom_order, desc(total_genomes))

  domain_df <- kingdom_df %>%
    group_by(domain, domain_order) %>%
    summarise(domain_genome_total = sum(total_genomes, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      y = (max(phyla_df$y) + min(phyla_df$y)) / 2,
      total_genomes = domain_genome_total,
      domain_color = unname(domain_colors[domain])
    ) %>%
    select(domain, domain_order, domain_color, total_genomes, y) %>%
    arrange(domain_order)

  root_total <- sum(domain_df$total_genomes, na.rm = TRUE)

  x_domain <- -0.25
  x_kingdom <- 3.00
  x_status <- 6.30
  x_phylum_label <- 7.35
  metric_shift <- if (show_domain) 2.50 else 6.75
  x_16s_bar_start <- 21.55 - metric_shift
  x_16s_bar_end <- 22.95 - metric_shift
  x_16s_text <- 23.15 - metric_shift
  x_mash_bar_start <- 25.70 - metric_shift
  x_mash_bar_end <- 27.10 - metric_shift
  x_mash_text <- 27.30 - metric_shift
  x_fastaai_bar_start <- 29.85 - metric_shift
  x_fastaai_bar_end <- 31.25 - metric_shift
  x_fastaai_text <- 31.45 - metric_shift
  x_fastani_bar_start <- 34.00 - metric_shift
  x_fastani_bar_end <- 35.40 - metric_shift
  x_fastani_text <- 35.60 - metric_shift
  percent_bar_width <- x_16s_bar_end - x_16s_bar_start

  phyla_df <- phyla_df %>%
    mutate(
      missing_fraction = pmin(pmax(percent_missing / 100, 0), 1),
      mash_fraction = pmin(pmax(mash_failed_percent / 100, 0), 1),
      fastaai_fraction = pmin(pmax(fastaai_failed_percent / 100, 0), 1),
      fastani_fraction = pmin(pmax(fastani_failed_percent / 100, 0), 1),
      x_16s_end = x_16s_bar_start + missing_fraction * percent_bar_width,
      x_mash_end = x_mash_bar_start + mash_fraction * percent_bar_width,
      x_fastaai_end = x_fastaai_bar_start + fastaai_fraction * percent_bar_width,
      x_fastani_end = x_fastani_bar_start + fastani_fraction * percent_bar_width
    )

  curves <- list()
  if (show_domain) for (i in seq_len(nrow(kingdom_df))) {
    current_domain_y <- domain_df$y[domain_df$domain == kingdom_df$domain[i]][1]
    curves[[length(curves) + 1]] <- curve_fun(
      x0 = x_domain, y0 = current_domain_y,
      x1 = x_kingdom, y1 = kingdom_df$y[i],
      id = paste0("domain_kingdom_", kingdom_df$domain[i], "_", kingdom_df$kingdom[i]),
      val = kingdom_df$total_genomes[i],
      col = kingdom_df$kingdom_color[i]
    )
  }

  for (i in seq_len(nrow(phyla_df))) {
    current_kingdom_y <- kingdom_df$y[
      kingdom_df$domain == phyla_df$domain[i] & kingdom_df$kingdom == phyla_df$kingdom[i]
    ][1]
    curves[[length(curves) + 1]] <- curve_fun(
      x0 = x_kingdom, y0 = current_kingdom_y,
      x1 = x_status, y1 = phyla_df$y[i],
      id = paste0("kingdom_phylum_", i),
      val = phyla_df$total_genomes[i],
      col = phyla_df$kingdom_color[i]
    )
  }
  curve_df <- bind_rows(curves)

  number_of_phyla <- nrow(phyla_df)
  phylum_text_size <- 4.45
  metric_text_size <- 4.12
  marker_size <- 4.05
  phyla_df <- phyla_df %>%
    mutate(phylum_size = if_else(nchar(phylum_label) > 52, 3.75, phylum_text_size))

  top_y <- max(phyla_df$y)
  bottom_y <- min(phyla_df$y)
  header_y <- top_y + 2.55
  legend_title_y <- top_y + 6.05
  legend_y_1 <- top_y + 4.80
  legend_y_2 <- legend_y_1
  legend_y_3 <- legend_y_1
  legend_point_x <- c(7.35, 16.50, 24.00)
  legend_text_x <- legend_point_x + 0.48
  legend_title_x <- x_phylum_label

  p <- ggplot() +
    geom_path(
      data = curve_df,
      aes(x = x, y = y, group = id, color = col, linewidth = val),
      alpha = 0.83, lineend = "round"
    ) +
    geom_label(
      data = if (show_domain) domain_df else domain_df[0, ],
      aes(x = x_domain, y = y, label = paste0(domain, "\n", comma(total_genomes)), fill = domain_color),
      color = "white", size = 4.60, fontface = "bold", label.size = 0,
      label.padding = unit(0.38, "lines"), lineheight = 0.90
    ) +
    geom_label(
      data = kingdom_df,
      aes(x = x_kingdom, y = y, label = kingdom_label, fill = kingdom_color),
      color = "white", size = 4.40, fontface = "bold", label.size = 0,
      label.padding = unit(0.30, "lines"), lineheight = 0.90
    ) +
    geom_point(
      data = phyla_df,
      aes(x = x_status, y = y, color = status_color),
      shape = 15, size = marker_size
    ) +
    geom_text(
      data = phyla_df,
      aes(x = x_phylum_label, y = y, label = phylum_label, size = phylum_size),
      color = "#17212B", hjust = 0, vjust = 0.5,
      fontface = "bold", family = "sans", lineheight = 1.0
    ) +
    geom_segment(
      data = phyla_df,
      aes(x = x_16s_bar_start - 0.45, xend = x_fastani_text + 1.35, y = y, yend = y),
      linewidth = 0.23, color = "#E1E5E9"
    ) +
    geom_segment(
      data = phyla_df %>% filter(!is.na(percent_missing)),
      aes(x = x_16s_bar_start, xend = x_16s_bar_end, y = y, yend = y),
      color = "#E4E8EC", linewidth = 3.2, lineend = "round"
    ) +
    geom_segment(
      data = phyla_df %>% filter(!is.na(percent_missing), percent_missing > 0),
      aes(x = x_16s_bar_start, xend = x_16s_end, y = y, yend = y, color = kingdom_color),
      linewidth = 4.1, lineend = "round"
    ) +
    geom_text(
      data = phyla_df,
      aes(x = x_16s_text, y = y, label = missing_label, color = kingdom_color),
      hjust = 0, size = metric_text_size, fontface = "bold"
    ) +
    geom_segment(
      data = phyla_df %>% filter(same_species_pairs > 0),
      aes(x = x_mash_bar_start, xend = x_mash_bar_end, y = y, yend = y),
      color = "#E4E8EC", linewidth = 3.2, lineend = "round"
    ) +
    geom_segment(
      data = phyla_df %>% filter(!is.na(mash_failed_percent), mash_failed_percent > 0),
      aes(x = x_mash_bar_start, xend = x_mash_end, y = y, yend = y, color = kingdom_color),
      linewidth = 4.1, lineend = "round"
    ) +
    geom_text(
      data = phyla_df,
      aes(x = x_mash_text, y = y, label = mash_label),
      color = "#202A36", hjust = 0, size = metric_text_size, fontface = "bold"
    ) +
    geom_segment(
      data = phyla_df %>% filter(same_species_pairs > 0),
      aes(x = x_fastaai_bar_start, xend = x_fastaai_bar_end, y = y, yend = y),
      color = "#E4E8EC", linewidth = 3.2, lineend = "round"
    ) +
    geom_segment(
      data = phyla_df %>% filter(!is.na(fastaai_failed_percent), fastaai_failed_percent > 0),
      aes(x = x_fastaai_bar_start, xend = x_fastaai_end, y = y, yend = y, color = kingdom_color),
      linewidth = 4.1, lineend = "round"
    ) +
    geom_text(
      data = phyla_df,
      aes(x = x_fastaai_text, y = y, label = fastaai_label),
      color = "#202A36", hjust = 0, size = metric_text_size, fontface = "bold"
    ) +
    geom_segment(
      data = phyla_df %>% filter(same_species_pairs > 0),
      aes(x = x_fastani_bar_start, xend = x_fastani_bar_end, y = y, yend = y),
      color = "#E4E8EC", linewidth = 3.2, lineend = "round"
    ) +
    geom_segment(
      data = phyla_df %>% filter(!is.na(fastani_failed_percent), fastani_failed_percent > 0),
      aes(x = x_fastani_bar_start, xend = x_fastani_end, y = y, yend = y, color = kingdom_color),
      linewidth = 4.1, lineend = "round"
    ) +
    geom_text(
      data = phyla_df,
      aes(x = x_fastani_text, y = y, label = fastani_label),
      color = "#202A36", hjust = 0, size = metric_text_size, fontface = "bold"
    ) +
    (if (show_domain) annotate(
      "text", x = x_domain, y = header_y, label = "Domain", size = 4.45, fontface = "bold"
    ) else NULL) +
    annotate("text", x = x_kingdom, y = header_y, label = "Kingdom", size = 4.45, fontface = "bold") +
    annotate("text", x = x_status, y = header_y, label = "Phylum\nstatus", size = 3.65, fontface = "bold", lineheight = 0.90) +
    annotate("text", x = x_phylum_label, y = header_y, label = "Phylum (genomes)", hjust = 0, size = 4.45, fontface = "bold") +
    annotate("text", x = mean(c(x_16s_bar_start, x_16s_bar_end)), y = header_y, label = "Missing 16S", hjust = 0.5, size = 4.30, fontface = "bold") +
    annotate("text", x = mean(c(x_mash_bar_start, x_mash_bar_end)), y = header_y, label = "Mash", hjust = 0.5, size = 4.30, fontface = "bold") +
    annotate("text", x = mean(c(x_fastaai_bar_start, x_fastaai_bar_end)), y = header_y, label = "FastAAI", hjust = 0.5, size = 4.30, fontface = "bold") +
    annotate("text", x = mean(c(x_fastani_bar_start, x_fastani_bar_end)), y = header_y, label = "FastANI", hjust = 0.5, size = 4.30, fontface = "bold") +
    annotate("text", x = legend_title_x, y = legend_title_y, label = "Phylum classification", hjust = 0, size = 4.05, fontface = "bold", color = "#111827") +
    annotate("point", x = legend_point_x[1], y = legend_y_1, shape = 15, size = 4.3, color = phylum_status_colors["Carl Woese phylum"]) +
    annotate("text", x = legend_text_x[1], y = legend_y_1, label = "Carl Woese phylum", hjust = 0, vjust = 0.5, size = 3.62, fontface = "bold", color = "#111827") +
    annotate("point", x = legend_point_x[2], y = legend_y_2, shape = 15, size = 4.3, color = phylum_status_colors["Additional official"]) +
    annotate("text", x = legend_text_x[2], y = legend_y_2, label = "Additional official phylum", hjust = 0, vjust = 0.5, size = 3.62, fontface = "bold", color = "#111827") +
    annotate("point", x = legend_point_x[3], y = legend_y_3, shape = 15, size = 4.3, color = phylum_status_colors["Unofficial / other"]) +
    annotate("text", x = legend_text_x[3], y = legend_y_3, label = "Unofficial / other phylum", hjust = 0, vjust = 0.5, size = 3.62, fontface = "bold", color = "#111827") +
    scale_color_identity() +
    scale_fill_identity() +
    scale_size_identity(guide = "none") +
    scale_linewidth(range = c(0.65, 9.0), trans = "sqrt", guide = "none") +
    coord_cartesian(
      xlim = if (show_domain) c(-2.15, 35.25) else c(0.55, 31.00),
      ylim = c(bottom_y - 3.0, legend_title_y + 1.15),
      clip = "off"
    ) +
    theme_void(base_size = 13) +
    theme(
      legend.position = "none",
      plot.margin = margin(22, 25, 40, 22),
      text = element_text(family = "sans", face = "bold")
    )

  p <- p + annotate(
    "text",
    x = x_status,
    y = bottom_y - 2.20,
    label = "Numbers in parentheses = genomes per phylum; N/A (similarity methods) = no eligible same-species pairs",
    hjust = 0,
    size = 3.35,
    color = "#354252",
    fontface = "plain"
  )

  print(p)
  out_width <- page_width
  out_height <- page_height

  ggsave(filename = paste0(output_prefix, ".png"), plot = p, width = out_width, height = out_height, dpi = png_dpi, bg = "white", limitsize = FALSE)
  ggsave(filename = paste0(output_prefix, ".pdf"), plot = p, width = out_width, height = out_height, bg = "white", limitsize = FALSE)
  ggsave(filename = paste0(output_prefix, ".svg"), plot = p, width = out_width, height = out_height, bg = "white", limitsize = FALSE)

  output_data <- phyla_df %>%
    select(
      domain, kingdom, phylum, phylum_status, total_genomes, same_species_pairs,
      total_genomes_16s, missing_16s, percent_missing, mash_failed_percent,
      fastaai_failed_percent, fastani_failed_percent, domain_color,
      kingdom_color, status_color
    )
  write_csv(output_data, paste0(output_prefix, "_plot_data.csv"))

  classification_summary <- output_data %>%
    group_by(domain, phylum_status) %>%
    summarise(
      number_of_phyla = n_distinct(phylum),
      number_of_genomes = sum(total_genomes, na.rm = TRUE),
      .groups = "drop"
    )
  write_csv(classification_summary, paste0(output_prefix, "_classification_summary.csv"))

  cat("\n============================================================\n")
  cat("DOMAIN FIGURE COMPLETED\n")
  cat("Genomes in ", domain_to_plot, ":", comma(root_total), "\n")
  cat("Total phyla:", number_of_phyla, "\n")
  cat("\nFiles created:\n")
  cat(paste0(output_prefix, ".png\n"))
  cat(paste0(output_prefix, ".pdf\n"))
  cat(paste0(output_prefix, ".svg\n"))
  cat(paste0(output_prefix, "_plot_data.csv\n"))
  cat(paste0(output_prefix, "_classification_summary.csv\n"))
  cat("============================================================\n\n")

  invisible(list(
    plot = p,
    phylum_data = output_data,
    kingdom_data = kingdom_df,
    domain_data = domain_df,
    classification_summary = classification_summary
  ))
}

# ============================================================
# 12. PRECISE KINGDOM SPLIT (ALL BACTERIAL PHYLA EXACTLY ONCE)
# ============================================================
bacteria_5B <- c("Bacillati", "Fusobacteriati", "Thermotogati", "Unclassified Bacteria")
bacteria_5C <- c("Pseudomonadati")

observed_kingdoms <- taxonomy_df %>%
  filter(domain == "Bacteria") %>%
  distinct(kingdom) %>%
  pull(kingdom)

if (!setequal(c(bacteria_5B, bacteria_5C), observed_kingdoms) ||
    length(intersect(bacteria_5B, bacteria_5C)) > 0) {
  stop(
    "Figure 5B/5C split does not cover the bacterial kingdoms exactly once. ",
    "Check taxonomy normalization and split vectors."
  )
}

split_check <- taxonomy_df %>%
  filter(domain == "Bacteria") %>%
  count(kingdom, phylum, name = "genomes") %>%
  mutate(
    figure = case_when(
      kingdom %in% bacteria_5B ~ "5B",
      kingdom %in% bacteria_5C ~ "5C",
      TRUE ~ NA_character_
    )
  )

stopifnot(
  !anyNA(split_check$figure),
  nrow(split_check) == nrow(distinct(split_check, phylum)),
  sum(split_check$genomes) == sum(taxonomy_df$domain == "Bacteria")
)

cat("\nFigure 5B phyla: ", sum(split_check$figure == "5B"),
    "; genomes: ", comma(sum(split_check$genomes[split_check$figure == "5B"])), "\n", sep = "")
cat("Figure 5C phyla: ", sum(split_check$figure == "5C"),
    "; genomes: ", comma(sum(split_check$genomes[split_check$figure == "5C"])), "\n", sep = "")

# ============================================================
# 13. GENERATE FIGURES
# ============================================================
archaea_result <- make_domain_taxonomy_plot(
  domain_to_plot = "Archaea",
  output_prefix = file.path(output_dir, "Figure_5A_Archaea"),
  show_domain = TRUE,
  page_width = 17.80,
  page_height = 12.60
)

bacteria_5B_result <- make_domain_taxonomy_plot(
  domain_to_plot = "Bacteria",
  output_prefix = file.path(output_dir, "Figure_5B_Bacteria"),
  kingdoms_to_plot = bacteria_5B,
  show_domain = FALSE,
  page_width = 16.54,
  page_height = 11.69
)

bacteria_5C_result <- make_domain_taxonomy_plot(
  domain_to_plot = "Bacteria",
  output_prefix = file.path(output_dir, "Figure_5C_Bacteria"),
  kingdoms_to_plot = bacteria_5C,
  show_domain = FALSE,
  page_width = 16.54,
  page_height = 15.20
)

cat("\nCompleted: Archaea, Bacteria 5B and Bacteria 5C.\n")
cat("Output files are in: ", normalizePath(output_dir), "\n", sep = "")
cat("Use PDF or SVG for sharp vector figures.\n")
