# ============================================================
# TAXONOMY DICTIONARY FIGURES (NO 16S, LARGER FONTS, SHIFTED LAYOUT)
# Recreates:
# FIGURE 1: Total -> Domain -> Kingdom
# FIGURE 2: Archaea -> Kingdom -> Phylum + lower ranks (extra shift for Archaea, reference removed)
# FIGURE 3: Bacteria -> Kingdom -> Phylum + lower ranks
# ============================================================

rm(list = ls())

library(tidyverse)
library(stringr)
library(ggplot2)
library(scales)
library(svglite)

# ============================================================
# INPUT FILES
# ============================================================

taxonomy_file <- "final_clean_taxonomy.csv"

# ============================================================
# OUTPUT DIRECTORY
# ============================================================

outdir <- "taxonomy_dictionary_figures_RECREATED"

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

# ============================================================
# COLORS
# ============================================================

domain_colors <- c(
  "Archaea"  = "#C44E52",
  "Bacteria" = "#4C72B0"
)

kingdom_colors <- c(
  "Methanobacteriati"     = "#D55E00",
  "Nanobdellati"          = "#E69F00",
  "Promethearchaeati"     = "#CC79A7",
  "Thermoproteati"        = "#F94144",
  "Unclassified Archaea"  = "#9C6644",
  "Bacillati"             = "#1B9E77",
  "Pseudomonadati"        = "#7570B3",
  "Fusobacteriati"        = "#00BFC4",
  "Thermotogati"          = "#E7298A",
  "Unclassified Bacteria" = "#66A61E"
)

# ============================================================
# HELPER FUNCTIONS
# ============================================================

clean_names_lower <- function(df) {
  names(df) <- names(df) |>
    str_trim() |>
    str_to_lower() |>
    str_replace_all("\\s+", "_")
  return(df)
}

curve_fun <- function(x0, y0, x1, y1, id, val, col) {
  t <- seq(0, 1, length.out = 180)
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

kingdom_order_fun <- function(x) {
  case_when(
    x == "Methanobacteriati" ~ 1,
    x == "Nanobdellati" ~ 2,
    x == "Promethearchaeati" ~ 3,
    x == "Thermoproteati" ~ 4,
    x == "Unclassified Archaea" ~ 5,
    x == "Bacillati" ~ 1,
    x == "Pseudomonadati" ~ 2,
    x == "Fusobacteriati" ~ 3,
    x == "Thermotogati" ~ 4,
    x == "Unclassified Bacteria" ~ 5,
    TRUE ~ 99
  )
}

# ============================================================
# READ TAXONOMY
# ============================================================

taxonomy_df <- read_csv(taxonomy_file, show_col_types = FALSE)
taxonomy_df <- clean_names_lower(taxonomy_df)

required_cols <- c(
  "accession", "domain", "kingdom", "phylum",
  "class", "order", "family", "genus", "species"
)

missing_cols <- setdiff(required_cols, names(taxonomy_df))
if (length(missing_cols) > 0) {
  stop("These required columns are missing: ", paste(missing_cols, collapse = ", "))
}

taxonomy_df <- taxonomy_df %>%
  mutate(
    accession = str_trim(as.character(accession)),
    domain    = str_trim(as.character(domain)),
    kingdom   = str_trim(as.character(kingdom)),
    phylum    = str_trim(as.character(phylum)),
    class     = str_trim(as.character(class)),
    order     = str_trim(as.character(order)),
    family    = str_trim(as.character(family)),
    genus     = str_trim(as.character(genus)),
    species   = str_trim(as.character(species)),
    kingdom = case_when(
      domain == "Bacteria" &
        kingdom %in% c("Bacteria incertae sedis",
                       "Unclassified_Deltaproteobacteria") ~ "Unclassified Bacteria",
      domain == "Archaea" &
        !kingdom %in% c("Methanobacteriati",
                        "Nanobdellati",
                        "Promethearchaeati",
                        "Thermoproteati") ~ "Unclassified Archaea",
      TRUE ~ kingdom
    )
  ) %>%
  filter(
    accession != "",
    domain != "",
    kingdom != "",
    phylum != ""
  ) %>%
  distinct(accession, .keep_all = TRUE)

cat("\nTotal genomes after cleaning:", nrow(taxonomy_df), "\n")

# ============================================================
# FIGURE 1: TOTAL -> DOMAIN -> KINGDOM (increased fonts)
# ============================================================

make_figure1_total_domain_kingdom <- function() {

  total_count <- n_distinct(taxonomy_df$accession)

  domain_df <- taxonomy_df %>%
    group_by(domain) %>%
    summarise(total_genomes = n_distinct(accession), .groups = "drop") %>%
    mutate(
      color = domain_colors[domain],
      y = case_when(
        domain == "Archaea"  ~ 6.7,
        domain == "Bacteria" ~ 3.2,
        TRUE ~ 0
      )
    )

  kingdom_df <- taxonomy_df %>%
    group_by(domain, kingdom) %>%
    summarise(total_genomes = n_distinct(accession), .groups = "drop") %>%
    mutate(
      color = kingdom_colors[kingdom],
      y = case_when(
        kingdom == "Methanobacteriati"     ~ 8.4,
        kingdom == "Nanobdellati"          ~ 7.0,
        kingdom == "Promethearchaeati"     ~ 5.5,
        kingdom == "Thermoproteati"        ~ 4.1,
        kingdom == "Bacillati"             ~ 2.8,
        kingdom == "Pseudomonadati"        ~ 1.1,
        kingdom == "Fusobacteriati"        ~ -0.4,
        kingdom == "Thermotogati"          ~ -1.9,
        kingdom == "Unclassified Bacteria" ~ -3.4,
        TRUE ~ 0
      )
    )

  total_y <- 4.9
  x_total   <- 0.0
  x_domain  <- 4.8
  x_kingdom <- 9.1

  curves <- list()
  for (i in seq_len(nrow(domain_df))) {
    curves[[length(curves) + 1]] <- curve_fun(
      x_total, total_y,
      x_domain, domain_df$y[i],
      paste0("total_domain_", i),
      domain_df$total_genomes[i],
      domain_df$color[i]
    )
  }
  for (i in seq_len(nrow(kingdom_df))) {
    dy <- domain_df$y[domain_df$domain == kingdom_df$domain[i]]
    curves[[length(curves) + 1]] <- curve_fun(
      x_domain, dy,
      x_kingdom, kingdom_df$y[i],
      paste0("domain_kingdom_", i),
      kingdom_df$total_genomes[i],
      kingdom_df$color[i]
    )
  }
  curve_df <- bind_rows(curves)

  p <- ggplot() +
    geom_path(
      data = curve_df,
      aes(x = x, y = y, group = id, color = col, linewidth = val),
      alpha = 0.32,
      lineend = "round"
    ) +
    geom_label(
      aes(
        x = x_total,
        y = total_y,
        label = paste0("Total type strains\n", comma(total_count))
      ),
      fill = "#4B5563",
      color = "white",
      size = 10.0,
      fontface = "bold",
      label.size = 0,
      label.padding = unit(0.48, "lines"),
      lineheight = 0.9
    ) +
    geom_label(
      data = domain_df,
      aes(
        x = x_domain,
        y = y,
        label = paste0(domain, "\n", comma(total_genomes)),
        fill = domain
      ),
      color = "white",
      size = 8.8,
      fontface = "bold",
      label.size = 0,
      label.padding = unit(0.42, "lines"),
      lineheight = 0.9
    ) +
    geom_label(
      data = kingdom_df,
      aes(
        x = x_kingdom,
        y = y,
        label = paste0(kingdom, "\n", comma(total_genomes))
      ),
      fill = "#8A8A8A",
      color = "white",
      size = 7.5,
      fontface = "bold",
      label.size = 0,
      label.padding = unit(0.36, "lines"),
      lineheight = 0.9
    ) +
    annotate(
      "text",
      x = 0.0,
      y = 10.2,
      label = "Type Strain Genomes –> Domains –> Kingdoms",
      hjust = 0,
      size = 10.0,
      fontface = "bold",
      color = "#111827"
    ) +
    scale_color_identity() +
    scale_fill_manual(values = c(domain_colors)) +
    scale_linewidth(range = c(0.2, 22), guide = "none") +
    coord_cartesian(
      xlim = c(-1.8, 10.8),
      ylim = c(-4.8, 11.0),
      clip = "off"
    ) +
    theme_void(base_size = 18) +
    theme(
      legend.position = "none",
      plot.margin = margin(50, 100, 50, 100),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      text = element_text(family = "Helvetica", face = "bold")
    )

  ggsave(file.path(outdir, "FIGURE_1_Total_Domain_Kingdom_RECREATED.png"),
         p, width = 16, height = 10, dpi = 500, bg = "white")
  ggsave(file.path(outdir, "FIGURE_1_Total_Domain_Kingdom_RECREATED.pdf"),
         p, width = 16, height = 10, bg = "white")
  ggsave(file.path(outdir, "FIGURE_1_Total_Domain_Kingdom_RECREATED.svg"),
         p, width = 16, height = 10, bg = "white")
}

# ============================================================
# FIGURE 2 / 3: DOMAIN -> KINGDOM -> PHYLUM + LOWER RANKS
# with extra shift for Archaea and removal of reference
# ============================================================

make_domain_phylum_plot <- function(selected_domain, output_prefix) {

  domain_tax <- taxonomy_df %>%
    filter(domain == selected_domain)

  if (nrow(domain_tax) == 0) {
    stop("No rows found for domain: ", selected_domain)
  }

  phyla_df <- domain_tax %>%
    group_by(domain, kingdom, phylum) %>%
    summarise(
      total_genomes = n_distinct(accession),
      class_count   = n_distinct(class[class != "" & !is.na(class)]),
      order_count   = n_distinct(order[order != "" & !is.na(order)]),
      family_count  = n_distinct(family[family != "" & !is.na(family)]),
      genus_count   = n_distinct(genus[genus != "" & !is.na(genus)]),
      species_count = n_distinct(species[species != "" & !is.na(species)]),
      .groups = "drop"
    ) %>%
    mutate(
      kingdom_order = kingdom_order_fun(kingdom),
      color = kingdom_colors[kingdom],
      color = ifelse(is.na(color), "#6B7280", color)
    ) %>%
    arrange(kingdom_order, desc(total_genomes), phylum) %>%
    mutate(
      y = rev(row_number()) * 1.22,
      genome_label = comma(total_genomes),
      # For the specific Archaea phylum, remove the parenthetical reference
      phylum_display = if_else(
        selected_domain == "Archaea" &
        str_detect(phylum, "Unclassified_Candidatus_Asgardarchaeia"),
        str_remove(phylum, " \\(ex Tamarit et al\\. 2024\\)"),
        phylum
      )
    )

  kingdom_df <- phyla_df %>%
    group_by(domain, kingdom, kingdom_order, color) %>%
    summarise(
      total_genomes = sum(total_genomes, na.rm = TRUE),
      y = mean(y),
      .groups = "drop"
    ) %>%
    arrange(kingdom_order)

  domain_count <- sum(kingdom_df$total_genomes)
  domain_y <- mean(kingdom_df$y)

  # Base x positions (same as before for Bacteria)
  x_domain  <- 0.0
  x_kingdom <- 3.1
  x_phylum  <- 6.0

  # Lower ranks base positions
  x_class   <- 10.05
  x_order   <- 10.75
  x_family  <- 11.45
  x_genus   <- 12.15
  x_species <- 12.85

  # Genome bar base positions
  x_genome_bar_start <- 14.30
  x_genome_bar_end   <- 18.30
  x_genome_text      <- 18.70

  # If Archaea, shift all lower-rank and genome bar elements further right
  if (selected_domain == "Archaea") {
    extra_shift <- 1.2   # adjust as needed
    x_class   <- x_class + extra_shift
    x_order   <- x_order + extra_shift
    x_family  <- x_family + extra_shift
    x_genus   <- x_genus + extra_shift
    x_species <- x_species + extra_shift
    x_genome_bar_start <- x_genome_bar_start + extra_shift
    x_genome_bar_end   <- x_genome_bar_end + extra_shift
    x_genome_text      <- x_genome_text + extra_shift
  }

  max_genome <- max(phyla_df$total_genomes, na.rm = TRUE)

  phyla_df <- phyla_df %>%
    mutate(
      genome_fraction = total_genomes / max_genome,
      genome_bar_end = x_genome_bar_start +
        genome_fraction * (x_genome_bar_end - x_genome_bar_start)
    )

  curves <- list()
  for (i in seq_len(nrow(kingdom_df))) {
    curves[[length(curves) + 1]] <- curve_fun(
      x_domain, domain_y,
      x_kingdom, kingdom_df$y[i],
      paste0("domain_kingdom_", i),
      kingdom_df$total_genomes[i],
      kingdom_df$color[i]
    )
  }
  for (i in seq_len(nrow(phyla_df))) {
    ky <- kingdom_df$y[kingdom_df$kingdom == phyla_df$kingdom[i]]
    curves[[length(curves) + 1]] <- curve_fun(
      x_kingdom, ky,
      x_phylum, phyla_df$y[i],
      paste0("kingdom_phylum_", i),
      phyla_df$total_genomes[i],
      phyla_df$color[i]
    )
  }
  curve_df <- bind_rows(curves)

  header_y <- max(phyla_df$y) + 1.65
  title_y  <- max(phyla_df$y) + 3.95
  axis_y   <- min(phyla_df$y) - 1.25
  axis_label_y <- min(phyla_df$y) - 1.95

  title_text <- paste0(
    selected_domain,
    " Type Strain Genomes Across Kingdoms, Phyla, and Lower Ranks"
  )

  p <- ggplot() +
    geom_path(
      data = curve_df,
      aes(x = x, y = y, group = id, color = col, linewidth = val),
      alpha = 0.28,
      lineend = "round"
    ) +
    geom_label(
      aes(
        x = x_domain,
        y = domain_y,
        label = paste0(selected_domain, "\n", comma(domain_count)),
        fill = selected_domain
      ),
      color = "white",
      size = 7.5,
      fontface = "bold",
      label.size = 0,
      label.padding = unit(0.35, "lines"),
      lineheight = 0.9
    ) +
    geom_label(
      data = kingdom_df,
      aes(
        x = x_kingdom,
        y = y,
        label = paste0(kingdom, "\n", comma(total_genomes)),
        fill = color
      ),
      color = "white",
      size = 5.8,
      fontface = "bold",
      label.size = 0,
      label.padding = unit(0.28, "lines"),
      lineheight = 0.9
    ) +
    geom_point(
      data = phyla_df,
      aes(x = x_phylum, y = y, color = color),
      size = 2.4
    ) +
    geom_text(
      data = phyla_df,
      aes(x = x_phylum + 0.55, y = y, label = phylum_display, color = color),
      hjust = 0,
      vjust = 0.5,
      size = 4.5,
      fontface = "bold",
      family = "Helvetica",
      lineheight = 0.9
    ) +
    # Lower-rank table background
    geom_tile(
      data = phyla_df,
      aes(x = x_class, y = y),
      width = 0.52, height = 0.38, fill = "#F3F4F6"
    ) +
    geom_tile(
      data = phyla_df,
      aes(x = x_order, y = y),
      width = 0.52, height = 0.38, fill = "#F3F4F6"
    ) +
    geom_tile(
      data = phyla_df,
      aes(x = x_family, y = y),
      width = 0.52, height = 0.38, fill = "#F3F4F6"
    ) +
    geom_tile(
      data = phyla_df,
      aes(x = x_genus, y = y),
      width = 0.52, height = 0.38, fill = "#F3F4F6"
    ) +
    geom_tile(
      data = phyla_df,
      aes(x = x_species, y = y),
      width = 0.52, height = 0.38, fill = "#F3F4F6"
    ) +
    geom_text(
      data = phyla_df,
      aes(x = x_class, y = y, label = comma(class_count), color = color),
      size = 4.0,
      fontface = "bold"
    ) +
    geom_text(
      data = phyla_df,
      aes(x = x_order, y = y, label = comma(order_count), color = color),
      size = 4.0,
      fontface = "bold"
    ) +
    geom_text(
      data = phyla_df,
      aes(x = x_family, y = y, label = comma(family_count), color = color),
      size = 4.0,
      fontface = "bold"
    ) +
    geom_text(
      data = phyla_df,
      aes(x = x_genus, y = y, label = comma(genus_count), color = color),
      size = 4.0,
      fontface = "bold"
    ) +
    geom_text(
      data = phyla_df,
      aes(x = x_species, y = y, label = comma(species_count), color = color),
      size = 4.0,
      fontface = "bold"
    ) +
    # Genome count bar
    geom_segment(
      data = phyla_df,
      aes(
        x = x_genome_bar_start,
        xend = x_genome_bar_end,
        y = y, yend = y
      ),
      color = "#E5E7EB",
      linewidth = 3.4,
      lineend = "round"
    ) +
    geom_segment(
      data = phyla_df,
      aes(
        x = x_genome_bar_start,
        xend = genome_bar_end,
        y = y, yend = y,
        color = color
      ),
      linewidth = 4.4,
      lineend = "round"
    ) +
    geom_text(
      data = phyla_df,
      aes(x = x_genome_text, y = y, label = genome_label, color = color),
      hjust = 0,
      size = 4.2,
      fontface = "bold"
    ) +
    # Headers
    annotate("text", x = x_domain, y = header_y, label = "Domain",
             size = 5.8, fontface = "bold") +
    annotate("text", x = x_kingdom, y = header_y, label = "Kingdom",
             size = 5.8, fontface = "bold") +
    annotate("text", x = x_phylum + 0.55, y = header_y, label = "Phylum",
             hjust = 0, size = 5.8, fontface = "bold") +
    annotate("text", x = mean(c(x_class, x_species)), y = header_y,
             label = "Lower-rank counts",
             size = 5.8, fontface = "bold") +
    annotate("text", x = x_class, y = header_y - 0.72,
             label = "Class", size = 4.2, fontface = "bold") +
    annotate("text", x = x_order, y = header_y - 0.72,
             label = "Order", size = 4.2, fontface = "bold") +
    annotate("text", x = x_family, y = header_y - 0.72,
             label = "Family", size = 4.2, fontface = "bold") +
    annotate("text", x = x_genus, y = header_y - 0.72,
             label = "Genus", size = 4.2, fontface = "bold") +
    annotate("text", x = x_species, y = header_y - 0.72,
             label = "Species", size = 4.2, fontface = "bold") +
    annotate("text", x = x_genome_bar_start, y = header_y,
             label = "Genome count", hjust = 0,
             size = 5.8, fontface = "bold") +
    # Axis for genome bar
    geom_segment(
      aes(
        x = x_genome_bar_start,
        xend = x_genome_bar_end,
        y = axis_y, yend = axis_y
      ),
      linewidth = 0.45,
      color = "black"
    ) +
    annotate("text", x = x_genome_bar_start, y = axis_label_y,
             label = "0", size = 4.8, hjust = 0) +
    annotate("text", x = x_genome_bar_end, y = axis_label_y,
             label = comma(max_genome), size = 4.8, hjust = 1) +
    annotate(
      "text",
      x = -0.45,
      y = title_y,
      label = title_text,
      hjust = 0,
      size = 7.5,
      fontface = "bold",
      color = "#111827"
    ) +
    scale_color_identity() +
    scale_fill_manual(values = c(domain_colors, kingdom_colors)) +
    scale_linewidth(range = c(0.08, 22), guide = "none") +
    coord_cartesian(
      xlim = c(-1.0, ifelse(selected_domain == "Archaea", 22.5, 21.0)),  # expand for Archaea
      ylim = c(min(phyla_df$y) - 2.5, max(phyla_df$y) + 4.5),
      clip = "off"
    ) +
    theme_void(base_size = 14) +
    theme(
      legend.position = "none",
      plot.margin = margin(45, ifelse(selected_domain == "Archaea", 180, 140), 45, 65),
      text = element_text(family = "Helvetica", face = "bold")
    )

  print(p)

  out_height <- ifelse(selected_domain == "Archaea", 18, 36)

  ggsave(
    file.path(outdir, paste0(output_prefix, ".png")),
    p, width = ifelse(selected_domain == "Archaea", 32, 30), height = out_height, dpi = 500, bg = "white"
  )
  ggsave(
    file.path(outdir, paste0(output_prefix, ".pdf")),
    p, width = ifelse(selected_domain == "Archaea", 32, 30), height = out_height, bg = "white"
  )
  ggsave(
    file.path(outdir, paste0(output_prefix, ".svg")),
    p, width = ifelse(selected_domain == "Archaea", 32, 30), height = out_height, bg = "white"
  )

  write_csv(
    phyla_df,
    file.path(outdir, paste0(output_prefix, "_plot_data.csv"))
  )
}

# ============================================================
# RUN ALL FIGURES
# ============================================================

make_figure1_total_domain_kingdom()

make_domain_phylum_plot(
  selected_domain = "Archaea",
  output_prefix = "FIGURE_2_Archaea_Kingdom_Phylum_LowerRanks_RECREATED"
)

make_domain_phylum_plot(
  selected_domain = "Bacteria",
  output_prefix = "FIGURE_3_Bacteria_Kingdom_Phylum_LowerRanks_RECREATED"
)

cat("\nDone! Recreated taxonomy dictionary figures saved in:\n", outdir, "\n")
