# Load required libraries
library(tidyverse)
library(stringr)
library(ggplot2)
library(scales)
library(svglite)

# Input file
taxonomy_file <- "final_clean_taxonomy.csv"

# Correctly formatted official_phyla vector (no inline comments)
official_phyla <- c(
  # Carl Woese / major historically recognized phyla
  "Pseudomonadota", "Campylobacterota", "Myxococcota", "Bdellovibrionota",
  "Actinomycetota", "Bacillota", "Mycoplasmatota", "Cyanobacteriota",
  "Spirochaetota", "Chlorobiota", "Bacteroidota", "Planctomycetota",
  "Chlamydiota", "Deinococcota", "Chloroflexota", "Methanobacteriota",
  "Thermoproteota",
  # Additional official LPSN phyla
  "Abditibacteriota", "Acidobacteriota", "Aquificota", "Armatimonadota",
  "Atribacterota", "Balneolota", "Caldisericota", "Calditrichota",
  "Chrysiogenota", "Coprothermobacterota", "Deferribacterota",
  "Elusimicrobiota", "Fibrobacterota", "Fidelibacterota", "Fusobacteriota",
  "Gemmatimonadota", "Kiritimatiellota", "Lentisphaerota", "Minisyncoccota",
  "Nitrospinota", "Nitrospirota", "Rhodothermota", "Synergistota",
  "Thermodesulfobacteriota", "Thermodesulfobiota", "Thermomicrobiota",
  "Thermosulfidibacterota", "Thermotogota", "Verrucomicrobiota",
  "Vulcanimicrobiota", "Nanobdellota", "Nitrososphaerota", "Thermoplasmatota",
  # Extra phyla retained from previous analysis
  "Bipolaricaulota", "Desulfobacterota", "Dictyoglomerota", "Dictyoglomota",
  "Halobacteriota", "Ignavibacteriota", "Poseidoniota", "Microcaldota",
  "Promethearchaeota"
)

# Woese group mapping
woese_mapping <- tribble(
  ~woese_group, ~phylum,
  "1. Purple Bacteria", "Pseudomonadota",
  "1. Purple Bacteria", "Campylobacterota",
  "1. Purple Bacteria", "Myxococcota",
  "1. Purple Bacteria", "Bdellovibrionota",
  "1. Purple Bacteria", "Thermodesulfobacteriota",
  "2. Gram-positive High G+C", "Actinomycetota",
  "3. Gram-positive Low G+C", "Bacillota",
  "3. Gram-positive Low G+C", "Mycoplasmatota",
  "4. Cyanobacteria", "Cyanobacteriota",
  "5. Spirochetes", "Spirochaetota",
  "6. Green Sulfur Bacteria", "Chlorobiota",
  "7. Bacteroides-Flavobacteria", "Bacteroidota",
  "8. Planctomyces", "Planctomycetota",
  "9. Chlamydiae", "Chlamydiota",
  "10. Radio-resistant Micrococci", "Deinococcota",
  "11. Green Non-sulfur Bacteria", "Chloroflexota",
  "12. Archaea", "Methanobacteriota",
  "12. Archaea", "Thermoproteota"
)

# Read and clean taxonomy
taxonomy_df <- read_csv(taxonomy_file, show_col_types = FALSE)
names(taxonomy_df) <- names(taxonomy_df) |> str_trim() |> str_to_lower()
taxonomy_df <- taxonomy_df %>%
  mutate(
    accession = str_trim(as.character(accession)),
    domain = str_trim(as.character(domain)),
    phylum = str_trim(as.character(phylum))
  ) %>%
  filter(accession != "", domain != "", phylum != "") %>%
  distinct(accession, .keep_all = TRUE)

# Count genomes per phylum
phylum_counts <- taxonomy_df %>%
  filter(phylum %in% official_phyla) %>%
  group_by(domain, phylum) %>%
  summarise(genome_count = n_distinct(accession), .groups = "drop")

# Split into Woese and new phyla
woese_df <- phylum_counts %>%
  inner_join(woese_mapping, by = "phylum")

new_df <- phylum_counts %>%
  filter(!phylum %in% woese_mapping$phylum) %>%
  mutate(new_group = ifelse(domain == "Archaea", "New Archaeal Phyla", "New Bacterial Phyla"))

# Helper function for bezier curves
curve_fun <- function(x0, y0, x1, y1, id, val, col) {
  t <- seq(0, 1, length.out = 180)
  tibble(
    x = (1 - t)^3 * x0 + 3 * (1 - t)^2 * t * (x0 + 0.65 * (x1 - x0)) +
      3 * (1 - t) * t^2 * (x1 - 0.65 * (x1 - x0)) + t^3 * x1,
    y = (1 - t)^3 * y0 + 3 * (1 - t)^2 * t * y0 +
      3 * (1 - t) * t^2 * y1 + t^3 * y1,
    id = id, val = val, col = col
  )
}

# ---------- Figure A: Woese groups ----------
woese_colors <- c(
  "1. Purple Bacteria" = "#7B6FD0",
  "2. Gram-positive High G+C" = "#55B7A5",
  "3. Gram-positive Low G+C" = "#9BCB73",
  "4. Cyanobacteria" = "#6BC7D6",
  "5. Spirochetes" = "#E6B85C",
  "6. Green Sulfur Bacteria" = "#B8D97A",
  "7. Bacteroides-Flavobacteria" = "#E7A36A",
  "8. Planctomyces" = "#E8A2C8",
  "9. Chlamydiae" = "#C79A6B",
  "10. Radio-resistant Micrococci" = "#E67AAE",
  "11. Green Non-sulfur Bacteria" = "#7DB7E8",
  "12. Archaea" = "#E58C8A"
)
woese_order <- names(woese_colors)

woese_plot_df <- woese_df %>%
  mutate(woese_group = factor(woese_group, levels = woese_order)) %>%
  arrange(woese_group, desc(genome_count), phylum) %>%
  mutate(
    y_phylum = rev(row_number()) * 1.4,
    color = woese_colors[as.character(woese_group)],
    phylum_label = paste0(phylum, " (", comma(genome_count), ")")
  )

woese_node_df <- woese_plot_df %>%
  group_by(woese_group, color) %>%
  summarise(
    total = sum(genome_count),
    y_group = mean(y_phylum),
    .groups = "drop"
  ) %>%
  mutate(label = paste0(as.character(woese_group), "\n(", comma(total), ")"))

woese_plot_df <- woese_plot_df %>%
  left_join(woese_node_df, by = c("woese_group", "color"))

x_left <- 0
x_right <- 6.8
x_label <- 7.25

curves <- list()
for (i in seq_len(nrow(woese_plot_df))) {
  curves[[i]] <- curve_fun(
    x_left, woese_plot_df$y_group[i],
    x_right, woese_plot_df$y_phylum[i],
    paste0("curve_", i),
    woese_plot_df$genome_count[i],
    woese_plot_df$color[i]
  )
}
curve_df <- bind_rows(curves)

p1 <- ggplot() +
  geom_path(
    data = curve_df,
    aes(x = x, y = y, group = id, color = col, linewidth = val),
    alpha = 0.28, lineend = "round"
  ) +
  geom_label(
    data = woese_node_df,
    aes(x = x_left, y = y_group, label = label, fill = color),
    size = 4.8, fontface = "bold", color = "white",
    label.size = 0, label.padding = unit(0.35, "lines"), lineheight = 0.9
  ) +
  geom_point(
    data = woese_plot_df,
    aes(x = x_right, y = y_phylum, color = color),
    size = 4
  ) +
  geom_text(
    data = woese_plot_df,
    aes(x = x_label, y = y_phylum, label = phylum_label),
    hjust = 0, vjust = 0.5, size = 4.7, fontface = "bold", color = "#111827"
  ) +
  annotate(
    "text", x = x_left, y = max(woese_plot_df$y_phylum) + 2.2,
    label = "Carl Woese Original Groups", size = 6.2, fontface = "bold"
  ) +
  annotate(
    "text", x = x_label, y = max(woese_plot_df$y_phylum) + 2.2,
    label = "Modern Phyla", hjust = 0, size = 6.2, fontface = "bold"
  ) +
  annotate(
    "text", x = -0.8, y = max(woese_plot_df$y_phylum) + 4.3,
    label = "A. Carl Woese Original 12 Groups", hjust = 0,
    size = 8.5, fontface = "bold", color = "#111827"
  ) +
  scale_color_identity() +
  scale_fill_identity() +
  scale_linewidth(range = c(0.5, 20), guide = "none") +
  coord_cartesian(
    xlim = c(-1.2, 11.5),
    ylim = c(min(woese_plot_df$y_phylum) - 2, max(woese_plot_df$y_phylum) + 5),
    clip = "off"
  ) +
  theme_void(base_size = 17) +
  theme(
    legend.position = "none",
    plot.margin = margin(40, 180, 40, 80),
    text = element_text(family = "Helvetica", face = "bold")
  )

print(p1)
ggsave("FIGURE_A_WOESE_ORIGINAL_GROUPS_SOFT_COLORS.png", p1, width = 18, height = 18, dpi = 500, bg = "white")
ggsave("FIGURE_A_WOESE_ORIGINAL_GROUPS_SOFT_COLORS.pdf", p1, width = 18, height = 18, bg = "white")
ggsave("FIGURE_A_WOESE_ORIGINAL_GROUPS_SOFT_COLORS.svg", p1, width = 18, height = 18, bg = "white")

# ---------- Figure B: New phyla only ----------
new_colors <- c(
  "New Bacterial Phyla" = "#A9B8D6",
  "New Archaeal Phyla" = "#C6B0E3"
)

new_plot_df <- new_df %>%
  mutate(new_group = factor(new_group, levels = c("New Bacterial Phyla", "New Archaeal Phyla"))) %>%
  arrange(new_group, phylum) %>%
  mutate(
    y_phylum = rev(row_number()) * 1.15,
    color = new_colors[as.character(new_group)],
    phylum_label = paste0(phylum, " (", comma(genome_count), ")")
  )

new_node_df <- new_plot_df %>%
  group_by(new_group, color) %>%
  summarise(
    total = sum(genome_count),
    y_group = mean(y_phylum),
    .groups = "drop"
  ) %>%
  mutate(label = paste0(as.character(new_group), "\n(", comma(total), ")"))

new_plot_df <- new_plot_df %>%
  left_join(new_node_df, by = c("new_group", "color"))

curves2 <- list()
for (i in seq_len(nrow(new_plot_df))) {
  curves2[[i]] <- curve_fun(
    x_left, new_plot_df$y_group[i],
    x_right, new_plot_df$y_phylum[i],
    paste0("curve_", i),
    new_plot_df$genome_count[i],
    new_plot_df$color[i]
  )
}
curve_df2 <- bind_rows(curves2)

p2 <- ggplot() +
  geom_path(
    data = curve_df2,
    aes(x = x, y = y, group = id, color = col, linewidth = val),
    alpha = 0.30, lineend = "round"
  ) +
  geom_label(
    data = new_node_df,
    aes(x = x_left, y = y_group, label = label, fill = color),
    size = 5.1, fontface = "bold", color = "white",
    label.size = 0, label.padding = unit(0.40, "lines"), lineheight = 0.9
  ) +
  geom_point(
    data = new_plot_df,
    aes(x = x_right, y = y_phylum, color = color),
    size = 3.4
  ) +
  geom_text(
    data = new_plot_df,
    aes(x = x_label, y = y_phylum, label = phylum_label),
    hjust = 0, vjust = 0.5, size = 4.1, fontface = "bold", color = "#111827"
  ) +
  annotate(
    "text", x = -0.8, y = max(new_plot_df$y_phylum) + 4.0,
    label = "B. New Phyla Outside the Original Woese Groups", hjust = 0,
    size = 8.2, fontface = "bold", color = "#111827"
  ) +
  scale_color_identity() +
  scale_fill_identity() +
  scale_linewidth(range = c(0.4, 16), guide = "none") +
  coord_cartesian(
    xlim = c(-1.2, 12),
    ylim = c(min(new_plot_df$y_phylum) - 2, max(new_plot_df$y_phylum) + 5),
    clip = "off"
  ) +
  theme_void(base_size = 17) +
  theme(
    legend.position = "none",
    plot.margin = margin(40, 180, 40, 80),
    text = element_text(family = "Helvetica", face = "bold")
  )

print(p2)
ggsave("FIGURE_B_NEW_PHYLA_ONLY_SOFT_COLORS.png", p2, width = 18, height = 26, dpi = 500, bg = "white")
ggsave("FIGURE_B_NEW_PHYLA_ONLY_SOFT_COLORS.pdf", p2, width = 18, height = 26, bg = "white")
ggsave("FIGURE_B_NEW_PHYLA_ONLY_SOFT_COLORS.svg", p2, width = 18, height = 26, bg = "white")

# Export data
write_csv(woese_plot_df, "FIGURE_A_WOESE_ORIGINAL_GROUPS_DATA.csv")
write_csv(new_plot_df, "FIGURE_B_NEW_PHYLA_ONLY_DATA.csv")

cat("\nDone! Two separate Sankey figures created with softer publication-style colors.\n")

library(tidyverse)
library(stringr)
library(ggplot2)
library(scales)
library(svglite)

taxonomy_file <- "final_clean_taxonomy.csv"

official_phyla <- c(
  "Pseudomonadota","Campylobacterota","Myxococcota","Bdellovibrionota",
  "Actinomycetota","Bacillota","Mycoplasmatota","Cyanobacteriota",
  "Spirochaetota","Chlorobiota","Bacteroidota","Planctomycetota",
  "Chlamydiota","Deinococcota","Chloroflexota","Methanobacteriota",
  "Thermoproteota","Thermodesulfobacteriota",
  "Abditibacteriota","Acidobacteriota","Aquificota","Armatimonadota",
  "Atribacterota","Balneolota","Caldisericota","Calditrichota",
  "Chrysiogenota","Coprothermobacterota","Deferribacterota",
  "Elusimicrobiota","Fibrobacterota","Fidelibacterota","Fusobacteriota",
  "Gemmatimonadota","Kiritimatiellota","Lentisphaerota","Minisyncoccota",
  "Nitrospinota","Nitrospirota","Rhodothermota","Synergistota",
  "Thermodesulfobiota","Thermomicrobiota","Thermosulfidibacterota",
  "Thermotogota","Verrucomicrobiota","Vulcanimicrobiota",
  "Nanobdellota","Nitrososphaerota","Thermoplasmatota",
  "Bipolaricaulota","Desulfobacterota","Dictyoglomerota","Dictyoglomota",
  "Halobacteriota","Ignavibacteriota","Poseidoniota","Microcaldota",
  "Promethearchaeota"
)

woese_mapping <- tribble(
  ~woese_group, ~phylum,
  "1. Purple Bacteria", "Pseudomonadota",
  "1. Purple Bacteria", "Campylobacterota",
  "1. Purple Bacteria", "Myxococcota",
  "1. Purple Bacteria", "Bdellovibrionota",
  "1. Purple Bacteria", "Thermodesulfobacteriota",
  "2. Gram-positive High G+C", "Actinomycetota",
  "3. Gram-positive Low G+C", "Bacillota",
  "3. Gram-positive Low G+C", "Mycoplasmatota",
  "4. Cyanobacteria", "Cyanobacteriota",
  "5. Spirochetes", "Spirochaetota",
  "6. Green Sulfur Bacteria", "Chlorobiota",
  "7. Bacteroides-Flavobacteria", "Bacteroidota",
  "8. Planctomyces", "Planctomycetota",
  "9. Chlamydiae", "Chlamydiota",
  "10. Radio-resistant Micrococci", "Deinococcota",
  "11. Green Non-sulfur Bacteria", "Chloroflexota",
  "12. Archaea", "Methanobacteriota",
  "12. Archaea", "Thermoproteota"
)

synonym_df <- tribble(
  ~phylum, ~synonym_label,
  "Cyanobacteriota", "Cyanobacteriophyta (synonym)",
  "Thermodesulfobacteriota", "Desulfobacterota (synonym)",
  "Methanobacteriota", "Halobacteriota (synonym)"
)

taxonomy_df <- read_csv(taxonomy_file, show_col_types = FALSE)

names(taxonomy_df) <- names(taxonomy_df) |>
  str_trim() |>
  str_to_lower()

taxonomy_df <- taxonomy_df %>%
  mutate(
    accession = str_trim(as.character(accession)),
    domain = str_trim(as.character(domain)),
    phylum = str_trim(as.character(phylum))
  ) %>%
  filter(accession != "", domain != "", phylum != "") %>%
  distinct(accession, .keep_all = TRUE)

phylum_counts <- taxonomy_df %>%
  filter(phylum %in% official_phyla) %>%
  group_by(domain, phylum) %>%
  summarise(genome_count = n_distinct(accession), .groups = "drop")

woese_df <- phylum_counts %>%
  inner_join(woese_mapping, by = "phylum")

new_df <- phylum_counts %>%
  filter(!phylum %in% woese_mapping$phylum) %>%
  mutate(
    new_group = ifelse(domain == "Archaea", "New Archaeal Phyla", "New Bacterial Phyla")
  )

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

woese_colors <- c(
  "1. Purple Bacteria" = "#7B6FD0",
  "2. Gram-positive High G+C" = "#55B7A5",
  "3. Gram-positive Low G+C" = "#9BCB73",
  "4. Cyanobacteria" = "#6BC7D6",
  "5. Spirochetes" = "#E6B85C",
  "6. Green Sulfur Bacteria" = "#B8D97A",
  "7. Bacteroides-Flavobacteria" = "#E7A36A",
  "8. Planctomyces" = "#E8A2C8",
  "9. Chlamydiae" = "#C79A6B",
  "10. Radio-resistant Micrococci" = "#E67AAE",
  "11. Green Non-sulfur Bacteria" = "#7DB7E8",
  "12. Archaea" = "#E58C8A"
)

woese_order <- names(woese_colors)

woese_plot_df <- woese_df %>%
  mutate(
    woese_group = factor(woese_group, levels = woese_order)
  ) %>%
  arrange(woese_group, desc(genome_count), phylum) %>%
  mutate(
    y_phylum = rev(row_number()) * 1.4,
    color = woese_colors[as.character(woese_group)],
    main_label = paste0(phylum, " (", comma(genome_count), ")")
  ) %>%
  left_join(synonym_df, by = "phylum")

woese_node_df <- woese_plot_df %>%
  group_by(woese_group, color) %>%
  summarise(
    total = sum(genome_count),
    y_group = mean(y_phylum),
    .groups = "drop"
  ) %>%
  mutate(
    label = paste0(as.character(woese_group), "\n(", comma(total), ")")
  )

woese_plot_df <- woese_plot_df %>%
  left_join(woese_node_df, by = c("woese_group", "color"))

x_left <- 0
x_right <- 6.8
x_label <- 7.25

curves <- list()

for (i in seq_len(nrow(woese_plot_df))) {
  curves[[i]] <- curve_fun(
    x_left,
    woese_plot_df$y_group[i],
    x_right,
    woese_plot_df$y_phylum[i],
    paste0("curve_", i),
    woese_plot_df$genome_count[i],
    woese_plot_df$color[i]
  )
}

curve_df <- bind_rows(curves)

p1 <- ggplot() +
  geom_path(
    data = curve_df,
    aes(x = x, y = y, group = id, color = col, linewidth = val),
    alpha = 0.28,
    lineend = "round"
  ) +
  geom_label(
    data = woese_node_df,
    aes(x = x_left, y = y_group, label = label, fill = color),
    size = 4.8,
    fontface = "bold",
    color = "white",
    label.size = 0,
    label.padding = unit(0.35, "lines"),
    lineheight = 0.9
  ) +
  geom_point(
    data = woese_plot_df,
    aes(x = x_right, y = y_phylum, color = color),
    size = 4
  ) +
  geom_text(
    data = woese_plot_df,
    aes(x = x_label, y = y_phylum + 0.14, label = main_label),
    hjust = 0,
    vjust = 0.5,
    size = 4.7,
    fontface = "bold",
    color = "#111827"
  ) +
  geom_text(
    data = woese_plot_df %>% filter(!is.na(synonym_label)),
    aes(x = x_label, y = y_phylum - 0.32, label = synonym_label),
    hjust = 0,
    vjust = 0.5,
    size = 3.8,
    fontface = "bold",
    color = "grey45"
  ) +
  annotate(
    "text",
    x = x_left,
    y = max(woese_plot_df$y_phylum) + 2.2,
    label = "Carl Woese Original Groups",
    size = 6.2,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = x_label,
    y = max(woese_plot_df$y_phylum) + 2.2,
    label = "Modern Phyla",
    hjust = 0,
    size = 6.2,
    fontface = "bold"
  ) +
  annotate(
    "text",
    x = -0.8,
    y = max(woese_plot_df$y_phylum) + 4.3,
    label = "A.Carl Woese Original 12 Groups −> Modern Phyla",
    hjust = 0,
    size = 8.5,
    fontface = "bold",
    color = "#111827"
  ) +
  scale_color_identity() +
  scale_fill_identity() +
  scale_linewidth(range = c(0.5, 20), guide = "none") +
  coord_cartesian(
    xlim = c(-1.2, 12.3),
    ylim = c(min(woese_plot_df$y_phylum) - 2, max(woese_plot_df$y_phylum) + 5),
    clip = "off"
  ) +
  theme_void(base_size = 17) +
  theme(
    legend.position = "none",
    plot.margin = margin(40, 220, 40, 80),
    text = element_text(family = "Helvetica", face = "bold")
  )

print(p1)

ggsave("FIGURE_A_WOESE_ORIGINAL_GROUPS_WITH_SYNONYMS.png", p1, width = 18, height = 18, dpi = 500, bg = "white")
ggsave("FIGURE_A_WOESE_ORIGINAL_GROUPS_WITH_SYNONYMS.pdf", p1, width = 18, height = 18, bg = "white")
ggsave("FIGURE_A_WOESE_ORIGINAL_GROUPS_WITH_SYNONYMS.svg", p1, width = 18, height = 18, bg = "white")

write_csv(woese_plot_df, "FIGURE_A_WOESE_ORIGINAL_GROUPS_WITH_SYNONYMS_DATA.csv")

cat("\nDone! Figure A created with synonym labels under accepted phyla.\n")