library(tidyverse)
library(dplyr)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(purrr) # Added for imap/reduce/walk

# ============================================================
# CONFIGURATION — add / remove lineages here only
# ============================================================
lineages <- list(
  DEEPOL1 = list(
    asv_list = "data/raw/DEEPOL1_EukBank.list",
    out_dir  = "data/edited/DEEPOL1"
  ),
  DEEPOL2 = list(
    asv_list = "data/raw/DEEPOL2_EukBank.list",
    out_dir  = "data/edited/DEEPOL2"
  ),
  DEEPOL3 = list(
    asv_list = "data/raw/DEEPOL3_EukBank.list",
    out_dir  = "data/edited/DEEPOL3"
  )
)

# Minimum number of samples an ASV must appear in to be retained
MIN_SAMPLES <- 2

# ============================================================
# METADATA PREPARATION (runs once, shared across all lineages)
# ============================================================
metadata_raw <- read.delim("data/raw/eukbank_18S_V4_samples.tsv",
                           header = TRUE, dec = ".")

# ============================================================
# OCEAN ASSIGNMENT
# ============================================================
assign_ocean_basin <- function(lat, lon) {
  lon <- ifelse(lon > 180, lon - 360, lon)
  
  case_when(
    lat >  66.5                                          ~ "Arctic Ocean",
    lat < -40                                            ~ "Southern Ocean",
    between(lat, 30, 46)  & between(lon, -6,  37)       ~ "Mediterranean Region",
    between(lat, 53, 66)  & between(lon, 12,  30)       ~ "Baltic Sea",
    between(lat, 36, 48)  & between(lon, 49,  55)       ~ NA_character_,
    between(lat, -40, 30) & between(lon,  20, 140)      ~ "Indian Ocean",
    lat >= 0  & (lon >= 100 | lon <= -70)               ~ "North Pacific Ocean",
    lat <  0  & (lon >= 140 | lon <= -70)               ~ "South Pacific Ocean",
    lat >= 0  & between(lon, -80, 20)                   ~ "North Atlantic Ocean",
    lat <  0  & between(lon, -70, 20)                   ~ "South Atlantic Ocean",
    TRUE                                                 ~ "Unassigned"
  )
}

marine_samples <- metadata_raw %>%
  filter(envplot %in% c("marine_sediment", "marine_water", "marine_organism", "marine_ice")) %>%
  filter(!is.na(latitude), !is.na(longitude)) %>%
  mutate(ocean = assign_ocean_basin(latitude, longitude))

metadata_raw <- metadata_raw %>%
  left_join(marine_samples %>% select(sample, ocean), by = "sample")

# ============================================================
# HABITAT & SIZE FRACTION PROCESSING
# ============================================================
metadata <- metadata_raw %>%
  mutate(habitat = case_when(
    envplot %in% c("marine_sediment", "marine_water", "marine_organism", "marine_ice") ~ "Marine",
    envplot %in% c("land_freshwater")                                    ~ "Freshwater",
    envplot %in% c("land_soil", "land_sediment",
                   "land_organism", "land_water")                        ~ "Terrestrial",
    envplot %in% c("none")                                               ~ "None"
  )) %>%
  mutate(
    size_fraction_lower_threshold = as.numeric(size_fraction_lower_threshold),
    size_fraction_upper_threshold = as.numeric(
      ifelse(is.na(size_fraction_upper_threshold), Inf, size_fraction_upper_threshold)
    ),
    size_fraction = case_when(
      size_fraction_lower_threshold <= 0.22 & size_fraction_upper_threshold == Inf ~ "Bulk_Total",
      size_fraction_lower_threshold >= 0.2  & size_fraction_upper_threshold <= 3   ~ "Pico",
      size_fraction_lower_threshold >= 3    & size_fraction_upper_threshold <= 20  ~ "Nano",
      size_fraction_lower_threshold >= 20   & size_fraction_upper_threshold <  200 ~ "Micro",
      size_fraction_lower_threshold >= 200  & size_fraction_upper_threshold <  Inf ~ "Meso",
      TRUE ~ "Non_Standard"
    )
  ) %>%
  mutate(salinity_cat = case_when(
    salinity <= 0.5                    ~ "freshwater",
    salinity > 0.5  & salinity <= 30  ~ "brackish",
    salinity > 30   & salinity <= 50  ~ "saline",
    salinity > 50                     ~ "briny",
    is.na(salinity)                   ~ "unknown"
  )) %>%
  mutate(
    ocean_layer_fine = case_when(
      depth <= 200                ~ "photic",
      depth > 200 & depth <= 1000 ~ "twilight",
      depth > 1000                ~ "aphotic",
      TRUE                        ~ NA_character_
    ),
    ocean_layer = case_when(
      depth <= 200 ~ "sunlit",
      depth > 200  ~ "dark",
      TRUE         ~ NA_character_
    )
  )

# ============================================================
# LOAD FULL EUKBANK COUNTS
# ============================================================
eukbank_counts <- read.delim("data/raw/eukbank_18S_V4_counts.tsv",
                             header = TRUE, dec = ".")

# ============================================================
# PER-LINEAGE PROCESSING LOOP
# ============================================================
for (lineage_name in names(lineages)) {
  cfg <- lineages[[lineage_name]]
  dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
  
  target_list <- read.table(cfg$asv_list, header = FALSE, col.names = "amplicon")
  
  lineage_counts <- eukbank_counts %>%
    filter(amplicon %in% target_list$amplicon)
  
  lineage_abundance <- lineage_counts %>%
    left_join(metadata %>% select(sample, total_nreads), by = "sample") %>%
    mutate(rel_abundance = nreads / total_nreads, rel_abundance_pct = rel_abundance * 100)
  
  out_abundance <- file.path(cfg$out_dir, paste0(lineage_name, "_abundance.csv"))
  write.csv(lineage_abundance, file = out_abundance, row.names = FALSE)
  
  lineage_filtered <- lineage_abundance %>%
    group_by(amplicon) %>%
    filter(sum(nreads > 0) > MIN_SAMPLES) %>%
    ungroup()
  
  out_filtered <- file.path(cfg$out_dir, paste0(lineage_name, "_filtered_abundance.csv"))
  write.csv(lineage_filtered, file = out_filtered, row.names = FALSE)
  
  common_samples   <- intersect(lineage_filtered$sample, metadata$sample)
  metadata_lineage <- metadata[metadata$sample %in% common_samples, ]
  
  out_meta <- file.path(cfg$out_dir, paste0(lineage_name, "_samples_edited.csv"))
  write.csv(metadata_lineage, file = out_meta, row.names = FALSE)
}

# ============================================================================
# NEW SECTION: BUILD PRESENCE/ABSENCE TABLE
# ============================================================================
cat("\nBuilding presence/absence table...\n")

pa_list <- imap(lineages, function(cfg, name) {
  # Path to the file we just created in the loop above
  path <- file.path(cfg$out_dir, paste0(name, "_abundance.csv"))
  
  if (!file.exists(path)) stop(sprintf("File not found: %s", path))
  
  read.csv(path, stringsAsFactors = FALSE) %>%
    filter(nreads > 0) %>%
    distinct(sample) %>%
    mutate(!!name := 1L)
})

# Start from all samples in metadata, left-join each lineage
pa_table <- metadata %>%
  select(sample) %>%
  reduce(pa_list,
         .init = .,
         function(acc, df) left_join(acc, df, by = "sample")) %>%
  mutate(across(all_of(names(lineages)), ~ replace_na(.x, 0L)))

cat(sprintf("  Presence/absence table: %d samples x %d lineages\n",
            nrow(pa_table), length(lineages)))

cat("  Detections per lineage:\n")
walk(names(lineages), function(nm) {
  cat(sprintf("    %s: %d samples (%.1f%%)\n",
              nm,
              sum(pa_table[[nm]]),
              100 * mean(pa_table[[nm]])))
})

# Join presence/absence back to full metadata for downstream analyses
meta_pa <- metadata %>%
  left_join(pa_table, by = "sample")

write.csv(meta_pa,
          file      = "data/edited/meta_with_presence_absence.csv",
          row.names = FALSE)

cat("\nAll lineages processed and Presence/Absence table saved.\n")