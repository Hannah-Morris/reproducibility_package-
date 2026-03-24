# Script: 05_core_microbiome.R
# Purpose: Generate global and species-level core microbiome heatmaps
# Input: physeq_rarefied_filtered (named list of phyloseq objects)
# Output: Global and per-species core microbiome plots at Family and Phylum level

library(phyloseq)
library(microbiome)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# =========================
# SETTINGS
# =========================
required_cols <- c(
  "Organism",
  "Host",
  "Habitat",
  "Study",
  "Project",
  "Region_amplified"
)

mean_threshold <- 0.002
prevalences <- seq(0.05, 1, 0.05)
core_cols <- c("yellow", "blue")

# =========================
# HELPER: drop zero-sum samples
# =========================
drop_zero_sum_samples <- function(ps) {
  keep <- sample_sums(ps) > 0
  prune_samples(keep, ps)
}

# =========================
# HELPER 1: Build long composition table
# =========================
build_composition_df <- function(phy_list, rank = "Family", mean_threshold = 0.002) {
  all_dfs <- list()
  
  for (nm in names(phy_list)) {
    cat("\n=== Processing", nm, "at rank", rank, "===\n")
    
    ps <- process_tax_level_safe(
      phy_list[[nm]],
      rank = rank,
      mean_threshold = mean_threshold
    )
    
    if (is.null(ps)) {
      cat("  -> skipped (NULL after processing)\n")
      next
    }
    
    df <- phyloseq::psmelt(ps)
    
    if (!rank %in% colnames(df)) {
      cat("  -> skipped (rank column missing)\n")
      next
    }
    
    keep_cols <- intersect(
      c("Sample", required_cols, rank, "Abundance"),
      colnames(df)
    )
    
    df_rank <- df[, keep_cols, drop = FALSE]
    colnames(df_rank)[colnames(df_rank) == rank] <- "Taxon"
    
    df_rank <- df_rank %>%
      filter(!is.na(Taxon), Taxon != "")
    
    df_rank$PRJNA <- nm
    all_dfs[[nm]] <- df_rank
  }
  
  if (length(all_dfs) == 0) {
    stop("No datasets survived processing at rank ", rank)
  }
  
  bind_rows(all_dfs)
}

# =========================
# HELPER 2: Build combined phyloseq
# =========================
build_phyloseq_from_comp <- function(comp_df, rank = "Family", required_cols = NULL) {
  
  comp_agg <- comp_df %>%
    group_by(Sample, Taxon) %>%
    summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop")
  
  otu_wide <- comp_agg %>%
    pivot_wider(
      names_from = Sample,
      values_from = Abundance,
      values_fill = 0
    )
  
  taxa_names <- otu_wide$Taxon
  otu_mat <- as.matrix(otu_wide[, -1, drop = FALSE])
  rownames(otu_mat) <- taxa_names
  
  otu_tab <- phyloseq::otu_table(otu_mat, taxa_are_rows = TRUE)
  
  meta_cols <- c("Sample", required_cols)
  meta_cols <- intersect(meta_cols, colnames(comp_df))
  
  meta_df <- comp_df %>%
    select(all_of(meta_cols)) %>%
    distinct() %>%
    group_by(Sample) %>%
    summarise(
      across(
        everything(),
        ~ {
          x <- .x[!is.na(.x)]
          if (length(x) == 0) NA_character_ else as.character(x[1])
        }
      ),
      .groups = "drop"
    )
  
  meta_df <- as.data.frame(meta_df)
  rownames(meta_df) <- meta_df$Sample
  meta_df$Sample <- NULL
  
  sd_tab <- phyloseq::sample_data(meta_df)
  
  tax_df <- data.frame(row.names = taxa_names)
  tax_df[[rank]] <- taxa_names
  tax_tab <- phyloseq::tax_table(as.matrix(tax_df))
  
  ps <- phyloseq::phyloseq(otu_tab, tax_tab, sd_tab)
  ps <- drop_zero_sum_samples(ps)
  ps <- microbiome::transform(ps, "compositional")
  
  return(ps)
}

# =========================
# HELPER 3: Plot global core heatmap
# =========================
make_global_core_plot <- function(ps, rank_label, detections_vec, prevalences_vec,
                                  colours_vec, min_prev) {
  
  ps <- drop_zero_sum_samples(ps)
  ps <- microbiome::transform(ps, "compositional")
  ps <- prune_taxa(taxa_sums(ps) > 0, ps)
  
  microbiome::plot_core(
    ps,
    plot.type      = "heatmap",
    prevalences    = prevalences_vec,
    colours        = colours_vec,
    detections     = detections_vec,
    min.prevalence = min_prev
  ) +
    ggtitle(paste0("Global Core Microbiome (", rank_label, ")")) +
    labs(
      subtitle = paste0("Minimum prevalence = ", min_prev, " (across all samples)"),
      x = "Detection Threshold\n(Relative Abundance)"
    )
}

# =========================
# HELPER 4: Plot species core heatmap
# =========================
try_core_plot <- function(ps, sp, rank_label,
                          detections_vec,
                          prevalences_vec,
                          colours_vec,
                          min_prev_vec) {
  
  ps <- drop_zero_sum_samples(ps)
  ps <- microbiome::transform(ps, "compositional")
  ps <- prune_taxa(taxa_sums(ps) > 0, ps)
  
  if (nsamples(ps) < 2 || ntaxa(ps) < 2) return(NULL)
  
  for (mp in min_prev_vec) {
    p <- try(
      microbiome::plot_core(
        ps,
        plot.type      = "heatmap",
        prevalences    = prevalences_vec,
        colours        = colours_vec,
        detections     = detections_vec,
        min.prevalence = mp
      ) +
        ggtitle(paste0("Core Microbiome (", rank_label, ")\n", sp)) +
        labs(
          subtitle = paste0("Minimum prevalence = ", mp),
          x = "Detection Threshold\n(Relative Abundance)"
        ),
      silent = TRUE
    )
    
    if (!inherits(p, "try-error")) return(p)
  }
  
  return(NULL)
}

# =========================
# HELPER 5: Run core workflow
# =========================
run_core_analysis <- function(phy_list,
                              rank_label = "Family",
                              mean_threshold = 0.002,
                              detections_global = c(0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2),
                              detections_species = c(0.0005, 0.001, 0.005, 0.01, 0.02, 0.05),
                              prevalences_vec = prevalences,
                              colours_vec = core_cols,
                              min_prev_global = 0.7,
                              min_prev_try = c(0.8, 0.7, 0.6, 0.5)) {
  
  comp_df <- build_composition_df(
    phy_list = phy_list,
    rank = rank_label,
    mean_threshold = mean_threshold
  )
  
  phy_core <- build_phyloseq_from_comp(
    comp_df = comp_df,
    rank = rank_label,
    required_cols = required_cols
  )
  
  phy_core <- drop_zero_sum_samples(phy_core)
  phy_core <- microbiome::transform(phy_core, "compositional")
  phy_core <- prune_taxa(taxa_sums(phy_core) > 0, phy_core)
  
  global_plot <- make_global_core_plot(
    ps = phy_core,
    rank_label = rank_label,
    detections_vec = detections_global,
    prevalences_vec = prevalences_vec,
    colours_vec = colours_vec,
    min_prev = min_prev_global
  )
  
  species_list <- sort(unique(as.character(sample_data(phy_core)$Organism)))
  species_plots <- list()
  
  for (sp in species_list) {
    ps_sp <- subset_samples(phy_core, Organism == sp)
    ps_sp <- prune_taxa(taxa_sums(ps_sp) > 0, ps_sp)
    
    p <- try_core_plot(
      ps = ps_sp,
      sp = sp,
      rank_label = rank_label,
      detections_vec = detections_species,
      prevalences_vec = prevalences_vec,
      colours_vec = colours_vec,
      min_prev_vec = min_prev_try
    )
    
    if (!is.null(p)) {
      species_plots[[sp]] <- p
    }
  }
  
  list(
    phy_core = phy_core,
    comp_df = comp_df,
    global_plot = global_plot,
    species_plots = species_plots
  )
}

# =========================
# 1. FAMILY CORE ANALYSIS
# =========================
family_core <- run_core_analysis(
  phy_list = physeq_rarefied_filtered,
  rank_label = "Family",
  mean_threshold = mean_threshold,
  detections_global = c(0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2),
  detections_species = c(0.0005, 0.001, 0.005, 0.01, 0.02, 0.05),
  min_prev_global = 0.7,
  min_prev_try = c(0.8, 0.7, 0.6, 0.5)
)

print(family_core$global_plot)

# Optional: show selected species plots
family_order <- c(
  "L.vannamei gut metagenome",
  "P.monodon gut metagenome",
  "P.clarkii gut metagenome",
  "P.leniusculus gut metagenome",
  "M.nipponese gut metagenome",
  "P.argus gut metagenome",
  "C.quadricarinatus gut metagenome",
  "C.cainii gut metagenome"
)

family_plots <- family_core$species_plots[family_order]
family_plots <- family_plots[!vapply(family_plots, is.null, logical(1))]

family_pages <- split(family_plots, ceiling(seq_along(family_plots) / 4))
for (i in seq_along(family_pages)) {
  print(wrap_plots(family_pages[[i]], ncol = 2))
}

# =========================
# 2. PHYLUM CORE ANALYSIS
# =========================
phylum_core <- run_core_analysis(
  phy_list = physeq_rarefied_filtered,
  rank_label = "Phylum",
  mean_threshold = mean_threshold,
  detections_global = c(0.01, 0.02, 0.05, 0.1, 0.2),
  detections_species = c(0.01, 0.02, 0.05, 0.1, 0.2),
  min_prev_global = 0.8,
  min_prev_try = c(0.8, 0.7, 0.6, 0.5)
)

print(phylum_core$global_plot)

phylum_order <- c(
  "L.vannamei gut metagenome",
  "P.monodon gut metagenome",
  "P.clarkii gut metagenome",
  "P.leniusculus gut metagenome",
  "M.nipponese gut metagenome",
  "P.argus gut metagenome",
  "C.quadricarinatus gut metagenome",
  "C.cainii gut metagenome"
)

phylum_plots <- phylum_core$species_plots[phylum_order]
phylum_plots <- phylum_plots[!vapply(phylum_plots, is.null, logical(1))]

phylum_pages <- split(phylum_plots, ceiling(seq_along(phylum_plots) / 4))
for (i in seq_along(phylum_pages)) {
  print(wrap_plots(phylum_pages[[i]], ncol = 2))
}