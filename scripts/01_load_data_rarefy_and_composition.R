# Script: 01_load_data_rarefy_and_composition.R
# Purpose: Load decapod gut microbiome datasets, rarefy reads, filter to bacterial taxa,
# and generate compositional plots across taxonomic ranks.
# Input files: OTU_<PRJNA>.csv, hiera_BLAST_<PRJNA>.csv, SraRunTable_<PRJNA>.csv
# Output objects: physeq_rarefied_filtered, comp_data_list, plots_list, plots_list_by_organism

library(tidyverse)
library(phyloseq)
library(microbiome)
library(vegan)
library(RColorBrewer)
library(patchwork)

# =========================
# 0. SETTINGS
# =========================

data_dir <- "~/Library/CloudStorage/OneDrive-UniversityofKent/all_data"
setwd(data_dir)

required_cols <- c(
  "Organism",
  "Host",
  "Study",
  "Project",
  "Region_amplified",
  "Habitat"
)

prjna_codes <- c(
  "422950","1087723","495902","291010","417739","1012318","780955",
  "731310","749331","717320","797514","577421","600476","682282",
  "553862","540737","831321","505066","609648","573062","354668",
  "549032","579035","689097","429671","1044169"
)

rarefy_depth <- 10000
rarefy_seed <- 123
mean_threshold <- 0.001
tax_levels <- c("Phylum", "Class", "Order", "Family", "Genus")

# =========================
# 1. LOAD AND RAREFY DATA
# =========================

physeq_list <- list()
physeq_rarefied_list <- list()

for (prjna in prjna_codes) {
  
  cat("\n===============================\n")
  cat("Processing PRJNA", prjna, "\n")
  cat("===============================\n")
  
  otu_file  <- paste0("OTU_", prjna, ".csv")
  tax_file  <- paste0("hiera_BLAST_", prjna, ".csv")
  meta_file <- paste0("SraRunTable_", prjna, ".csv")
  
  otu_df  <- read.csv(otu_file, row.names = 1, check.names = FALSE)
  tax_df  <- read.csv(tax_file, row.names = 1, check.names = FALSE)
  meta_df <- read.csv(meta_file, check.names = FALSE, stringsAsFactors = FALSE)
  
  stopifnot("Run" %in% colnames(meta_df))
  
  meta_df$Run <- trimws(meta_df$Run)
  meta_df <- meta_df[!is.na(meta_df$Run) & meta_df$Run != "", , drop = FALSE]
  
  if (anyDuplicated(meta_df$Run)) {
    dup <- unique(meta_df$Run[duplicated(meta_df$Run)])
    stop("Duplicate Run IDs in ", meta_file, ": ", paste(dup, collapse = ", "))
  }
  
  rownames(meta_df) <- meta_df$Run
  
  missing_cols <- setdiff(required_cols, colnames(meta_df))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required metadata columns in ", meta_file, ": ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  common_otus <- intersect(rownames(otu_df), rownames(tax_df))
  
  dropped_otu <- setdiff(rownames(otu_df), common_otus)
  dropped_tax <- setdiff(rownames(tax_df), common_otus)
  
  if (length(dropped_otu) > 0) {
    cat("Dropping", length(dropped_otu), "OTUs not found in taxonomy for PRJNA", prjna, "\n")
    otu_df <- otu_df[common_otus, , drop = FALSE]
  }
  
  if (length(dropped_tax) > 0) {
    cat("Dropping", length(dropped_tax), "taxonomy rows not found in OTU table for PRJNA", prjna, "\n")
    tax_df <- tax_df[common_otus, , drop = FALSE]
  }
  
  otu_tab  <- otu_table(as.matrix(otu_df), taxa_are_rows = TRUE)
  tax_tab  <- tax_table(as.matrix(tax_df))
  meta_tab <- sample_data(meta_df)
  
  colnames(otu_tab) <- rownames(meta_df)
  
  physeq <- phyloseq(otu_tab, tax_tab, meta_tab)
  physeq_list[[prjna]] <- physeq
  
  physeq_rarefied <- rarefy_even_depth(
    physeq,
    sample.size = rarefy_depth,
    rngseed = rarefy_seed,
    verbose = FALSE
  )
  
  physeq_rarefied_list[[prjna]] <- physeq_rarefied
  
  cat("Taxa after rarefaction:", ntaxa(physeq_rarefied), "\n")
}

physeq_rarefied_filtered <- Filter(
  function(x) ntaxa(x) > 0 && nsamples(x) > 0,
  physeq_rarefied_list
)

# =========================
# 2. HELPER FUNCTIONS
# =========================

process_tax_level_safe <- function(phy, rank = "Phylum", mean_threshold = 0.002) {
  
  if (ntaxa(phy) == 0 || nsamples(phy) == 0) return(NULL)
  
  tax_df <- as.data.frame(as.matrix(tax_table(phy)))
  
  if (!"Domain" %in% colnames(tax_df)) {
    cat("No Domain column found; skipping dataset.\n")
    return(NULL)
  }
  
  keep_taxa <- rownames(tax_df)[tax_df$Domain == "Bacteria"]
  keep_taxa <- keep_taxa[!is.na(keep_taxa)]
  
  if (length(keep_taxa) < 1) {
    cat("No bacterial taxa left; skipping dataset.\n")
    return(NULL)
  }
  
  phy <- prune_taxa(keep_taxa, phy)
  
  if (ntaxa(phy) < 1) {
    cat("Dataset empty after bacteria filter; skipping.\n")
    return(NULL)
  }
  
  relab <- transform_sample_counts(phy, function(x) x / sum(x))
  otu_mat <- as.matrix(otu_table(relab))
  otu_mat[is.na(otu_mat)] <- 0
  
  if (nrow(otu_mat) < 1) {
    cat("OTU table empty after relative abundance transform; skipping.\n")
    return(NULL)
  }
  
  keep <- apply(otu_mat, 1, function(x) mean(x, na.rm = TRUE) > mean_threshold)
  keep[is.na(keep)] <- FALSE
  
  n_keep <- sum(keep)
  cat("Taxa kept after abundance filter:", n_keep, "/", nrow(otu_mat), "\n")
  
  if (n_keep < 1) {
    cat("All taxa removed by abundance filter; skipping dataset.\n")
    return(NULL)
  }
  
  phy_filt <- prune_taxa(names(keep)[keep], phy)
  
  if (ntaxa(phy_filt) < 1) {
    cat("Dataset empty after prune_taxa; skipping.\n")
    return(NULL)
  }
  
  agg <- tryCatch(
    microbiome::aggregate_taxa(phy_filt, level = rank),
    error = function(e) {
      cat("aggregate_taxa failed:", e$message, "\n")
      NULL
    }
  )
  
  if (is.null(agg) || ntaxa(agg) < 1) return(NULL)
  
  comp <- tryCatch(
    microbiome::transform(agg, "compositional"),
    error = function(e) {
      cat("transform('compositional') failed:", e$message, "\n")
      NULL
    }
  )
  
  if (is.null(comp) || ntaxa(comp) < 1) return(NULL)
  
  return(comp)
}

build_composition_df <- function(phy_list, rank = "Phylum", mean_threshold = 0.001) {
  
  all_dfs <- list()
  
  for (nm in names(phy_list)) {
    cat("\n=== Processing", nm, "at rank", rank, "===\n")
    
    ps <- process_tax_level_safe(phy_list[[nm]], rank = rank, mean_threshold = mean_threshold)
    if (is.null(ps)) next
    
    df <- phyloseq::psmelt(ps)
    
    if (!all(required_cols %in% colnames(df))) {
      cat("Skipping", nm, "- missing required metadata columns in melted data.\n")
      next
    }
    
    if (!rank %in% colnames(df)) {
      cat("Skipping", nm, "- rank", rank, "not found in psmelt output.\n")
      next
    }
    
    df_rank <- df %>%
      dplyr::select(
        Sample,
        Organism,
        Host,
        Study,
        Project,
        Region_amplified,
        Habitat,
        Taxon = dplyr::all_of(rank),
        Abundance
      ) %>%
      filter(!is.na(Taxon), Taxon != "")
    
    all_dfs[[nm]] <- df_rank
  }
  
  if (length(all_dfs) == 0) {
    stop("No datasets survived processing at rank ", rank)
  }
  
  bind_rows(all_dfs, .id = "PRJNA")
}

make_tax_palette <- function(taxa_names) {
  uniq <- unique(taxa_names)
  n <- length(uniq)
  
  if (n <= 12) {
    cols <- brewer.pal(max(n, 3), "Set3")[seq_len(n)]
  } else {
    cols <- colorRampPalette(brewer.pal(12, "Set3"))(n)
  }
  
  names(cols) <- uniq
  cols
}

# =========================
# 3. COMPOSITION PLOTS BY TAXONOMIC RANK
# =========================

plots_list <- list()
comp_data_list <- list()

for (rank in tax_levels) {
  
  cat("\n===============================\n")
  cat("Processing rank:", rank, "\n")
  cat("===============================\n")
  
  comp_df <- build_composition_df(
    phy_list = physeq_rarefied_filtered,
    rank = rank,
    mean_threshold = mean_threshold
  )
  
  comp_data_list[[rank]] <- comp_df
  
  comp_org <- comp_df %>%
    group_by(Organism, Taxon) %>%
    summarise(mean_abund = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
    group_by(Organism) %>%
    mutate(Abundance = mean_abund / sum(mean_abund)) %>%
    ungroup()
  
  pal <- make_tax_palette(comp_org$Taxon)
  
  p <- ggplot(comp_org, aes(x = Organism, y = Abundance, fill = Taxon)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = pal, drop = FALSE) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      panel.grid.major.x = element_blank()
    ) +
    ggtitle(paste("0.1% Composition at", rank)) +
    ylab("Mean relative abundance") +
    xlab("Organism")
  
  plots_list[[rank]] <- p
}

# View plots as needed:
# plots_list[["Phylum"]]
# plots_list[["Class"]]
# plots_list[["Order"]]
# plots_list[["Family"]]
# plots_list[["Genus"]]

# =========================
# 4. COMPOSITION PLOTS BY ORGANISM WITH REGION STRIP
# =========================

rank_to_plot <- "Phylum"

comp_df <- comp_data_list[[rank_to_plot]]
if (is.null(comp_df)) stop("comp_data_list does not contain this rank")

plot_df <- comp_df %>%
  mutate(
    PRJNA = factor(PRJNA),
    Region_amplified = factor(Region_amplified)
  ) %>%
  group_by(Organism, PRJNA, Region_amplified, Taxon) %>%
  summarise(Abundance = mean(Abundance), .groups = "drop") %>%
  group_by(Organism, PRJNA) %>%
  mutate(Abundance = Abundance / sum(Abundance)) %>%
  ungroup()

tax_pal <- make_tax_palette(plot_df$Taxon)

region_pal <- brewer.pal(max(3, length(unique(plot_df$Region_amplified))), "Set2")
names(region_pal) <- unique(plot_df$Region_amplified)

organisms <- unique(plot_df$Organism)
plots_list_by_organism <- list()

for (org in organisms) {
  
  df_sub <- plot_df %>% filter(Organism == org)
  
  p_main <- ggplot(df_sub, aes(x = PRJNA, y = Abundance, fill = Taxon)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = tax_pal, drop = FALSE) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank(),
      plot.margin = margin(5, 5, 0, 5)
    ) +
    labs(
      title = paste("Taxonomic Composition at", rank_to_plot, "\nOrganism:", org),
      x = NULL,
      y = "Relative abundance"
    )
  
  region_df <- df_sub %>%
    distinct(PRJNA, Region_amplified)
  
  p_region <- ggplot(region_df, aes(x = PRJNA, y = 1, fill = Region_amplified)) +
    geom_tile() +
    scale_fill_manual(values = region_pal, drop = FALSE) +
    theme_minimal(base_size = 12) +
    theme(
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(0, 5, 5, 5),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10)
    ) +
    guides(fill = guide_legend(title = "Region amplified"))
  
  combined <- p_main / p_region + plot_layout(heights = c(4, 0.5))
  
  plots_list_by_organism[[org]] <- combined
  
  # Save individual plots if needed:
  # ggsave(
  #   filename = paste0("composition_", rank_to_plot, "_", gsub(" ", "_", org), "_with_region.png"),
  #   plot = combined,
  #   width = 10,
  #   height = 7,
  #   dpi = 300
  # )
}

# View individual plots as needed:
# plots_list_by_organism[["L.vannamei gut metagenome"]]

# Create multi-panel pages manually if needed, for example:
# wrap_plots(plots_list_by_organism, ncol = 2)