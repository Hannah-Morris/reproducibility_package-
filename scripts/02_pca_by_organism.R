# Script: 02_pca_by_organism.R
# Purpose: Run PCA separately for each organism and plot samples coloured by PRJNA project
# Input: physeq_rarefied_filtered (named list of phyloseq objects)
# Output: Faceted PCA plot by organism

library(phyloseq)
library(tidyverse)
library(purrr)
library(Polychrome)
library(colorspace)

## ----------------------------
## 1. Label each phyloseq object with its PRJNA name
## ----------------------------
# physeq_rarefied_filtered must be a named list, with names = PRJNA codes

phy_list_labeled <- imap(physeq_rarefied_filtered, ~{
  ps <- .x
  proj <- .y
  
  if (is.null(sample_data(ps, errorIfNULL = FALSE))) {
    sample_data(ps) <- data.frame(row.names = sample_names(ps))
  }
  
  sample_data(ps)$Project <- proj
  ps
})

## ----------------------------
## 2. Global palette for all PRJNA groups
## ----------------------------
all_projects <- sort(unique(unlist(
  map(phy_list_labeled, ~ as.character(sample_data(.x)$Project))
)))

proj_levels <- all_projects

pal_vec_use <- colorspace::qualitative_hcl(
  n = length(proj_levels),
  palette = "Dynamic"
)
names(pal_vec_use) <- proj_levels

## ----------------------------
## 3. Helper: union OTU matrix across projects within an organism
## ----------------------------
make_otu_union_matrix <- function(ps_list) {
  mats <- map(ps_list, ~{
    otu <- as(otu_table(.x), "matrix")
    if (taxa_are_rows(.x)) otu <- t(otu)
    otu
  })
  
  all_taxa <- sort(unique(unlist(map(mats, colnames))))
  
  mats_aligned <- map(mats, ~{
    miss <- setdiff(all_taxa, colnames(.x))
    if (length(miss) > 0) {
      add <- matrix(
        0,
        nrow = nrow(.x),
        ncol = length(miss),
        dimnames = list(rownames(.x), miss)
      )
      .x <- cbind(.x, add)
    }
    .x[, all_taxa, drop = FALSE]
  })
  
  do.call(rbind, mats_aligned)
}

## ----------------------------
## 4. Run PCA separately for each organism
## ----------------------------
all_organisms <- unique(unlist(
  map(phy_list_labeled, ~ as.character(sample_data(.x)$Organism))
))
all_organisms <- all_organisms[!is.na(all_organisms) & all_organisms != ""]

pca_all_df <- map_dfr(
  all_organisms,
  function(org) {
    
    # Select phyloseq objects containing this organism
    ps_sub <- keep(phy_list_labeled, ~ any(as.character(sample_data(.x)$Organism) == org))
    if (length(ps_sub) == 0) return(NULL)
    
    # Keep only this organism within each object
    ps_sub <- map(ps_sub, ~ prune_samples(as.character(sample_data(.x)$Organism) == org, .x))
    
    # Drop empty objects
    ps_sub <- keep(ps_sub, ~ nsamples(.x) > 1 && ntaxa(.x) > 1)
    if (length(ps_sub) == 0) return(NULL)
    
    # Build union OTU matrix
    otu_union <- make_otu_union_matrix(ps_sub)
    
    # Stabilise counts
    otu_union <- log1p(otu_union)
    
    # Remove zero-variance taxa
    keep_cols <- apply(otu_union, 2, var) > 0
    otu_union <- otu_union[, keep_cols, drop = FALSE]
    if (ncol(otu_union) < 2) return(NULL)
    
    # PCA
    pca <- prcomp(otu_union, center = TRUE, scale. = TRUE)
    
    # Variance explained
    pca_var <- pca$sdev^2
    pca_var_exp <- pca_var / sum(pca_var)
    pc1_pct <- round(100 * pca_var_exp[1], 1)
    pc2_pct <- round(100 * pca_var_exp[2], 1)
    
    scores <- as.data.frame(pca$x) %>%
      tibble::rownames_to_column("SampleID") %>%
      mutate(
        Organism = org,
        pc1_pct = pc1_pct,
        pc2_pct = pc2_pct
      )
    
    # Metadata
    meta <- map_dfr(ps_sub, ~{
      md <- data.frame(sample_data(.x), stringsAsFactors = FALSE)
      md$SampleID <- rownames(md)
      md %>% mutate(across(everything(), as.character))
    }) %>%
      distinct(SampleID, .keep_all = TRUE) %>%
      select(-any_of("Organism"))
    
    left_join(scores, meta, by = "SampleID")
  }
)

pca_all_df$Project <- factor(as.character(pca_all_df$Project), levels = all_projects)
pca_all_df <- pca_all_df %>% mutate(Organism_facet = Organism)

## ----------------------------
## 5. Per-facet variance labels
## ----------------------------
facet_labels <- pca_all_df %>%
  distinct(Organism_facet, pc1_pct, pc2_pct) %>%
  group_by(Organism_facet) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(label = paste0("PC1 = ", pc1_pct, "%\nPC2 = ", pc2_pct, "%"))

## ----------------------------
## 6. Only draw ellipses where enough samples exist
## ----------------------------
ellipse_df <- pca_all_df %>%
  count(Organism_facet, Project) %>%
  filter(n >= 3) %>%
  select(Organism_facet, Project) %>%
  inner_join(pca_all_df, by = c("Organism_facet", "Project"))

## ----------------------------
## 7. Plot
## ----------------------------
p_pca <- ggplot(pca_all_df, aes(PC1, PC2, colour = Project)) +
  geom_point(size = 2.5, alpha = 0.85) +
  stat_ellipse(
    data = ellipse_df,
    aes(group = Project),
    type = "t",
    level = 0.68,
    linewidth = 0.6
  ) +
  facet_wrap(~ Organism_facet, scales = "free") +
  geom_text(
    data = facet_labels,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.1, vjust = 1.1,
    size = 3
  ) +
  scale_colour_manual(values = pal_vec_use, drop = TRUE) +
  theme_minimal() +
  labs(
    title = "PCA of microbial communities by organism",
    subtitle = "Separate PCA per organism; colours indicate PRJNA project",
    colour = "PRJNA Project",
    x = "PC1",
    y = "PC2"
  )

print(p_pca)

# Optional save
# ggsave("PCA_by_organism.png", p_pca, width = 12, height = 8, dpi = 300)