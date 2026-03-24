# Script: 04_beta_diversity_pcoa.R
# Purpose: Generate PCoA plots across multiple beta diversity distance metrics
# Input: phy_all (merged phyloseq object with Organism in sample_data)
# Output: Combined PCoA figure with PERMANOVA summaries

library(phyloseq)
library(vegan)
library(ggplot2)
library(patchwork)

## ----------------------------
## 1. Input settings
## ----------------------------
ps <- phy_all
group_var <- "Organism"
methods <- c("bray", "jaccard_pa", "canberra", "ruzicka")

## ----------------------------
## 2. Helper function
## ----------------------------
make_pcoa_plot <- function(ps, method, group_var = "Organism") {
  
  # Distance matrix
  if (method == "jaccard_pa") {
    otu <- as(otu_table(ps), "matrix")
    if (taxa_are_rows(ps)) otu <- t(otu)
    dist_mat <- vegan::vegdist(otu, method = "jaccard", binary = TRUE)
    method_title <- "Jaccard (presence/absence)"
    
  } else if (method == "ruzicka") {
    otu <- as(otu_table(ps), "matrix")
    if (taxa_are_rows(ps)) otu <- t(otu)
    dist_mat <- vegan::vegdist(otu, method = "jaccard", binary = FALSE)
    method_title <- "Ružička (quantitative Jaccard)"
    
  } else {
    dist_mat <- phyloseq::distance(ps, method = method)
    method_title <- toupper(method)
  }
  
  # Align metadata to distance labels
  samples <- attr(dist_mat, "Labels")
  meta <- data.frame(sample_data(ps))[samples, , drop = FALSE]
  meta$Group <- meta[[group_var]]
  
  # Ordination
  ord <- ordinate(ps, method = "PCoA", distance = dist_mat)
  
  # PERMANOVA
  ad <- vegan::adonis2(dist_mat ~ Group, data = meta, permutations = 999)
  r2 <- round(ad$R2[1], 3)
  p  <- signif(ad$`Pr(>F)`[1], 3)
  
  # Plot
  p_plot <- plot_ordination(ps, ord, color = group_var) +
    geom_point(size = 2.5, alpha = 0.8) +
    theme_minimal() +
    labs(
      title    = paste("PCoA –", method_title),
      subtitle = paste0("PERMANOVA: R² = ", r2, ", p = ", p),
      color    = group_var
    ) +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold")
    )
  
  list(plot = p_plot, dist = dist_mat, permanova = ad)
}

## ----------------------------
## 3. Generate plots
## ----------------------------
plots <- list()
dist_mats <- list()
permanova_results <- list()

for (m in methods) {
  res <- make_pcoa_plot(ps, m, group_var = group_var)
  plots[[m]] <- res$plot
  dist_mats[[m]] <- res$dist
  permanova_results[[m]] <- res$permanova
}

## ----------------------------
## 4. Combine figure
## ----------------------------
combined_pcoa_plot <- wrap_plots(plots, ncol = 2, guides = "collect") &
  theme(legend.position = "right")

print(combined_pcoa_plot)

# Optional save
# ggsave("beta_diversity_pcoa_combined.png", combined_pcoa_plot, width = 14, height = 10, dpi = 300)

## ----------------------------
## 5. Print PERMANOVA results
## ----------------------------
permanova_results