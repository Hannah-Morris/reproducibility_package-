# Script: 03_alpha_diversity.R
# Purpose: Calculate alpha diversity metrics across host species and generate combined violin plots
# Input: phy_all (merged phyloseq object with Organism in sample_data)
# Output: Shannon, Simpson, Chao1, and Observed diversity plots combined into a 2x2 figure

library(phyloseq)
library(dplyr)
library(tibble)
library(ggpubr)
library(rstatix)
library(patchwork)

## ----------------------------
## 1. Estimate alpha diversity
## ----------------------------
alpha_diversity <- estimate_richness(phy_all)

alpha_diversity$Sample <- rownames(alpha_diversity)

## ----------------------------
## 2. Build metadata + statframe
## ----------------------------
metadata <- sample_data(phy_all) %>%
  data.frame(stringsAsFactors = FALSE)
metadata$Sample <- rownames(metadata)

statframe <- alpha_diversity %>%
  dplyr::select(Sample, Observed, Chao1, Simpson, Shannon) %>%
  left_join(
    metadata %>% dplyr::select(Sample, Organism),
    by = "Sample"
  ) %>%
  rename(
    organisms    = Organism,
    observed     = Observed,
    chaoscore    = Chao1,
    simpsonscore = Simpson,
    shannonscore = Shannon
  )

## ----------------------------
## 3. Helper: choose global test + pairwise test
## ----------------------------
run_alpha_tests <- function(df, response, group_var = "organisms") {
  response_vec <- df[[response]]
  shapiro_p <- shapiro_test(response_vec)$p.value
  
  formula_txt <- paste(response, "~", group_var)
  fml <- as.formula(formula_txt)
  
  if (shapiro_p < 0.05) {
    global_res <- df %>% kruskal_test(fml)
    pairwise_res <- df %>% dunn_test(fml, p.adjust.method = "bonferroni")
  } else {
    global_res <- df %>% anova_test(fml)
    pairwise_res <- df %>% tukey_hsd(fml)
  }
  
  list(
    shapiro_p = shapiro_p,
    global_res = global_res,
    pairwise_res = pairwise_res
  )
}

## ----------------------------
## 4. Run tests for each metric
## ----------------------------
res_shannon  <- run_alpha_tests(statframe, "shannonscore")
res_simpson  <- run_alpha_tests(statframe, "simpsonscore")
res_chao     <- run_alpha_tests(statframe, "chaoscore")
res_observed <- run_alpha_tests(statframe, "observed")

## Add y-positions for p-value annotations
pwc_shannon  <- res_shannon$pairwise_res  %>% add_xy_position(x = "organisms")
pwc_simpson  <- res_simpson$pairwise_res  %>% add_xy_position(x = "organisms")
pwc_chao     <- res_chao$pairwise_res     %>% add_xy_position(x = "organisms")
pwc_observed <- res_observed$pairwise_res %>% add_xy_position(x = "organisms")

## ----------------------------
## 5. Helper: build alpha diversity plot
## ----------------------------
make_alpha_plot <- function(df, yvar, title_txt, pwc_res, global_res) {
  p <- ggviolin(
    df,
    x = "organisms",
    y = yvar,
    title = title_txt,
    fill = "organisms",
    palette = "Paired"
  ) +
    stat_pvalue_manual(
      pwc_res,
      hide.ns = TRUE,
      step.increase = 0.1
    ) +
    labs(caption = get_pwc_label(pwc_res)) +
    theme(
      plot.title   = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title   = element_text(size = 14, face = "bold"),
      axis.text.x  = element_text(size = 0),
      axis.text.y  = element_text(size = 12),
      legend.text  = element_text(size = 12),
      legend.title = element_text(size = 12, face = "bold"),
      legend.position = "right"
    )
  
  p <- p + labs(subtitle = get_test_label(global_res, detailed = FALSE))
  return(p)
}

## ----------------------------
## 6. Generate plots
## ----------------------------
statplot_shannon <- make_alpha_plot(
  statframe, "shannonscore", "Shannon Diversity of Species",
  pwc_shannon, res_shannon$global_res
)

statplot_simpson <- make_alpha_plot(
  statframe, "simpsonscore", "Simpson Diversity of Species",
  pwc_simpson, res_simpson$global_res
)

statplot_chao <- make_alpha_plot(
  statframe, "chaoscore", "Chao1 Diversity of Species",
  pwc_chao, res_chao$global_res
)

statplot_observed <- make_alpha_plot(
  statframe, "observed", "Observed Diversity of Species",
  pwc_observed, res_observed$global_res
)

## ----------------------------
## 7. Combined figure
## ----------------------------
alpha_div_plot <- (statplot_shannon | statplot_simpson) /
  (statplot_chao | statplot_observed) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

print(alpha_div_plot)

# Optional save
# ggsave("alpha_diversity_combined.png", alpha_div_plot, width = 14, height = 10, dpi = 300)