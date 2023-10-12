#This script employs the chooseR iterative subsampling procedure on the myeloid cluster
#chooseR-identified clusters will be annotated in a later script

source("chooseR_helpers/pipeline.R")
library(Seurat)
library(ggplot2)
`%>%` <- magrittr::`%>%`

obj <- readRDS("data_objects/05_clean_monomacs_sub.RDS")

dir.create("chooseR/monomacs")

npcs <- 30
resolutions <- c(0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 2, 4, 6, 8, 12)
assay <- "SCT"
reduction <- "pca"
results_path <- "chooseR/monomacs/"

# Run pipeline
for (res in resolutions) {
  message(paste0("Clustering ", res, "..."))
  message("\tFinding ground truth...")
  
  # "Truths" will be stored at glue::glue("{reduction}.{assay}_res.{res}")
  obj <- find_clusters(
    obj,
    reduction = reduction,
    assay = assay,
    resolution = res
  )
  clusters <- obj[[glue::glue("{reduction}.{assay}_res.{res}")]]
  
  # Now perform iterative, sub-sampled clusters
  results <- multiple_cluster(
    obj,
    n = 100,
    size = 0.8,
    npcs = npcs,
    res = res,
    reduction = reduction,
    assay = assay
  )
  
  # Now calculate the co-clustering frequencies
  message(paste0("Tallying ", res, "..."))
  # This is the more time efficient vectorisation
  # However, it exhausts vector memory for (nearly) all datasets
  # matches <- purrr::map(columns, find_matches, df = results)
  # matches <- purrr::reduce(matches, `+`)
  columns <- colnames(dplyr::select(results, -cell))
  mtchs <- matrix(0, nrow = dim(results)[1], ncol = dim(results)[1])
  i <- 1 # Counter
  for (col in columns) {
    message(paste0("\tRound ", i, "..."))
    mtchs <- Reduce("+", list(
      mtchs,
      find_matches(col, df = results)
    ))
    i <- i + 1
  }
  
  message(paste0("Scoring ", res, "..."))
  mtchs <- dplyr::mutate_all(
    dplyr::as_tibble(mtchs),
    function(x) dplyr::if_else(Re(x) > 0, percent_match(x), 0)
  )
  
  # Now calculate silhouette scores
  message(paste0("Silhouette ", res, "..."))
  sil <- cluster::silhouette(
    x = as.numeric(as.character(unlist(clusters))),
    dmatrix = (1 - as.matrix(mtchs))
  )
  saveRDS(sil, paste0(results_path, "silhouette_", res, ".rds"))
  
  # Finally, calculate grouped metrics
  message(paste0("Grouping ", res, "..."))
  grp <- group_scores(mtchs, unlist(clusters))
  saveRDS(grp, paste0(results_path, "frequency_grouped_", res, ".rds"))
  sil <- group_sil(sil, res)
  saveRDS(sil, paste0(results_path, "silhouette_grouped_", res, ".rds"))
}

# Save original data, with ground truth labels
saveRDS(obj, paste0(results_path, "clustered_data.rds"))

# Create silhouette plot
# Read in scores and calculate CIs
scores <- purrr::map(
  paste0(results_path, "silhouette_grouped_", resolutions, ".rds"),
  readRDS
)
scores <- dplyr::bind_rows(scores) %>%
  dplyr::group_by(res) %>%
  dplyr::mutate("n_clusters" = dplyr::n()) %>%
  dplyr::ungroup()
meds <- scores %>%
  dplyr::group_by(res) %>%
  dplyr::summarise(
    "boot" = list(boot_median(avg_sil)),
    "n_clusters" = mean(n_clusters)
  ) %>%
  tidyr::unnest_wider(boot)

writexl::write_xlsx(meds, paste0(results_path, "median_ci.xlsx"))

# Find thresholds
threshold <- max(meds$low_med)
choice <- as.character(
  meds %>%
    dplyr::filter(med >= threshold) %>%
    dplyr::arrange(n_clusters) %>%
    tail(n = 1) %>%
    dplyr::pull(res)
)

# And plot!
ggplot(meds, aes(factor(res), med)) +
  geom_crossbar(
    aes(ymin = low_med, ymax = high_med),
    fill = "grey",
    size = 0.25
  ) +
  #geom_hline(aes(yintercept = threshold), colour = "blue") +
  #geom_vline(aes(xintercept = choice), colour = "red") +
  geom_jitter(
    data = scores,
    aes(factor(res), avg_sil),
    size = 0.35,
    width = 0.15
  ) +
  scale_x_discrete("Resolution") +
  scale_y_continuous(
    "Silhouette Score",
    expand = c(0, 0),
    limits = c(-1, 1),
    breaks = seq(-1, 1, 0.25),
    oob = scales::squish
  ) +
  cowplot::theme_minimal_hgrid() +
  theme(
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7),
    axis.line.x = element_line(colour = "black"),
    axis.line.y = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black"),
  )

ggsave(
  filename = paste0(results_path, "silhouette_distribution_plot.png"),
  dpi = 300,
  height = 3.5,
  width = 3.5,
  units = "in"
)
