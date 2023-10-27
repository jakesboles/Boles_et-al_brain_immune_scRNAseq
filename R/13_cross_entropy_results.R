library(tidyverse)
library(rlist)
library(Seurat)
library(gt)

classes <- c("blymph", "cd45neg", "dendritic_cells", "granulocytes", 
             "macrophages", "microglia", "monocytes", "natural_killer", "tcells")

tbl <- list()

for (i in seq_along(classes)){
  tbl[[i]] <- read.table(paste0("cross_entropy_test/02_tsne_diff/", classes[i], "/umap_ce_diff_result_group.txt"),
                         sep = ",", skip = 6, nrows = 1)
}

tbl <- list.rbind(tbl)

tbl <- tbl %>%
  mutate(V1 = str_remove_all(V1, "D = ") %>%
           as.numeric(),
         V2 = str_remove_all(V2, "p-value = ") %>%
           str_remove_all("p-value < ") %>%
           as.numeric()) %>%
  dplyr::rename("D" = "V1", "p" = "V2") %>%
  mutate(cell = classes)

(gt_tbl <- tbl %>%
    mutate(p = round(p, digits = 3)) %>%
    mutate(p = if_else(p < 0.001, "< 0.001",
                       paste0(p))) %>%
    mutate(cell = case_when(cell == "blymph" ~ "B-lymphocytes",
                            cell == "cd45neg" ~ "CD45- cells",
                            cell == "dendritic_cells" ~ "Dendritic cells",
                            cell == "granulocytes" ~ "Granulocytes",
                            cell == "macrophages" ~ "Macrophages",
                            cell == "microglia" ~ "Microglia",
                            cell == "monocytes" ~ "Monocytes",
                            cell == "natural_killer" ~ "Natural killer cells",
                            cell == "tcells" ~ "T-cells")) %>%
    arrange(cell) %>%
    gt(groupname_col = "cell",
       row_group_as_column = T,
       process_md = T) %>%
    fmt_number(decimals = 3) %>%
    cols_label(
      D = "Kullback-Leibler divergence",
      p = "Holm-adjusted {{*p*-value}}"
    ) %>%
    tab_style(
      style = list(
        cell_text(align = "center",
                  v_align = "middle")),
      locations = list(
        cells_body(),
        cells_column_labels(),
        cells_row_groups(),
        cells_stub())) %>%
    tab_style(
      style = list(
        cell_text(weight = "bold"),
        cell_fill(color = "red",
                  alpha = 0.5)),
      locations = cells_body(
        columns = c(D, p),
        rows = p < 0.05)
    ) 
)

gtsave(gt_tbl,
       filename = "tabular_output/cross_entropy_results.html")