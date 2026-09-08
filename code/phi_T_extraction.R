#### Script to extract proteomics data for the T sector
#### Author: Diana Szeliova

library(dplyr)

source("~/offset/code/read_proteomics.R")

# -----------------------------------------------------------------------
# Gene annotation (Wu 2023, Supplementary Table 3)
# -----------------------------------------------------------------------
eftu <- c("tufA", "tufB")

trna_syn <- c(
  "alaS", "argS", "asnS", "aspS", "cysS", "glnS", "gltX", "glyQ", "glyS", "hisS",
  "ileS", "leuS", "lysS", "lysU", "metG", "pheS", "pheT", "proS", "serS", "thrS",
  "trpS", "tyrS", "valS"
)


# -----------------------------------------------------------------------
# Helper: sum tRNA synthetases and EF-Tu mass fractions, and combine with growth rates
# into one data frame.
# -----------------------------------------------------------------------
summarise_dataset <- function(proteomics, mu, sample_ids, dataset_name) {
  missing_ids <- setdiff(sample_ids, colnames(proteomics))
  if (length(missing_ids) > 0) {
    stop(
      "summarise_dataset (", dataset_name, "): sample IDs not found in proteomics columns: ",
      paste(missing_ids, collapse = ", ")
    )
  }
  
  gene_sets <- list(
    sum_trnas = trna_syn,
    sum_eftu = eftu
  )
  for (set_name in names(gene_sets)) {
    n_found <- sum(rownames(proteomics) %in% gene_sets[[set_name]])
    if (n_found == 0) {
      warning(
        "summarise_dataset (", dataset_name, "): no genes matched for ", set_name,
        " - check gene naming in this dataset"
      )
    }
  }
  
  trna_rows <- proteomics[rownames(proteomics) %in% trna_syn, sample_ids, drop = FALSE]
  eftu_rows <- proteomics[rownames(proteomics) %in% eftu, sample_ids, drop = FALSE]

  data.frame(
    mu = mu,
    id = sample_ids,
    dataset = dataset_name,
    sum_trnas = colSums(trna_rows, na.rm = TRUE),
    sum_eftu = colSums(eftu_rows, na.rm = TRUE)
  ) %>%
    mutate(
      phi_T = sum_trnas + sum_eftu,
    )
}

# -----------------------------------------------------------------------
# Mori 2021 - carbon-limited chemostat samples (EV3 / EV9)
# -----------------------------------------------------------------------
growth_rates_mori1 <- summarise_dataset(
  mori_proteomics1, c_lim_samples$Growth.rate..1.h., sample_ids1, "mori2021"
)
growth_rates_mori1$id <- c_lim_samples$Sample.ID

# -----------------------------------------------------------------------
# Mori 2021 - minimal media samples (EV2 / EV8)
# -----------------------------------------------------------------------
growth_rates_mori2 <- summarise_dataset(
  mori_proteomics2, mori_mu2, sample_ids2, "mori2021"
)
growth_rates_mori2$id <- min_media_samples$Short.Description

# -----------------------------------------------------------------------
# Wu 2023
# -----------------------------------------------------------------------
growth_rates_wu <- summarise_dataset(
  wu_proteomics, wu_metadata$Growth.rate..1.h., rownames(wu_metadata), "wu2023"
)


# -----------------------------------------------------------------------
# Combine all datasets
# -----------------------------------------------------------------------
growth_rates <- rbind(
  growth_rates_mori1, growth_rates_mori2, growth_rates_wu
)

# -----------------------------------------------------------------------
# Plots - points colored by dataset
# -----------------------------------------------------------------------
dataset_colors <- c(mori2021 = "#7570b3", wu2023 = "#d95f02")

point_colors <- dataset_colors[growth_rates$dataset]

plot_by_dataset <- function(y, main) {
  plot(growth_rates$mu, y, col = point_colors, pch = 16,
       ylim = c(0, max(y, na.rm = TRUE) * 1.05),
       xlab = "growth rate (1/h)", ylab = "mass fraction", main = main)
  legend("topleft", legend = names(dataset_colors), col = dataset_colors, pch = 16)
}

plot_by_dataset(growth_rates$sum_trnas, "tRNA synthetases")
plot_by_dataset(growth_rates$sum_eftu, "EF-Tu")
plot_by_dataset(growth_rates$phi_T, "T sector")

# -----------------------------------------------------------------------
# Save results
# -----------------------------------------------------------------------
write.csv(growth_rates, "~/offset/phi_T_data.csv", row.names = FALSE)
