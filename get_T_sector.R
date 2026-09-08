#### Script to extract proteomics data for the T sector
#### Author: Diana Szeliova

library(dplyr)

data_dir <- "~/offset/proteomics_data"

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
mori_metadata1 <- read.csv(file.path(data_dir, "EV3-Samples-2.csv"))
c_lim_samples <- mori_metadata1[8:22, ]

mori_proteomics1 <- read.csv(file.path(data_dir, "EV9-AbsoluteMassFractions-2.csv"))
mori_proteomics1 <- mori_proteomics1[mori_proteomics1$Gene.name != "", ]
rownames(mori_proteomics1) <- mori_proteomics1$Gene.name

# read.csv converts "-" to "." in column names, so sample IDs from the
# metadata file need the same conversion before they're used to subset
# proteomics columns 
sample_ids1 <- gsub("-", ".", c_lim_samples$Sample.ID)

growth_rates_mori1 <- summarise_dataset(
  mori_proteomics1, c_lim_samples$Growth.rate..1.h., sample_ids1, "mori2021"
)
growth_rates_mori1$id <- c_lim_samples$Sample.ID

# -----------------------------------------------------------------------
# Mori 2021 - minimal media samples (EV2 / EV8)
# -----------------------------------------------------------------------
mori_metadata2 <- read.csv(file.path(data_dir, "EV2-Samples-1.csv"))
target_rows <- c(
  "Acetate", "Carbon 46", "Carbon 61", "Carbon 37",
  "Carbon 85", "Carbon 50", "Carbon 54", "Carbon 49"
)
min_media_samples <- mori_metadata2[mori_metadata2$Short.Description %in% target_rows, ]
mori_mu2 <- log(2) / (min_media_samples$Doubling.time..min. / 60)
sample_ids2 <- gsub("-", ".", min_media_samples$Sample.ID)

mori_proteomics2 <- read.csv(file.path(data_dir, "EV8-AbsoluteMassFractions-1.csv"))
mori_proteomics2 <- mori_proteomics2[mori_proteomics2$Gene.name != "", ]
rownames(mori_proteomics2) <- mori_proteomics2$Gene.name

growth_rates_mori2 <- summarise_dataset(
  mori_proteomics2, mori_mu2, sample_ids2, "mori2021"
)
growth_rates_mori2$id <- min_media_samples$Short.Description

# -----------------------------------------------------------------------
# Wu 2023 - steady-state samples only (M1-M4, N5-N8, P1, P8-P10)
# -----------------------------------------------------------------------
wu_metadata <- read.csv(file.path(data_dir, "wu2023_metadata.csv"), row.names = 1)
wu_metadata <- wu_metadata[!is.na(wu_metadata$Growth.rate..1.h.), ]
steady_state_ids <- c(
  "M1", "M2", "M3", "M4", "N5", "N6", "N7", "N8",
  "P1", "P8", "P9", "P10"
)
wu_metadata <- wu_metadata[steady_state_ids, ]

wu_proteomics1 <- read.csv(file.path(data_dir, "wu2023_group1.csv"))
wu_proteomics2 <- read.csv(file.path(data_dir, "wu2023_group2.csv"))
wu_proteomics3 <- read.csv(file.path(data_dir, "wu2023_group3.csv"))

# unlike the other datasets, missing Gene.name here is coded as 0, not ""
drop_no_locus <- function(df) df[!is.na(df$Gene.name) & df$Gene.name != 0, ]
wu_list <- lapply(list(wu_proteomics1, wu_proteomics2, wu_proteomics3), drop_no_locus)

wu_proteomics <- Reduce(function(x, y) merge(x, y, by = "Gene.name", all = TRUE), wu_list)
rownames(wu_proteomics) <- wu_proteomics$Gene.name

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
