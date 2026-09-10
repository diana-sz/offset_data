#### Shared raw-data loading and normalization for proteomics datasets
#### used by phi_Rb_plot.R and phi_T_extraction.R
#### Author: Diana Szeliova

data_dir <- "~/offset/proteomics_data"

# normalize a proteomics table so each numeric (sample) column sums to 1
normalize_to_sum1 <- function(df) {
  if (is.matrix(df)) {
    return(sweep(df, 2, colSums(df, na.rm = TRUE), "/"))
  }
  numeric_cols <- sapply(df, is.numeric)
  df[numeric_cols] <- sweep(df[numeric_cols], 2, colSums(df[numeric_cols], na.rm = TRUE), "/")
  return(df)
}

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
wu_proteomics <- normalize_to_sum1(wu_proteomics)

# -----------------------------------------------------------------------
# Mori 2021 - carbon-limited chemostat samples (EV3 / EV9)
# -----------------------------------------------------------------------
mori_metadata1 <- read.csv(file.path(data_dir, "EV3-Samples-2.csv"))
c_lim_samples <- mori_metadata1[8:22, ]
# read.csv converts "-" to "." in column names, so sample IDs from the
# metadata file need the same conversion before they're used to subset
# proteomics columns
sample_ids1 <- gsub("-", ".", c_lim_samples$Sample.ID)

mori_proteomics1 <- read.csv(file.path(data_dir, "EV9-AbsoluteMassFractions-2.csv"))
mori_proteomics1 <- mori_proteomics1[mori_proteomics1$Gene.name != "", ]
rownames(mori_proteomics1) <- mori_proteomics1$Gene.name
mori_proteomics1 <- normalize_to_sum1(mori_proteomics1)

# -----------------------------------------------------------------------
# Mori 2021 - minimal media samples (EV2 / EV8)
# -----------------------------------------------------------------------
mori_metadata2 <- read.csv(file.path(data_dir, "EV2-Samples-1.csv"))
target_rows <- c(
  "Acetate", "Carbon 46", "Carbon 61", "Carbon 37",
  "Carbon 85", "Carbon 50", "Carbon 54", "Carbon 49"
)
min_media_samples <- mori_metadata2[mori_metadata2$Short.Description %in% target_rows, ]
mori_lambda2 <- log(2) / (min_media_samples$Doubling.time..min. / 60)
sample_ids2 <- gsub("-", ".", min_media_samples$Sample.ID)

mori_proteomics2 <- read.csv(file.path(data_dir, "EV8-AbsoluteMassFractions-1.csv"))
mori_proteomics2 <- mori_proteomics2[mori_proteomics2$Gene.name != "", ]
rownames(mori_proteomics2) <- mori_proteomics2$Gene.name
mori_proteomics2 <- normalize_to_sum1(mori_proteomics2)
