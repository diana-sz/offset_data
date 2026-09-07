#### Script to extract proteomics data for the R sector
#### Author: Diana Szeliova
#### Date: 1.9.2026

library(dplyr)

data_dir <- "~/offset/proteomics_data"

# -----------------------------------------------------------------------
# Gene annotation (Wu 2023, Supplementary Table 3)
# -----------------------------------------------------------------------
affiliated_trans <- c(
  "arfA", "arfB", "efp", "frr", "fusA", "infA", "infB", "infC", "lepA", "prfA",
  "prfB", "prfC", "tsf", "tufA", "tufB"
)

eftu <- c("tufA", "tufB") # subset of affiliated_trans

trna_syn <- c(
  "alaS", "argS", "asnS", "aspS", "cysS", "glnS", "gltX", "glyQ", "glyS", "hisS",
  "ileS", "leuS", "lysS", "lysU", "metG", "pheS", "pheT", "proS", "serS", "thrS",
  "trpS", "tyrS", "valS"
)

ribosome <- c("rplA","rplB","rplC","rplD","rplE","rplF","rplI","rplJ","rplK","rplL",
              "rplM","rplN","rplO","rplP","rplQ","rplR","rplS","rplT","rplU","rplV",
              "rplW","rplX","rplY","rpmA","rpmB","rpmC","rpmD","rpmE","rpmF","rpmG",
              "rpmH","rpmI","rpmJ","rpsA","rpsB","rpsC","rpsD","rpsE","rpsF","rpsG",
              "rpsH","rpsI","rpsJ","rpsK","rpsL","rpsM","rpsN","rpsO","rpsP","rpsQ",
              "rpsR","rpsS","rpsT","rpsU","sra")

# -----------------------------------------------------------------------
# Helper: sum tRNA synthetase, EF-Tu, affiliated translational protein, and
# ribosomal protein abundances per sample, and combine with growth rates
# into one data frame. Two proteome-sector definitions are computed:
#   sum_all_klumpp - tRNA synthetases + EF-Tu (= tufA + tufB) + tsf
#   sum_all_wu     - tRNA synthetases + all affiliated translational proteins
# sum_rb is the ribosomal protein sum, returned separately.
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
    sum_eftu = eftu,
    sum_affil = affiliated_trans,
    sum_rb = ribosome
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
  affil_rows <- proteomics[rownames(proteomics) %in% affiliated_trans, sample_ids, drop = FALSE]
  rb_rows <- proteomics[rownames(proteomics) %in% ribosome, sample_ids, drop = FALSE]
  
  data.frame(
    mu = mu,
    id = sample_ids,
    dataset = dataset_name,
    sum_trnas = colSums(trna_rows, na.rm = TRUE),
    sum_eftu = colSums(eftu_rows, na.rm = TRUE),
    sum_affil = colSums(affil_rows, na.rm = TRUE),
    sum_rb = colSums(rb_rows, na.rm = TRUE)
  ) %>%
    mutate(
      sum_all_klumpp = sum_trnas + sum_eftu,
      sum_all_wu = sum_trnas + sum_affil
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

plot_by_dataset(growth_rates$sum_all_klumpp, "tRNAsyn + EF-Tu")
plot_by_dataset(growth_rates$sum_all_wu, "Wu definition (tRNAsyn + affil. trans.)")
plot_by_dataset(growth_rates$sum_trnas, "tRNA synthetases")
plot_by_dataset(growth_rates$sum_eftu, "EF-Tu + tsf")
plot_by_dataset(growth_rates$sum_rb, "Ribosomal proteins")

# -----------------------------------------------------------------------
# Save results
# -----------------------------------------------------------------------
write.csv(growth_rates, file.path(data_dir, "R_data.csv"), row.names = FALSE)
