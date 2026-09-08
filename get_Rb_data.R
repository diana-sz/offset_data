#### Script to merge data for ribosomal proteome mass fraction (phi )
#### Author: Diana Szeliova

library(RColorBrewer)
library(scales)
require(dplyr)

data_dir <- "~/offset/proteomics_data"

# ribosomal genes from EcoCyc 
ribo_genes <- c("b4202", "b0169", "b2606", "b3407", "b3165", "b3301", "b3299", 
                "b3342", "b1089", "b3986", "b3230", "b3312", "b3305", "b4200", 
                "b3309", "b3320", "b3186", "b3316", "b3314", "b3983", "b2609", 
                "b3294", "b3298", "b3636", "b3231", "b3321", "b3302", "b4203",
                "b3341", "b2185", "b3319", "b3185", "b0023", "b3296", "b3315",
                "b3984", "b4506", "b3311", "b3304", "b3307", "b3703", "b3310", 
                "b3297", "b3936", "b3985", "b3306", "b3637", "b3308", "b3313",
                "b3065", "b3303", "b3318", "b3317", "b0296", "b1717", "b0911", 
                "b1716")

chure <- read.csv(file.path(data_dir, "Chure_2023.csv"), skip=8)

# data digitized from their figure
wu <- read.csv(file.path(data_dir, "Wu_2023_ext_fig_1.csv"))

zhu <- data.frame(mu=c(1.88,1.26,0.97,0.69,0.41), 
                  phi_Rb=c(0.484,0.364,0.294,0.227,0.172)*0.5, # convert RP ratio to phi_Rb
                  name = "Zhu et al. 2025", source = "Zhu et al. 2025")


#### Mori data #################################################################
metadata <- read.csv(file.path(data_dir, "EV3-Samples-2.csv"))
c_lim_samples <- metadata[8:22,]
mori_growth_rates <- c_lim_samples$Growth.rate..1.h.

mori_proteomics <- read.csv(file.path(data_dir, "EV9-AbsoluteMassFractions-2.csv"))
mori_proteomics <- mori_proteomics[mori_proteomics$Gene.locus != "", ]
rownames(mori_proteomics) <- mori_proteomics$Gene.locus
mori_proteomics <- mori_proteomics[, c_lim_samples$Sample.ID]


#### Mori data #################################################################
metadata <- read.csv(file.path(data_dir, "EV2-Samples-1.csv"))

mori_proteomics_b <- read.csv(file.path(data_dir, "EV8-AbsoluteMassFractions-1.csv"))
mori_proteomics_b <- mori_proteomics_b[mori_proteomics_b$Gene.locus != "", ]
rownames(mori_proteomics_b) <- mori_proteomics_b$Gene.locus

target_rows <- c("Acetate", "Carbon 46", "Carbon 61", "Carbon 37",
                 "Carbon 85", "Carbon 50", "Carbon 54", "Carbon 49")             
min_media_samples <- metadata[metadata$Short.Description %in% target_rows,]
mori_growth_rates2 <- log(2)/(min_media_samples$Doubling.time..min./60)
sample_ids <- gsub("-", ".", min_media_samples$Sample.ID)
mori_proteomics2 <- mori_proteomics_b[, sample_ids]


get_sector_fraction_gene_list <- function(target_genes, proteomics){
  target_rows <- proteomics[target_genes, ]
  target_sums <- colSums(target_rows, na.rm = TRUE)
  fractions <- target_sums/colSums(proteomics, na.rm=TRUE)
  names(fractions) <- colnames(proteomics)
  return(fractions)
}


ribosomes_mori1 <- get_sector_fraction_gene_list(ribo_genes, mori_proteomics)
ribosomes_mori2 <- get_sector_fraction_gene_list(ribo_genes, mori_proteomics2)

mori1 <- data.frame(mu = mori_growth_rates, phi_Rb = ribosomes_mori1,
                    name = "Mori et al. 2021", source = "Mori et al. 2021")
mori2 <- data.frame(mu = mori_growth_rates2, phi_Rb = ribosomes_mori2,
                    name = "Mori et al. 2021", source = "Mori et al. 2021")


chure$name <- "Chure & Cremer 2023"
colnames(chure)[colnames(chure) == "growth_rate_hr"] <- "mu"
colnames(chure)[colnames(chure) == "mass_fraction"]  <- "phi_Rb"

wu$phi_Rb <- wu$phi_R/100
wu$source <- "Wu et al. 2023"
wu$name <- "Wu et al. 2023"

# combine all data
## Keep only required columns before binding
keep_cols <- c("mu", "phi_Rb", "name", "source")

chure  <- chure[,  keep_cols]
wu     <- wu[,     keep_cols]
mori1  <- mori1[,  keep_cols]
mori2  <- mori2[,  keep_cols]


all_data <- rbind(chure, wu, mori1, mori2, zhu) 

# exlude Si
all_data <- all_data[all_data$source != "Si et al., 2017", ]


write.csv(all_data, "~/offset/phi_Rb_data.csv")

# define colors and shapes
colors <- brewer.pal(length(unique(all_data$name)), "Dark2")
color_map <- setNames(colors, unique(all_data$name))
all_data$color <- color_map[all_data$name]

shape_map <- setNames(
  seq_along(unique(all_data$source))+1,
  unique(all_data$source)
)
all_data$shape <- shape_map[all_data$source]

plot(phi_Rb~mu, data=all_data,
     col = all_data$color,
     pch = all_data$shape,
     xlab = NA, ylab = NA, axes = FALSE,
     ylim = c(0,0.28),
     xlim = c(0, 2.2),
     xaxs = "i", yaxs = "i",
     cex=1)

axis(1, padj = -0.25)
axis(2, las = 2, at = seq(0,0.25,0.05), labels = c("0%", "5%", "10%", "15%", "20%", "25%"), hadj = 0.75)
box()

mtext(bquote("Growth rate" ~ lambda ~ "[" * h^-1 * "]"), 
      side = 1, outer = FALSE, line = 2.5, cex = 1.1)
mtext(expression("Ribosome proteome fraction " * Phi[R] ),
      side = 2, outer = FALSE, line = 2.3, cex = 1.1)

leg <- c("Chure & Cremer 2023", "Wu et al. 2023", "Mori et al. 2021", "Zhu et al. 2025")
shape_map["Chure & Cremer 2023"] <- 18
legend("topleft", legend = leg, col = color_map[leg], pch = shape_map[leg], cex = 0.8, pt.cex = 1.1)

leg2 <- sort(unique(chure$source))
legend("bottomright", legend = leg2, 
       col = color_map["Chure & Cremer 2023"], 
       pch = shape_map[leg2], cex = 0.5, ncol = 1,
       title = "Sources in Chure & Cremer 2023")



