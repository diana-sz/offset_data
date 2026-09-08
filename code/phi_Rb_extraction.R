#### Script to merge data for ribosomal proteome mass fraction (phi_Rb)
#### Author: Diana Szeliova

library(RColorBrewer)
library(scales)
require(dplyr)

source("~/offset/code/read_proteomics.R")

get_sector_fraction_gene_list <- function(target_genes, proteomics){
  target_rows <- proteomics[target_genes, ]
  fractions <- colSums(target_rows)
  names(fractions) <- colnames(proteomics)
  return(fractions)
}

ribo_genes <- c("rplA","rplB","rplC","rplD","rplE","rplF","rplI","rplJ","rplK","rplL",
              "rplM","rplN","rplO","rplP","rplQ","rplR","rplS","rplT","rplU","rplV",
              "rplW","rplX","rplY","rpmA","rpmB","rpmC","rpmD","rpmE","rpmF","rpmG",
              "rpmH","rpmI","rpmJ","rpsA","rpsB","rpsC","rpsD","rpsE","rpsF","rpsG",
              "rpsH","rpsI","rpsJ","rpsK","rpsL","rpsM","rpsN","rpsO","rpsP","rpsQ",
              "rpsR","rpsS","rpsT","rpsU","sra")


# Chure 2023 data ##############################################################
# Note: one data point from Bremer 2008 (LB condition) was added
chure <- read.csv(file.path(data_dir, "Chure_2023.csv"), skip=8)

# Wu 2023 data #################################################################
wu_proteomics_ss <- wu_proteomics[, steady_state_ids]
ribosomes_wu <- get_sector_fraction_gene_list(ribo_genes, wu_proteomics_ss)

wu <- data.frame(mu = wu_metadata$Growth.rate..1.h., phi_Rb = ribosomes_wu,
                name = "Wu et al. 2023", source = "Wu et al. 2023")

#### Zhu data ##################################################################
zhu <- data.frame(mu=c(1.88,1.26,0.97,0.69,0.41), 
                  phi_Rb=c(0.484,0.364,0.294,0.227,0.172)*0.5, # convert RP ratio to phi_Rb
                  name = "Zhu et al. 2025", source = "Zhu et al. 2025")


#### Mori data - carbon-limited chemostat #####################################
mori_proteomics1_ss <- mori_proteomics1[, sample_ids1]
ribosomes_mori1 <- get_sector_fraction_gene_list(ribo_genes, mori_proteomics1_ss)

mori1 <- data.frame(mu = c_lim_samples$Growth.rate..1.h., phi_Rb = ribosomes_mori1,
                    name = "Mori et al. 2021", source = "Mori et al. 2021")

#### Mori data - minimal media #################################################
mori_proteomics2_ss <- mori_proteomics2[, sample_ids2]
ribosomes_mori2 <- get_sector_fraction_gene_list(ribo_genes, mori_proteomics2_ss)

mori2 <- data.frame(mu = mori_mu2, phi_Rb = ribosomes_mori2,
                    name = "Mori et al. 2021", source = "Mori et al. 2021")


chure$name <- "Chure & Cremer 2023"
colnames(chure)[colnames(chure) == "growth_rate_hr"] <- "mu"
colnames(chure)[colnames(chure) == "mass_fraction"]  <- "phi_Rb"


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


#### Plotting ##################################################################
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

print(lm(phi_Rb~mu, data=all_data))
