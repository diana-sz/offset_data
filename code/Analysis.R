require('rstudioapi') 

# Sets working directory to source file location ###############################

directory <- dirname(getActiveDocumentContext()$path)

setwd(directory) 

# reads Rb protein data ########################################################

Rbdata   <- read.table("../Rb_protein_data.csv", sep = ',', header = TRUE)

lambda_Rb  <- as.numeric(Rbdata[-1,1])

phi_Rb <- as.numeric(Rbdata[-1,2])

reg_Rb <- lm(phi_Rb ~ lambda_Rb)

plot(lambda_Rb, phi_Rb)
abline(reg_Rb)

summary(reg_Rb)

# reads T protein data #########################################################

phitdata <- as.matrix(read.table("../T_protein_data.csv", sep = ',', header = TRUE)) 

lambda <- as.numeric(phitdata[-1,1])

phit   <- as.numeric(phitdata[-1,6])

reg_T <- lm(phit ~ lambda)

reg_T2 <-lm(phit[lambda>0.5] ~ lambda[lambda>0.5])

plot(lambda, phit)
abline(lm(phit ~ lambda))


summary(reg_T)



