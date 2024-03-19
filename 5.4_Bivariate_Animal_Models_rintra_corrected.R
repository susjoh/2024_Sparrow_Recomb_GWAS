#
# Animal Models - Bivariate
# EW & SEJ
# Feb 2022
#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 0. Set up Workspace              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

library(dplyr)
library(asreml)
library(ggplot2)
library(tidyr)
source("r/ASReml4.EstEffects.R")
source("r/makeGRM.R")

pedigree <- read.table("clean_data_for_John/Helgeland_Pedigree_Abridged.txt", header = T)
ainv <- ainverse(pedigree)

recsumm <- read.table("clean_data_for_John/3_Full_Recombination_Phenotypes_QCed.txt", header = T)

load("clean_data_for_John/70K_200K_maf_geno_mind_v5.GRM.RData")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Structure & Run Animal Models       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

str(recsumm)

recsumm_hold <- recsumm
recsumm <- recsumm_hold
recsumm <- subset(recsumm, select = c(parent, meiosis, sex, Mother, matriline, Total_Coverage,
                                      yapp_CO_count_QCed, yapp_CO_count_no_micro, 
                                      yapp_CO_count_macro, intra_shuff, intra_shuff_macro, 
                                      intra_shuff_no_micro, intra_shuff_gene, intra_shuff_gene_macro, 
                                      intra_shuff_gene_no_micro))

recsumm$parent <- factor(recsumm$parent)
recsumm$Mother <- factor(recsumm$Mother)
recsumm$matriline <- factor(recsumm$matriline)
recsumm$sex <- factor(recsumm$sex)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Structure & Run Animal Models       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~ Which response variables to run?

modeltab <- expand.grid(RespVar = c("intra_shuff_no_micro"),
                        GRM = c("Pedigree", "GRM"))

#~~ Create master results table

resfixed  <- NULL
resrandom <- NULL
resVp     <- NULL
resvparam <- list()

#~~ RUN MODELS

for(i in 1:nrow(modeltab)){
  
  x <- recsumm
  x <- x[,c(as.character(modeltab$RespVar[i]), "parent", "sex", "Total_Coverage", "matriline", "yapp_CO_count_no_micro")]
  x <- na.omit(x) %>% droplevels
  
  
  eval(parse(text = paste0("x$", as.character(modeltab$RespVar[i]), "_f <- ifelse(x$sex == \"F\", recsumm$", as.character(modeltab$RespVar[i]), ", NA)")))
  eval(parse(text = paste0("x$", as.character(modeltab$RespVar[i]), "_m <- ifelse(x$sex == \"M\", recsumm$", as.character(modeltab$RespVar[i]), ", NA)")))
  
  
  x$Total_Coverage2 <- x$Total_Coverage^2
  if(modeltab$GRM[i] == "Pedigree"){
    x <- subset(x, parent %in% attr(ainv, "rowNames"))
  } else {
    x <- subset(x, parent %in% attr(grminv, "rowNames"))
  }
  
  x1 <- as.character(modeltab$RespVar[i])
  
  eval(parse(text = paste0("model1 <- asreml(fixed  = cbind(", x1, "_m, ", x1, "_f) ~ trait + trait:yapp_CO_count_no_micro + trait:Total_Coverage + trait:Total_Coverage2,
                     random = ~ corgh(trait):vm(parent, ", ifelse(modeltab$GRM[i] == "Pedigree", "ainv", "grminv"),") + idh(trait):ide(parent),
                     residual   = ~ units:us(trait, init = NA),
                     data = x,
                     workspace=\"1000mb\",
                     maxit = 20)")))
  
  resvparam[[i]] <- model1$vparameters
  
  # In this case, V1 is the actual correlation!
  
  xfixed  <- summary.asreml(model1, coef = T)$coef.fixed %>% data.frame
  xrandom <- summary.asreml(model1, coef = T)$varcomp
  
  xfixed$Effect  <- row.names(xfixed)
  xrandom$Effect <- row.names(xrandom)
  
  xf <- subset(x, sex == "F")
  xm <- subset(x, sex == "M")
  
  xVp <- data.frame(loglik = model1$loglik,
                    nedf   = model1$nedf,
                    N      = nrow(x),
                    Nids   = length(unique(x$id)),
                    N_f    = nrow(xf),
                    Nids_f = length(unique(xf$id)),
                    N_m    = nrow(xm),
                    Nids_m = length(unique(xm$id)))
  
  x <- data.frame(RespVar = modeltab$RespVar[i],
                  GRM = modeltab$GRM[i])
  
  xfixed  <- cbind(x, xfixed)
  xrandom <- cbind(x, xrandom)
  xVp     <- cbind(x, xVp)
  
  
  resfixed  <- rbind(resfixed, xfixed)
  resrandom <- rbind(resrandom, xrandom)
  resVp     <- rbind(resVp, xVp)
  
  write.table(resfixed,  "results/2_Bivar_Fixed_Effects_rintracor_temp.txt", row.names = F, sep = "\t", quote = F)
  write.table(resrandom, "results/2_Bivar_Random_Effects_rintracor_temp.txt", row.names = F, sep = "\t", quote = F)
  write.table(resVp,     "results/2_Bivar_Info_rintracor_temp.txt", row.names = F, sep = "\t", quote = F)
  
  
  rm(model1, xfixed, xrandom, x, x1, xm, xf, xVp)
  
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Parse Animal Model Results          #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


write.table(resfixed,  "results/2_Bivar_Fixed_Effects_rintracor.txt", row.names = F, sep = "\t", quote = F)
write.table(resrandom, "results/2_Bivar_Random_Effects_rintracor.txt", row.names = F, sep = "\t", quote = F)
write.table(resVp,     "results/2_Bivar_Info_rintracor.txt", row.names = F, sep = "\t", quote = F)
save(resvparam, file = "results/2_Bivar_Vparams_rintracor.RData")

