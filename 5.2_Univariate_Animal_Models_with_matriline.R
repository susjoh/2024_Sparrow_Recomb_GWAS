#
# Animal Models
# SEJ
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

# system("gcta64.exe --bfile clean_data_for_John\\70K_200K_maf_geno_mind_v5 --autosome --autosome-num 29 --make-grm-gz --out clean_data_for_John\\70K_200K_maf_geno_mind_v5.GRM")
# system("gcta64.exe --grm-gz clean_data_for_John\\70K_200K_maf_geno_mind_v5.GRM --grm-adj 0 --make-grm-gz --out clean_data_for_John\\70K_200K_maf_geno_mind_v5.GRM.adj")
# 
# grm.auto <- read.table("clean_data_for_John/70K_200K_maf_geno_mind_v5.GRM.adj.grm.gz")  
# ids.auto <- read.table("clean_data_for_John/70K_200K_maf_geno_mind_v5.GRM.adj.grm.id")  
# 
# grminv <- makeGRM(grm.auto, ids.auto, id.vector = recsumm$parent) 
# dim(grminv)
# 
# attr(grminv, which = "INVERSE") <- TRUE
# 
# save(grminv, file = "clean_data_for_John/70K_200K_maf_geno_mind_v5.GRM.RData")
# rm(grm.auto, ids.auto, grminv)

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

#~~ Which response variables to run?

modeltab <- expand.grid(RespVar = c("yapp_CO_count_no_micro", "intra_shuff_no_micro"),
                        Sex = c("M", "F"),
                        GRM = c("Pedigree"))

#~~ Create master results table

resfixed  <- NULL
resrandom <- NULL
respin    <- NULL
resVp     <- NULL
reswald   <- NULL

#~~ RUN MODELS

for(i in 1:nrow(modeltab)){
  x <- subset(recsumm, sex == modeltab$Sex[i])
  x <- x[,c(as.character(modeltab$RespVar[i]), "parent", "sex", "Total_Coverage", "matriline")]
  x$Total_Coverage2 <- x$Total_Coverage^2
  if(modeltab$GRM[i] == "Pedigree"){
    x <- subset(x, parent %in% attr(ainv, "rowNames"))
  } else {
    x <- subset(x, parent %in% attr(grminv, "rowNames"))
  }
  
  x <- na.omit(x) %>% droplevels
  
  eval(parse(text = paste0("model1 <- asreml(fixed = ", modeltab$RespVar[i], " ~ Total_Coverage + Total_Coverage2,
                      random = ~ vm(parent, ", ifelse(modeltab$GRM[i] == "Pedigree", "ainv", "grminv"),") + ide(parent) + matriline,
                      data = x,
                      na.action = na.method(x = \"omit\", y = \"omit\"),
                      residual = ~idv(units), workspace=\"500mb\")")))
  
  xfixed  <- summary.asreml(model1, coef = T)$coef.fixed %>% data.frame
  xrandom <- summary.asreml(model1, coef = T)$varcomp
  xpin    <- asreml4pin(model1)
  xVp     <- asreml4estVp(model1)
  xwald   <- wald.asreml(model1) %>% data.frame
  
  xfixed$Effect  <- row.names(xfixed)
  xrandom$Effect <- row.names(xrandom)
  xwald$Effect   <- row.names(xwald)
  
  xVp$loglik <- model1$loglik
  xVp$nedf   <- model1$nedf
  xVp$N      <- nrow(x)
  xVp$Nids   <- length(unique(x$id))
  
  x <- data.frame(RespVar = modeltab$RespVar[i], Sex = modeltab$Sex[i], GRM = modeltab$GRM[i])
  
  xfixed  <- cbind(x, xfixed)
  xrandom <- cbind(x, xrandom)
  xpin    <- cbind(x, xpin)
  xVp     <- cbind(x, xVp)
  xwald   <- cbind(x, xwald)
  
  resfixed  <- rbind(resfixed, xfixed)
  resrandom <- rbind(resrandom, xrandom)
  respin    <- rbind(respin, xpin)
  resVp     <- rbind(resVp, xVp)
  reswald   <- rbind(reswald, xwald)
  
  
  rm(model1, xfixed, xrandom, xpin, xVp, x, xwald)
  
}



#~~ Add shuffling with CO count as fixed effect


modeltab <- expand.grid(RespVar = c("intra_shuff_no_micro"),
                        Sex = c("M", "F"),
                        GRM = c("Pedigree"))


for(i in 1:nrow(modeltab)){
  x <- subset(recsumm, sex == modeltab$Sex[i])
  x <- x[,c(as.character(modeltab$RespVar[i]), "parent", "sex", "Total_Coverage", "yapp_CO_count_no_micro", "matriline")]
  
  if(modeltab$GRM[i] == "Pedigree"){
    x <- subset(x, parent %in% attr(ainv, "rowNames"))
  } else {
    x <- subset(x, parent %in% attr(grminv, "rowNames"))
  }
  x$Total_Coverage2 <- x$Total_Coverage^2
  
  x <- na.omit(x) %>% droplevels
  
  eval(parse(text = paste0("model1 <- asreml(fixed = ", modeltab$RespVar[i], " ~ Total_Coverage + Total_Coverage2 + yapp_CO_count_no_micro,
                      random = ~ vm(parent, ", ifelse(modeltab$GRM[i] == "Pedigree", "ainv", "grminv"),") + ide(parent) + matriline,
                      data = x,
                      na.action = na.method(x = \"omit\", y = \"omit\"),
                      residual = ~idv(units), workspace=\"500mb\")")))
  
  xfixed  <- summary.asreml(model1, coef = T)$coef.fixed %>% data.frame
  xrandom <- summary.asreml(model1, coef = T)$varcomp
  xpin    <- asreml4pin(model1)
  xVp     <- asreml4estVp(model1)
  xwald   <- wald.asreml(model1) %>% data.frame
  
  xfixed$Effect  <- row.names(xfixed)
  xrandom$Effect <- row.names(xrandom)
  xwald$Effect   <- row.names(xwald)
  
  xVp$loglik <- model1$loglik
  xVp$nedf   <- model1$nedf
  xVp$N      <- nrow(x)
  xVp$Nids   <- length(unique(x$id))
  
  x <- data.frame(RespVar = paste0(modeltab$RespVar[i], "_cofix"), Sex = modeltab$Sex[i], GRM = modeltab$GRM[i])
  
  xfixed  <- cbind(x, xfixed)
  xrandom <- cbind(x, xrandom)
  xpin    <- cbind(x, xpin)
  xVp     <- cbind(x, xVp)
  xwald   <- cbind(x, xwald)
  
  resfixed  <- rbind(resfixed, xfixed)
  resrandom <- rbind(resrandom, xrandom)
  respin    <- rbind(respin, xpin)
  resVp     <- rbind(resVp, xVp)
  reswald   <- rbind(reswald, xwald)
  
  
  rm(model1, xfixed, xrandom, xpin, xVp, x, xwald)
  
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Parse Animal Model Results          #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

write.table(resfixed,  "results/5_Univar_Fixed_Effects_matriline.txt", row.names = F, sep = "\t", quote = F)
write.table(resrandom, "results/5_Univar_Random_Effects_matriline.txt", row.names = F, sep = "\t", quote = F)
write.table(respin,    "results/5_Univar_Prop_of_Variance_matriline.txt", row.names = F, sep = "\t", quote = F)
write.table(resVp,     "results/5_Univar_Vp_and_Info_matriline.txt", row.names = F, sep = "\t", quote = F)
write.table(reswald,   "results/5_Univar_Wald_Fixed_Effects_matriline.txt", row.names = F, sep = "\t", quote = F)

#~~ Quick visualisation

respin$ModelName <- paste(respin$RespVar, respin$Sex, respin$GRM, sep = "_")

ggplot(respin, aes(ModelName, Estimate, fill = Effect)) +
  geom_bar(stat = "identity", colour = "grey50") +
  theme_bw() +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  theme(axis.text.x  = element_text (size = 8, angle = 270),
        axis.text.y  = element_text (size = 11),
        strip.text.x = element_text (size = 11),
        strip.text.y = element_text (size = 11),
        axis.title.y = element_text (size = 11),
        axis.title.x = element_text (size = 11))


