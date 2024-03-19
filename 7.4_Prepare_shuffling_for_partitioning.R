

library(dplyr)
library(asreml)

load("reg_h2/6_Data_for_reg_h2_TRANS.RData")

#~~ Deal with the shuffling

shufftab <- read.table("results/3_Full_Shuffling_Summary_by_Chr.txt", header = T)

shufftab_w <- subset(shufftab, chrom == h)
shufftab_w$intra_shuff_no_micro <- shufftab_w$IS_no_micro
shufftab_wo <- subset(shufftab, chrom != h)


shuffsumm <- shufftab_wo %>%
  group_by(meiosis) %>%
  summarise(intra_shuff_no_micro = sum(IS_no_micro))


snptab <- read.table("70k_data/70K_200K_maf_geno_mind_v5.bim")

modeltab <- subset(modeltab, RespVar == "intra_shuff_no_micro")

writeLines(paste(1, as.character(unique(recsumm$parent))), "chr_h2_rintra/idlist.txt")

#~~ Create master results table

resfixed  <- NULL
resrandom <- NULL
respin    <- NULL
resVp     <- NULL
reswald   <- NULL

system(paste0("gcta-win-1.94.1.exe --bfile 70k_data/70K_200K_maf_geno_mind_v5 --autosome --autosome-num 33 --chr ", h, " --keep chr_h2_rintra/idlist.txt --make-grm-gz --out chr_h2_rintra/chr_", h))
system(paste0("gcta-win-1.94.1.exe --grm-gz chr_h2_rintra/chr_", h, " --grm-adj 0 --make-grm-gz --out chr_h2_rintra/chr_", h, ".adj"))

grm.auto <- read.table(paste0("chr_h2_rintra/chr_", h, ".adj.grm.gz"))
ids.auto <- read.table(paste0("chr_h2_rintra/chr_", h, ".adj.grm.id"))

grmreg <- makeGRM(grm.auto, ids.auto, id.vector = recsumm$parent)

system(paste0("rm chr_h2_rintra/chr_", h, ".grm.id"))
system(paste0("rm chr_h2_rintra/chr_", h, ".grm.gz"))
system(paste0("rm chr_h2_rintra/chr_", h, ".adj.grm.id"))
system(paste0("rm chr_h2_rintra/chr_", h, ".adj.grm.gz"))

rm(grm.auto, ids.auto)

writeLines(subset(snptab, V1 == h)$V2, paste0("chr_h2_rintra/chr_", h, "_snps.txt"))

system(paste0("gcta-win-1.94.1.exe --bfile 70k_data/70K_200K_maf_geno_mind_v5 --autosome --autosome-num 33 --exclude chr_h2_rintra/chr_", h, "_snps.txt --keep chr_h2_rintra/idlist.txt --make-grm-gz --out chr_h2_rintra/chr_", h))
system(paste0("gcta-win-1.94.1.exe --grm-gz chr_h2_rintra/chr_", h, " --grm-adj 0 --make-grm-gz --out chr_h2_rintra/chr_", h, ".adj"))

grm.auto <- read.table(paste0("chr_h2_rintra/chr_", h, ".adj.grm.gz"))
ids.auto <- read.table(paste0("chr_h2_rintra/chr_", h, ".adj.grm.id"))

grminv <- makeGRM(grm.auto, ids.auto, id.vector = recsumm$parent)

system(paste0("rm chr_h2_rintra/chr_", h, ".grm.id"))
system(paste0("rm chr_h2_rintra/chr_", h, ".grm.gz"))
system(paste0("rm chr_h2_rintra/chr_", h, ".adj.grm.id"))
system(paste0("rm chr_h2_rintra/chr_", h, ".adj.grm.gz"))

rm(grm.auto, ids.auto)
gc()

attr(grmreg, which = "INVERSE") <- TRUE
attr(grminv, which = "INVERSE") <- TRUE

#~~ Calculate trans r intra

y <- subset(shufftab, chrom != h)

y <- y %>% group_by(meiosis) %>% summarise(trans_intra_shuff_no_micro = sum(IS_no_micro))

recsumm <- left_join(recsumm, y)

#~~ RUN MODELS

for(i in 1:nrow(modeltab)){
  
  x <- recsumm[,c(as.character(modeltab$RespVar[i]), "parent", "sex", "Total_Coverage", "trans_intra_shuff_no_micro")]
  x$Total_Coverage2 <- x$Total_Coverage^2
  
  x$parent2 <- x$parent
  
  x <- subset(x, sex == modeltab$Sex[i])
  
  x <- na.omit(x) %>% droplevels
  
  eval(parse(text = paste0("model1 <- asreml(fixed = ", modeltab$RespVar[i], " ~ Total_Coverage + Total_Coverage2,
                      random = ~ vm(parent, grminv) + vm(parent2, grmreg) + ide(parent),
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
  
  x1 <- data.frame(RespVar = modeltab$RespVar[i], Sex = modeltab$Sex[i], Chromosome = h, Model = "cistrans")
  
  xfixed  <- cbind(x1, xfixed)
  xrandom <- cbind(x1, xrandom)
  xpin    <- cbind(x1, xpin)
  xVp     <- cbind(x1, xVp)
  xwald   <- cbind(x1, xwald)
  
  resfixed  <- rbind(resfixed, xfixed)
  resrandom <- rbind(resrandom, xrandom)
  respin    <- rbind(respin, xpin)
  resVp     <- rbind(resVp, xVp)
  reswald   <- rbind(reswald, xwald)
  
  eval(parse(text = paste0("model1 <- asreml(fixed = trans_intra_shuff_no_micro ~ Total_Coverage + ", ifelse(modeltab$Sex[i] == "Both", "sex + ", ""), "Total_Coverage2,
                      random = ~ vm(parent, grminv) + vm(parent2, grmreg) + ide(parent),
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
  
  x1 <- data.frame(RespVar = modeltab$RespVar[i], Sex = modeltab$Sex[i], Chromosome = h, Model = "trans")
  
  xfixed  <- cbind(x1, xfixed)
  xrandom <- cbind(x1, xrandom)
  xpin    <- cbind(x1, xpin)
  xVp     <- cbind(x1, xVp)
  xwald   <- cbind(x1, xwald)
  
  resfixed  <- rbind(resfixed, xfixed)
  resrandom <- rbind(resrandom, xrandom)
  respin    <- rbind(respin, xpin)
  resVp     <- rbind(resVp, xVp)
  reswald   <- rbind(reswald, xwald)
  
  rm(model1, xfixed, xrandom, xpin, xVp, x, xwald)
  
  save(resfixed, resrandom, respin, resVp, reswald, file = paste0("chr_h2_rintra/chr_", h, "_intra.RData"))
}