

library(dplyr)
library(asreml)

load("../reg_h2/6_Data_for_reg_h2_TRANS.RData")

snptab <- read.table("../70k_data/70K_200K_maf_geno_mind_v5.bim")

#~~ Create master results table

system(paste0("./gcta-1.94.1 --bfile ../70k_data/70K_200K_maf_geno_mind_v5 --autosome --autosome-num 33 --chr ", h, " --keep idlist.txt --make-grm-gz --out chr_", h))
system(paste0("./gcta-1.94.1 --grm-gz chr_", h, " --grm-adj 0 --make-grm-gz --out chr_", h, ".adj"))

grm.auto <- read.table(paste0("chr_", h, ".adj.grm.gz"))
ids.auto <- read.table(paste0("chr_", h, ".adj.grm.id"))

grmreg <- makeGRM(grm.auto, ids.auto, id.vector = recsumm$parent)
attr(grmreg, which = "INVERSE") <- TRUE

save(grmreg, file = paste0("chr_", h, "_Region.RData"))

system(paste0("rm chr_", h, ".grm.id"))
system(paste0("rm chr_", h, ".grm.gz"))
system(paste0("rm chr_", h, ".adj.grm.id"))
system(paste0("rm chr_", h, ".adj.grm.gz"))

rm(grm.auto, ids.auto)

writeLines(subset(snptab, V1 == h)$V2, paste0("chr_", h, "_snps.txt"))

system(paste0("./gcta-1.94.1 --bfile ../70k_data/70K_200K_maf_geno_mind_v5 --autosome --autosome-num 33 --exclude chr_", h, "_snps.txt --keep idlist.txt --make-grm-gz --out chr_", h))
system(paste0("./gcta-1.94.1 --grm-gz chr_", h, " --grm-adj 0 --make-grm-gz --out chr_", h, ".adj"))

grm.auto <- read.table(paste0("chr_", h, ".adj.grm.gz"))
ids.auto <- read.table(paste0("chr_", h, ".adj.grm.id"))

grminv <- makeGRM(grm.auto, ids.auto, id.vector = recsumm$parent)
attr(grminv, which = "INVERSE") <- TRUE

save(grminv, file = paste0("chr_", h, "_Genome.RData"))

system(paste0("rm chr_", h, ".grm.id"))
system(paste0("rm chr_", h, ".grm.gz"))
system(paste0("rm chr_", h, ".adj.grm.id"))
system(paste0("rm chr_", h, ".adj.grm.gz"))

rm(grm.auto, ids.auto)
gc()

