
library(dplyr)
library(ggplot2)
library(GenABEL)
library(RepeatABEL)

source("r/rGLSadj.R")

load("clean_data_for_John/70K_200K_maf_geno_mind_v5.GRM.RData")
load("clean_data_for_John/70K_200K_maf_geno_mind.RData")

recsumm <- read.table("clean_data_for_John/3_Full_Recombination_Phenotypes_QCed.txt", header = T, stringsAsFactors = F)
recsumm$id <- recsumm$parent

recsumm <- subset(recsumm, select = c(id, parent, meiosis, sex, Total_Coverage,
                                      yapp_CO_count_QCed, yapp_CO_count_no_micro, 
                                      yapp_CO_count_macro, intra_shuff, intra_shuff_macro, 
                                      intra_shuff_no_micro, intra_shuff_gene, intra_shuff_gene_macro, 
                                      intra_shuff_gene_no_micro))

recsumm$parent <- factor(recsumm$parent)
recsumm$sex <- factor(recsumm$sex)

grminv2 <- grminv
grminv2[lower.tri(grminv2)] = t(grminv2)[lower.tri(grminv2)]


firstRun <- FALSE

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. Run the GWASs                                     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if(firstRun){
  
  # Female COs
  
  x <- subset(recsumm, sex == "F") %>% droplevels %>% na.omit
  
  COf_prefit <- preFitModel(fixed = yapp_CO_count_no_micro ~ Total_Coverage,
                            id.name = "id",
                            genabel.data = sparabel70,
                            phenotype.data = x, 
                            corStruc = list(id = list("GRM")),
                            GRM = grminv2)
  
  COf_gwas <- rGLSadj(formula.FixedEffects=yapp_CO_count_no_micro~Total_Coverage,
                      genabel.data=sparabel70,
                      phenotype.data=x,
                      id="id",
                      V = COf_prefit$V,
                      GRM=grminv2)
  
  COf_gwas <- process_rGLSadj_results(COf_gwas, sparabel70)
  
  save(COf_gwas, file = "results/6_COf_gwas.RData")
  rm(COf_gwas, COf_prefit)
  gc()
  
  # Male COs
  
  x <- subset(recsumm, sex == "M") %>% droplevels %>% na.omit
  
  COm_prefit <- preFitModel(fixed = yapp_CO_count_no_micro ~ 1+Total_Coverage,
                            id.name = "id",
                            genabel.data = sparabel70,
                            phenotype.data = x, 
                            corStruc = list(id = list("GRM")),
                            GRM = grminv2)
  
  COm_gwas <- rGLSadj(formula.FixedEffects=yapp_CO_count_no_micro~Total_Coverage,
                      genabel.data=sparabel70,
                      phenotype.data=x,
                      id="id",
                      V = COm_prefit$V,
                      GRM=grminv2)
  
  COm_gwas <- process_rGLSadj_results(COm_gwas, sparabel70)
  
  save(COm_gwas, file = "results/6_COm_gwas.RData")
  
  rm(COm_gwas, COm_prefit)
  gc()
  
  # Female shuffling
  
  x <- subset(recsumm, sex == "F") %>% droplevels %>% na.omit
  
  ISfuncor_prefit <- preFitModel(fixed = intra_shuff_no_micro ~ 1+Total_Coverage,
                            id.name = "id",
                            genabel.data = sparabel70,
                            phenotype.data = x, 
                            corStruc = list(id = list("GRM")),
                            GRM = grminv2)
  
  ISfuncor_gwas <- rGLSadj(formula.FixedEffects=intra_shuff_no_micro~Total_Coverage,
                      genabel.data=sparabel70,
                      phenotype.data=x,
                      id="id",
                      V = ISfuncor_prefit$V,
                      GRM=grminv2)
  
  ISfuncor_gwas <- process_rGLSadj_results(ISfuncor_gwas, sparabel70)
  
  save(ISfuncor_gwas, file = "results/6_ISfuncor_gwas.RData")
  #rm(ISfuncor_gwas, ISfuncor_prefit)
  gc()
  
  # Male shuffling
  
  x <- subset(recsumm, sex == "M") %>% droplevels %>% na.omit
  
  ISmuncor_prefit <- preFitModel(fixed = intra_shuff_no_micro ~ 1+Total_Coverage,
                            id.name = "id",
                            genabel.data = sparabel70,
                            phenotype.data = x, 
                            corStruc = list(id = list("GRM")),
                            GRM = grminv2)
  
  ISmuncor_gwas <- rGLSadj(formula.FixedEffects=intra_shuff_no_micro~Total_Coverage,
                      genabel.data=sparabel70,
                      phenotype.data=x,
                      id="id",
                      V = ISmuncor_prefit$V,
                      GRM=grminv2)
  
  ISmuncor_gwas <- process_rGLSadj_results(ISmuncor_gwas, sparabel70)
  
  save(ISmuncor_gwas, file = "results/6_ISmuncor_gwas.RData")
  
  #rm(ISmuncor_gwas, ISmuncor_gwas)
  gc()
  
  # Female shuffling corrected
  
  x <- subset(recsumm, sex == "F") %>% droplevels %>% na.omit
  
  ISf_prefit <- preFitModel(fixed = intra_shuff_no_micro ~ 1+Total_Coverage + yapp_CO_count_no_micro,
                            id.name = "id",
                            genabel.data = sparabel70,
                            phenotype.data = x, 
                            corStruc = list(id = list("GRM")),
                            GRM = grminv2)
  
  ISf_gwas <- rGLSadj(formula.FixedEffects=intra_shuff_no_micro~Total_Coverage+ yapp_CO_count_no_micro,
                      genabel.data=sparabel70,
                      phenotype.data=x,
                      id="id",
                      V = ISf_prefit$V,
                      GRM=grminv2)
  
  ISf_gwas <- process_rGLSadj_results(ISf_gwas, sparabel70)
  
  save(ISf_gwas, file = "results/6_ISf_gwas_cofix.RData")
  rm(ISf_gwas, ISf_prefit)
  gc()
  
  # Male shuffling corrected
  
  x <- subset(recsumm, sex == "M") %>% droplevels %>% na.omit
  
  ISm_prefit <- preFitModel(fixed = intra_shuff_no_micro ~ 1+Total_Coverage+ yapp_CO_count_no_micro,
                            id.name = "id",
                            genabel.data = sparabel70,
                            phenotype.data = x, 
                            corStruc = list(id = list("GRM")),
                            GRM = grminv2)
  
  ISm_gwas <- rGLSadj(formula.FixedEffects=intra_shuff_no_micro~Total_Coverage+ yapp_CO_count_no_micro,
                      genabel.data=sparabel70,
                      phenotype.data=x,
                      id="id",
                      V = ISm_prefit$V,
                      GRM=grminv2)
  
  ISm_gwas <- process_rGLSadj_results(ISm_gwas, sparabel70)
  
  save(ISm_gwas, file = "results/6_ISm_gwas_cofix.RData")
  
  rm(ISm_gwas, ISm_gwas)
  gc()
  
  
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. Compile the results                              #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

load("results/6_COf_gwas.RData")
load("results/6_COm_gwas.RData")
load("results/6_ISfuncor_gwas.RData")
load("results/6_ISmuncor_gwas.RData")
load("results/6_ISm_gwas_cofix.RData")
load("results/6_ISf_gwas_cofix.RData")

gwas.results <- rbind(cbind(COf_gwas, Sex = "F", Model = "ACC"),
                      cbind(COm_gwas, Sex = "M", Model = "ACC"),
                      cbind(ISfuncor_gwas, Sex = "F", Model = "IntraR"),
                      cbind(ISmuncor_gwas, Sex = "M", Model = "IntraR"),
                      cbind(ISf_gwas, Sex = "F", Model = "IntraR_cofix"),
                      cbind(ISm_gwas, Sex = "M", Model = "IntraR_cofix"))

#~~ Get MAF information

idvec <- unique(recsumm$id) %>% as.character
idvec <- idvec[which(idvec %in% idnames(sparabel70))]

snpinfo <- summary.snp.data(gtdata(sparabel70[idvec,]))
str(snpinfo)
snpinfo$Chromosome <- as.numeric(as.character(snpinfo$Chromosome))

snpinfo$SNP.Name <- row.names(snpinfo)

snpinfo <- arrange(snpinfo, Chromosome, Position)
snpinfo$Diff <- c(1, diff(snpinfo$Position))
snpinfo$Diff[which(snpinfo$Diff < 0)] <- 10000
snpinfo$Cumu <- cumsum(snpinfo$Diff)


snpinfo <- subset(snpinfo, select = c(SNP.Name, Q.2, Cumu))
head(snpinfo)

gwas.results <- left_join(gwas.results, snpinfo)
table(gwas.results$Chromosome)
gwas.results$Chromosome <- as.numeric(as.character(gwas.results$Chromosome))


plot(-log10(Pc1df) ~ Cumu, col = (Chromosome %% 2)+1, data = subset(gwas.results, Sex == "F" & Model == "ACC"))
plot(-log10(Pc1df) ~ Cumu, col = (Chromosome %% 2)+1, data = subset(gwas.results, Sex == "M" & Model == "ACC"))
plot(-log10(Pc1df) ~ Cumu, col = (Chromosome %% 2)+1, data = subset(gwas.results, Sex == "F" & Model == "IntraR"))
plot(-log10(Pc1df) ~ Cumu, col = (Chromosome %% 2)+1, data = subset(gwas.results, Sex == "M" & Model == "IntraR"))
plot(-log10(Pc1df) ~ Cumu, col = (Chromosome %% 2)+1, data = subset(gwas.results, Sex == "F" & Model == "IntraR_cofix"))
plot(-log10(Pc1df) ~ Cumu, col = (Chromosome %% 2)+1, data = subset(gwas.results, Sex == "M" & Model == "IntraR_cofix"))


write.table(gwas.results, file = "results/6_GWAS_Results.txt", row.names = F, sep = "\t", quote = F)


COf_gwas <- arrange(COf_gwas, Pc1df)
head(COf_gwas)

spargenes <- read.table("prev_data/20221109_sparrow_genes_formatted.txt",  header = T, sep = "\t")

genehold <- NULL

for(i in 1:100){
  genehold <- rbind(genehold,
                    subset(spargenes, Chr == as.character(COf_gwas$Chromosome[i]) &
                             Start < (COf_gwas$Position[i]+2e5) &
                             Stop > (COf_gwas$Position[i]-2e5))) %>% unique
  
  
}

writeLines(unique(as.character(genehold$Gene)), "Female_GWAS.txt")



genehold <- NULL

for(i in 1:100){
  genehold <- rbind(genehold,
                    subset(spargenes, Chr == as.character(COm_gwas$Chromosome[i]) &
                             Start < (COm_gwas$Position[i]+1e5) &
                             Stop > (COm_gwas$Position[i]-1e5))) %>% unique
  
  
}

writeLines(unique(as.character(genehold$Gene)), "Male_GWAS.txt")
