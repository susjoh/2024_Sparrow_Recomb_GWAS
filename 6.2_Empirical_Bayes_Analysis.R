

library(ashr)
library(ggplot2)
library(dplyr)

gwasres <- read.table("results/6_GWAS_Results.txt", header = T, sep = "\t")

gwasres <- subset(gwasres, Model != "IntraR_cofix")

modeltab <- subset(gwasres, select = c(Sex, Model)) %>% unique

fullres <- list()

for(i in 1:nrow(modeltab)){
  
  x <- subset(gwasres, Sex == modeltab$Sex[i] & Model == modeltab$Model[i])
  
  measure_ash <- ash(x$effB,
                     x$se_effB,
                     mixcompdist = "uniform",
                     method= "fdr",
                     optmethod="mixEM", 
                     control = list(maxiter = 10000),
                     outputlevel= 3)
  
  ash_results <- measure_ash$result
  ash_results <- cbind(x, ash_results)
  
  fullres[[i]] <- ash_results
  rm(ash_results, measure_ash)
  
}

fullres <- bind_rows(fullres)

sigres <- subset(fullres,  lfsr < 0.05)
table(sigres$Model, sigres$Sex)

sigres1 <- subset(fullres,  lfsr < 0.01)
table(sigres1$Model, sigres1$Sex)

write.table(sigres, "results/6_Significant_SNPs_GWAS_EB.txt", row.names = F, sep = "\t", quote = F)
write.table(fullres, "results/6_GWAS_Results_EB.txt", row.names = F, sep = "\t", quote = F)

