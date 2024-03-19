
library(dplyr)
library(ggplot2)

#~~ Parse the chromosome partitioning

x <- dir("chr_h2_rintra")
x <- x[grep("_rintra.RData", x)]

fullfixed <- list()
fullpin <- list()
fullrandom <- list()
fullVp <- list()
fullwald <- list()

for(i in x){
  
  vecno <- length(fullfixed)+1
  load(paste0("chr_h2_rintra/", i))
  fullfixed[[vecno]] <- resfixed
  fullpin[[vecno]] <- respin
  fullrandom[[vecno]] <- resrandom
  fullVp[[vecno]] <- resVp
  fullwald[[vecno]] <- reswald
  
  rm(resfixed, respin, resrandom, resVp, reswald)
  gc()
}

fullfixed <- bind_rows(fullfixed)
fullpin <- bind_rows(fullpin)
fullrandom <- bind_rows(fullrandom)
fullVp <- bind_rows(fullVp)
fullwald <- bind_rows(fullwald)

write.table(fullfixed, "results/8_Chr_h2_fixed_effects_rintra.txt", row.names = F, sep = "\t", quote = F)
write.table(fullpin, "results/8_Chr_h2_prop_of_variance_rintra.txt", row.names = F, sep = "\t", quote = F)
write.table(fullrandom, "results/8_Chr_h2_random_effects_rintra.txt", row.names = F, sep = "\t", quote = F)
write.table(fullVp, "results/8_Chr_h2_Vp_rintra.txt", row.names = F, sep = "\t", quote = F)
write.table(fullwald, "results/8_Chr_h2_wald_stats_rintra.txt", row.names = F, sep = "\t", quote = F)

