library(tidyverse)

dir("chr_h2_perm/")

x <- dir("chr_h2_perm/")
x <- x[grep("RData", x)]
x <- x[-grep("temp", x)]

fullpin <- NULL

for(i in x){
  load(paste0("chr_h2_perm/", i))
  fullpin <- rbind(fullpin, cbind(respin, File = i))
}

fullpin$File <- gsub("chr_", "", fullpin$File)
fullpin$File <- gsub(".RData", "", fullpin$File, fixed = T)
fullpin$File <- gsub("_perm_iter", "", fullpin$File, fixed = T)
fullpin <- separate(fullpin, File, c("Chr", "Iteration"), sep = "_", remove = T)
fullpin$Chr <- NULL

fullpin <- unique(fullpin)

snptab <- read.table("70k_data/70K_200K_maf_geno_mind.bim")
head(snptab)
chrtab <- snptab %>% group_by(V1) %>% summarise(Length = max(V4)/1e6)
names(chrtab)[1] <- "Chromosome"

x <- subset(fullpin, Effect == "vm(parent2, grmreg)")

x <- left_join(x, chrtab)

head(x)

x$Chromosome[which(x$Chromosome == 29)] <- "1A"

# x_hold <- x
# 
# x <- x_hold
# 
# x <- subset(x, Estimate > 0.001)

modeltab <- subset(x, select = c(RespVar, Sex, Model, Iteration)) %>% unique

fullslopes <- NULL

for(i in 1:nrow(modeltab)){
  
  fit1 <- summary(lm(Estimate ~ Length, data = filter(x, RespVar == modeltab$RespVar[i] &
                                                        Sex == modeltab$Sex[i] &
                                                        Model == modeltab$Model[i] &
                                                        Iteration == modeltab$Iteration[i],
                                                        Chromosome %in% c(1:20, "1A"))))
  fullslopes <- rbind(fullslopes, 
                      cbind(RespVar = modeltab$RespVar[i],
                            Sex = modeltab$Sex[i],
                            Model = modeltab$Model[i],
                            Iteration = modeltab$Iteration[i],
                            fit1$coefficients,
                            adjR2 = fit1$adj.r.squared,
                            Chromosomes = "no_micro"))
}

for(i in 1:nrow(modeltab)){
  
  fit1 <- summary(lm(Estimate ~ Length, data = filter(x, RespVar == modeltab$RespVar[i] &
                                                        Sex == modeltab$Sex[i] &
                                                        Model == modeltab$Model[i] &
                                                        Iteration == modeltab$Iteration[i])))
  fullslopes <- rbind(fullslopes, 
                      cbind(RespVar = modeltab$RespVar[i],
                            Sex = modeltab$Sex[i],
                            Model = modeltab$Model[i],
                            Iteration = modeltab$Iteration[i],
                            fit1$coefficients,
                            adjR2 = fit1$adj.r.squared,
                            Chromosomes = "all"))
}

fullslopes <- data.frame(fullslopes)
names(fullslopes) <- c("RespVar", "Sex", "Model", "Iteration", "Estimate", "SE", "t", "P", "AdjR2", "Analysis")
fullslopes$Effect <- c("Intercept", "Length")

for(i in 5:9) fullslopes[,i] <- as.numeric(fullslopes[,i])

fullslopes$RespVar <- ifelse(fullslopes$RespVar == "yapp_CO_no_micro", "ACC", "r_intra")

test <- subset(fullslopes, P < 0.05 & Effect == "Length")

write.table(fullslopes, "results/8_Chr_h2_regression_permutation_results.txt", row.names = F, sep = "\t", quote = F)

x$RespVar <- ifelse(x$RespVar == "yapp_CO_no_micro", "ACC", "r_intra")

x$Effect <- NULL

write.table(x, "results/8_Chr_h2_FULL_permutation_results.txt", row.names = F, sep = "\t", quote = F)

#~~~~~~~~~~ HOW TO PLOT?

truex <- read.table("results/7_Chr_h2_FULL_results.txt", header = T, sep = "\t")
truex$Iteration <- 0

x$Sex <- ifelse(x$Sex == "F", "Females", "Males")

ggplot(x, aes(Length,  Estimate, group = Iteration)) +
  geom_hline(yintercept = 0) +
  stat_smooth(method = "lm") +
  stat_smooth(data = truex, aes(Length,  Estimate), method = "lm", col = "red") +
  geom_errorbar(data = truex, aes(Length,  Estimate, ymin = Estimate - SE, ymax = Estimate + SE), alpha = 0.5) +
  geom_point(data = truex, aes(Length,  Estimate), size = 2, alpha = 0.7) +
  facet_grid(RespVar~Sex, scales = "fixed") +
  theme_bw() +
  theme(strip.text = element_text(size = 10))+
  labs(x = "Chromosome Length (Mb)", y = "Proportion of Phenotypic Variance")
ggsave("figs/8_Chr_Partitioning_FULL_landscape.png", width = 8, height = 5)

