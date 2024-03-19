library(tidyverse)


fullpin <- rbind(read.table("results/8_Chr_h2_prop_of_variance.txt", header = T, sep = "\t"),
                 read.table("results/8_Chr_h2_prop_of_variance_rintra.txt", header = T, sep = "\t"))
fullpin <- unique(fullpin)
fullpin <- subset(fullpin)


snptab <- read.table("70k_data/70K_200K_maf_geno_mind.bim")
head(snptab)
chrtab <- snptab %>% group_by(V1) %>% summarise(Length = max(V4)/1e6)
names(chrtab)[1] <- "Chromosome"

x <- subset(fullpin, Effect == "vm(parent2, grmreg)")

x <- left_join(x, chrtab)

head(x)

x$Sex <- ifelse(x$Sex == "M", "Males", "Females")
x$Model <- ifelse(x$Model == "cistrans", "All Chromosomes (cis & trans)", "Excluding Chromosome (trans)")
x$Chromosome[which(x$Chromosome == 29)] <- "1A"

x_hold <- x

x <- x_hold

x <- subset(x, Estimate > 0.001)

modeltab <- subset(x, select = c(RespVar, Sex, Model)) %>% unique

fullslopes <- NULL

for(i in 1:nrow(modeltab)){
  
  fit1 <- summary(lm(Estimate ~ Length, data = filter(x, RespVar == modeltab$RespVar[i] &
                                                        Sex == modeltab$Sex[i] &
                                                        Model == modeltab$Model[i] &
                                                        Chromosome %in% c(1:20, "1A"))))
  fullslopes <- rbind(fullslopes, 
                      cbind(RespVar = modeltab$RespVar[i],
                            Sex = modeltab$Sex[i],
                            Model = modeltab$Model[i],
                            fit1$coefficients,
                            adjR2 = fit1$adj.r.squared,
                            Chromosomes = "no_micro"))
}

for(i in 1:nrow(modeltab)){
  
  fit1 <- summary(lm(Estimate ~ Length, data = filter(x, RespVar == modeltab$RespVar[i] &
                                                        Sex == modeltab$Sex[i] &
                                                        Model == modeltab$Model[i])))
  fullslopes <- rbind(fullslopes, 
                      cbind(RespVar = modeltab$RespVar[i],
                            Sex = modeltab$Sex[i],
                            Model = modeltab$Model[i],
                            fit1$coefficients,
                            adjR2 = fit1$adj.r.squared,
                            Chromosomes = "all"))
}

fullslopes <- data.frame(fullslopes)
names(fullslopes) <- c("RespVar", "Sex", "Model", "Estimate", "SE", "t", "P", "AdjR2", "Analysis")
fullslopes$Effect <- c("Intercept", "Length")

for(i in 4:8) fullslopes[,i] <- as.numeric(fullslopes[,i])

fullslopes$RespVar <- ifelse(fullslopes$RespVar == "yapp_CO_no_micro", "ACC", "r_intra")

test <- subset(fullslopes, P < 0.05 & Effect == "Length")

write.table(fullslopes, "results/7_Chr_h2_regression_results.txt", row.names = F, sep = "\t", quote = F)

x$RespVar <- ifelse(x$RespVar == "yapp_CO_no_micro", "ACC", "r_intra")

x$Effect <- NULL
head(x)

write.table(x, "results/7_Chr_h2_FULL_results.txt", row.names = F, sep = "\t", quote = F)

#~~~~~~~~~~ HOW TO PLOT?



x$RespVar2 <- ifelse(x$RespVar == "ACC", "Autosomal Crossover Count (ACC)", "Intrachromosomal Shuffling (r_intra)")

x$RespModel <- paste0(x$RespVar, ": ", x$Model)
x$RespSex <- paste0(x$RespVar, ": ", x$Sex)

ggplot(x, aes(Length,  Estimate)) +
  geom_hline(yintercept = 0) +
  stat_smooth(method = "lm") +
  geom_errorbar(aes(ymin = Estimate - SE, ymax = Estimate + SE), alpha = 0.5) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text(aes(label = Chromosome), size = 2, col = "white") +
  facet_grid(Model~RespSex, scales = "fixed") +
  theme_bw() +
  theme(strip.text = element_text(size = 10))+
  labs(x = "Chromosome Length (Mb)", y = "Proportion of Phenotypic Variance")
ggsave("figs/7_Chr_Partitioning_FULL_landscape.png", width = 10, height = 6)

fullslopes$RespSex <- paste0(fullslopes$RespVar, ": ", fullslopes$Sex)
fullslopes$Slope <- paste0("slope = ", round(fullslopes$Estimate, 5))
fullslopes$Intercept <- paste0("intercept = ", round(fullslopes$Estimate, 5))

ggplot() +
  geom_hline(yintercept = 0) +
  stat_smooth(data = x, aes(Length,  Estimate), method = "lm") +
  geom_errorbar(data = x, aes(Length,  Estimate, ymin = Estimate - SE, ymax = Estimate + SE), alpha = 0.5) +
  geom_point(data = x, aes(Length,  Estimate), size = 2, alpha = 0.5) +
  geom_text(data = subset(fullslopes, Effect == "Length"),
            aes(x = 10, y = 0.13, label = Slope), size = 3.5, hjust = 0) +
  geom_text(data = subset(fullslopes, Effect == "Intercept"),
            aes(x = 10, y = 0.11, label = Intercept), size = 3.5, hjust = 0) +
  facet_grid(RespSex ~ Model, scales = "free_y") +
  coord_cartesian(ylim = c(-0.01, 0.14)) +
  theme_bw() +
  labs(x = "Chromosome Length (Mb)", y = "Proportion of Phenotypic Variance")
ggsave("figs/7_Chr_Partitioning_FULL_stats.png", width = 6, height = 10)

#~~ Make a better table for the paper

fullpin <- rbind(read.table("results/8_Chr_h2_prop_of_variance.txt", header = T, sep = "\t"),
                 read.table("results/8_Chr_h2_prop_of_variance_rintra.txt", header = T, sep = "\t"))
fullpin$RespVar <- ifelse(fullpin$RespVar == "yapp_CO_no_micro", "ACC", "r_intra")

fullpin$Effect2 <- NA
fullpin$Effect2[which(fullpin$Effect == "vm(parent, grminv)")] <- "Additive Genetic"
fullpin$Effect2[which(fullpin$Effect == "vm(parent2, grmreg)")] <- "Additive Genetic (Chromosome)"
fullpin$Effect2[which(fullpin$Effect == "ide(parent)")] <- "Permanent Environment"
fullpin$Effect2[which(fullpin$Effect == "units!units")] <- "Residual"

fullpin$Effect <- fullpin$Effect2
fullpin$Effect2 <- NULL

fullpin <- unique(fullpin)
fullpin <- left_join(fullpin, chrtab)
head(fullpin)

write.table(fullpin, "results/7_Chr_h2_FULL_results_for_paper.txt", row.names = F, sep = "\t", quote = F)
