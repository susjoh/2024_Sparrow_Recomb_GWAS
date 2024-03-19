
library(dplyr)
library(ggplot2)

#~~ Parse the chromosome partitioning

x <- dir("chr_h2")
x <- x[grep(".RData", x)]

fullfixed <- list()
fullpin <- list()
fullrandom <- list()
fullVp <- list()
fullwald <- list()

for(i in x){
  
  vecno <- length(fullfixed)+1
  load(paste0("chr_h2/", i))
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

write.table(fullfixed, "results/8_Chr_h2_fixed_effects.txt", row.names = F, sep = "\t", quote = F)
write.table(fullpin, "results/8_Chr_h2_prop_of_variance.txt", row.names = F, sep = "\t", quote = F)
write.table(fullrandom, "results/8_Chr_h2_random_effects.txt", row.names = F, sep = "\t", quote = F)
write.table(fullVp, "results/8_Chr_h2_Vp.txt", row.names = F, sep = "\t", quote = F)
write.table(fullwald, "results/8_Chr_h2_wald_stats.txt", row.names = F, sep = "\t", quote = F)


fullfixed <- read.table("results/8_Chr_h2_fixed_effects.txt", header = T, sep = "\t")
fullpin <- read.table("results/8_Chr_h2_prop_of_variance.txt", header = T, sep = "\t")
fullrandom <- read.table("results/8_Chr_h2_random_effects.txt", header = T, sep = "\t")
fullVp <- read.table("results/8_Chr_h2_Vp.txt", header = T, sep = "\t")
fullwald <- read.table("results/8_Chr_h2_wald_stats.txt", header = T, sep = "\t")


snptab <- read.table("70k_data/70K_200K_maf_geno_mind.bim")
head(snptab)
chrtab <- snptab %>% group_by(V1) %>% summarise(Length = max(V4)/1e6)
names(chrtab)[1] <- "Chromosome"

x <- subset(fullpin, Effect == "vm(parent2, grmreg)")
x <- left_join(x, chrtab)

x1 <- summary(lm(Estimate ~ Length, data = subset(x, Model == "cistrans" & Sex == "F")))
x2 <- summary(lm(Estimate ~ Length, data = subset(x, Model == "trans" & Sex == "F")))
x3 <- summary(lm(Estimate ~ Length, data = subset(x, Model == "cistrans" & Sex == "M")))
x4 <- summary(lm(Estimate ~ Length, data = subset(x, Model == "trans" & Sex == "M")))

# y = 

x$Sex <- ifelse(x$Sex == "M", "Males", "Females")
x$Model <- ifelse(x$Model == "cistrans", "All Chromosomes (cis & trans)", "Excluding Chromosome (trans)")
x$Chromosome[which(x$Chromosome == 29)] <- "1A"



ggplot(x, aes(Length,  Estimate)) +
  stat_smooth(method = "lm") +
  geom_errorbar(aes(ymin = Estimate - SE, ymax = Estimate + SE)) +
  geom_point(size = 2) +
  #geom_text(aes(label = Chromosome), size = 2, col = "white") +
  facet_grid(Model ~ Sex, scales = "free_y") +
  theme_bw() +
  labs(x = "Chromosome Length (Mb)", y = "Proportion of Phenotypic Variance")

ggplot(subset(x, Chromosome %in% c(1:20, "1A")), aes(Length,  Estimate)) +
  stat_smooth(method = "lm") +
  geom_errorbar(aes(ymin = Estimate - SE, ymax = Estimate + SE)) +
  geom_point(size = 4) +
  geom_text(aes(label = Chromosome), size = 2, col = "white") +
  facet_grid(Sex ~ Model, scales = "free_y") +
  theme_bw() +
  labs(x = "Chromosome Length (Mb)", y = "Heritability")
ggsave("figs/8_Chromosome_Partitioning.png", width = 6, height = 4)



x <- subset(fullpin, Effect == "vm(parent2, grmreg)" & Chromosome %in% c(1:20, 29))
x <- left_join(x, chrtab)
summary(lm(Estimate ~ Length, data = subset(x, Model == "cistrans" & Sex == "F")))
summary(lm(Estimate ~ Length, data = subset(x, Model == "trans" & Sex == "F")))
summary(lm(Estimate ~ Length, data = subset(x, Model == "cistrans" & Sex == "M")))
summary(lm(Estimate ~ Length, data = subset(x, Model == "trans" & Sex == "M")))


ggplot(subset(x, Chromosome %in% c(1:20, "1A")& Model == "Excluding SNP Chromosome (trans)"), aes(Length,  Estimate)) +
  #stat_smooth(method = "lm") +
  geom_errorbar(aes(ymin = Estimate - SE, ymax = Estimate + SE), alpha = 0) +
  geom_point(size = 4, alpha = 0) +
  geom_text(aes(label = Chromosome), size = 2, col = "white",alpha = 0) +
  facet_wrap( ~ Sex, scales = "free_y") +
  theme_bw() +
  labs(x = "Chromosome Length (Mb)", y = "Heritability")
ggsave("figs/8_Chromosome_Partitioning_blank.png", width = 5, height = 3)

ggplot(subset(x, Chromosome %in% c(1:20, "1A")& Model == "Excluding SNP Chromosome (trans)"), aes(Length,  Estimate)) +
  stat_smooth(method = "lm") +
  geom_errorbar(aes(ymin = Estimate - SE, ymax = Estimate + SE)) +
  geom_point(size = 4) +
  geom_text(aes(label = Chromosome), size = 2, col = "white") +
  facet_wrap( ~ Sex, scales = "free_y") +
  theme_bw() +
  labs(x = "Chromosome Length (Mb)", y = "Heritability")
ggsave("figs/8_Chromosome_Partitioning_trans.png", width = 5, height = 3)








x <- subset(fullpin, Effect == "vm(parent2, grmreg)")
x <- left_join(x, chrtab)

x$Sex <- ifelse(x$Sex == "M", "B. Males", "A. Females")
x$Chromosome[which(x$Chromosome == 29)] <- "1A"


ggplot(subset(x, Model == "trans"), aes(Length,  Estimate)) +
  stat_smooth(method = "lm") +
  geom_errorbar(aes(ymin = Estimate - SE, ymax = Estimate + SE)) +
  geom_point(size = 2) +
  #geom_text(aes(label = Chromosome), size = 2, col = "white") +
  geom_hline(yintercept = 0) +
  facet_wrap(~Sex, scales = "free") +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(hjust = 0, size = 12)) +
  scale_y_continuous(breaks = seq(0, 0.12, 0.02)) +
  coord_cartesian(ylim = c(-0.012, 0.11)) +
  labs(x = "Chromosome Length (Mb)", y = "Proportion of phenotypic variance")
ggsave("figs/7_Chromosome_Partitioning_for_paper.png", width = 8, height = 4)

summary(lm(Estimate ~ Length, data = subset(x, Model == "trans" & Sex == "Females")))
summary(lm(Estimate ~ Length, data = subset(x, Model == "trans" & Sex == "Males")))


ggsave("figs/8_Chromosome_Partitioning_blank.png", width = 5, height = 3)

ggplot(subset(x, Chromosome %in% c(1:20, "1A")& Model == "Excluding SNP Chromosome (trans)"), aes(Length,  Estimate)) +
  stat_smooth(method = "lm") +
  geom_errorbar(aes(ymin = Estimate - SE, ymax = Estimate + SE)) +
  geom_point(size = 4) +
  geom_text(aes(label = Chromosome), size = 2, col = "white") +
  facet_wrap( ~ Sex, scales = "free_y") +
  theme_bw() +
  labs(x = "Chromosome Length (Mb)", y = "Heritability")
ggsave("figs/8_Chromosome_Partitioning_trans.png", width = 5, height = 3)

