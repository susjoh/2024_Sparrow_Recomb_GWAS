#
# Statistics and sanity check
# SEJ
#

library(dplyr)
library(ggplot2)

load("prev_data/1_Cleaned_Sparrow_Yapp_Data.RData")

rm(ped, recchr, rectab, linkmap)
recsumm_old <- recsumm

load("results/2_Full_QC_Sparrow_Yapp_Data.RData")

recsumm <- read.table("results/3_Full_Recombination_Phenotypes_QCed.txt", header = T)
recsumm <- subset(recsumm, !is.na(yapp_CO_count_QCed))

fam180 <- read.table("prev_data/Pdo_200k_n3960_21032017.fam")


ggplot(recsumm, aes(yapp_kindepth, yapp_CO_count_QCed)) +
  geom_jitter(alpha = 0.01, width = 0.2, height = 0) +
  stat_smooth()

ggplot(recsumm, aes(factor(yapp_kindepth), yapp_CO_count_QCed)) +
  geom_boxplot() +
  stat_smooth()

ggplot(recsumm, aes(factor(ParentCount), yapp_CO_count_QCed)) +
  geom_boxplot(notch = T)

ggplot(recsumm, aes(factor(yapp_kindepth), yapp_CO_count_QCed, col = factor(ParentCount))) +
  geom_boxplot() +
  stat_smooth()

recsumm <- add_count(recsumm, parent, name = "meiosis_count")

ggplot(recsumm, aes(meiosis_count, yapp_CO_count_QCed)) +
  geom_point(alpha = 0.01) +
  stat_smooth(method = "lm")

ggplot(recsumm, aes(F, yapp_CO_count_QCed)) +
  geom_point(alpha = 0.01) +
  stat_smooth(method = "lm")

recsumm$GenoFlag <- ifelse(recsumm$parent %in% fam180$V2, "180k", "70k")
ggplot(recsumm, aes(factor(GenoFlag), yapp_CO_count_QCed)) +
  geom_boxplot() +
  stat_smooth()
table(recsumm$GenoFlag)

recsumm$GenoFlag_offspring <- ifelse(recsumm$offspring %in% fam180$V2, "180k", "70k")
ggplot(recsumm, aes(factor(GenoFlag_offspring), yapp_CO_count_QCed)) +
  geom_boxplot() +
  stat_smooth()

#~~~~~~~ Make a figure of the #COs per chromosome

str(recchr_v2)
recchr_v2$chrom[which(recchr_v2$chrom == 29)] <- "1A"
recchr_v2$chrom <- factor(recchr_v2$chrom, levels = c(1, "1A", 2:15, 17:25, 27, 28))
recchr_v2$n <- factor(recchr_v2$n, levels = 8:0)


ggplot(recchr_v2, aes(chrom, fill = n)) +
  geom_bar() +
  geom_hline(yintercept = nrow(recsumm)/2, colour = "white", linetype = "dashed") +
  theme_bw() +
  scale_fill_brewer(palette = "GnBu") +
  labs(x = "Chromosome", y = "Proportion of Gametes", fill = "CO Count") +
  scale_y_continuous(breaks = seq(0, nrow(recsumm), nrow(recsumm)/5), labels = seq(0, 1, 0.2))
ggsave("figs/4_Distribution_of_CO_counts_per_Chromosome.png", width = 6, height = 4)

snppos <- read.table("70k_data/70K_200K_maf_geno_mind_v5.bim")
snppos <- subset(snppos, V1 != 26)
head(snppos)
snppos$V1[which(snppos$V1 == 29)] <- "1A"
snppos$V1 <- factor(snppos$V1, levels = c(1, "1A", 2:15, 17:25, 27, 28))

x <- snppos[sample(nrow(snppos), 1000),]

ggplot(snppos, aes(V1, V4/1e6)) +
  geom_point(shape = 3, alpha = 0.1) +
  theme_bw() +
  labs(x = "Chromosome", y = "Chromosome Position (Mb)")
ggsave("figs/4_Distribution_of_SNPs_per_Chromosome.png", width = 6, height = 4)

#~~ Make figures of CO count per chromosome

recsumm$sex2 <- ifelse(recsumm$sex == "M", "Male", "Female")
ggplot(recsumm, aes(sex2, yapp_CO_count_no_micro, fill = sex)) +
  geom_boxplot(outlier.alpha = 0.2, width = 0.6) +
  theme_bw() +
  labs(x = "Sex", y = "Autosomal Crossover Count") +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.position = "none")
ggsave("figs/4_ACC_by_Sex.png", width = 2.5, height = 2.5)

ggplot(recsumm, aes(sex2, yapp_CO_count_no_micro, col = sex)) +
  geom_jitter(width = 0.2, alpha = 0.1) +
  geom_boxplot(colour = "black", alpha = 0, outlier.alpha = 0.2, width = 0.5) +
  theme_bw() +
  labs(x = "Sex", y = "Autosomal Crossover Count") +
  scale_colour_brewer(palette = "Set1") +
  theme(legend.position = "none")


ggplot(recsumm, aes(sex2, intra_shuff_no_micro, fill = sex)) +
  geom_boxplot(outlier.alpha = 0.2) +
  theme_bw() +
  labs(x = "Sex", y = "Autosomal Crossover Count") +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.position = "none")  
ggsave("figs/4_Shuff_by_Sex.png", width = 2.5, height = 2.5)


#~~ How many matrilines?


x <- recsumm %>% subset(sex == "F") %>% group_by(matriline) %>% summarise(MatCount = n()) %>% arrange(-MatCount)
x$CumMat <- cumsum(x$MatCount)
x$Order <- 1:nrow(x)

ggplot(x, aes(Order, CumMat)) +
  geom_line() +
  geom_hline(yintercept = 0.9*6409, colour = "red", linetype = "dashed") +
  #geom_vline(xintercept = 150, colour = "red", linetype = "dashed") +
  theme_bw()+
  labs(x = "Matriline (Ordered from most to least common)", y ="Cumulative number of gametes")
